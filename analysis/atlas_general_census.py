"""Census by exact rational arithmetic at random parameter points.

The coefficients A,B,C,D are polynomials in the rate constants.  Substituting
random rationals for every constant except the labeled one turns each symbolic
zero-test into an exact rational zero-test, which is enormously cheaper than
cancelling the full expression.  A polynomial that vanishes at several
independent random points is zero with overwhelming probability
(Schwartz-Zippel), and the classification is required to agree across draws;
disagreement is reported rather than resolved silently.
"""
import sys, json, time, collections
sys.path.insert(0,'/home/dong/Workspace/WritePaper/MDPI_Entropy_QuantumB/github-package/analysis')
import sympy as sp
from random import Random
from multiprocessing import Pool
import os
from atlas_general import enumerate_networks, build

x = sp.Symbol("x", positive=True)


def fp_first(m, subs):
    n = m.states
    A = sp.zeros(n, n); b = sp.zeros(n, 1)
    for i in range(1, n + 1):
        out = [(j, sy) for (a, j, sy) in m.edges if a == i]
        A[i-1, i-1] = sum(sy.subs(subs) for (_, sy) in out)
        for (j, sy) in out:
            if j == -1: b[i-1] += sy.subs(subs)
            elif j > 0:  A[i-1, j-1] -= sy.subs(subs)
    return A.LUsolve(b)[0]


def classify_network(arg):
    n, internal, diss, prod, seed = arg
    m0, labels = build(n, internal, diss, prod, 0)
    syms = [sy for (_, _, sy) in labels]
    rng = Random(seed)
    out, disagree = [], 0
    for e, lab_sym in enumerate(syms):
        votes = set()
        for trial in range(3):
            vals = {s: sp.Rational(rng.randint(2, 97), rng.randint(2, 97))
                    for s in syms if s is not lab_sym}
            kT = sp.Rational(rng.randint(2, 97), rng.randint(2, 97))
            f = fp_first(m0, vals)                       # f as a function of lab_sym
            num = sp.cancel(f.subs(lab_sym, x * kT))
            den = sp.cancel(f.subs(lab_sym, kT))
            K = sp.cancel(sp.together(num / den))
            pn, pd = sp.fraction(K)
            Pn, Pd = sp.Poly(sp.expand(pn), x), sp.Poly(sp.expand(pd), x)
            if Pn.degree() > 1 or Pd.degree() > 1:
                votes.add("not Moebius"); continue
            A = Pn.coeff_monomial(x); B = Pn.coeff_monomial(1)
            C = Pd.coeff_monomial(x); D = Pd.coeff_monomial(1)
            if A * D - B * C == 0: votes.add("blind"); continue
            votes.add("inverse (A=0)" if A == 0 else
                      "identity" if B == 0 and C == 0 else
                      "concave (B=0)" if B == 0 else
                      "convex (C=0)" if C == 0 else "mixed")
        if len(votes) > 1: disagree += 1
        out.append(sorted(votes)[0])
    return out, disagree


if __name__ == "__main__":
    nmax = int(sys.argv[1])
    jobs = []
    for n in range(1, nmax + 1):
        for k, (i, d, p) in enumerate(enumerate_networks(n)):
            jobs.append((n, i, d, p, 1000 * n + k))
    print(f"{len(jobs)} networks", flush=True)
    t0 = time.time(); cen = collections.Counter(); dis = 0
    with Pool(processes=max(1, os.cpu_count() - 2)) as pool:
        for k, (res, dd) in enumerate(pool.imap_unordered(classify_network, jobs, chunksize=2)):
            cen.update(res); dis += dd
            if (k + 1) % 100 == 0:
                print(f"  {k+1}/{len(jobs)}  {sum(cen.values())} edges  {time.time()-t0:.0f}s", flush=True)
    print()
    for kk, v in cen.most_common(): print(f"   {kk:16s} {v:6d}")
    print(f"   {'TOTAL':16s} {sum(cen.values()):6d}")
    print(f"   draws disagreeing: {dis}")
    print(f"   elapsed {time.time()-t0:.0f}s")
