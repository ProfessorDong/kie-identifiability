"""Atlas over general directed networks, not only chains.

network_atlas.py enumerates chain-shaped mechanisms with optional features.
That family is easy to describe but it is a restricted one, and the classification
it supports should not be quoted as though it covered every steady-state
mechanism.  This module removes the restriction for small networks: it
enumerates every directed graph on the bound states, keeps the ones that are
well posed as competitive experiments, and classifies each.

A network is well posed here when free enzyme can be reached (something
dissociates), product can be reached from the complex, and every bound state is
reachable from the complex, so that no state is inert.  Graphs are deduplicated
under relabeling of states 2..n; state 1 is distinguished as the complex formed
on binding.

The second question this settles is the one the chain atlas flagged as open:
isotope-sensitive steps in parallel, reached from a common state.

Run from analysis/:   python atlas_general.py [nmax]
"""
from __future__ import annotations

import itertools
from itertools import permutations

import numpy as np
import sympy as sp

import masses as M
from network_atlas import (Mech, blind_by_topology, classify, first_passage,
                           profile_geometry)

G = M.gamma_sc("C")
FREE, PROD = 0, -1


def _wellposed(n, internal, diss, prod):
    """Free enzyme reachable, product reachable from state 1, no inert state."""
    if not diss or not prod:
        return False
    adj = {i: set() for i in range(1, n + 1)}
    for (i, j) in internal:
        adj[i].add(j)
    # reachability from state 1 over internal edges
    seen, stack = {1}, [1]
    while stack:
        i = stack.pop()
        for j in adj[i]:
            if j not in seen:
                seen.add(j)
                stack.append(j)
    if seen != set(range(1, n + 1)):
        return False
    if not (set(prod) & seen):
        return False
    if not (set(diss) & seen):
        return False
    # first passage is defined only if every state can reach an absorbing one:
    # a bound state from which neither product nor free enzyme is reachable
    # traps the walk and leaves the linear system singular
    absorbing = set(diss) | set(prod)
    for start in range(1, n + 1):
        s2, st = {start}, [start]
        ok = start in absorbing
        while st and not ok:
            i = st.pop()
            for j in adj[i]:
                if j in absorbing:
                    ok = True; break
                if j not in s2:
                    s2.add(j); st.append(j)
        if not ok:
            return False
    return True


def _canon(n, internal, diss, prod):
    """Canonical form under relabeling of states 2..n, keeping state 1 fixed."""
    best = None
    for perm in permutations(range(2, n + 1)):
        m = {1: 1}
        m.update({old: new for old, new in zip(range(2, n + 1), perm)})
        key = (tuple(sorted((m[i], m[j]) for (i, j) in internal)),
               tuple(sorted(m[i] for i in diss)),
               tuple(sorted(m[i] for i in prod)))
        if best is None or key < best:
            best = key
    return best


def enumerate_networks(n):
    """Every well-posed directed network on n bound states, up to relabeling."""
    pairs = [(i, j) for i in range(1, n + 1) for j in range(1, n + 1) if i != j]
    states = list(range(1, n + 1))
    seen, out = set(), []
    for nint in range(len(pairs) + 1):
        for internal in itertools.combinations(pairs, nint):
            for dmask in range(1, 1 << n):
                diss = [s for k, s in enumerate(states) if dmask >> k & 1]
                for pmask in range(1, 1 << n):
                    prod = [s for k, s in enumerate(states) if pmask >> k & 1]
                    if not _wellposed(n, internal, diss, prod):
                        continue
                    key = _canon(n, internal, diss, prod)
                    if key in seen:
                        continue
                    seen.add(key)
                    out.append((internal, tuple(diss), tuple(prod)))
    return out


def build(n, internal, diss, prod, iso_edge):
    """A Mech with the isotope label on `iso_edge`, given as an edge index."""
    edges, labels = [], []
    for (i, j) in internal:
        s = sp.Symbol(f"k_{i}{j}", positive=True)
        edges.append((i, j, s)); labels.append(("int", (i, j), s))
    for i in diss:
        s = sp.Symbol(f"koff_{i}", positive=True)
        edges.append((i, FREE, s)); labels.append(("off", i, s))
    for i in prod:
        s = sp.Symbol(f"kcat_{i}", positive=True)
        edges.append((i, PROD, s)); labels.append(("cat", i, s))
    return Mech(name=f"n={n} {internal}|off{diss}|cat{prod}", states=n,
                edges=edges, iso=[labels[iso_edge][2]]), labels


def moebius_of(m):
    """(A,B,C,D) or None if V/K is not Moebius in the labeled constant."""
    k = m.iso[0]
    x, kT = sp.symbols("x k_T", positive=True)
    f = first_passage(m)
    K = sp.cancel(sp.together(sp.simplify(f.subs(k, x * kT) / f.subs(k, kT))))
    num, den = sp.fraction(K)
    pn, pd = sp.Poly(sp.expand(num), x), sp.Poly(sp.expand(den), x)
    if pn.degree() > 1 or pd.degree() > 1:
        return None
    return (pn.coeff_monomial(x), pn.coeff_monomial(1),
            pd.coeff_monomial(x), pd.coeff_monomial(1))


def draws_for(co, tag, k=6, seed=7):
    rng = np.random.default_rng(seed + len(tag))
    syms = sorted(set().union(*(z.free_symbols for z in co)), key=str)
    return [{s: float(v) for s, v in zip(syms, rng.lognormal(0, 1.0, len(syms)))}
            for _ in range(k)]


def run(nmax=3, verbose=False):
    """Classify every labeled edge of every well-posed network up to nmax states."""
    import collections
    cen = collections.Counter()
    notmoeb, mism = [], 0
    for n in range(1, nmax + 1):
        for (internal, diss, prod) in enumerate_networks(n):
            m0, labels = build(n, internal, diss, prod, 0)
            for e in range(len(labels)):
                m, _ = build(n, internal, diss, prod, e)
                if blind_by_topology(m):
                    cen["blind"] += 1
                    continue
                co = moebius_of(m)
                if co is None:
                    notmoeb.append(m.name)
                    cen["not Moebius"] += 1
                    continue
                A, B, C, D = (sp.simplify(z) for z in co)
                if sp.simplify(A * D - B * C) == 0:
                    cen["blind"] += 1          # K constant; not a curvature class
                    continue
                Bz, Cz, Az = B == 0, C == 0, A == 0
                cls = ("inverse (A=0)" if Az else
                       "identity" if Bz and Cz else "concave (B=0)" if Bz
                       else "convex (C=0)" if Cz else "mixed")
                cen[cls] += 1
                for sub in draws_for(co, m.name, k=3):
                    _, geo, _ = classify(co, sub)
                    pg, _ = profile_geometry(co, sub)
                    if pg != geo:
                        mism += 1
                        if verbose:
                            print(f"  MISMATCH {m.name} {geo} vs {pg}")
    return cen, notmoeb, mism


def parallel_case(ntrial=400, seed=11):
    """Two isotope-sensitive steps reached from a common state.

    This is the configuration the chain atlas flags as outside the Moebius
    classification.  Here the two branches carry separate intrinsic pairs, which
    is the general situation; the question is whether the offset can change sign
    along the mass-scaling ray.
    """
    rng = np.random.default_rng(seed)
    x = sp.Symbol("x", positive=True)
    koff, kb, kc1, kc2, ko2, ko3, kT = sp.symbols(
        "k_off k_b k_c1 k_c2 k_o2 k_o3 k_T", positive=True)
    m = Mech("parallel", 3,
             [(1, 0, koff), (1, 2, kb), (1, 3, sp.Symbol("k_d", positive=True)),
              (2, -1, kc1), (2, 0, ko2), (3, -1, kc2), (3, 0, ko3)],
             [kc1, kc2])
    f = first_passage(m)
    sw = 0
    for _ in range(ntrial):
        sub = {s: float(v) for s, v in
               zip([koff, kb, sp.Symbol("k_d", positive=True), ko2, ko3, kT],
                   rng.lognormal(0, 1.1, 6))}
        Kf = sp.lambdify(x, sp.simplify(
            f.subs({kc1: x * kT, kc2: x * kT}).subs(sub)
            / f.subs({kc1: kT, kc2: kT}).subs(sub)), "numpy")
        xd = np.geomspace(1 + 1e-6, 1e4, 3000)
        F = np.log(Kf(xd ** G)) - G * np.log(Kf(xd))
        F = F[np.isfinite(F)]
        if F.size and F.max() > 1e-12 and F.min() < -1e-12:
            sw += 1
    return ntrial, sw


if __name__ == "__main__":
    import sys
    nmax = int(sys.argv[1]) if len(sys.argv) > 1 else 3
    cen, nm, mism = run(nmax)
    print(f"general directed networks, up to {nmax} bound states")
    tot = sum(cen.values())
    for k, v in cen.most_common():
        print(f"   {k:16s} {v:5d}")
    print(f"   {'TOTAL':16s} {tot:5d}")
    print(f"\n   not Moebius in a single labeled constant: {len(nm)}")
    print(f"   classification vs numerical profile mismatches: {mism}")
    n, sw = parallel_case()
    print(f"\n   parallel sensitive steps: {sw} of {n} draws change sign")
