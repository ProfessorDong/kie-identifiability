"""Atlas of kinetic networks and the identification geometry they produce.

For a competitive experiment the specificity constant is

    V/K = k_on * f,

with f the probability that a complex formed from free enzyme reaches product
before dissociating.  f is a first-passage probability on the bound states, so
it is a ratio of sums of path products, and it is AFFINE in any rate constant
that labels exactly one edge.  Referencing to tritium therefore gives a Moebius
observation map whenever a single edge carries the isotope-sensitive constant,

    K(x) = (A x + B)/(C x + D),      K(1) = 1,

and the geometry of the identified set for the offset F follows from the signs
of B and C (Proposition S5 of the supplement):

    B = 0, C > 0   concave      identified set opens ABOVE
    C = 0, B > 0   convex       identified set opens BELOW
    B, C > 0       mixed        switches at x* = sqrt(BD/(AC))
    B = C = 0      linear       no masking

This module derives f symbolically for an enumerated family of mechanisms,
extracts the coefficients, classifies each network, and checks every
classification against direct numerical profiling of the identified set.
Nothing here is asserted that is not verified against the profile.
"""
from __future__ import annotations

import itertools
from dataclasses import dataclass, field

import numpy as np
import sympy as sp

import masses as M

FLOOR = 1e-12          # noise floor of ln K_HT - gamma ln K_DT near x = 1
GSC = M.gamma_sc("C")
F0 = M.offset_F0("C")


# ---------------------------------------------------------------- mechanisms
@dataclass
class Mech:
    """A mechanism, as first-passage structure on the bound states.

    states   : number of bound intermediates, labelled 1..n
    edges    : (i, j, sym) transitions between bound states; j == 0 means
               dissociation back to free enzyme, j == -1 means release of
               product.  State 1 is the complex formed on binding.
    iso      : symbols carrying the isotope-sensitive rate constant
    name     : label for the atlas table
    """
    name: str
    states: int
    edges: list
    iso: list
    note: str = ""
    subs: dict = field(default_factory=dict)


def first_passage(m: Mech):
    """P(reach product before free enzyme | start in state 1), symbolically."""
    n = m.states
    f = sp.symbols(f"f1:{n + 1}")                      # f[i] for state i+1
    eqs = []
    for i in range(1, n + 1):
        out = [(j, s) for (a, j, s) in m.edges if a == i]
        tot = sum(s for (_, s) in out)
        rhs = 0
        for (j, s) in out:
            if j == -1:                                 # product: success
                rhs += s * 1
            elif j == 0:                                # free enzyme: failure
                rhs += s * 0
            else:
                rhs += s * f[j - 1]
        eqs.append(sp.Eq(f[i - 1] * tot, rhs))
    sol = sp.solve(eqs, list(f), dict=True)
    if not sol:
        raise ValueError(f"no steady state for {m.name}")
    return sp.simplify(sol[0][f[0]])


def moebius(m: Mech):
    """(A, B, C, D) of K(x) = (V/K)(x k_T)/(V/K)(k_T), or None if not Moebius."""
    if len(m.iso) != 1:
        return None                                     # >1 sensitive edge
    k = m.iso[0]
    x, kT = sp.symbols("x k_T", positive=True)
    f = first_passage(m)
    K = sp.simplify(f.subs(k, x * kT) / f.subs(k, kT))
    num, den = sp.fraction(sp.cancel(sp.together(K)))
    pn, pd = sp.Poly(sp.expand(num), x), sp.Poly(sp.expand(den), x)
    if pn.degree() > 1 or pd.degree() > 1:
        return None
    A, B = pn.coeff_monomial(x), pn.coeff_monomial(1)
    C, D = pd.coeff_monomial(x), pd.coeff_monomial(1)
    return tuple(sp.simplify(z) for z in (A, B, C, D))


def curvature_point(coeffs, subs):
    """The inflection x* = sqrt(BD/(AC)) of ln K in log coordinates.

    Writing h(u) = ln K(e^u), p = B/A and q = D/C, one finds
    h''(u) proportional to (p - q)(1 - p q e^{-2u}).  Normalization K(1) = 1
    gives A + B = C + D, and a normal (increasing) effect gives AD - BC > 0,
    which is exactly q > p.  So the first factor is negative always, and

        h is CONVEX for x < x*,  CONCAVE for x > x*,   x* = sqrt(BD/(AC)).

    Intrinsic effects are normal, so only x > 1 is physical.  A switch point
    at or below 1 therefore lies outside the accessible domain and the network
    behaves as a purely concave one.
    """
    A, B, C, D = (float(sp.N(z.subs(subs))) for z in coeffs)
    if A == 0:
        return np.inf          # x* = sqrt(BD/(AC)) diverges as A -> 0
    if C == 0:
        return np.inf
    if B == 0:
        return 0.0
    return float(np.sqrt(B * D / (A * C)))


def classify(coeffs, subs):
    """Curvature class and identified-set geometry at a given parameter point.

    Returns (class, geometry, x*).  The geometry is the direction of the
    identified set over the PHYSICAL domain of intrinsic effects, x > 1, which
    is what an experiment can constrain; a curvature switch below x = 1 is not
    reachable and does not make the network ambiguous.
    """
    A, B, C, D = (float(sp.N(z.subs(subs))) for z in coeffs)
    tol = 1e-12
    det = A * D - B * C
    # K is constant exactly when the determinant vanishes.  A = 0 alone is not
    # constancy: K = B/(Cx+D) still decreases in x, which is an inverse observed
    # effect, and its log is concave, so the set still opens above.  Labeling a
    # dissociation step produces exactly this, which is why a family that labels
    # only forward steps never meets it.
    if abs(det) < tol:
        return "blind", "unrestricted", np.nan
    if abs(B) < tol and abs(C) < tol:
        return "identity", "no masking", np.nan
    xs = curvature_point(coeffs, subs)
    # h(u) = ln K(e^u) is concave exactly where sign(det) * (x - x*) > 0.
    if det > 0:
        if xs <= 1.0:
            return ("concave" if B == 0 else "concave (effectively)"), "opens above", xs
        if np.isinf(xs):
            return "convex", "opens below", xs
        return "mixed", "switches", xs
    # det < 0: the network turns a normal intrinsic effect into an INVERSE
    # observed one, and the curvature regions exchange places.
    if xs <= 1.0:
        return "inverse, convex", "opens below", xs
    if np.isinf(xs):
        return "inverse, concave", "opens above", xs
    return "inverse, mixed", "switches", xs


def profile_geometry(coeffs, subs, xd_hi=1e6, n=4000):
    """Direct numerical check of the identified-set direction over x > 1.

    Sweeps intrinsic pairs along the mass-scaling ray and records the sign of
    the offset F = ln K_HT - gamma ln K_DT.  The sweep starts away from x = 1,
    where F vanishes identically and floating point noise would otherwise
    manufacture a spurious sign change, and judges signs relative to the
    largest magnitude seen rather than against an absolute floor.
    """
    A, B, C, D = (float(sp.N(z.subs(subs))) for z in coeffs)
    if abs(A * D - B * C) < 1e-12:
        return "unrestricted", 0.0
    K = lambda x: (A * x + B) / (C * x + D)
    xd = np.geomspace(1.0 + 1e-6, xd_hi, n)
    # resolve the neighborhood of the curvature switch, which can be a very
    # thin interval just above x = 1 and would otherwise fall between grid
    # points; a switch that is real but confined there is still a switch, and
    # the returned magnitude is what says whether it could ever be measured
    xs = curvature_point(coeffs, subs)
    if np.isfinite(xs) and xs > 1.0:
        lo, hi = 1.0 + (xs - 1.0) * 1e-4, xs * 4.0
        xd = np.unique(np.concatenate([xd, np.geomspace(lo, hi, n // 2)]))
    F = np.log(K(xd ** GSC)) - GSC * np.log(K(xd))
    F = F[np.isfinite(F)]
    if F.size < 10:
        return "undetermined", 0.0
    scale = np.max(np.abs(F))
    if scale < 1e-14:
        return "no masking", 0.0
    # F is a difference of logarithms and cancels to zero along x = 1, so a
    # minority branch smaller than the double-precision noise of that
    # cancellation is not a sign change.  The floor is absolute and reported,
    # not tuned: below it the two branches are indistinguishable from zero.
    up, dn = float(F.max()), float(F.min())
    if up <= FLOOR:
        return "opens above", up
    if dn >= -FLOOR:
        return "opens below", dn
    return "switches", float(min(up, -dn))          # minority-branch magnitude


# ---------------------------------------------------------------- enumeration
def _sym(t):
    return sp.symbols(t, positive=True)


def enumerate_family(nmax=3):
    """Systematic family: chains of 1..nmax bound states with optional features.

    Every mechanism has a complex formed on binding (state 1) that can
    dissociate, a chain of forward steps to the last state, and product release
    from the last state.  The optional features are the ones that are argued
    in the literature to matter for isotope effects: reversibility of any
    forward step, an isotope-blind parallel route to product, and a second
    exit from the last state.  The isotope-sensitive constant is placed on
    each forward step in turn, so the atlas covers early and late chemistry.
    """
    out = []
    for n in range(1, nmax + 1):
        koff = _sym("k_off")
        fwd = [_sym(f"k_f{i}") for i in range(1, n)] + [_sym("k_cat")]
        rev = [_sym(f"k_r{i}") for i in range(1, n)]
        nrev = len(rev)
        for rmask in range(1 << nrev):                  # which steps reverse
            for byp in (0, 1):                          # blind route to product
                for leak in (0, 1):                     # second exit at the end
                    if leak and n == 1:
                        continue                        # duplicates k_off
                    for iso_at in range(len(fwd)):      # where the label sits
                        E = [(1, 0, koff)]
                        for i in range(n - 1):
                            E.append((i + 1, i + 2, fwd[i]))
                            if rmask >> i & 1:
                                E.append((i + 2, i + 1, rev[i]))
                        E.append((n, -1, fwd[-1]))
                        if byp:
                            E.append((1, -1, _sym("k_byp")))
                        if leak:
                            E.append((n, 0, _sym("k_off2")))
                        feats = []
                        if rmask:
                            feats.append("reversible " +
                                         ",".join(str(i + 1) for i in range(nrev)
                                                  if rmask >> i & 1))
                        if byp:
                            feats.append("blind bypass")
                        if leak:
                            feats.append("late exit")
                        pos = "chemistry" if iso_at == len(fwd) - 1 else f"step {iso_at + 1}"
                        out.append(Mech(
                            name=f"n={n}; iso on {pos}" +
                                 ("; " + "; ".join(feats) if feats else ""),
                            states=n, edges=E, iso=[fwd[iso_at]],
                            note="; ".join(feats)))
    return out


def draws(coeffs, tag="", k=6, seed=20260810):
    """Random positive values for every free symbol appearing in the coefficients.

    Drawing from the coefficients rather than the edge list matters: after the
    labelled constant is rewritten as x*k_T, the reference constant k_T survives
    in A..D and would otherwise be left unsubstituted.
    """
    rng = np.random.default_rng(seed + len(tag))
    syms = sorted(set().union(*(z.free_symbols for z in coeffs)), key=str)
    return [{s: float(v) for s, v in zip(syms, rng.lognormal(0, 1.1, len(syms)))}
            for _ in range(k)]


def structural_class(coeffs):
    """Topology-only class, from SYMBOLIC vanishing of the coefficients.

    This is the part fixed by the network alone.  "conditional" marks the
    networks where topology leaves both directions open and the outcome is
    decided by an inequality among the rate constants, x* > 1 or x* <= 1.
    """
    A, B, C, D = (sp.simplify(z) for z in coeffs)
    if sp.simplify(A) == 0 or sp.simplify(A * D - B * C) == 0:
        return "blind", "unrestricted"
    zb, zc = sp.simplify(B) == 0, sp.simplify(C) == 0
    if zb and not zc:
        return "concave", "opens above"          # B=0 forces det=AD>0
    if zc and not zb:
        return "convex", "opens below"           # C=0 forces det=AD>0
    if zb and zc:
        return "identity", "no masking"
    return "conditional", "set by x* vs 1"


NDRAW = 12


def atlas(nmax=3, verbose=True):
    """Build the atlas and check every entry against numerical profiling."""
    rows, flagged, bad = [], [], 0
    for m in enumerate_family(nmax):
        co = moebius(m)
        if co is None:
            flagged.append((m.name, "not Moebius in the labelled constant"))
            continue
        cls, geo = structural_class(co)
        agree, seen = True, set()
        for sub in draws(co, m.name, k=NDRAW):
            ncls, ngeo, xs = classify(co, sub)
            pgeo, _ = profile_geometry(co, sub)
            seen.add(ngeo)
            # the numerical direction must match the direct profile, always;
            # the topology-only class must match unless it is conditional
            ok = (pgeo == ngeo) and (cls == "conditional" or ncls.startswith(cls))
            if cls == "conditional":
                ok = ok and pgeo in ("opens above", "opens below", "switches")
            if not ok:
                agree = False
                bad += 1
                if verbose:
                    print(f"  MISMATCH {m.name}: struct={cls}/{geo} "
                          f"num={ncls}/{ngeo} x*={xs:.4g} profile={pgeo}")
        rows.append(dict(name=m.name, states=m.states, cls=cls, geometry=geo,
                         realized=sorted(seen), verified=agree, note=m.note))
    return rows, flagged, bad


def reaches_free_enzyme(m: Mech, start: int):
    """Can the free-enzyme state be reached from `start` along bound-state edges?"""
    seen, stack = set(), [start]
    while stack:
        i = stack.pop()
        if i in seen:
            continue
        seen.add(i)
        for (a, j, _) in m.edges:
            if a != i:
                continue
            if j == 0:
                return True
            if j > 0:
                stack.append(j)
    return False


def all_first_passage(m: Mech):
    """First-passage success probability from every bound state, symbolically."""
    n = m.states
    f = sp.symbols(f"f1:{n + 1}")
    eqs = []
    for i in range(1, n + 1):
        out = [(j, sy) for (a, j, sy) in m.edges if a == i]
        tot = sum(sy for (_, sy) in out)
        rhs = sum(sy * (1 if j == -1 else 0 if j == 0 else f[j - 1])
                  for (j, sy) in out)
        eqs.append(sp.Eq(f[i - 1] * tot, rhs))
    # A linear solve rather than sp.solve: the first-passage equations are
    # linear in f, and the general solver can return a dict missing variables
    # it has eliminated, which silently breaks the caller on some graphs.
    A = sp.zeros(n, n)
    b = sp.zeros(n, 1)
    for i in range(1, n + 1):
        out = [(j, sy) for (a, j, sy) in m.edges if a == i]
        A[i - 1, i - 1] = sum(sy for (_, sy) in out)
        for (j, sy) in out:
            if j == -1:
                b[i - 1] += sy
            elif j > 0:
                A[i - 1, j - 1] -= sy
    sol = A.LUsolve(b)
    val = {0: sp.Integer(0), -1: sp.Integer(1)}
    val.update({i + 1: sp.simplify(sol[i]) for i in range(n)})
    return val


def blind_by_topology(m: Mech):
    """A labelled edge carries no competitive information unless the branch it
    controls leads somewhere with a different fate.

    Rate constants enter a first-passage probability only through the way they
    split flux at a branch point.  The label on an edge i -> j is therefore
    invisible unless state i has another exit whose own success probability
    differs from that of j.  Two familiar situations are the special cases: a
    step that is the sole exit from its state is taken with probability one no
    matter how fast it is, and any step downstream of the last irreversible
    commitment leads to product with certainty along every route, so all its
    competing exits share the same fate.
    """
    val = all_first_passage(m)
    for (a, j, sym) in m.edges:
        if sym not in m.iso:
            continue
        others = [jj for (aa, jj, sy) in m.edges if aa == a and sy is not sym]
        if any(sp.simplify(val[jj] - val[j]) != 0 for jj in others):
            return False
    return True


def two_label_family(nmax=3):
    """Mechanisms with TWO isotope-sensitive edges: reversible chemistry.

    A chemical step that is reversible is isotope sensitive in both directions,
    and the two are tied by the equilibrium isotope effect rather than being
    equal.  Labelling the forward constant x and the reverse x/eie covers this;
    eie = 1 is the special case of a purely kinetic label.
    """
    out = []
    for n in range(2, nmax + 1):
        koff = _sym("k_off")
        fwd = [_sym(f"k_f{i}") for i in range(1, n)] + [_sym("k_cat")]
        rev = [_sym(f"k_r{i}") for i in range(1, n)]
        for i in range(n - 1):                          # chemistry at step i+1
            for byp in (0, 1):
                E = [(1, 0, koff)]
                for j in range(n - 1):
                    E.append((j + 1, j + 2, fwd[j]))
                    if j == i:
                        E.append((j + 2, j + 1, rev[j]))
                E.append((n, -1, fwd[-1]))
                if byp:
                    E.append((1, -1, _sym("k_byp")))
                out.append(Mech(
                    name=f"n={n}; reversible chemistry at step {i + 1}"
                         + ("; blind bypass" if byp else ""),
                    states=n, edges=E, iso=[fwd[i], rev[i]],
                    note="two labelled edges" + ("; blind bypass" if byp else "")))
    return out


def moebius_two(m: Mech, eie):
    """(A,B,C,D) when the forward label is x and the reverse label is x/eie."""
    kf, kr = m.iso
    x, kT = sp.symbols("x k_T", positive=True)
    f = first_passage(m)
    sub_x = {kf: x * kT, kr: x * kT / eie}
    sub_1 = {kf: kT, kr: kT / eie}
    K = sp.simplify(f.subs(sub_x) / f.subs(sub_1))
    num, den = sp.fraction(sp.cancel(sp.together(K)))
    pn, pd = sp.Poly(sp.expand(num), x), sp.Poly(sp.expand(den), x)
    if pn.degree() > 1 or pd.degree() > 1:
        return None
    A, B = pn.coeff_monomial(x), pn.coeff_monomial(1)
    C, D = pd.coeff_monomial(x), pd.coeff_monomial(1)
    return tuple(sp.simplify(z) for z in (A, B, C, D))


# ---------------------------------------------------------------- reporting
def verify(nmax=4, ndraw=12, verbose=True):
    """Full check.  Returns the number of failures; zero is the required result.

    Three independent statements are checked for every mechanism:
      1. the branch-point criterion for blindness agrees with the algebra;
      2. B vanishes identically exactly when there is no isotope-blind route
         to product;
      3. the direction predicted from (sign of AD - BC, x* against 1) agrees
         with direct numerical profiling of the offset, at every draw.
    """
    fails = 0
    fam = enumerate_family(nmax) + two_label_family(nmax)
    for m in fam:
        two = len(m.iso) == 2
        co = moebius_two(m, sp.Rational(6, 5)) if two else moebius(m)
        if co is None:
            if verbose:
                print(f"  FLAGGED (not Moebius): {m.name}")
            continue
        blind = blind_by_topology(m)
        alg_blind = (sp.simplify(co[0]) == 0 or
                     sp.simplify(co[0] * co[3] - co[1] * co[2]) == 0)
        if not two and blind != alg_blind:
            fails += 1
            print(f"  FAIL blindness: {m.name}")
        if not blind and not two:
            if (sp.simplify(co[1]) == 0) == ("bypass" in m.note):
                fails += 1
                print(f"  FAIL bypass<->B: {m.name}")
        if not blind:
            for sub in draws(co, m.name, k=ndraw):
                _, geo, xs = classify(co, sub)
                pgeo, mag = profile_geometry(co, sub)
                if pgeo != geo:
                    fails += 1
                    print(f"  FAIL direction: {m.name} x*={xs:.6g} "
                          f"predicted={geo} profile={pgeo} |F|={abs(mag):.2e}")
    for k in (2, 3, 4, 5):
        if not series_factorizes(k):
            fails += 1
            print(f"  FAIL series factorization at n={k}")
    return fails


def report(nmax=3):
    """Print the atlas as it appears in the supplement."""
    fam = enumerate_family(nmax)
    print(f"{'mechanism':58s} {'class':14s} {'identified set':16s}")
    print("-" * 90)
    for m in fam:
        co = moebius(m)
        if blind_by_topology(m):
            cls, geo = "blind", "unrestricted"
        elif sp.simplify(co[1]) == 0:
            cls, geo = "concave", "opens above"
        else:
            cls, geo = "conditional", "above if x* <= 1"
        print(f"{m.name:58s} {cls:14s} {geo:16s}")


def series_chain(n):
    """Chain of n isotope-sensitive steps, each partitioning against its own exit."""
    kf = [sp.Symbol(f"k_f{i}", positive=True) for i in range(n)]
    ko = [sp.Symbol(f"k_o{i}", positive=True) for i in range(n)]
    E = []
    for i in range(n):
        E.append((i + 1, 0, ko[i]))
        E.append((i + 1, (i + 2) if i < n - 1 else -1, kf[i]))
    return Mech(f"chain of {n} sensitive steps", n, E, kf,
                note="multiple labels, series"), kf, ko


def series_factorizes(n):
    """K for a series chain is the PRODUCT of n series maps.

    The consequence matters more than the algebra: ln K is then a sum of
    functions each concave in log coordinates, so ln K is concave, and the
    identified set opens above exactly as in the single-step case.  Such a
    mechanism is not Moebius, so Proposition S5 does not reach it, yet its
    geometry is the same.
    """
    m, kf, ko = series_chain(n)
    x, kT = sp.symbols("x k_T", positive=True)
    f = first_passage(m)
    K = sp.cancel(sp.together(f.subs({s: x * kT for s in kf})
                              / f.subs({s: kT for s in kf})))
    prod = sp.prod([x * (kT + ko[i]) / (kT * x + ko[i]) for i in range(n)])
    return sp.simplify(K - prod) == 0


def check_manuscript(path=None):
    """Assert that the counts quoted in the supplement match this enumeration.

    Scoped to the atlas subsection, so a number that happens to appear
    elsewhere in the document cannot satisfy the check by accident.  This is
    the guard against the failure mode where the enumeration changes and the
    prose keeps the old census.
    """
    from pathlib import Path
    if path is None:
        path = ("/home/dong/Workspace/WritePaper/MDPI_Entropy_QuantumB/"
                "v3_quantum/manuscript/si/si_body.tex")
    txt = Path(path).read_text()
    a = txt.index(r"\subsection{An atlas of mechanisms}")
    b = txt.index(r"\subsection{A bypass contracts")
    sec = txt[a:b]
    main = Path(path).parent.parent / "pnas_manuscript.tex"
    mtxt = main.read_text()
    ma = mtxt.index("We tested this by enumeration")
    sec += mtxt[ma:ma + 1400]
    fam = enumerate_family(4)
    nblind = sum(1 for m in fam if blind_by_topology(m))
    nbyp = sum(1 for m in fam if not blind_by_topology(m) and "bypass" in m.note)
    want = {"total": len(fam), "blind": nblind,
            "informative": len(fam) - nblind, "with bypass": nbyp,
            "without bypass": len(fam) - nblind - nbyp}
    bad = 0
    for lab, v in want.items():
        if f"${v}$" not in sec:
            print(f"  FAIL manuscript census: {lab} = {v} not in the atlas section")
            bad += 1
    return bad


if __name__ == "__main__":
    import sys as _s
    print("network atlas: verification")
    n = verify(nmax=int(_s.argv[1]) if len(_s.argv) > 1 else 3)
    n += check_manuscript()
    print(f"\n{'FAIL' if n else 'PASS'}: {n} failures")
    print()
    report()
    _s.exit(1 if n else 0)


def sensitive_branch_points(m: Mech, scale):
    """States where the isotope label changes how flux is SPLIT.

    A rate constant enters a first-passage probability only through branching
    ratios.  At a state whose exits all carry the same power of the label, the
    ratios are label-free and the state contributes nothing; only states with
    mixed scaling do.  `scale` maps each edge symbol to its exponent of x.

    One such state is SUFFICIENT for the observation map to be Moebius, but it
    is not necessary: a reversible pair produces two such states and the map is
    Moebius anyway, because flux returning through the loop cancels the
    quadratic term.  The atlas therefore does not rely on this count; it
    computes the degree of every mechanism directly.
    """
    out = []
    val = all_first_passage(m)
    for i in range(1, m.states + 1):
        exits = [(j, sy) for (a, j, sy) in m.edges if a == i]
        if len(exits) < 2:
            continue
        if len({scale.get(sy, 0) for (_, sy) in exits}) < 2:
            continue                                    # ratios are label-free
        fates = {sp.simplify(val[j]) for (j, _) in exits}
        if len(fates) > 1:                              # and the fates differ
            out.append(i)
    return out
