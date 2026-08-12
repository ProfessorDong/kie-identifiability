"""The identified set under a reversible chemical step.

The main text inverts the observation map of the minimal irreversible scheme.
Hydride transfer in DHFR and formate oxidation in FDH are reversible, so the
relevant question is what survives when the chemical step can run backwards.
This module derives that map from the kinetics rather than by analogy, and finds
that reversibility opens a bounded WINDOW of equilibrium isotope effects inside
which the identified set is the whole real line and the experiment determines
nothing.  Outside the window, on either side, the set is bounded below.  The
window is empty whenever the observed offset is positive.

THE MAP.  For

    E + S  <=(k1,k2)=>  ES  <=(k3,k4)=>  EP  --(k5)-->  E + P

with k1, k2, k5 isotope independent, V/K = k1 k3 k5/(k2 k4 + k2 k5 + k3 k5).
Taking the ratio for isotopologue i against the tritium reference R, writing
x_i = k3i/k3R for the intrinsic effect and using k4i = k3i/Keq_i,

    K_iR = x_i (1 + C_f + C_r) / (1 + C_f x_i + C_r x_i / E_i),

    C_f = k3R/k2,   C_r = k4R/k5,   E_i = Keq_i/Keq_R.

Both commitments are defined on the reference isotopologue and are therefore
shared across the H/T and D/T measurements, which is the property the analysis
needs.  Setting C_r = 0 and C_f = 1/c recovers Eq. (1) of the main text exactly.

THE INVERSION.  With S = 1 + C_f + C_r and G_i = C_f + C_r/E_i,

    x_i = K_i / D_i,     D_i = S - K_i G_i = 1 + C_f(1-K_i) + C_r(1 - K_i/E_i),

so the offset takes the additively separable form

    F = F_obs + gamma ln D_D - ln D_H,

which is what makes the structure transparent: F diverges to -infinity exactly
when D_D can be driven to zero while D_H stays positive and smaller than
D_D K_H/K_D (the constraint x_H > x_D).

THE VACUITY WINDOW.  Both D vanish simultaneously only at the solution of a 2x2
linear system.  With a = K_H - 1, b = K_D - 1, u = K_H/E_H - 1, v = K_D/E_D - 1,

    C_f* = (v-u)/(a v - b u),      C_r* = (a-b)/(a v - b u),   den = a v - b u.

Since a > b under normal ordering, C_r* >= 0 needs den > 0, and C_f* >= 0 then
needs v >= u.  BOTH conditions matter.  Tying the equilibrium isotope effects by
mass scaling, E_H = E_D^gamma, the second gives

    E_D >= (K_HT / K_DT)^{1/(gamma-1)}  ==  E_D* ,

but that alone is not sufficient: den -> -(a-b) < 0 as E_D -> infinity, so den
changes sign again at some E_D** and the vacuity region is the bounded interval

    [E_D*, E_D**],   E_D** the root above E_D* of  a(K_DT/E - 1) = b(K_HT/E^g - 1),

which is transcendental in E for non-integer gamma and has no closed form.
Outside that window the set is bounded below again.

At E_D = E_D* one has v = u and hence den = u(a-b), so the window is nonempty
exactly when u > 0 there.  That reduces to K_DT^gamma > K_HT, i.e. to
F_obs < 0.  When the OBSERVED OFFSET IS POSITIVE the window is empty and no
equilibrium isotope effect can make the identified set vacuous, which is why the
yeast ADH conclusion of decisive_case.py is unconditional in E_D.

Checked against the optimizer and on 3000 random admissible pairs.

Note the contrast that appears here is the RATIO K_HT/K_DT, not the difference
contrast L_H = (K_HT-1)/(K_DT-1) that governs the irreversible case.
Reversibility is a different question and answers to a different statistic.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

import masses as M
from partial_id import F_min_exact

GSC = M.gamma_sc("C")
F0 = M.offset_F0("C")


# ------------------------------------------------------------------ the map
def observed(x, cf, cr, eie):
    """K_obs from the intrinsic effect under the reversible scheme."""
    return x * (1.0 + cf + cr) / (1.0 + cf * x + cr * x / eie)


def D_factor(K, cf, cr, eie):
    """D = 1 + C_f(1-K) + C_r(1 - K/E); the inversion is x = K/D."""
    return 1.0 + cf * (1.0 - K) + cr * (1.0 - K / eie)


def invert(K, cf, cr, eie):
    D = D_factor(K, cf, cr, eie)
    return np.where(D > 0, K / np.maximum(D, 1e-300), np.inf)


def eie_threshold(kht, kdt, gamma=GSC):
    """Lower edge E_D* of the vacuity window: the root of v = u.

    This is NECESSARY for vacuity but not sufficient. See vacuity_window.
    """
    return (kht / kdt) ** (1.0 / (gamma - 1.0))


def _den(E, kht, kdt, gamma=GSC):
    """den = a v - b u, whose positivity is the second vacuity condition."""
    a, b = kht - 1.0, kdt - 1.0
    return a * (kdt / E - 1.0) - b * (kht / E ** gamma - 1.0)


def vacuity_window(kht, kdt, gamma=GSC, emax=1e6):
    """Closed interval of E_DT on which the identified set is unbounded below.

    Vacuity requires the corner (C_f*, C_r*) to lie in the physical quadrant,
    which needs BOTH v >= u (giving E_DT >= E_D*) AND den = a v - b u > 0.
    The second condition fails again at large E_DT, because den -> -(a-b) < 0,
    so the window is a bounded interval [E_D*, E_D**] rather than a half-line.

    At E_DT = E_D* one has v = u and hence den = u(a-b), so the window is
    nonempty exactly when u > 0 there, which reduces to K_DT^gamma > K_HT,
    i.e. to F_obs < 0.  When the observed offset is POSITIVE no equilibrium
    isotope effect makes the set vacuous.

    Returns (E_lo, E_hi), or None when the window is empty.
    """
    from scipy.optimize import brentq
    f_obs = np.log(kht) - gamma * np.log(kdt)
    if f_obs >= 0.0:
        return None
    lo = eie_threshold(kht, kdt, gamma)
    if _den(lo, kht, kdt, gamma) <= 0.0:
        return None
    hi = brentq(_den, lo * (1 + 1e-12), emax, args=(kht, kdt, gamma))
    return lo, hi


# --------------------------------------------------------- the identified set
def F_reversible(kht, kdt, cf, cr, edt, gamma=GSC):
    """Offset implied by a commitment pair and an equilibrium isotope effect.

    Returns nan where the pair is inadmissible (x_H > x_D > 1 is required).
    """
    eht = edt ** gamma                       # EIEs tied by mass scaling
    DH = D_factor(kht, cf, cr, eht)
    DD = D_factor(kdt, cf, cr, edt)
    if DH <= 0 or DD <= 0:
        return np.nan
    xh, xd = kht / DH, kdt / DD
    if not (xh > xd > 1.0):
        return np.nan
    return np.log(xh) - gamma * np.log(xd)


def F_min_reversible(kht, kdt, edt, n=900, cmax=None, gamma=GSC):
    """Numerical infimum of F over (C_f, C_r) >= 0 at fixed E_DT.

    Returns (F_min, C_f, C_r).  A value of -inf signals the vacuous case; the
    grid is placed relative to the analytic corner so the divergence is actually
    resolved rather than missed between grid points.
    """
    a, b = kht - 1.0, kdt - 1.0
    eht = edt ** gamma
    u, v = kht / eht - 1.0, kdt / edt - 1.0
    den = a * v - b * u
    grids = [np.concatenate([[0.0], np.logspace(-4, 3, n)])]
    if abs(den) > 1e-14:                      # refine around the corner
        cfs, crs = (v - u) / den, (a - b) / den
        if np.isfinite(cfs) and np.isfinite(crs) and cfs >= 0 and crs >= 0:
            eps = np.logspace(-9, -0.3, 60)
            grids.append(np.unique(np.concatenate([
                cfs * (1 - eps), cfs * (1 + eps), [cfs]])))
            grids.append(np.unique(np.concatenate([
                crs * (1 - eps), crs * (1 + eps), [crs]])))
    cf_grid = np.unique(np.concatenate([g[g >= 0] for g in grids]))
    CF, CR = np.meshgrid(cf_grid, cf_grid, indexing="ij")
    with np.errstate(divide="ignore", invalid="ignore"):
        DH = 1.0 + CF * (1.0 - kht) + CR * (1.0 - kht / eht)
        DD = 1.0 + CF * (1.0 - kdt) + CR * (1.0 - kdt / edt)
        xh, xd = kht / DH, kdt / DD
        F = np.log(xh) - gamma * np.log(xd)
    ok = (DH > 0) & (DD > 0) & (xh > xd) & (xd > 1.0) & np.isfinite(F)
    if not ok.any():
        return np.nan, np.nan, np.nan
    F = np.where(ok, F, np.inf)
    k = np.unravel_index(np.argmin(F), F.shape)
    return float(F[k]), float(CF[k]), float(CR[k])


def series_bound_reversible(g, edt, gamma=GSC):
    """Max over temperatures of the reversible endpoint, as in the main text."""
    vals = [F_min_reversible(r.K_HT, r.K_DT, edt, gamma=gamma)[0]
            for _, r in g.iterrows()]
    vals = [v for v in vals if np.isfinite(v)]
    return max(vals) if vals else -np.inf


# =============================================================== verification
def _verify():
    print("=" * 78)
    print("1. THE MAP REDUCES TO THE IRREVERSIBLE CASE AT C_r = 0")
    print("=" * 78)
    worst = 0.0
    for kht, kdt in ((5.04, 1.65), (2.362, 1.49), (6.44, 1.92)):
        one = F_min_exact(kht, kdt)[0]
        rev, cf, cr = F_min_reversible(kht, kdt, edt=1.0)
        worst = max(worst, abs(one - rev))
        print(f"  K=({kht:5.2f},{kdt:5.2f})  irreversible {one:+.6f}   "
              f"reversible at E=1 {rev:+.6f}   diff {abs(one-rev):.2e}")
    print(f"  worst {worst:.2e}   {'PASS' if worst < 1e-4 else 'FAIL'}")
    print("  (at E = 1 the two commitments enter only through C_f + C_r,")
    print("   so the reversible family collapses onto the irreversible one.)")

    print()
    print("=" * 78)
    print("2. THE VACUITY THRESHOLD, ANALYTIC vs NUMERICAL")
    print("=" * 78)
    print(f"{'K_HT':>7}{'K_DT':>7}{'F_obs':>9}{'window [E_D*, E_D**]':>24}"
          f"{'below':>9}{'inside':>9}{'above':>9}")
    for kht, kdt in ((5.04, 1.65), (2.362, 1.49), (6.44, 1.92), (1.93, 1.29),
                     (7.13, 1.73)):
        fobs = np.log(kht) - GSC * np.log(kdt)
        win = vacuity_window(kht, kdt)
        if win is None:
            print(f"{kht:7.2f}{kdt:7.2f}{fobs:+9.4f}{'empty':>24}"
                  f"{F_min_reversible(kht,kdt,2.0)[0]:9.3f}{'--':>9}{'--':>9}")
            continue
        lo, hi = win
        b = F_min_reversible(kht, kdt, lo * 0.98)[0]
        m = F_min_reversible(kht, kdt, 0.5 * (lo + hi))[0]
        a_ = F_min_reversible(kht, kdt, hi * 1.02)[0]
        print(f"{kht:7.2f}{kdt:7.2f}{fobs:+9.4f}"
              f"{'['+format(lo,'.4f')+', '+format(hi,'.4f')+']':>24}"
              f"{b:9.3f}{m:9.1f}{a_:9.3f}")
    print("  The window is a bounded INTERVAL, not a half-line: the set is")
    print("  bounded below on both sides of it. It is empty when F_obs > 0.")

    print()
    print("=" * 78)
    print("3. HOW THE SET WIDENS BELOW THE THRESHOLD")
    print("=" * 78)
    kht, kdt = 5.04, 1.65
    es = eie_threshold(kht, kdt)
    print(f"  ecDHFR W133F at 25 C: K_HT={kht}, K_DT={kdt}, E_D* = {es:.4f}")
    print(f"{'E_DT':>8}{'E_DT/E_D*':>12}{'F_min':>12}{'C_f':>10}{'C_r':>10}")
    for e in (1.0, 1.05, 1.1, 1.2, 1.3, 1.4, 0.95 * es):
        f, cf, cr = F_min_reversible(kht, kdt, e)
        print(f"{e:8.3f}{e/es:12.3f}{f:12.4f}{cf:10.4f}{cr:10.4f}")
    print("  The endpoint falls monotonically as the equilibrium isotope effect")
    print("  approaches its threshold, so ignoring reversibility is optimistic.")


def _benchmark():
    print()
    print("=" * 78)
    print("4. THRESHOLD ACROSS THE BENCHMARK")
    print("=" * 78)
    d = pd.read_csv("../data/trinomial_benchmark.csv")
    d["E_star"] = eie_threshold(d.K_HT.values, d.K_DT.values)
    wins = [vacuity_window(r.K_HT, r.K_DT) for _, r in d.iterrows()]
    d["E_lo"] = [w[0] if w else np.nan for w in wins]
    d["E_hi"] = [w[1] if w else np.nan for w in wins]
    n_empty = int(d.E_lo.isna().sum())
    print(f"  73 records: lower edge E_D* ranges {d.E_star.min():.3f} to "
          f"{d.E_star.max():.3f}, median {d.E_star.median():.3f}")
    print(f"  upper edge E_D** median {d.E_hi.median():.3f}; "
          f"median window width {np.nanmedian(d.E_hi - d.E_lo):.3f}")
    print(f"  {n_empty} of 73 records have an EMPTY window (F_obs > 0), so no")
    print(f"  equilibrium isotope effect can make them vacuous.")
    print()
    print(f"{'family':10s}{'variant':26s}{'step':9s}"
          f"{'min E_D*':>10}{'median E_D*':>13}")
    for key, g in d.groupby(["family", "variant", "step"]):
        print(f"{key[0]:10s}{key[1]:26s}{key[2]:9s}"
              f"{g.E_star.min():10.3f}{g.E_star.median():13.3f}")
    d[["family","variant","step","T_C","K_HT","K_DT","E_star","E_lo","E_hi"]].to_csv(
        "../results/reversible_thresholds.csv", index=False)
    print()
    for e in (1.0, 1.1, 1.2, 1.3, 1.5):
        n = int(((d.E_lo <= e) & (e <= d.E_hi)).sum())
        print(f"  at E_DT = {e:.1f}: {n:2d} of 73 records are inside their window")

    print()
    print("=" * 78)
    print("5. SERIES BOUNDS UNDER REVERSIBILITY")
    print("=" * 78)
    print(f"{'series':40s}{'E=1.0':>10}{'E=1.1':>10}{'E=1.2':>10}"
          f"{'min E_D*':>11}")
    rows = []
    for key, g in d.groupby(["family", "variant", "step"]):
        vals = [series_bound_reversible(g, e) for e in (1.0, 1.1, 1.2)]
        es = g.E_star.min()
        txt = "".join(f"{v:10.4f}" if np.isfinite(v) else f"{'vacuous':>10}"
                      for v in vals)
        print(f"{key[0]+' '+key[1]+' ('+key[2]+')':40s}{txt}{es:11.3f}")
        rows.append(dict(family=key[0], variant=key[1], step=key[2],
                         F_min_E1p0=vals[0], F_min_E1p1=vals[1],
                         F_min_E1p2=vals[2], E_star_min=es))
    res = pd.DataFrame(rows)
    n_ref = {e: int((res[c] > F0).sum()) for e, c in
             ((1.0, 'F_min_E1p0'), (1.1, 'F_min_E1p1'), (1.2, 'F_min_E1p2'))}
    print()
    for e, n in n_ref.items():
        print(f"  at E_DT = {e}: {n} of 16 series exclude the gated model")
    print("  Reversibility can only lower the endpoint, so the irreversible")
    print("  bounds of the main text are the optimistic case.")
    res.to_csv("../results/reversible_bounds.csv", index=False)
    print("\n[written] ../results/reversible_{thresholds,bounds}.csv")


if __name__ == "__main__":
    _verify()
    _benchmark()
