"""The identified set under a reversible chemical step.

The main text inverts the observation map of the minimal irreversible scheme.
Hydride transfer in DHFR and formate oxidation in FDH are reversible, so the
relevant question is what survives when the chemical step can run backwards.
This module derives that map from the kinetics rather than by analogy, and shows
that reversibility introduces a sharp threshold: below it the identified set is
still a bounded half-line, above it the set is the whole real line and the
experiment determines nothing about the offset.

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

THE THRESHOLD.  Both D vanish simultaneously only at the solution of a 2x2
linear system.  With a = K_H - 1, b = K_D - 1, u = K_H/E_H - 1, v = K_D/E_D - 1,

    C_f* = (v-u)/(a v - b u),      C_r* = (a-b)/(a v - b u).

Since a > b under normal ordering, that corner lies in the physical quadrant
C_f, C_r >= 0 exactly when v >= u.  Tying the equilibrium isotope effects by
mass scaling, E_H = E_D^gamma, this reduces to

    E_D >= (K_HT / K_DT)^{1/(gamma-1)}  ==  E_D*,

the vacuity threshold.  Below it the set is bounded below; at or above it the
set is all of R.  Note the contrast that appears here is the RATIO K_HT/K_DT,
not the difference contrast L_H = (K_HT-1)/(K_DT-1) that governs the
irreversible case.  Reversibility is a different question and answers to a
different statistic.
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
    """E_D* above which the identified set is the whole real line."""
    return (kht / kdt) ** (1.0 / (gamma - 1.0))


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
    print(f"{'K_HT':>7}{'K_DT':>7}{'E_D* (analytic)':>17}"
          f"{'F_min at 0.9 E*':>17}{'F_min at 1.1 E*':>17}")
    for kht, kdt in ((5.04, 1.65), (2.362, 1.49), (6.44, 1.92), (1.93, 1.29)):
        es = eie_threshold(kht, kdt)
        lo, _, _ = F_min_reversible(kht, kdt, 0.9 * es)
        hi, _, _ = F_min_reversible(kht, kdt, 1.1 * es)
        print(f"{kht:7.2f}{kdt:7.2f}{es:17.4f}{lo:17.3f}{hi:17.3f}")
    print("  Below the threshold F_min is finite; above it the search runs away,")
    print("  which is the numerical signature of an unbounded identified set.")

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
    print(f"  73 records: E_D* ranges {d.E_star.min():.3f} to "
          f"{d.E_star.max():.3f}, median {d.E_star.median():.3f}")
    print()
    print(f"{'family':10s}{'variant':26s}{'step':9s}"
          f"{'min E_D*':>10}{'median E_D*':>13}")
    for key, g in d.groupby(["family", "variant", "step"]):
        print(f"{key[0]:10s}{key[1]:26s}{key[2]:9s}"
              f"{g.E_star.min():10.3f}{g.E_star.median():13.3f}")
    print()
    for e in (1.0, 1.1, 1.2, 1.3):
        n = int((d.E_star <= e).sum())
        print(f"  at E_DT = {e:.1f}: {n:2d} of 73 records would be vacuous")

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
    d[["family", "variant", "step", "T_C", "K_HT", "K_DT", "E_star"]].to_csv(
        "../results/reversible_thresholds.csv", index=False)
    print("\n[written] ../results/reversible_{thresholds,bounds}.csv")


if __name__ == "__main__":
    _verify()
    _benchmark()
