"""Profiles whose optimum sits on the edge of the F grid.

offset_analysis.py profiles the offset on |F| <= 1.5.  For four series the
smallest chi2 on that grid is at F = +1.5, so the tabulated "optimum" is a
property of the grid, not of the fit.  This script continues each of those
profiles outward, warm-started, to F = 40, and records

    chi2_inf   the infimum of the profile beyond the edge,
    F_lo       the 95% end recomputed against chi2_inf rather than the grid
               minimum (a lower infimum can only raise this end),
    finite     whether the profile has a finite minimum at all.

In every case found here chi2 falls monotonically and flattens as F grows: the
data favor ever larger intrinsic offsets hidden by ever heavier masking, the
same one-sidedness as the half-line of the main text.  There is then no point
estimate, the profile set is [F_lo, +inf), and goodness of fit is judged at
chi2_inf.

Run from analysis/ after offset_analysis.py:   python profile_edges.py
"""
from __future__ import annotations

import numpy as np
import pandas as pd

import offset_analysis as O

EDGE = 1.5
FAR = np.r_[np.arange(1.6, 8.01, 0.1), [10.0, 15.0, 20.0, 30.0, 40.0]]
KEY = ["family", "system", "variant", "step"]


def extend(T, kht, sh, kdt, sd):
    """Continue the profile from the grid edge outward; return (chi2_inf, F_lo, finite)."""
    r = O.fit_at(T, kht, sh, kdt, sd, EDGE, n_rand=400)
    edge_sol, warm = r.x, r.x
    vals = [(EDGE, 2 * r.cost)]
    for F in FAR:
        rr = O.fit_at(T, kht, sh, kdt, sd, F, n_rand=20, warm=warm)
        vals.append((F, 2 * rr.cost))
        warm = rr.x
    vals = np.array(vals)
    chi2_inf = float(vals[:, 1].min())
    # a finite optimum would show chi2 rising again before the far end
    finite = bool(vals[-1, 1] > chi2_inf + 1e-3)
    # walk inward from the edge, then take the 95% end against the infimum
    # over both walks: the smallest F whose chi2 is within 3.841 of it
    warm = edge_sol
    for F in np.arange(EDGE - 0.01, -1.5001, -0.01):
        rr = O.fit_at(T, kht, sh, kdt, sd, F, n_rand=6, warm=warm)
        warm = rr.x
        vals = np.vstack([vals, [F, 2 * rr.cost]])
        if 2 * rr.cost > chi2_inf + 3.841:
            break
    inside = vals[vals[:, 1] <= chi2_inf + 3.841, 0]
    F_lo = round(float(inside.min()), 2)
    return chi2_inf, F_lo, finite


def main():
    O.RNG = np.random.default_rng(20261012)
    prof = pd.read_csv("../results/offset_profiles.csv")
    d = pd.read_csv("../data/trinomial_benchmark.csv")
    rows = []
    for _, p in prof.iterrows():
        if abs(p.F_hat - EDGE) > 1e-9:
            continue
        g = d[(d[KEY] == p[KEY]).all(axis=1)].sort_values("T_C")
        T = g.T_C.to_numpy() + 273.15
        chi2_inf, F_lo, finite = extend(T, g.K_HT.to_numpy(), g.K_HT_se.to_numpy(),
                                        g.K_DT.to_numpy(), g.K_DT_se.to_numpy())
        rows.append(dict(family=p.family, system=p.system, variant=p.variant, step=p.step,
                         chi2_grid=p.chi2_min, chi2_inf=chi2_inf, F_lo_grid=p.F_lo,
                         F_lo=F_lo, finite_optimum=finite))
        print(f"{p.family} {p.variant} ({p.step}): chi2 {p.chi2_min:.3f} at the edge, "
              f"infimum {chi2_inf:.3f}, finite optimum {finite}; "
              f"95% set [{F_lo:+.2f}, +inf) (grid gave {p.F_lo:+.2f})")
    pd.DataFrame(rows).to_csv("../results/profile_edges.csv", index=False)


if __name__ == "__main__":
    main()
