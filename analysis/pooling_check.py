"""Is a common offset within a temperature series consistent with the errors?

The series statistic used for the bounds is a maximum over temperatures, which
is inefficient. An inverse-variance weighted mean would be far more powerful if
the offset were constant within a series. This script tests that premise, and
reports the numbers quoted in the supplement, so the claim is reproducible
rather than asserted.

For each series the offset per record is F = ln K_HT - gamma ln K_DT with

    sigma_F^2 = (s_HT/K_HT)^2 + gamma^2 (s_DT/K_DT)^2
                - 2 gamma rho (s_HT/K_HT)(s_DT/K_DT),

the two effects sharing a tritium reference so their correlation is unknown.
rho = 0 is the neutral choice; rho = -1 inflates the errors most and is
therefore the case most favourable to pooling.

Run from analysis/:   python pooling_check.py
"""
from __future__ import annotations

import numpy as np
import pandas as pd

import masses as M

G = M.gamma_sc("C")
KEY = ["family", "system", "variant", "step"]


def series_chi2(d, rho):
    """chi2 per degree of freedom for a constant offset, series by series."""
    out = []
    for key, g in d.groupby(KEY):
        F = np.log(g.K_HT) - G * np.log(g.K_DT)
        a, b = g.K_HT_se / g.K_HT, g.K_DT_se / g.K_DT
        sig = np.sqrt(a ** 2 + (G * b) ** 2 - 2 * G * rho * a * b)
        w = 1.0 / sig ** 2
        Fbar = float((F * w).sum() / w.sum())
        chi2 = float((((F - Fbar) / sig) ** 2).sum())
        dof = len(g) - 1
        if dof > 0:
            out.append((f"{key[0]} {key[2]} ({key[3]})", len(g), chi2 / dof))
    return pd.DataFrame(out, columns=["series", "n_T", "chi2_dof"])


def main():
    d = pd.read_csv("../data/trinomial_benchmark.csv")
    print(f"constant-offset fit within each series, {d.groupby(KEY).ngroups} series")
    for rho, note in ((0.0, "independent errors"),
                      (-1.0, "maximally adverse, most favourable to pooling")):
        r = series_chi2(d, rho)
        print(f"\n  rho = {rho:+.1f}  ({note})")
        print(f"    median chi2/dof {r.chi2_dof.median():.2f}   "
              f"max {r.chi2_dof.max():.1f}  ({r.loc[r.chi2_dof.idxmax(), 'series']})")
        print(f"    series with chi2/dof > 1: "
              f"{int((r.chi2_dof > 1).sum())} of {len(r)}")
    r0 = series_chi2(d, 0.0)
    print("\n  worst five at rho = 0:")
    print(r0.sort_values("chi2_dof", ascending=False).head(5).to_string(index=False))


if __name__ == "__main__":
    main()
