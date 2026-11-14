"""How each analysis unit stands against the semiclassical reference, F = 0.

The census is reported against two hypotheses.  The accounting against the
ground-channel tunneling model (units whose largest endpoint lies below F0,
units limited by precision, nominal exclusions) is made by
bounds_uncertainty.py and drawn by export_fig_asym.py.  This module makes the
same accounting against the reference relation itself, which is the comparison
the Swain-Schaad exponent has always been used for:

  below     the observed exponent lies below gamma_SC at every admissible
            temperature.  The allowed range of the intrinsic offset then
            contains F = 0 (Proposition S2), so no gain in precision at those
            values can show an intrinsic departure from the reference.  The
            count depends on the mass convention of the reference; the count
            for the bare-mass exponent is given beside it.
  above     an observed exponent exceeds gamma_SC.  By Theorem 1 that record
            lies on the open branch and its endpoint is F_obs > 0.
  above_95  of those, the one-sided 95% bound at rho = -1 stays above zero.
  above_adj of those, the bound adjusted across all analysis units stays above
            zero.

Bounds are the exact rho = -1 quantiles of bounds_uncertainty.py (Bonferroni
over the temperatures of a series; across units, over units as well).

Run from analysis/:   python reference_counts.py
"""
from __future__ import annotations

import numpy as np
import pandas as pd

import bounds_uncertainty as BU
import corpus
import masses as M
from partial_id import F_min_exact

GSC = M.gamma_sc("C")
F0 = M.offset_F0("C")
# the exponent for bare isotope masses (the familiar 3.26), for comparison
G_BARE = ((M.M_H_ATOMIC ** -0.5 - M.M_T_ATOMIC ** -0.5)
          / (M.M_D_ATOMIC ** -0.5 - M.M_T_ATOMIC ** -0.5))


def unit_table() -> pd.DataFrame:
    """One row per analysis unit: largest endpoint, largest observed exponent,
    per-unit bound and bound adjusted across all units (rho = -1)."""
    units = corpus.counts()["units"]
    rows = []
    s = corpus.series()
    for (fam, var, step), g in s.groupby(["family", "variant", "step"], sort=False):
        rows.append((f"{fam} {var} ({step})", g))
    for _, r in corpus.single_condition().iterrows():
        rows.append((f"{r.family} {r.variant}, {r.T_C:g} C ({r.source})", pd.DataFrame([r])))
    # labeling design of each form: "labeled" (mixed-label: the D/T substrate
    # and its tritiated reference carry deuterium at the non-transferred
    # position) or "shared" (one tritiated molecule serves both experiments)
    design = {(d.family, d.variant, d.step): d.design
              for d in pd.read_csv(corpus.DATA / "tracer_design.csv").itertuples()}
    out = []
    for label, g in rows:
        adm = g[(g.K_HT > g.K_DT) & (g.K_DT > 1)]
        point = max(F_min_exact(kh, kd)[0] for kh, kd in zip(adm.K_HT, adm.K_DT))
        gmax = float(np.max(np.log(adm.K_HT) / np.log(adm.K_DT)))
        key = (g.family.iloc[0], g.variant.iloc[0], g.step.iloc[0])
        out.append(dict(unit=label, design=design[key], n_T=len(g), point=point, gamma_obs_max=gmax,
                        lcb=BU.exact_series_bound(g, units=1),
                        lcb_all_units=BU.exact_series_bound(g, units=units)))
    return pd.DataFrame(out)


def counts(t: pd.DataFrame | None = None) -> dict:
    """The accounting against F = 0, with the accounting against F0 beside it."""
    t = unit_table() if t is None else t
    above = t.point > 0
    # the two descriptions of "above" must agree: an exponent above gamma_SC
    # puts the record on the open branch, where the endpoint is F_obs > 0
    if not (above == (t.gamma_obs_max > GSC)).all():
        raise RuntimeError("endpoint sign and observed exponent disagree for a unit")
    return dict(units=len(t),
                below=int((~above).sum()),
                below_bare=int((t.gamma_obs_max < G_BARE).sum()),
                above=int(above.sum()),
                above_95=int((t.lcb > 0).sum()),
                above_adj=int((t.lcb_all_units > 0).sum()),
                above_adj_units=list(t.unit[t.lcb_all_units > 0]),
                above_mixed=int((above & (t.design == "labeled")).sum()),
                model_excluded_mixed=int(((t.lcb > F0) & (t.design == "labeled")).sum()),
                model_uninformative=int((t.point < F0).sum()),
                model_precision_limited=int(((t.point >= F0) & (t.lcb <= F0)).sum()),
                model_excluded=int((t.lcb > F0).sum()),
                model_excluded_adj=int((t.lcb_all_units > F0).sum()))


def main():
    t = unit_table()
    c = counts(t)
    t.sort_values("point").to_csv("../results/reference_counts.csv", index=False,
                                  float_format="%.6f")
    print(f"{c['units']} analysis units against the semiclassical reference F = 0")
    print(f"  observed below the reference at every admissible temperature: {c['below']}")
    print(f"  observed above it: {c['above']}")
    print(f"  (with the bare-mass exponent {G_BARE:.5f} as the reference: "
          f"{c['below_bare']} below, {c['units'] - c['below_bare']} above)")
    print(f"    one-sided 95% bound above zero (rho = -1): {c['above_95']}")
    print(f"    still above zero after adjustment across all units: {c['above_adj']}"
          f"  {c['above_adj_units']}")
    print(f"  mixed-label designs among the units above the reference: {c['above_mixed']} of "
          f"{c['above']}; among the nominal exclusions of the model: {c['model_excluded_mixed']} "
          f"of {c['model_excluded']}")
    print(f"against the ground-channel model, F0 = {F0:.4f}: "
          f"{c['model_uninformative']} / {c['model_precision_limited']} / "
          f"{c['model_excluded']}; {c['model_excluded_adj']} after adjustment")
    print("wrote ../results/reference_counts.csv")


if __name__ == "__main__":
    main()
