"""The reference asymmetry r, computed from the sources' secondary isotope effects.

In the primary protocol the H/T and D/T tritium references differ at the
non-transferred position, so their commitments stand in a ratio r = c_D/c_H.
Where a source reports secondary H/T and D/T effects on those references, r is
inferred as their ratio (under the separability condition of the supplement):

    r = (secondary k_H/k_T) / (secondary k_D/k_T)

with first-order error propagation.  That propagation reproduces the values the
supplement quotes, r = 1.3107 +- 0.0164 for yeast ADH and 1.137 +- 0.019 for
bovine serum amine oxidase.

Until 2026-09-14 every r in the manuscript was a hand-typed constant, and none of
the amine oxidase secondary effects was in any data file.  Monoamine oxidase B at
pH 6.1 was consequently tabulated as 1.14, the value shared by bovine serum amine
oxidase and by monoamine oxidase B at pH 7.5; its own data give 1.13 to 1.21.
Every r stated in the text is now checked against this module by the audit.

Run from analysis/:   python reference_asymmetry.py
"""
from __future__ import annotations

import math
from pathlib import Path

import pandas as pd

DATA = Path(__file__).resolve().parent.parent / "data"


def secondary_effects() -> pd.DataFrame:
    """All reported secondary H/T and D/T effects, one row per measured condition."""
    a = pd.read_csv(DATA / "secondary_adh.csv")[
        ["family", "variant", "T_C", "K_HT", "K_HT_se", "K_DT", "K_DT_se", "source"]]
    b = pd.read_csv(DATA / "secondary_reference.csv")[
        ["family", "variant", "T_C", "K_HT", "K_HT_se", "K_DT", "K_DT_se", "source"]]
    d = pd.concat([a, b], ignore_index=True)
    dup = d.duplicated(["family", "variant", "T_C"], keep=False)
    if dup.any():
        raise ValueError("a condition appears in both secondary-effect files:\n"
                         + d[dup].to_string())
    return d


def r_of(row) -> tuple[float, float]:
    """r and its propagated standard error for one measured condition."""
    r = row.K_HT / row.K_DT
    return r, r * math.hypot(row.K_HT_se / row.K_HT, row.K_DT_se / row.K_DT)


def table() -> pd.DataFrame:
    d = secondary_effects().copy()
    rs = d.apply(r_of, axis=1, result_type="expand")
    d["r"], d["r_se"] = rs[0], rs[1]
    return d.sort_values(["family", "variant", "T_C"]).reset_index(drop=True)


def r_at(family: str, variant: str, T_C: float) -> float:
    """r for a single measured condition; raises if it was not measured."""
    d = table()
    hit = d[(d.family == family) & (d.variant == variant) & (d.T_C == T_C)]
    if len(hit) != 1:
        raise KeyError(f"no secondary effects for {family} {variant} at {T_C} C")
    return float(hit.r.iloc[0])


def r_range(family: str, variant: str) -> tuple[float, float]:
    """Smallest and largest r over the temperatures measured for one form."""
    d = table()
    x = d[(d.family == family) & (d.variant == variant)]
    if x.empty:
        raise KeyError(f"no secondary effects for {family} {variant}")
    return float(x.r.min()), float(x.r.max())


def tracer_design() -> pd.DataFrame:
    """How each form's two tritium references were made, read from the sources.

    'shared': one tritiated isotopologue serves both the H/T and the D/T
    experiment and only the bulk substrate changes, so the two references are
    the same molecule, c_D = c_H and r = 1 identically.  'labeled': the protocol
    places deuterium at the non-transferred position of the D/T reference, so
    the references differ and r must be inferred from secondary effects.
    """
    return pd.read_csv(DATA / "tracer_design.csv")


def r_design(family: str, variant: str) -> float:
    """r for a form: exactly 1 for a shared tracer, else raises (use r_range)."""
    d = tracer_design()
    hit = d[(d.family == family) & (d.variant == variant)]
    if hit.empty:
        raise KeyError(f"no tracer design recorded for {family} {variant}")
    if set(hit.design) != {"shared"}:
        raise ValueError(f"{family} {variant}: references differ; r is inferred")
    return 1.0


def main():
    d = table()
    print(f"{'family':6s} {'variant':14s} {'T/C':>5} {'2nd H/T':>8} {'2nd D/T':>8} "
          f"{'r':>7} {'se':>6}  source")
    for _, x in d.iterrows():
        print(f"{x.family:6s} {x.variant:14s} {x.T_C:5.1f} {x.K_HT:8.3f} {x.K_DT:8.3f} "
              f"{x.r:7.4f} {x.r_se:6.4f}  {x.source}")
    t = tracer_design()
    print(f"\ntracer design over {len(t)} forms: "
          + ", ".join(f"{k} {v}" for k, v in t.groupby("design").size().items()))
    labeled = set(map(tuple, t[t.design == "labeled"][["family", "variant"]].values))
    measured = set(map(tuple, d[["family", "variant"]].drop_duplicates().values))
    print(f"forms with secondary effects that are also 'labeled': "
          f"{len(measured & labeled)} of {len(measured)}")
    lo = d.loc[d.r.idxmin()]
    print(f"\nforms with secondary effects: {d.groupby(['family', 'variant']).ngroups}")
    print(f"r over all measured conditions: {d.r.min():.4f} to {d.r.max():.4f}")
    print(f"smallest: {lo.family} {lo.variant} at {lo.T_C:g} C, r = {lo.r:.4f} +- "
          f"{lo.r_se:.4f}, {(lo.r - 1) / lo.r_se:.1f} standard errors above unity")


if __name__ == "__main__":
    main()
