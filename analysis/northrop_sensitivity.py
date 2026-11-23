"""How intrinsic effects obtained by Northrop's method depend on the exponent it assumes.

Northrop's method closes the two equations of the minimal scheme (one tritiated
tracer for both comparisons) by imposing ln x_H = gamma ln x_D on the intrinsic
pair and solving for the commitment.  For a pair observed below gamma a
solution exists for every gamma above the observed exponent, so the pair cannot
test the gamma it is given.  This module asks how much the answer depends on
that choice:

  magnitude    the relative change in the derived intrinsic H/T effect when
               gamma changes by 0.1 near the semiclassical value (every record
               of the dataset with a solution at gamma_SC - 0.05 and at
               gamma_SC + 0.05);
  temperature  the isotope effect on the activation energy, Ea(T) - Ea(H), and
               the prefactor ratio A_H/A_T, fitted (Arrhenius, least squares in
               ln x_H against 1/T) to the derived intrinsic H/T effects of each
               temperature series with at least three temperatures, at a
               constant gamma of 3.26, gamma_SC and 3.50.  A series enters
               only if the method has a solution at every temperature for all
               three values.

  worked pair  the reductase pair of the main text's worked example
               (K_HT = 3.06, K_DT = 1.51) at gamma_SC and 0.1 either side, with
               the commitments the method implies;
  tracer ratio for mixed-label records, the change in the derived intrinsic
               H/T effect when the shared-tracer form (r = 1) is replaced by
               the ordered form with r set to the lower bound R that the
               reported secondary effects at the same condition give.

Writes results/northrop_sensitivity.txt and results/northrop_sensitivity.csv.
Run from analysis/:  python3 northrop_sensitivity.py
"""
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import brentq

import masses as M
import reference_counts as RC

GSC = M.gamma_sc("C")
GAMMAS = (3.26, GSC, 3.50)            # constant exponents compared over temperature
STEP = 0.1                            # change in the exponent for the magnitude test
R_KCAL = 8.314462618 / 4184.0         # gas constant, kcal mol^-1 K^-1
RESULTS = Path(__file__).resolve().parent.parent / "results"


def intrinsic_ht(k_ht, k_dt, gamma):
    """Intrinsic H/T effect returned by Northrop's method at exponent gamma.

    Shared-tracer map x = K / (1 - (K - 1) C) for each comparison, with the
    commitment C chosen so that ln x_H / ln x_D = gamma.  Returns nan when the
    pair is not normally ordered or when no positive commitment exists (the
    observed exponent is already at or above gamma)."""
    if not (k_ht > k_dt > 1):
        return np.nan

    def x(k, c):
        return k / (1 - (k - 1) * c)

    def f(c):
        return np.log(x(k_ht, c)) / np.log(x(k_dt, c)) - gamma

    lo, hi = 1e-14, (1 - 1e-12) / (k_ht - 1)
    if f(lo) >= 0:
        return np.nan
    return x(k_ht, brentq(f, lo, hi, xtol=1e-15))


def intrinsic_ht_ordered(k_ht, k_dt, gamma, r):
    """Intrinsic H/T effect returned by Northrop's method when the two tracers
    differ: reciprocal commitment c for the H/T tracer and r c for the D/T
    tracer.  The root in c is bracketed on a log grid; nan when the pair is not
    normally ordered or has no solution.  Raises if the root is not unique."""
    if not (k_ht > k_dt > 1):
        return np.nan

    def g(c):
        x_h = k_ht * c / (1 + c - k_ht)
        x_d = k_dt * r * c / (1 + r * c - k_dt)
        return np.log(x_h) / np.log(x_d) - gamma

    cs = (k_ht - 1) * (1 + np.logspace(-12, 8, 4001))
    v = np.array([g(c) for c in cs])
    idx = np.where(np.sign(v[:-1]) != np.sign(v[1:]))[0]
    if len(idx) == 0:
        return np.nan
    if len(idx) > 1:
        raise ValueError(f"several roots for ({k_ht}, {k_dt}) at r = {r}")
    c = brentq(g, cs[idx[0]], cs[idx[0] + 1], xtol=1e-14)
    return k_ht * c / (1 + c - k_ht)


def commitments_at(k_ht, k_dt, gamma):
    """(C_T, C_H, C_D) implied by Northrop's method at exponent gamma for a
    shared tracer: C_T from x_H = K/(1 - (K - 1) C_T), C_L = x_L C_T."""
    x_h = intrinsic_ht(k_ht, k_dt, gamma)
    c_t = (1 - k_ht / x_h) / (k_ht - 1)
    x_d = k_dt / (1 - (k_dt - 1) * c_t)
    return float(c_t), float(x_h * c_t), float(x_d * c_t)


def worked_pair(k_ht=3.06, k_dt=1.51, step=STEP, center=GSC):
    """The worked reductase pair: intrinsic H/T effect at center - step, center
    and center + step, and the commitments implied at center."""
    lo, mid, hi = (intrinsic_ht(k_ht, k_dt, center + d) for d in (-step, 0.0, step))
    c_t, c_h, c_d = commitments_at(k_ht, k_dt, center)
    # the source applied the method at gamma = 3.34 and reported C_f = 1.52 at 5 C
    c_h_334 = commitments_at(k_ht, k_dt, 3.34)[1]
    return dict(x_lo=lo, x_mid=mid, x_hi=hi, C_T=c_t, C_H=c_h, C_D=c_d, C_H_334=c_h_334)


def tracer_ratio_bias(gamma=GSC):
    """Percentage change in the intrinsic H/T effect from replacing r = 1 by the
    lower bound R (secondary H/T over secondary D/T effect at the same
    condition) for every mixed-label record with such secondary effects.
    Returns a table with one row per record that has both solutions, and the
    number of mixed-label records with secondary effects."""
    import reference_asymmetry as RA
    sec = RA.table()
    design = RA.tracer_design()
    labeled = set(map(tuple, design[design.design == "labeled"][["family", "variant"]].values))
    rows, n = [], 0
    for unit, df in RC.reference_groups().items():
        for r in df.itertuples():
            if (r.family, r.variant) not in labeled:
                continue
            hit = sec[(sec.family == r.family) & (sec.variant == r.variant) & np.isclose(sec.T_C, r.T_C)]
            if len(hit) != 1:
                continue
            n += 1
            R = float(hit.r.iloc[0])
            a = intrinsic_ht(r.K_HT, r.K_DT, gamma)
            b = intrinsic_ht_ordered(r.K_HT, r.K_DT, gamma, R)
            if np.isfinite(a) and np.isfinite(b):
                rows.append(dict(unit=unit, T_C=r.T_C, R=R, x_shared=a, x_ordered=b, pct=100 * (b / a - 1)))
    return pd.DataFrame(rows), n


def magnitude_sensitivity(step=STEP, center=GSC):
    """Percentage change in the derived intrinsic H/T effect for a change of
    `step` in the assumed exponent, centered on `center`, for every record with
    a solution at both ends.  Returns (values, number of records examined)."""
    groups = RC.reference_groups()
    vals, n = [], 0
    for df in groups.values():
        for r in df.itertuples():
            n += 1
            a = intrinsic_ht(r.K_HT, r.K_DT, center - step / 2)
            b = intrinsic_ht(r.K_HT, r.K_DT, center + step / 2)
            if np.isfinite(a) and np.isfinite(b):
                vals.append(100 * (b / a - 1))
    return np.array(vals), n


def temperature_sensitivity(gammas=GAMMAS, min_temperatures=3):
    """Arrhenius parameters of the derived intrinsic H/T effects of each series
    at each constant exponent in `gammas`.  Returns (table, number of series
    with at least `min_temperatures` temperatures)."""
    groups = RC.reference_groups()
    rows, eligible = [], 0
    for unit, df in groups.items():
        if len(df) < min_temperatures:
            continue
        eligible += 1
        inv_t = 1 / (df.T_C.to_numpy(float) + 273.15)
        fit = {}
        for g in gammas:
            xh = np.array([intrinsic_ht(a, b, g) for a, b in zip(df.K_HT, df.K_DT)])
            if not np.isfinite(xh).all():
                fit = None
                break
            slope, intercept = np.polyfit(inv_t, np.log(xh), 1)
            fit[g] = (R_KCAL * slope, np.exp(intercept))
        if fit is None:
            continue
        row = dict(unit=unit, n_T=len(df))
        for g, (dea, ahat) in fit.items():
            tag = f"{g:.2f}" if abs(g - GSC) > 1e-9 else "sc"
            row[f"dEa_{tag}"] = dea
            row[f"AH_AT_{tag}"] = ahat
        rows.append(row)
    t = pd.DataFrame(rows)
    t["dEa_change"] = t["dEa_3.50"] - t["dEa_3.26"]
    t["AH_AT_factor"] = t["AH_AT_3.50"] / t["AH_AT_3.26"]
    return t, eligible


def summary():
    vals, n_records = magnitude_sensitivity()
    t, eligible = temperature_sensitivity()
    bias, n_mixed = tracer_ratio_bias()
    wp = worked_pair()
    return dict(
        worked_lo=wp["x_lo"], worked_mid=wp["x_mid"], worked_hi=wp["x_hi"],
        worked_C_T=wp["C_T"], worked_C_H=wp["C_H"], worked_C_D=wp["C_D"],
        n_mixed=n_mixed, n_mixed_solved=len(bias),
        bias_median=float(bias.pct.median()), bias_min=float(bias.pct.min()), bias_max=float(bias.pct.max()),
        n_mixed_units=int(bias.unit.nunique()),
        n_records=n_records, n_solved=len(vals),
        pct_median=float(np.median(vals)), pct_lo=float(vals.min()), pct_hi=float(vals.max()),
        n_series=len(t), n_eligible=eligible,
        dEa_change_median=float(t.dEa_change.median()),
        dEa_change_max=float(t.dEa_change.abs().max()),
        factor_lo=float(t.AH_AT_factor.min()), factor_hi=float(t.AH_AT_factor.max()),
        factor_median=float(t.AH_AT_factor.median()),
    ), t


def _label(unit):
    """Series name as the supplement's tables print it, with the species named
    where one family spans two (corpus.display_family)."""
    import re
    import corpus
    fam, var, step = re.match(r"(\S+) (.*) \((\w+)\)$", unit).groups()
    name = f"{corpus.display_family(fam, var)} {var}" + ("" if step == "hydride" else f" ({step})")
    return name.replace("wild type", "WT").replace("%", "\\%").replace("Mg2+", "Mg$^{2+}$")


def latex_table(t):
    """Rows of the supplement's table: each series, its number of temperatures,
    the isotope effect on the activation energy at the three constant exponents,
    and the factor on A_H/A_T from 3.26 to 3.50."""
    lines = ["% Generated by analysis/northrop_sensitivity.py -- do not edit; edit the generator."]
    for _, r in t.sort_values("unit").iterrows():
        lines.append(f"{_label(r['unit'])} & {int(r['n_T'])} & ${r['dEa_3.26']:+.2f}$ & ${r['dEa_sc']:+.2f}$ & "
                     f"${r['dEa_3.50']:+.2f}$ & ${r['AH_AT_factor']:.2f}$\\\\")
    return "\n".join(lines) + "\n"


def main():
    s, t = summary()
    lines = [
        "Northrop's method: dependence of the derived intrinsic H/T effect on the assumed exponent",
        f"  gamma_SC = {GSC:.5f}; shared-tracer map; Arrhenius fits in ln x_H against 1/T",
        "",
        f"Magnitude: {s['n_solved']} of {s['n_records']} records have a solution at gamma_SC -/+ {STEP / 2:.2f}.",
        f"  change in x_H for a change of {STEP} in gamma: median {s['pct_median']:.1f}%, "
        f"range {s['pct_lo']:.1f}% to {s['pct_hi']:.1f}%",
        "",
        f"Temperature dependence: {s['n_series']} of the {s['n_eligible']} series with three or more temperatures "
        "have a solution at every temperature for gamma = 3.26, gamma_SC and 3.50.",
        f"  change in Ea(T) - Ea(H) from gamma = 3.26 to 3.50: median {s['dEa_change_median']:+.2f} kcal/mol, "
        f"largest {s['dEa_change_max']:.2f} kcal/mol",
        f"  factor on A_H/A_T from gamma = 3.26 to 3.50: median {s['factor_median']:.2f}, "
        f"range {s['factor_lo']:.2f} to {s['factor_hi']:.2f}",
        "",
        f"Worked pair (3.06, 1.51): x_H = {s['worked_lo']:.3f}, {s['worked_mid']:.3f}, {s['worked_hi']:.3f} "
        f"at gamma_SC - {STEP}, gamma_SC, gamma_SC + {STEP}; at gamma_SC C_T = {s['worked_C_T']:.4f}, "
        f"C_H = {s['worked_C_H']:.3f}, C_D = {s['worked_C_D']:.3f}",
        f"Tracer ratio: {s['n_mixed_solved']} of {s['n_mixed']} mixed-label records with secondary effects "
        f"({s['n_mixed_units']} units) have both solutions; r = R instead of 1 changes x_H by "
        f"median {s['bias_median']:+.1f}%, range {s['bias_min']:+.1f}% to {s['bias_max']:+.1f}%",
        "",
        t.to_string(index=False, float_format=lambda v: f"{v:+.3f}"),
    ]
    (RESULTS / "northrop_sensitivity.txt").write_text("\n".join(lines) + "\n")
    t.to_csv(RESULTS / "northrop_sensitivity.csv", index=False, float_format="%.6g")
    (RESULTS / "sm_table_northrop.tex").write_text(latex_table(t))
    print("\n".join(lines[:15]))
    print("[written] results/northrop_sensitivity.txt, northrop_sensitivity.csv, sm_table_northrop.tex")


if __name__ == "__main__":
    main()
