"""The published cofactor titration: a control experiment and a test of the scheme.

Tsai and Klinman, Biochemistry 40, 2303–2311 (2001), Table 5, DOI
10.1021/bi002075l. All four concentrations used the F93W horse liver enzyme,
benzyl alcohol, 3 C, 40% methanol, 100 mM cacodylate, 200 mM semicarbazide,
pH* 7.0. "Error" retains the source's description of replicate deviations.
The primary measurements are apparent alcohol V/K effects at fixed NAD+.

The source used the titration to examine the order of substrate binding and
made its other measurements at 10 mM NAD+.  The 10 mM arm is the 3 C record of
the temperature series and is counted only there; the other three
concentrations are not analysis units of the dataset.  This module (i) tests
whether one intrinsic pair can reproduce the four conditions under the
minimal scheme (scheme_test, scheme_fit: it cannot at the two subsaturating
concentrations; random_order_check: a random order of binding changes only the
commitment and would raise both effects at low cofactor), and (ii) repeats the counts with the three additional
concentrations treated as units of their own, under the same working kinetic,
labeling and log-error assumptions.  Run from analysis/ or as a script.
"""
import numpy as np
import pandas as pd

import bounds_uncertainty as bounds
import corpus
from partial_id import F_min_exact
import reference_counts as reference
from yadh_robustness import pchloro_series, t_tail_bound

RESULTS = corpus.DATA.parent / "results"
SOURCE = "TsaiKlinman2001"
DOI = "10.1021/bi002075l"
VARIANT = "F93W (40% MeOH)"
SERIES = f"LADH {VARIANT} (hydride)"
COMMERCIAL = "FDH WT (commercial) (hydride)"
EFFECTS = ["K_HT", "K_HT_se", "K_DT", "K_DT_se"]
SOURCE_EFFECTS = ["primary_HT", "primary_HT_error", "primary_DT", "primary_DT_error"]
SECONDARY = ["secondary_HT", "secondary_HT_error", "secondary_DT", "secondary_DT_error"]
FIT_RATIOS = (1.0, 1.2, 2.0)      # tracer ratios at which one intrinsic pair is fitted

# Cofactor concentrations and Michaelis constants as the sources state them (mM).
KM_NAD_F93W = 0.75    # Tsai and Klinman 2001, p. 2309: F93W in 40% methanol
YEAST_NAD = 10.0      # Cha et al. 1989, Table 2 legend; Rucker et al. 1992, Experimental Procedures
YEAST_KM_NAD = 5.0    # Cha et al. 1989, Table 2 legend (their ref. 39)


def cofactor_saturation(tests):
    """Where the scheme test fails, in units of the cofactor's Michaelis constant,
    and where the yeast pairs were measured in the same units.

    `tests` is scheme_test(...).  Returns the largest concentration at which the
    observed D/T effect falls short of the scheme's minimum by more than two
    standard errors, that concentration over K_M of the titrated enzyme, and the
    yeast cofactor concentration over its own K_M.
    """
    failing = tests[tests.shortfall_z > 2].NAD_mM
    top = float(failing.max())
    return dict(fails_up_to_mM=top, fails_up_to_KM=top / KM_NAD_F93W,
                yeast_over_KM=YEAST_NAD / YEAST_KM_NAD)


def cofactor_conditions(groups):
    """Read Table 5 and verify its common arm against both parts of Table 4."""
    t = pd.read_csv(corpus.DATA / "tsai2001_nad_titration.csv")
    if t.NAD_mM.duplicated().any() or t.NAD_mM.eq(10).sum() != 1:
        raise ValueError("expected distinct NAD+ concentrations and exactly one common 10 mM arm")
    known = groups[SERIES].query("T_C == 3")
    sec = pd.read_csv(corpus.DATA / "secondary_reference.csv")
    sec = sec[(sec.variant == VARIANT) & (sec.source == SOURCE) & (sec.T_C == 3)]
    common = t[t.NAD_mM == 10]
    for original, columns in [(known, SOURCE_EFFECTS), (sec, SECONDARY)]:
        if not np.array_equal(original[EFFECTS].to_numpy(), common[columns].to_numpy()):
            raise ValueError("Tsai Tables 4 and 5 disagree at the common arm")
    t = t.rename(columns=dict(zip(SOURCE_EFFECTS, EFFECTS)))
    return t.assign(family="LADH", variant=VARIANT, step="hydride", T_C=3,
                    source=SOURCE, source_DOI=DOI, source_table="Table 5")


def unit_name(mM):
    return f"LADH {VARIANT}, {mM:g} mM NAD+, 3 C ({SOURCE})"


def ledger(groups, challenge_labels):
    """One row per unique measurement, linked to the analysis unit and source."""
    records = []
    purity = pd.read_csv(corpus.DATA / "deuterium_purity.csv")
    dois = dict(zip(purity.source, purity.doi))
    for unit, group in groups.items():
        g = group.copy()
        g["unit"] = unit
        g["stratum"] = "NAD titration" if unit in challenge_labels else "dataset"
        if "source_DOI" not in g.columns:
            g["source_DOI"] = g.source.map(dois)
        if "source_table" not in g.columns:
            g["source_table"] = g.note
        if unit == SERIES:
            g["NAD_mM"] = 10.0
            g.loc[g.T_C == 3, "source_table"] = "Tables 4 and 5 (same determination)"
        records.append(g)
    fields = ["unit", "stratum", "family", "variant", "step", "T_C", "NAD_mM",
              *EFFECTS, "source_DOI", "source_table"]
    data = pd.concat(records, ignore_index=True)[fields]
    if data["source_DOI"].isna().any() or data.duplicated(["unit", "T_C"]).any():
        raise ValueError("a benchmark record has missing provenance or a duplicate unit and temperature")
    data["normal_order_at_center"] = (data.K_HT > data.K_DT) & (data.K_DT > 1)
    return data


def assay_bounds(t, groups):
    """Row-level bounds with each row's actual multiplicity family membership."""
    out = []
    for r in t.itertuples():
        is_common = r.NAD_mM == 10
        unit = SERIES if is_common else unit_name(r.NAD_mM)
        values = (r.K_HT, r.K_HT_se, r.K_DT, r.K_DT_se)
        record_alpha = 0.05 / (len(groups) * len(groups[unit]))
        out.append(dict(NAD_mM=r.NAD_mM, unit=unit, shared_with_table4=is_common,
                        K_HT=r.K_HT, K_HT_se=r.K_HT_se, K_DT=r.K_DT, K_DT_se=r.K_DT_se,
                        secondary_HT=r.secondary_HT, secondary_HT_error=r.secondary_HT_error,
                        secondary_DT=r.secondary_DT, secondary_DT_error=r.secondary_DT_error,
                        secondary_ratio=r.secondary_HT / r.secondary_DT,
                        gamma_observed=np.log(r.K_HT) / np.log(r.K_DT),
                        endpoint=F_min_exact(r.K_HT, r.K_DT)[0],
                        lcb_record=bounds.exact_bound(*values),
                        alpha_adjusted=record_alpha,
                        lcb_record_adjusted=bounds.exact_bound(*values, alpha=record_alpha)))
    return pd.DataFrame(out)


def scheme_test(t):
    """Can one intrinsic pair reproduce all four cofactor concentrations?

    The four conditions share the enzyme, the substrate, the temperature and
    the solvent, and therefore one intrinsic pair (x_H, x_D).  In the minimal
    scheme with ordered tracers,

        K_H = x_H (1 + c) / (x_H + c),   K_D = x_D (1 + r c) / (x_D + r c),   r >= 1,

    only the reciprocal commitment c can differ from one concentration to the
    next.  Three consequences are exact: c > K_H - 1 (since K_H < 1 + c);
    x_D >= K*, the largest observed D/T effect (an observed effect never
    exceeds the intrinsic one); and K_D increases with x_D, with c and with r.
    Every condition therefore requires

        K_D > K* K_H / (K* + K_H - 1).

    Returns that minimum for each concentration, the shortfall of the observed
    D/T effect, and the shortfall in standard errors (delta method, the three
    reported errors taken as independent standard errors).  A row whose D/T
    effect equals K* satisfies the inequality identically, because the right
    side is below K* whenever K_H > 1, and so carries no test: here the 10 mM
    and the 20 mM rows.
    """
    star = t.loc[t.K_DT.idxmax()]
    out = []
    for r in t.itertuples():
        den = star.K_DT + r.K_HT - 1
        need = star.K_DT * r.K_HT / den
        if r.K_DT == star.K_DT:
            out.append(dict(NAD_mM=r.NAD_mM, DT_minimum=need, DT_shortfall=need - r.K_DT,
                            shortfall_se=np.nan, shortfall_z=np.nan))
            continue
        d_h = star.K_DT * (star.K_DT - 1) / den ** 2
        d_star = r.K_HT * (r.K_HT - 1) / den ** 2
        se = np.sqrt(r.K_DT_se ** 2 + (d_h * r.K_HT_se) ** 2 + (d_star * star.K_DT_se) ** 2)
        out.append(dict(NAD_mM=r.NAD_mM, DT_minimum=need, DT_shortfall=need - r.K_DT,
                        shortfall_se=se, shortfall_z=(need - r.K_DT) / se))
    return pd.DataFrame(out)


def saturating_fall(t, low=0.02, high=10.0):
    """How far each primary effect falls between the saturating and the lowest
    cofactor concentration, with the difference in standard errors (the two
    reported errors taken as independent standard errors)."""
    a, b = t[t.NAD_mM == high].iloc[0], t[t.NAD_mM == low].iloc[0]
    out = {}
    for k in ("K_HT", "K_DT"):
        d = a[k] - b[k]
        out[k] = (float(d), float(d / np.hypot(a[k + "_se"], b[k + "_se"])))
    return out


def random_order_effect(k_h, k_t, cofactor, k, alcohol=1.0, second=0.5):
    """Competitive isotope effect on V/K of the alcohol when alcohol (A) and
    cofactor (B) add in random order through one ternary complex, from the
    steady state of the whole scheme; no closed form is used.  The two
    isotopologues are present at concentrations alcohol and second * alcohol:
    the ratio of their fluxes per unit concentration does not depend on
    either, so no tracer limit is needed (a vanishing tracer only costs
    floating-point precision).

    k = (k1, k2, k3, k4, k5, k6, k7, k8):
        E + A <-> EA    k1, k2        E + B <-> EB     k3, k4
        EA + B <-> EAB  k5, k6        EB + A <-> EAB   k7, k8
        EAB -> E + products           k_h (bulk isotopologue), k_t (tracer)
    """
    k1, k2, k3, k4, k5, k6, k7, k8 = k
    e, eb, ea_h, eab_h, ea_t, eab_t = range(6)
    q = np.zeros((6, 6))
    q[e, eb], q[eb, e] = k3 * cofactor, k4
    tracer = second * alcohol
    for ea, eab, k_l, conc in ((ea_h, eab_h, k_h, alcohol), (ea_t, eab_t, k_t, tracer)):
        q[e, ea], q[ea, e] = k1 * conc, k2
        q[eb, eab], q[eab, eb] = k7 * conc, k8
        q[ea, eab], q[eab, ea] = k5 * cofactor, k6
        q[eab, e] = k_l
    gen = q - np.diag(q.sum(axis=1))
    lhs = np.vstack([np.ones(6), gen.T[1:]])        # normalization replaces the balance of free enzyme
    prob = np.linalg.solve(lhs, np.r_[1.0, np.zeros(5)])
    return (k_h * prob[eab_h] / alcohol) / (k_t * prob[eab_t] / tracer)


def random_order_escape(cofactor, k):
    """Rate constant with which bound alcohol escapes from the ternary complex
    in the random scheme of random_order_effect: directly (k8), or after the
    cofactor leaves (k6) and the binary complex releases the alcohol (k2)
    before the cofactor returns (k5 [B])."""
    k1, k2, k3, k4, k5, k6, k7, k8 = k
    return k8 + k6 * k2 / (k2 + k5 * cofactor)


def random_order_check(n=400, seed=20261120):
    """A random order of addition changes only the commitment, and lowering
    the cofactor concentration can only lower it.

    For n sets of rate constants obeying detailed balance around the binding
    cycle, the steady-state effect equals x (1 + c)/(x + c), the minimal
    scheme's map, with c = random_order_escape / k_T, and it does not rise
    with the cofactor concentration (c falls as the cofactor rises).  Returns
    the largest relative deviation from that form and the largest rise of the
    effect between successive cofactor concentrations (zero or negative when
    the effect never rises).
    """
    rng = np.random.default_rng(seed)
    grid = np.logspace(-3, 3, 13)
    worst, rise = 0.0, -np.inf
    for _ in range(n):
        k1, k2, k3, k4, k5, k6, k7 = np.exp(rng.uniform(-3, 3, 7))
        k8 = k3 * k7 * k6 * k2 / (k1 * k5 * k4)       # detailed balance of the binding cycle
        k = (k1, k2, k3, k4, k5, k6, k7, k8)
        k_t = np.exp(rng.uniform(-3, 3))
        x = 1 + np.exp(rng.uniform(-2, 3))
        alcohol = np.exp(rng.uniform(-3, 3))
        vals = []
        for b in grid:
            got = random_order_effect(x * k_t, k_t, b, k, alcohol=alcohol)
            c = random_order_escape(b, k) / k_t
            worst = max(worst, abs(got / (x * (1 + c) / (x + c)) - 1))
            vals.append(got)
        rise = max(rise, float(np.max(np.diff(vals))))
    return float(worst), float(rise)


def scheme_fit(t, r=1.0):
    """Smallest chi-square of one intrinsic pair fitted to all four conditions.

    The fit has one intrinsic pair, a reciprocal commitment for each
    concentration and the stated tracer ratio r (c_D = r c_H).  For a given
    intrinsic pair each commitment is fitted on its own, so the search is over
    two parameters and is deterministic.  At r = 1 the scale of the pair is
    free (pooling_check.ladh_shared_pair_fit), which leaves five identifiable
    parameters for the eight effects; the minimum itself is unique.
    """
    from scipy.optimize import minimize, minimize_scalar
    rows = t[EFFECTS].to_numpy()
    grid = np.linspace(-8.0, 14.0, 89)

    def one(x_h, x_d, row):
        k_h, s_h, k_d, s_d = row

        def f(log_c):
            c = np.exp(log_c)
            return (((x_h * (1 + c) / (x_h + c) - k_h) / s_h) ** 2
                    + ((x_d * (1 + r * c) / (x_d + r * c) - k_d) / s_d) ** 2)
        i = int(np.argmin(f(grid)))
        lo, hi = grid[max(i - 1, 0)], grid[min(i + 1, len(grid) - 1)]
        return minimize_scalar(f, bounds=(lo, hi), method="bounded",
                               options=dict(xatol=1e-12)).fun

    def total(q):
        if np.any(np.abs(q) > 30):
            return np.inf
        return sum(one(1 + np.exp(q[0]), 1 + np.exp(q[1]), row) for row in rows)

    return float(min(minimize(total, start, method="Nelder-Mead",
                              options=dict(xatol=1e-10, fatol=1e-12)).fun
                     for start in ([2.0, -0.4], [3.0, 0.5])))


def latex_table(assays):
    """Rows of the main-text cofactor table: the measured pair, its exponent,
    the smallest D/T effect the minimal scheme allows at that H/T effect
    (scheme_test) and the measured D/T effect minus that minimum, in standard
    errors (blank where the row attains the largest D/T effect and carries no
    test).  The row shared with the temperature series carries a dagger."""
    lines = ["% Generated by analysis/cofactor_benchmark.py from Tsai and Klinman Table 5."]
    for r in assays.itertuples():
        # Source D/T was reported as 1.2(0.1) at 1 mM, two decimals elsewhere.
        decimals = 1 if r.NAD_mM == 1 else 2
        common = r"$^{\dagger}$" if r.shared_with_table4 else ""
        diff = "" if np.isnan(r.shortfall_z) else f"${-r.shortfall_z:+.1f}$"
        fields = [f"{r.NAD_mM:g}{common}", f"${r.K_HT:.1f}\\pm{r.K_HT_se:.1f}$",
                  f"${r.K_DT:.{decimals}f}\\pm{r.K_DT_se:.{decimals}f}$", f"{r.gamma_observed:.2f}",
                  f"${r.DT_minimum:.2f}$", diff]
        lines.append(" & ".join(fields) + r"\\")
    lines.append(r"\bottomrule")
    return "\n".join(lines) + "\n"


def main():
    baseline = reference.reference_groups()
    raw = cofactor_conditions(baseline)
    challenges = {unit_name(m): g for m, g in raw[raw.NAD_mM != 10].groupby("NAD_mM")}
    all_groups = {**baseline, **challenges}
    data = ledger(all_groups, challenges)
    data.to_csv(RESULTS / "matched_measurements.csv", index=False)
    assays = assay_bounds(raw, all_groups).merge(scheme_test(raw), on="NAD_mM")
    assays.to_csv(RESULTS / "nad_condition_bounds.csv", index=False, float_format="%.9g")
    (RESULTS / "sm_table_nad.tex").write_text(latex_table(assays))
    fits = pd.DataFrame([dict(tracer_ratio=r, chi2=scheme_fit(raw, r), effects=2 * len(raw))
                         for r in FIT_RATIOS])
    fits.to_csv(RESULTS / "nad_scheme_fit.csv", index=False, float_format="%.6g")
    summary = []
    for name, groups in [("dataset", baseline), ("all", all_groups),
                         ("dataset_without_commercial_FDH", {k: v for k, v in baseline.items() if k != COMMERCIAL}),
                         ("all_without_commercial_FDH", {k: v for k, v in all_groups.items() if k != COMMERCIAL})]:
        per_unit = reference.unit_table(groups)
        c = reference.counts(per_unit)
        used = data[data.unit.isin(groups)]
        summary.append(dict(screen=name, records=len(used), admissible=int(used.normal_order_at_center.sum()), **c))
        if name == "all":
            per_unit["stratum"] = np.where(per_unit.unit.isin(challenges), "NAD titration", "dataset")
            per_unit.to_csv(RESULTS / "all_condition_units.csv", index=False, float_format="%.9g")
    pd.DataFrame(summary).drop(columns="above_adj_units").to_csv(
        RESULTS / "condition_screen_sensitivity.csv", index=False)
    n_units = len(all_groups)
    yeast_benzyl = t_tail_bound(df=10, units=n_units)
    yeast_chloro = pchloro_series(n_units)["t10"]
    print(f"Adjusted yeast bounds with Student t(10) tails and {n_units} units: "
          f"benzyl {yeast_benzyl:+.6f}, p-chlorobenzyl {yeast_chloro:+.6f}")
    print(assays[["NAD_mM", "gamma_observed", "endpoint", "lcb_record", "lcb_record_adjusted"]].to_string(index=False))
    print("Test of the minimal scheme across the four concentrations (one intrinsic pair):")
    print(assays[["NAD_mM", "K_DT", "DT_minimum", "DT_shortfall", "shortfall_se", "shortfall_z"]].to_string(index=False))
    print(fits.to_string(index=False))
    fall = saturating_fall(raw)
    print(f"From 10 to 0.02 mM NAD+ the H/T effect falls by {fall['K_HT'][0]:.2f} ({fall['K_HT'][1]:.1f} s.e.) "
          f"and the D/T effect by {fall['K_DT'][0]:.2f} ({fall['K_DT'][1]:.1f} s.e.).")
    worst, rise = random_order_check()
    print(f"Random order of addition, one ternary complex: steady state equals x(1 + c)/(x + c) to {worst:.1e}; "
          f"largest rise of the effect with cofactor concentration {rise:.1e} (never rises).")
    sat = cofactor_saturation(assays)
    print(f"The scheme fails up to {sat['fails_up_to_mM']:g} mM NAD+, {sat['fails_up_to_KM']:.2f} times the Michaelis constant "
          f"of the titrated enzyme; the yeast pairs were measured at {sat['yeast_over_KM']:.1f} times theirs.")
    print(pd.DataFrame(summary).drop(columns="above_adj_units").to_string(index=False))


if __name__ == "__main__":
    main()
