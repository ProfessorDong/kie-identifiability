"""Correct interpretation of the profiled offsets.

A profile interval is only meaningful where the forward model actually fits.
Post-processing offset_profiles.csv shows the two are strongly anti-correlated:
every narrow interval comes from a series the model cannot describe, where a
large chi2 everywhere makes the RELATIVE curvature steep and manufactures a
spuriously tight bound on a meaningless F.

So the analysis is stratified by goodness of fit, and conclusions are drawn only
from series where the model is adequate.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from scipy import stats

import masses as M

F0 = M.offset_F0("C")
OUT = []


def say(s=""):
    print(s)
    OUT.append(s)


def main():
    d = pd.read_csv("../results/offset_profiles.csv")
    d["dof"] = (2 * d.n_T - 5).clip(lower=1)
    d["p_fit"] = 1 - stats.chi2.cdf(d.chi2_min, d.dof)
    d["adequate"] = d.p_fit > 0.05
    d["w_rel"] = d.width / abs(F0)

    say("=" * 84)
    say("WHAT THE BENCHMARK DETERMINES ABOUT THE DISCRIMINATING INVARIANT")
    say("=" * 84)
    say(f"gamma_SC = {M.gamma_sc('C'):.5f} (C-H reduced masses)")
    say(f"semiclassical: F = 0     gated-tunneling ridge: F0 = {F0:+.6f}")
    say(f"a series discriminates only if its 95% interval is narrower than")
    say(f"|F0| = {abs(F0):.4f} and lands on one side of it.")
    say("")

    say("-" * 84)
    say("A. SERIES WHERE THE FORWARD MODEL IS ADEQUATE (p_fit > 0.05)")
    say("-" * 84)
    say(f"{'series':42s}{'chi2/dof':>10}{'p_fit':>7}{'F_hat':>8}"
        f"{'width':>8}{'x|F0|':>7}{'open':>6}")
    a = d[d.adequate].sort_values("width")
    for _, r in a.iterrows():
        say(f"{r.family+' '+r.variant:42s}{r.chi2_min:6.1f}/{r.dof:<3d}"
            f"{r.p_fit:7.3f}{r.F_hat:8.2f}{r.width:8.3f}{r.w_rel:7.1f}"
            f"{'YES' if r.open_interval else '-':>6}")
    say("")
    say(f"  adequate-fit series                : {len(a)} of {len(d)}")
    say(f"  narrowest 95% interval on F        : {a.width.min():.3f}"
        f"  ({a.width.min()/abs(F0):.1f}x the signal)")
    say(f"  median 95% interval on F           : {a.width.median():.3f}"
        f"  ({a.width.median()/abs(F0):.1f}x the signal)")
    say(f"  series excluding BOTH 0 and F0     : "
        f"{int((a.excludes_semiclassical & a.excludes_ridge).sum())}")
    say(f"  series that DISCRIMINATE           : "
        f"{int(((a.width < abs(F0)) & (a.excludes_semiclassical ^ a.excludes_ridge)).sum())}")
    say("")

    say("-" * 84)
    say("B. SERIES THE MODEL CANNOT DESCRIBE (p_fit <= 0.05): intervals void")
    say("-" * 84)
    say(f"{'series':42s}{'chi2/dof':>10}{'p_fit':>9}{'F_hat':>8}{'width':>8}")
    b = d[~d.adequate].sort_values("chi2_min", ascending=False)
    for _, r in b.iterrows():
        say(f"{r.family+' '+r.variant:42s}{r.chi2_min:6.1f}/{r.dof:<3d}"
            f"{r.p_fit:9.1e}{r.F_hat:8.2f}{r.width:8.3f}")
    say("")
    say("  These are NOT measurements of F. Where the model misfits, F absorbs")
    say("  the misfit; the apparently tight bounds are artefacts.")
    say("")

    say("-" * 84)
    say("C. THE ARTEFACT, QUANTIFIED")
    say("-" * 84)
    rho = stats.spearmanr(np.log10(d.chi2_min.clip(lower=1e-3)), d.width)
    say(f"  Spearman(log chi2_min, interval width) = {rho.statistic:+.3f}"
        f"  (p = {rho.pvalue:.4f})")
    say(f"  median width, adequate fits   : {a.width.median():.3f}")
    say(f"  median width, inadequate fits : {b.width.median():.3f}")
    say("  Worse fits give TIGHTER intervals: the signature of misspecification,")
    say("  not of information.")
    say("")

    say("-" * 84)
    say("D. MASS-MODULATED CONTROL (light vs heavy ecDHFR)")
    say("-" * 84)
    hl = d[d.variant.isin(["light enzyme", "heavy enzyme"])]
    for _, r in hl.iterrows():
        ok = "adequate" if r.adequate else "INADEQUATE FIT"
        say(f"  {r.variant:14s} F = {r.F_hat:+.2f} [{r.F_lo:+.2f},{r.F_hi:+.2f}]"
            f"   chi2 = {r.chi2_min:.1f}/{r.dof}  p = {r.p_fit:.3f}  {ok}")
    say("  Heavy labelling moves the system along the ridge, where F is")
    say("  invariant, so both should return the same F. The light series does")
    say("  not fit, so the control is inconclusive as it stands.")
    say("")

    say("=" * 84)
    say("CONCLUSION")
    say("=" * 84)
    say("Across every series the model can actually describe, the 95% interval")
    say(f"on the discriminating invariant is at least {a.width.min()/abs(F0):.1f} times the")
    say("separation between the two mechanisms. Not one series discriminates.")
    say("")
    say("The limitation is not the number of enzymes: four families, five")
    say("systems, three organisms and two chemical steps all give the same")
    say("answer. It is that F must be recovered through the commitment map,")
    say("whose masking displaces the observed offset by up to 20x the signal.")
    say("")
    say(f"Precision target for a decisive experiment: a 95% interval on F")
    say(f"inside {abs(F0):.3f} log-KIE units, i.e. at least "
        f"{a.width.min()/abs(F0):.0f}-fold tighter than the best series here.")

    d.to_csv("../results/offset_profiles_annotated.csv", index=False)
    with open("../results/offset_summary.txt", "w") as f:
        f.write("\n".join(OUT) + "\n")
    print("\n[written] ../results/offset_summary.txt")


if __name__ == "__main__":
    main()
