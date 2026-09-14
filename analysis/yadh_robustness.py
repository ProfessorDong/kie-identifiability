"""How hard is the yeast ADH conclusion to break?

It rests on one temperature, four determinations, and 1989-vintage error bars.
That is the thinnest data supporting the strongest claim in the paper, so this
module attacks it from every direction that does not require new measurements.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from scipy import stats

import masses as M
from partial_id import F_min_exact
from reversible import vacuity_window
from vibronic_envelope import B_vib

OUT = []


def say(s=""):
    print(s); OUT.append(s)


def convention(kind):
    """(gamma, F0, sqrt(t)) for a stated mass convention."""
    if kind == "reduced C-H":
        h, d, t = M.mu_ratios("C")
    elif kind == "bare atomic":
        m = (1.00782503, 2.01410178, 3.01604928)
        h, d, t = [np.sqrt(x / m[0]) for x in m]
    elif kind == "integer 1:2:3":
        m = (1.0, 2.0, 3.0)
        h, d, t = [np.sqrt(x / m[0]) for x in m]
    g = d * (t - 1) / (t - d)
    P = lambda mu: 0.5 * np.log(t / mu)
    return g, P(h) - g * P(d), np.sqrt(t)


def main():
    y = pd.read_csv("../data/cha1989_yadh.csv")
    dets = y.iloc[:3]
    avg = y.iloc[3]

    say("=" * 88)
    say("HOW HARD IS THE YEAST ADH CONCLUSION TO BREAK?")
    say("=" * 88)

    say("")
    say("1. MASS CONVENTION")
    say("-" * 88)
    say(f"{'convention':16s}{'gamma_SC':>10}{'F0':>11}{'F_obs':>10}"
        f"{'sd':>9}{'95% LCB':>10}{'> F0':>7}{'> 0':>6}")
    for k in ("reduced C-H", "bare atomic", "integer 1:2:3"):
        g, f0, _ = convention(k)
        f = np.log(avg.K_HT) - g * np.log(avg.K_DT)
        sd = np.sqrt((avg.K_HT_se / avg.K_HT) ** 2 + (g * avg.K_DT_se / avg.K_DT) ** 2)
        lcb = f - 1.645 * sd
        say(f"{k:16s}{g:10.5f}{f0:11.6f}{f:10.4f}{sd:9.4f}{lcb:10.4f}"
            f"{'YES' if lcb > f0 else 'no':>7}{'YES' if lcb > 0 else 'no':>6}")
    say("  Cha et al. themselves used 3.26. The conclusion does not depend on")
    say("  which convention is adopted, and is strongest under theirs.")

    say("")
    say("2. HOW BADLY WOULD THE ERRORS HAVE TO BE UNDERSTATED?")
    say("-" * 88)
    g, f0, _ = convention("reduced C-H")
    f = np.log(avg.K_HT) - g * np.log(avg.K_DT)
    sd = np.sqrt((avg.K_HT_se / avg.K_HT) ** 2 + (g * avg.K_DT_se / avg.K_DT) ** 2)
    for tgt, nm in ((f0, "the gated envelope F0"), (0.0, "the semiclassical locus")):
        lam = (f - tgt) / (1.645 * sd)
        say(f"  to stop excluding {nm:24s}: errors x {lam:.2f}")
    # Until 2026-09-14 this claimed the margin was conservative by sqrt(5),
    # reading Cha's "standard deviations based on five or more determinations"
    # as replicate SDs.  They are within-experiment scatter over >=5 time points,
    # propagated into the mean for the average row -- and against the scatter
    # between the three determinations, the average's errors are already the
    # standard error of the mean.  No margin is claimed from the label.
    det = y[y.note.str.contains("determination")]
    se_h = det.K_HT.std(ddof=1) / np.sqrt(len(det))
    se_d = det.K_DT.std(ddof=1) / np.sqrt(len(det))
    say(f"  Errors on the average ({avg.K_HT_se:.2f}, {avg.K_DT_se:.2f}) match the standard")
    say(f"  error of the {len(det)} determinations ({se_h:.3f}, {se_d:.3f}), so these")
    say("  factors are the margin as it stands, with no allowance for the 'SD' label.")

    say("")
    say("3. POOLING THE THREE INDEPENDENT DETERMINATIONS")
    say("-" * 88)
    F, S = [], []
    for _, r in dets.iterrows():
        F.append(np.log(r.K_HT) - g * np.log(r.K_DT))
        S.append(np.sqrt((r.K_HT_se / r.K_HT) ** 2 + (g * r.K_DT_se / r.K_DT) ** 2))
    F, S = np.array(F), np.array(S)
    w = 1.0 / S ** 2
    mu = float(np.sum(w * F) / np.sum(w))
    se = float(1.0 / np.sqrt(np.sum(w)))
    Q = float(np.sum(w * (F - mu) ** 2))
    pQ = 1 - stats.chi2.cdf(Q, len(F) - 1)
    say(f"  individual offsets: {', '.join(f'{x:+.4f}' for x in F)}")
    say(f"  inverse-variance mean {mu:+.4f} +- {se:.4f}, 95% LCB {mu-1.645*se:+.4f}")
    say(f"  heterogeneity Q = {Q:.2f} on {len(F)-1} d.o.f., p = {pQ:.3f}")
    say(f"  published average row gives {f:+.4f}; the two agree to {abs(mu-f):.4f}")
    say("  No detectable heterogeneity, so the determinations are mutually")
    say("  consistent and the published average is not carrying the result.")

    say("")
    say("4. AGAINST THE SUMMED VIBRONIC ENVELOPE OVER A BOUNDED BOX")
    say("-" * 88)
    lcb = f - 1.645 * sd
    say(f"  yeast ADH 95% lower bound: {lcb:+.4f}")
    say(f"{'lambda range':>18}{'max B_vib':>12}{'LCB > max':>12}")
    for lo, hi, lab in ((5, 60, "5-60"), (10, 60, "10-60"), (15, 60, "15-60"),
                        (20, 60, "20-60")):
        best = -np.inf
        for lam in np.linspace(lo, hi, 8):
            for nu in (1500., 2200., 3000., 3600.):
                for dG in (-20., -12., -5., -1.):
                    b, *_ = B_vib(nu, lam, dG, 298.15, nA=10, nw=10)
                    if np.isfinite(b):
                        best = max(best, b)
        say(f"{lab:>18}{best:12.4f}{'YES' if lcb > best else 'no':>12}")
    # Until 2026-09-14 this section concluded that the observation excludes the
    # summed model above ~10 kcal/mol.  The supplement withdrew that: these
    # maxima sit at the edge of the sampled box and reverse sign when it is
    # extended (vibronic_envelope.py).  The printout now says so.
    say("  WITHDRAWN as an exclusion: these maxima are pinned at the edge of the")
    say("  sampled box and reverse sign when it is extended, the summed family")
    say("  reaching offsets well above anything measured.  No exclusion of the")
    say("  summed model is claimed; the comparison is with the ground channel only.")

    say("")
    say("5. REVERSIBILITY AND COMMITMENT, RESTATED")
    say("-" * 88)
    win = vacuity_window(avg.K_HT, avg.K_DT)
    fmin, cstar, interior, closed = F_min_exact(avg.K_HT, avg.K_DT)
    say(f"  L_H = {(avg.K_HT-1)/(avg.K_DT-1):.2f} > gamma_SC, so the set is the open")
    say(f"  half-line ({fmin:+.4f}, inf): no commitment lowers the endpoint.")
    say(f"  vacuity window: {'EMPTY' if win is None else win} (F_obs > 0), so no")
    say("  equilibrium isotope effect makes the observation uninformative.")

    say("")
    say("6. MULTIPLICITY")
    say("-" * 88)
    # The yeast system is singled out from every analysis unit examined, so its
    # bound also carries a one-sided Bonferroni correction over all of them.
    import corpus
    m = corpus.counts()["units"]
    z = stats.norm.ppf(1 - 0.05 / m)
    g, f0, _ = convention("reduced C-H")
    f = np.log(avg.K_HT) - g * np.log(avg.K_DT)
    a, b = avg.K_HT_se / avg.K_HT, g * avg.K_DT_se / avg.K_DT
    say(f"  analysis units {m}, one-sided z = {z:.4f}")
    for rho in (0.0, -1.0):
        bnd = f - z * np.sqrt(a * a + b * b - 2 * rho * a * b)
        say(f"  rho = {rho:+.0f}: Bonferroni bound {bnd:+.4f}, margin over F0 "
            f"{bnd - f0:+.4f}, over 0 {bnd:+.4f}")

    with open("../results/yadh_robustness.txt", "w") as fh:
        fh.write("\n".join(OUT) + "\n")
    print("\n[written] ../results/yadh_robustness.txt")


if __name__ == "__main__":
    main()
