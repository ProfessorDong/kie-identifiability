"""Conditional reanalysis of Cha, Murray and Klinman's yeast ADH measurements.

The original H/T and D/T experiments use different tritiated preparations.
The primary-pair test assumes ordered tracer commitments, a shared primary
chemical effect, isotope-insensitive binding and one irreversible transfer step.
The reversible example is a separate, hypothetical shared-tracer calculation.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

import masses as M
from partial_id import F_min_exact
from reversible import eie_threshold, F_min_reversible, vacuity_window

GSC = M.gamma_sc("C")
F0 = M.offset_F0("C")
OUT = []


def say(s=""):
    print(s); OUT.append(s)


def main():
    d = pd.read_csv("../data/cha1989_yadh.csv")
    say("=" * 90)
    say("CONDITIONAL PRIMARY-PAIR TEST: YEAST ADH PRIMARY EFFECTS (Cha, Murray & Klinman 1989)")
    say("=" * 90)
    say(f"gamma_SC = {GSC:.5f}   F0 = {F0:+.6f}")
    say("Errors as published. Cha et al. call them standard deviations of >=5")
    say("determinations (time points within one experiment); the errors on the")
    say("average are consistent with an estimate of the standard")
    say("error. Bounds use normal log errors as a working model.")
    say("")
    say(f"{'record':28s}{'K_HT':>13}{'K_DT':>13}{'gam_obs':>9}{'L_H':>7}"
        f"{'F_min':>9}{'sd':>8}{'95% LCB':>9}{'>F0':>6}{'>0':>5}")
    rows = []
    for _, r in d.iterrows():
        fmin, cstar, interior, closed = F_min_exact(r.K_HT, r.K_DT)
        sd = np.sqrt((r.K_HT_se / r.K_HT) ** 2 + (GSC * r.K_DT_se / r.K_DT) ** 2)
        lcb = fmin - 1.645 * sd
        say(f"{r.note:28s}{r.K_HT:8.2f}+-{r.K_HT_se:<4.2f}"
            f"{r.K_DT:8.2f}+-{r.K_DT_se:<4.2f}"
            f"{np.log(r.K_HT)/np.log(r.K_DT):9.3f}"
            f"{(r.K_HT-1)/(r.K_DT-1):7.2f}{fmin:9.4f}{sd:8.4f}{lcb:9.4f}"
            f"{'YES' if lcb > F0 else 'no':>6}{'YES' if lcb > 0 else 'no':>5}")
        rows.append(dict(note=r.note, F_min=fmin, sd=sd, lcb=lcb,
                         interior=interior, closed=closed,
                         excl_gated=bool(lcb > F0), excl_semiclassical=bool(lcb > 0)))
    res = pd.DataFrame(rows)
    say("")
    say(f"  All {len(res)} records have L_H > gamma_SC, so every identified set is the")
    say("  OPEN half-line (F_obs, inf): the endpoint is the commitment-free value and")
    say("  ordered forward commitments alone cannot lower it under the assumed tracer model.")
    say(f"  {int(res.excl_gated.sum())}/{len(res)} exclude the ground-channel gated model;"
        f" {int(res.excl_semiclassical.sum())}/{len(res)} also exclude the")
    say("  semiclassical locus F = 0.")

    avg = d.iloc[-1]
    say("")
    say("-" * 90)
    say("ROBUSTNESS OF THE REPORTED AVERAGE")
    say("-" * 90)
    fmin = F_min_exact(avg.K_HT, avg.K_DT)[0]
    sd = np.sqrt((avg.K_HT_se / avg.K_HT) ** 2 + (GSC * avg.K_DT_se / avg.K_DT) ** 2)
    es = eie_threshold(avg.K_HT, avg.K_DT)
    win = vacuity_window(avg.K_HT, avg.K_DT)
    say(f"  F_min = {fmin:+.4f}, sd = {sd:.4f}, 95% LCB = {fmin-1.645*sd:+.4f}")
    say(f"  Shared-tracer reversible model: algebraic crossing E_D* = {es:.3f}, vacuity window = "
        f"{'EMPTY' if win is None else win}")
    say("  For a shared tracer and mass-scaled equilibrium isotope effects, a positive")
    say("  observed offset bounds the intrinsic offset at every E_D (SI Note 13).")
    say(f"{'E_DT':>10}{'F_min':>11}{'95% LCB':>11}{'excludes gated?':>18}")
    for e in (1.0, 1.1, 1.2, 1.4, 1.6, 2.0, 5.0, 10.0):
        f, _, _ = F_min_reversible(avg.K_HT, avg.K_DT, e)
        say(f"{e:10.2f}{f:11.4f}{f-1.645*sd:11.4f}"
            f"{'YES' if f-1.645*sd > F0 else 'no':>18}")
    say("  The measured mixed-tracer pair requires additional information linking")
    say("  forward and reverse commitments before this reversible result can apply.")
    say("")
    say("  Correlation acts in both directions: rho > 0 shrinks sd(F), rho < 0")
    say("  inflates it.  The paper reports the maximally adverse rho = -1.")
    say(f"{'rho':>10}{'sd(F)':>11}{'95% LCB':>11}")
    sH, sD = avg.K_HT_se/avg.K_HT, avg.K_DT_se/avg.K_DT
    for rho in (-1.0, 0.0, 0.5, 0.9):
        s = np.sqrt(sH**2 + GSC**2*sD**2 - 2*GSC*rho*sH*sD)
        say(f"{rho:10.2f}{s:11.4f}{fmin-1.645*s:11.4f}")

    res.to_csv("../results/decisive_case.csv", index=False)
    with open("../results/decisive_case.txt", "w") as f:
        f.write("\n".join(OUT) + "\n")
    print("\n[written] ../results/decisive_case.{csv,txt}")


if __name__ == "__main__":
    main()
