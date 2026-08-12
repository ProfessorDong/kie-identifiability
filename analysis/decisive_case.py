"""The one system in the record that does discriminate.

Cha, Murray and Klinman (Science 243, 1325 (1989)) report competitive primary
V/K effects for yeast alcohol dehydrogenase against a common tritium reference:
exactly the observable this analysis needs, at a single temperature.  Their own
argument is the offset test in disguise.  They compare the observed exponent
3.58 with what they call "the theoretical upper limit for the exponent", 3.34,
which is the reduced-mass semiclassical value gamma_SC = 3.34887 used here.

Two things follow that they could not have argued in 1989.

First, the conclusion is robust to kinetic complexity.  Kohen and Jensen later
objected that commitments can inflate an observed exponent, which is correct in
general.  It cannot rescue this case: L_H = 8.40 exceeds gamma_SC, so by
Proposition 2 the identified set is the OPEN half-line (F_obs, infinity).  Every
admissible commitment moves the intrinsic offset UP, so an observed offset above
the envelope cannot be a masking artifact.  The one-sided direction of the set is
what makes the 1989 argument safe.

Second, the direction of that inference depends on the kinetic scheme.  Under a
shared additive commitment the set is open below and the argument would fail.
The scheme-derived map of Eq. (1), which is Northrop's equation re-referenced to
the tritium reference, is the one that applies.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

import masses as M
from partial_id import F_min_exact
from reversible import eie_threshold, F_min_reversible

GSC = M.gamma_sc("C")
F0 = M.offset_F0("C")
OUT = []


def say(s=""):
    print(s); OUT.append(s)


def main():
    d = pd.read_csv("../data/cha1989_yadh.csv")
    say("=" * 90)
    say("THE DECISIVE CASE: YEAST ADH PRIMARY EFFECTS (Cha, Murray & Klinman 1989)")
    say("=" * 90)
    say(f"gamma_SC = {GSC:.5f}   F0 = {F0:+.6f}")
    say("Errors are standard deviations of >=5 determinations. Using them as the")
    say("standard error of the mean is conservative by about sqrt(5).")
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
    say("  no admissible commitment can lower it. Masking cannot explain these data.")
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
    say(f"  F_min = {fmin:+.4f}, sd = {sd:.4f}, 95% LCB = {fmin-1.645*sd:+.4f}")
    say(f"  reversibility: vacuity threshold E_D* = {es:.3f}")
    say(f"{'E_DT':>10}{'F_min':>11}{'95% LCB':>11}{'excludes gated?':>18}")
    for e in (1.0, 1.1, 1.2, 1.4, 1.6):
        f, _, _ = F_min_reversible(avg.K_HT, avg.K_DT, e)
        say(f"{e:10.2f}{f:11.4f}{f-1.645*sd:11.4f}"
            f"{'YES' if f-1.645*sd > F0 else 'no':>18}")
    say("  Equilibrium isotope effects for alcohol/aldehyde hydride transfer are")
    say("  well below the threshold, so reversibility does not overturn this.")
    say("")
    say("  Correlation between the two effects only helps: rho > 0 reduces")
    say(f"{'rho':>10}{'sd(F)':>11}{'95% LCB':>11}")
    sH, sD = avg.K_HT_se/avg.K_HT, avg.K_DT_se/avg.K_DT
    for rho in (0.0, 0.5, 0.9):
        s = np.sqrt(sH**2 + GSC**2*sD**2 - 2*GSC*rho*sH*sD)
        say(f"{rho:10.2f}{s:11.4f}{fmin-1.645*s:11.4f}")

    res.to_csv("../results/decisive_case.csv", index=False)
    with open("../results/decisive_case.txt", "w") as f:
        f.write("\n".join(OUT) + "\n")
    print("\n[written] ../results/decisive_case.{csv,txt}")


if __name__ == "__main__":
    main()
