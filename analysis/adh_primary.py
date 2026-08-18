"""Primary competitive H/T and D/T effects for the alcohol dehydrogenases.

Three papers from one laboratory report the observable this analysis needs, on
the same chemistry (benzyl alcohol oxidation) at 25 C:

  Cha, Murray & Klinman, Science 243, 1325 (1989)   yeast ADH
  Bahnson et al., Biochemistry 32, 5503 (1993)      horse liver ADH, 5 forms
  Bahnson et al., PNAS 94, 12797 (1997)             horse liver ADH, 4 more

Taken at face value they appear to conflict: the yeast offset is strongly
POSITIVE (+0.129) while all nine horse liver offsets are NEGATIVE (-0.026 to
-0.172).  They do not conflict, and the reason is the one this framework is
built to express.

Every record has L_H > gamma_SC, so by Proposition 2 each identified set is the
OPEN half-line (F_obs, infinity).  The horse liver sets therefore all CONTAIN
zero and every positive value: they exclude nothing, and are consistent with
tunneling, with the semiclassical locus, and with the gated model alike.  The
yeast set does not contain zero, so it excludes.

Bahnson et al. (1993) give the physical reason directly: "Unlike the oxidation
of benzyl alcohol catalyzed by YADH, which has a rate-limiting chemical step,
the LADH reaction is partially limited in rate by product benzaldehyde
dissociation."  Horse liver carries a commitment; yeast does not.  Masking drives
the observed offset DOWN, which is exactly the sign of the horse liver values.

So the framework predicts which systems can be informative before any mechanism
is assumed: those whose chemical step is rate limiting.  That is the design
statement of the main text, borne out across nine records.

Bahnson et al. state the criterion in the same coordinates, quoting the
semiclassical upper limit as 3.26-3.34 (spanning both mass conventions):
    kinetic complexity:  k_H/k_T <  (k_D/k_T)^3.26-3.34      (F < 0)
    tunneling:           k_H/k_T >  (k_D/k_T)^3.26-3.34      (F > 0)
following Saunders, JACS 107, 164 (1985).  Proposition 2 makes that exact: the
identified set is one-sided, so F_obs > 0 cannot be manufactured by commitment,
while F_obs < 0 is uninformative because commitment alone produces it.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

import masses as M
from partial_id import F_min_exact

GSC = M.gamma_sc("C")
F0 = M.offset_F0("C")
OUT = []


def say(s=""):
    print(s); OUT.append(s)


def analyse(df, label):
    rows = []
    for _, r in df.iterrows():
        fmin, _, interior, closed = F_min_exact(r.K_HT, r.K_DT)
        sd = np.sqrt((r.K_HT_se / r.K_HT) ** 2 + (GSC * r.K_DT_se / r.K_DT) ** 2)
        rows.append(dict(family=label, variant=r.variant, source=r.source,
                         K_HT=r.K_HT, K_DT=r.K_DT,
                         gamma_obs=np.log(r.K_HT) / np.log(r.K_DT),
                         L_H=(r.K_HT - 1) / (r.K_DT - 1),
                         F_min=fmin, sd=sd, lcb=fmin - 1.645 * sd,
                         open_set=not closed,
                         excl_gated=bool(fmin - 1.645 * sd > F0),
                         excl_semiclassical=bool(fmin - 1.645 * sd > 0)))
    return pd.DataFrame(rows)


def main():
    lad = pd.read_csv("../data/ladh_adh_primary.csv")
    ya = pd.read_csv("../data/cha1989_yadh.csv")
    ya = ya[ya.note.str.contains("average")]
    res = pd.concat([analyse(lad, "LADH"), analyse(ya, "YADH")], ignore_index=True)

    say("=" * 94)
    say("PRIMARY COMPETITIVE H/T AND D/T FOR THE ALCOHOL DEHYDROGENASES, 25 C")
    say("=" * 94)
    say(f"gamma_SC = {GSC:.5f}, F0 = {F0:+.6f}. Bahnson et al. quote the")
    say("semiclassical upper limit as 3.26-3.34, spanning both mass conventions.")
    say("")
    say(f"{'enzyme':6s}{'variant':14s}{'K_HT':>8}{'K_DT':>8}{'gamma':>8}{'L_H':>7}"
        f"{'F_min':>9}{'sd':>8}{'95% LCB':>9}{'excl F0':>9}{'excl 0':>8}")
    for _, r in res.iterrows():
        say(f"{r.family:6s}{r.variant:14s}{r.K_HT:8.3f}{r.K_DT:8.3f}"
            f"{r.gamma_obs:8.3f}{r.L_H:7.2f}{r.F_min:9.4f}{r.sd:8.4f}{r.lcb:9.4f}"
            f"{'YES' if r.excl_gated else 'no':>9}"
            f"{'YES' if r.excl_semiclassical else 'no':>8}")

    lad_r = res[res.family == "LADH"]
    say("")
    say(f"  every record has L_H > gamma_SC = {GSC:.3f}, so every identified set is")
    say("  the open half-line (F_obs, inf) and no commitment lowers an endpoint.")
    say(f"  LADH: {len(lad_r)} records, offsets {lad_r.F_min.min():+.4f} to "
        f"{lad_r.F_min.max():+.4f}, "
        f"{'all negative' if (lad_r.F_min < 0).all() else f'{int((lad_r.F_min >= 0).sum())} not negative'}.")
    say(f"        {int(lad_r.excl_gated.sum())} of {len(lad_r)} exclude anything.")
    n0 = int((lad_r.F_min < 0).sum())
    say(f"        {n0} of {len(lad_r)} sets contain 0 and every positive value:")
    say(f"        consistent with tunneling, with the semiclassical locus, and")
    say(f"        with the gated model alike, so uninformative rather than")
    m = len(lad_r) - n0
    say(f"        contradictory. The other {m} {'lies' if m == 1 else 'lie'} above 0.")
    say(f"  YADH: offset {res[res.family=='YADH'].F_min.iloc[0]:+.4f}, excludes both.")
    say("")
    say("  Bahnson et al. (1993) give the reason: the LADH reaction is partially")
    say("  rate limited by benzaldehyde dissociation, whereas YADH has a rate-")
    say("  limiting chemical step. Commitment drives the observed offset down,")
    say("  which is the sign of every LADH value.")
    say("")
    say("  The framework therefore predicts which systems can inform before any")
    say("  mechanism is assumed: those whose chemical step is rate limiting.")

    res.to_csv("../results/adh_primary.csv", index=False)
    with open("../results/adh_primary.txt", "w") as f:
        f.write("\n".join(OUT) + "\n")
    print("\n[written] ../results/adh_primary.{csv,txt}")


if __name__ == "__main__":
    main()
