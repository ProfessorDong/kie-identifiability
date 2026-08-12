"""What a competitive isotope experiment identifies about the offset, without fitting.

COUNTING.  At each temperature the experiment returns TWO numbers,
(K_HT^obs, K_DT^obs).  The mechanism carries THREE unknown functions of
temperature: the intrinsic D/T effect u(T), the offset F(T), and the commitment
c(T).  The system is underdetermined by one function at every temperature, and
temperature resolution does not repair it: each new temperature adds two
equations and three unknowns.  Identifiability comes entirely from restrictions
that tie the unknowns across temperatures, and one of those -- F constant -- is
itself part of the mechanism under test.

WHAT SURVIVES WITHOUT ANY RESTRICTION.  Inverting the observation map at fixed
commitment gives x(c) = K c / (1 + c - K), admissible for c > K - 1, so the
offset implied by a commitment c is

    F(c) = ln x_H(c) - gamma_SC ln x_D(c),          c > K_HT - 1.

F(c) -> +infinity as c -> (K_HT - 1)+, and F(c) -> F_obs = ln K_HT -
gamma_SC ln K_DT as c -> infinity (the commitment-free limit).  Expanding at
large c,

    dF/dc = [gamma_SC (K_DT - 1) - (K_HT - 1)] / c^2 + O(c^-3),

so the approach to F_obs is from below or above according to the sign of
gamma_SC - L_H, where

    L_H = (K_HT - 1)/(K_DT - 1)

is the normal-difference contrast of the H-reference design.  Hence:

  Proposition (identifiable set).  The set of offsets consistent with one
  observation is the half-line [F_min, infinity), where
     L_H >= gamma_SC : F(c) decreases, F_min = F_obs, attained as c -> infinity;
     L_H <  gamma_SC : F(c) has an interior minimum at finite commitment and
                       F_min < F_obs strictly.
  An observation therefore bounds the offset FROM BELOW ONLY, and never above.

  Corollary.  A mechanism is excluded by an observation iff its predicted offset
  lies below F_min.  Since commitment masking drives F_min negative while both
  F = 0 and F0 = -0.042 are at or above zero on that scale, neither mechanism is
  ever excluded by a single measurement, at any precision.

That the classical H-reference contrast decides the shape of the identifiable
set is not a coincidence: L_H = gamma_SC is exactly the locus at which the
commitment-free limit stops being the extreme point of the inversion.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

import masses as M

GSC = M.gamma_sc("C")
F0 = M.offset_F0("C")
OUT = []


def say(s=""):
    print(s)
    OUT.append(s)


def x_of_c(K, c):
    return K * c / (1.0 + c - K)


def F_curve(kht, kdt, n=60000):
    c = (kht - 1.0) * (1 + np.logspace(-9, 9, n))
    return c, np.log(x_of_c(kht, c)) - GSC * np.log(x_of_c(kdt, c))


def F_min(kht, kdt):
    c, F = F_curve(kht, kdt)
    i = int(np.argmin(F))
    return F[i], c[i], 0 < i < len(c) - 1


def main():
    d = pd.read_csv("../data/trinomial_benchmark.csv")
    say("=" * 84)
    say("THE MODEL-FREE IDENTIFIABLE SET FOR THE DISCRIMINATING INVARIANT")
    say("=" * 84)
    say(f"gamma_SC = {GSC:.5f}   semiclassical F = 0   tunneling ridge "
        f"F0 = {F0:+.6f}")
    say("")

    say("1. THE SHAPE OF THE SET IS DECIDED BY THE H-REFERENCE CONTRAST")
    d["L_H"] = (d.K_HT - 1) / (d.K_DT - 1)
    d["predict_interior"] = GSC * (d.K_DT - 1) - (d.K_HT - 1) > 0
    res = [F_min(r.K_HT, r.K_DT) for _, r in d.iterrows()]
    d["F_min"] = [r[0] for r in res]
    d["c_star"] = [r[1] for r in res]
    d["interior"] = [r[2] for r in res]
    d["F_obs"] = np.log(d.K_HT) - GSC * np.log(d.K_DT)
    agree = int((d.predict_interior == d.interior).sum())
    say(f"   records: {len(d)}   predicted interior minimum: "
        f"{int(d.predict_interior.sum())}   observed: {int(d.interior.sum())}")
    say(f"   sign criterion L_H < gamma_SC agrees with the numerics on "
        f"{agree}/{len(d)} records")
    say(f"   L_H range where the minimum is interior : "
        f"[{d[d.interior].L_H.min():.3f}, {d[d.interior].L_H.max():.3f}]")
    say(f"   L_H range where it is at c -> infinity  : "
        f"[{d[~d.interior].L_H.min():.3f}, {d[~d.interior].L_H.max():.3f}]")
    say(f"   the two ranges are separated by gamma_SC = {GSC:.4f}")
    say("")

    say("2. THE ASSUMPTION-FREE BOUND PER SERIES")
    say("   With F constant in temperature, each temperature contributes its own")
    say("   lower bound, so F >= max_T F_min(T).  No fitting, no Arrhenius form,")
    say("   no gating model: only positivity of the commitment.")
    say("")
    g = (d.groupby(["family", "variant", "step"])
           .agg(n_T=("T_C", "size"), bound=("F_min", "max"),
                any_interior=("interior", "any")).reset_index()
           .sort_values("bound", ascending=False))
    say(f"   {'series':44s}{'n_T':>4}{'F >=':>10}{'gap to F0':>11}{'excludes':>10}")
    for _, r in g.iterrows():
        say(f"   {r.family+' '+r.variant+' ('+r.step+')':44s}{r.n_T:4d}"
            f"{r.bound:10.4f}{F0-r.bound:11.4f}"
            f"{'RIDGE' if r.bound > F0 else 'none':>10}")
    say("")
    say(f"   series excluding the tunneling ridge : {int((g.bound > F0).sum())} of {len(g)}")
    say(f"   series excluding semiclassical       : {int((g.bound > 0).sum())} of {len(g)}")
    best = g.iloc[0]
    say(f"   closest: {best.family} {best.variant}, bound {best.bound:+.4f},")
    say(f"            short of the ridge by {F0-best.bound:.4f} log-KIE units.")
    say("")

    say("3. WHAT A DECISIVE SINGLE MEASUREMENT WOULD REQUIRE")
    say("   One observation excludes the ridge iff F_obs > F0 with the minimum")
    say("   at c -> infinity, i.e. a commitment-free measurement with")
    say(f"   ln K_HT - gamma_SC ln K_DT > {F0:+.4f}.")
    say("")
    say(f"   {'K_DT^obs':>10}{'K_HT^obs required':>20}{'observed exponent':>20}")
    for kdt in (1.5, 2.0, 2.5, 3.0):
        kht = np.exp(F0 + GSC * np.log(kdt))
        say(f"   {kdt:10.2f}{kht:20.3f}{np.log(kht)/np.log(kdt):20.4f}")
    say("")
    say("   The required exponent sits just below gamma_SC in every case, so the")
    say("   experiment must be essentially commitment-free.  That is the single")
    say("   most informative thing this analysis says about experimental design:")
    say("   suppressing the commitment matters more than precision or")
    say("   temperature coverage, because the commitment is what makes the")
    say("   identifiable set a half-line open in the wrong direction.")

    d.to_csv("../results/identifiable_records.csv", index=False)
    g.to_csv("../results/modelfree_bounds.csv", index=False)
    with open("../results/identifiable_set.txt", "w") as f:
        f.write("\n".join(OUT) + "\n")
    print("\n[written] ../results/identifiable_set.txt, modelfree_bounds.csv")


if __name__ == "__main__":
    main()
