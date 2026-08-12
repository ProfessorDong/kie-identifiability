"""What an experiment must achieve to make the two sets disjoint.

The empirical result of the main text is null, and this module asks whether that
is a property of the systems measured or of the measurement itself.  It is the
latter: the discriminating signal is |F0| = 0.042 in log-KIE units, and the
precision of published competitive KIEs puts the sampling standard deviation of
F within a factor of two of that.  The experiment is underpowered by
construction, and the requirement can be stated exactly.

FOUR RESULTS.

1. Required precision.  Testing at one-sided level alpha with power 1-beta
   against a system on the semiclassical locus (F = 0) and with the commitment
   suppressed, the offset must be resolved to

       sd(F)  <=  |F0| / (z_{1-alpha} + z_{1-beta}) .

2. The shared reference is an asset, not a nuisance.  Because
   F = ln K_HT - gamma ln K_DT with both effects measured against the same
   tritium,

       Var(F) = s_H^2 + gamma^2 s_D^2 - 2 gamma rho s_H s_D ,

   which DECREASES in rho.  If the reference error is common mode, writing
   s_L^2 = s_ref^2 + s_L,ind^2 and Cov = s_ref^2 gives

       Var(F) = (gamma-1)^2 s_ref^2 + s_H,ind^2 + gamma^2 s_D,ind^2 ,

   so the shared component enters with weight (gamma-1)^2 = 5.52 instead of
   1 + gamma^2 = 12.21.  Measuring H/T and D/T in one triple-label mixture
   therefore buys a factor of 2.2 in that variance component over measuring them
   separately.  Sources do not report the correlation, so this is currently
   discarded.

3. The commitment costs more than the precision does.  The gap between the
   observed offset and the endpoint of the identified set is a deterministic
   function of the observed pair, and across the benchmark it is one to two
   orders of magnitude larger than |F0|.

4. An independent determination of the commitment closes the set, and the
   accuracy it needs is modest compared with the precision the KIEs need.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from scipy import stats

import masses as M
from partial_id import F_min_exact

GSC = M.gamma_sc("C")
F0 = M.offset_F0("C")
OUT = []


def say(s=""):
    print(s)
    OUT.append(s)


# ------------------------------------------------------------ 1. precision
def sd_required(alpha=0.05, power=0.80, signal=abs(F0)):
    """Largest sd(F) admitting the stated power against F_true = 0."""
    return signal / (stats.norm.ppf(1 - alpha) + stats.norm.ppf(power))


def var_F(s_H, s_D, rho, gamma=GSC):
    return s_H ** 2 + gamma ** 2 * s_D ** 2 - 2 * gamma * rho * s_H * s_D


def s_required(rho, alpha=0.05, power=0.80, gamma=GSC):
    """Equal relative precision on both KIEs meeting the requirement."""
    return sd_required(alpha, power) / np.sqrt(1 + gamma ** 2 - 2 * gamma * rho)


def main():
    d = pd.read_csv("../data/trinomial_benchmark.csv")
    s_H = (d.K_HT_se / d.K_HT).values
    s_D = (d.K_DT_se / d.K_DT).values

    say("=" * 84)
    say("WHAT THE EXPERIMENT MUST ACHIEVE")
    say("=" * 84)
    say(f"gamma_SC = {GSC:.5f},  discriminating signal |F0| = {abs(F0):.6f}")
    say("")

    say("-" * 84)
    say("1. REQUIRED PRECISION ON THE OFFSET")
    say("-" * 84)
    say(f"{'power':>8}{'sd(F) required':>18}{'vs benchmark median':>22}")
    med = np.median(np.sqrt(var_F(s_H, s_D, 0.0)))
    for p in (0.50, 0.80, 0.90, 0.95):
        r = sd_required(power=p)
        say(f"{p:8.2f}{r:18.5f}{med/r:20.2f}x")
    say(f"  benchmark median sd(F) at rho = 0 is {med:.5f}")
    say(f"  so reaching 80% power needs a "
        f"{med/sd_required(power=0.80):.2f}-fold reduction in sd,")
    say(f"  i.e. about {(med/sd_required(power=0.80))**2:.1f}x more replicates "
        f"at fixed per-assay precision.")

    say("")
    say("-" * 84)
    say("2. THE SHARED TRITIUM REFERENCE REDUCES THE VARIANCE")
    say("-" * 84)
    say(f"{'rho':>8}{'sqrt(1+g^2-2 g rho)':>22}{'required s (%)':>18}"
        f"{'gain vs rho=0':>16}")
    base = s_required(0.0)
    for rho in (0.0, 0.3, 0.5, 0.7, 0.9, 0.99):
        s = s_required(rho)
        say(f"{rho:8.2f}{np.sqrt(1+GSC**2-2*GSC*rho):22.4f}"
            f"{100*s:18.3f}{s/base:16.2f}x")
    say(f"  benchmark median relative precision: K_HT {100*np.median(s_H):.2f}%,"
        f" K_DT {100*np.median(s_D):.2f}%")
    say("")
    say("  Common-mode decomposition, s_L^2 = s_ref^2 + s_L,ind^2:")
    say(f"    shared reference enters Var(F) with weight (gamma-1)^2 = "
        f"{(GSC-1)**2:.3f}")
    say(f"    independent errors enter with weight 1 + gamma^2 = "
        f"{1+GSC**2:.3f}")
    say(f"    ratio {(1+GSC**2)/(GSC-1)**2:.2f}, the factor bought by measuring")
    say("    H/T and D/T in one triple-label mixture rather than separately.")
    say("")
    say("  At fixed s_D the variance is minimized at s_H = gamma rho s_D, where")
    say("  Var(F) = gamma^2 s_D^2 (1 - rho^2): the reference error cancels")
    say("  exactly in the limit rho -> 1. Reporting the covariance would let a")
    say("  reader exploit this; none of the sources do.")

    say("")
    say("-" * 84)
    say("3. THE COMMITMENT COSTS MORE THAN THE PRECISION")
    say("-" * 84)
    gaps = []
    for _, r in d.iterrows():
        f_obs = np.log(r.K_HT) - GSC * np.log(r.K_DT)
        f_min = F_min_exact(r.K_HT, r.K_DT)[0]
        gaps.append(f_obs - f_min)
    gaps = np.array(gaps)
    say(f"  F_obs - F_min over 73 records: median {np.median(gaps):.4f}, "
        f"max {gaps.max():.4f}")
    say(f"  in units of |F0|: median {np.median(gaps)/abs(F0):.1f}, "
        f"max {gaps.max()/abs(F0):.1f}")
    say(f"  {int((gaps==0).sum())} records have zero gap (L_H >= gamma, so the")
    say("  endpoint is the commitment-free value and nothing is lost).")
    say("  Where the gap is nonzero it dwarfs the signal, so suppressing the")
    say("  commitment matters more than improving precision.")

    say("")
    say("-" * 84)
    say("4. HOW WELL AN INDEPENDENT COMMITMENT MEASUREMENT MUST BE KNOWN")
    say("-" * 84)
    say("  Knowing c point identifies F, with dF/dc from Proposition 2.")
    say("  Requiring the induced sd(F) to stay below half the signal:")
    say(f"{'K_HT':>7}{'K_DT':>7}{'c assumed':>11}{'|dF/dc|':>11}"
        f"{'max sd(c)/c':>14}")
    for kht, kdt in ((5.04, 1.65), (6.44, 1.92), (2.36, 1.49)):
        a, b = kht - 1, kdt - 1
        for c in (2 * a, 10 * a):
            dFdc = (c * (GSC * b - a) - a * b * (GSC - 1)) / (c * (c - a) * (c - b))
            tol = (abs(F0) / 2) / abs(dFdc) / c if dFdc != 0 else np.inf
            say(f"{kht:7.2f}{kdt:7.2f}{c:11.2f}{abs(dFdc):11.4f}{100*tol:13.1f}%")
    say("  A commitment known to a few tens of percent already contributes less")
    say("  than the KIE precision does, so pre-steady-state or isotope-trapping")
    say("  determination of c is the cheaper route to identification.")

    say("")
    say("-" * 84)
    say("5. TEMPERATURES DO NOT SUBSTITUTE FOR PRECISION")
    say("-" * 84)
    say("  The series statistic is a maximum over temperatures and does not")
    say("  average down: adding temperatures raises the selection bias rather")
    say("  than shrinking the standard error. Pooling across temperatures needs")
    say("  F constant in T, which is part of what is under test. Under that")
    say("  assumption sd falls as 1/sqrt(n_T):")
    say(f"{'n_T':>6}{'sd(F) pooled':>16}{'power at rho=0':>18}")
    for n in (1, 3, 5, 10, 20):
        sd = med / np.sqrt(n)
        pw = stats.norm.cdf(abs(F0) / sd - stats.norm.ppf(0.95))
        say(f"{n:6d}{sd:16.5f}{pw:18.3f}")
    say("  Reaching 80% power by temperatures alone needs "
        f"{int(np.ceil((med/sd_required(power=0.80))**2))} of them at present "
        "precision.")

    with open("../results/design_power.txt", "w") as f:
        f.write("\n".join(OUT) + "\n")
    print("\n[written] ../results/design_power.txt")


if __name__ == "__main__":
    main()
