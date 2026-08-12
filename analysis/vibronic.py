"""Does the separation bound survive the full vibronic sum?

The main-text theorem is proved for the ground vibronic channel.  A referee will
ask whether excited channels destroy it.  This module answers that in two parts,
one exact and one numerical.

PART 1 (exact): how bounds behave under channel summation.
The total rate is a sum over channels, k_L = sum_j k_L^(j).  Writing
t_j = k_T^(j)/k_T for the tritium channel weights and
Y_j = k_D^(j)/k_T^(j), Z_j = k_H^(j)/k_T^(j), the summed effects are
    K_DT = E_t[Y],  K_HT = E_t[Z],
so the summed offset is  F = ln E[Z] - gamma ln E[Y].
If every channel obeys F^(j) = ln Z_j - gamma ln Y_j <= F*, then
Z_j <= e^{F*} Y_j^gamma, hence
    F <= F* + ( ln E[Y^gamma] - gamma ln E[Y] )  ==  F* + Delta.
By Jensen, since y -> y^gamma is CONVEX for gamma > 1, Delta >= 0, with equality
iff Y is constant across channels.  The inequality therefore runs the WRONG WAY:
summation can only push the offset UP.  A bound proved channel-by-channel does
not transfer to the sum for free; it transfers with a penalty Delta that is
exactly the dispersion of the channel D/T ratios.

Writing Y = E[Y] e^Z with E[e^Z] = 1, Delta = ln E[e^{gamma Z}], and for small
dispersion Delta ~ (1/2) gamma (gamma - 1) Var(Z).  With gamma = 3.349 this is
3.93 Var(Z), so the bound stays informative only while
    sd( ln (channel D/T ratio) )  <<  sqrt(|F0| / 3.93) = 0.103 .

PART 2 (numerical): are the individual excited channels even bounded by F0?
For harmonic wells the 0 -> n Franck-Condon factor is |S^{0n}|^2 =
(S^n/n!) e^{-S} with Huang-Rhys factor S_L = kappa sqrt(mu_L) R^2, which we
average over the same Gaussian gating distribution.  Excited channels also carry
isotope-dependent vibronic energy gaps n*hbar*omega_L with
omega_L ~ mu_L^{-1/2}, entering the Marcus factor.  Both are included below.
"""
from __future__ import annotations

import math

import numpy as np
from scipy import integrate

import masses as M

GSC = M.gamma_sc("C")
F0 = M.offset_F0("C")
MU = dict(zip("HDT", M.mu_ratios("C")))          # sqrt-mass ratios
KB_KCAL = 1.98720425e-3


# ---------------------------------------------------------------- overlaps
def gated_overlap(n, mu, A, w):
    """<|S^{0n}|^2> over Gaussian gating, in reduced groups (A, w).

    With kappa*sqrt(mu_H) = 1 by choice of units, R0^2 = A and sigma^2 = w/2,
    the Huang-Rhys factor is S = mu * R^2 and
        |S^{0n}(R)|^2 = (S^n / n!) exp(-S).
    """
    R0, sig = np.sqrt(A), np.sqrt(w / 2.0)
    lo, hi = R0 - 12 * sig, R0 + 12 * sig

    def f(R):
        S = mu * R ** 2
        return (np.exp(-(R - R0) ** 2 / (2 * sig ** 2)) / np.sqrt(2 * np.pi) / sig
                * np.exp(n * np.log(np.maximum(S, 1e-300)) - S)
                / math.factorial(n))

    val, _ = integrate.quad(f, lo, hi, limit=300)
    return max(val, 1e-300)


def channel_offset(n, A, w):
    """F for a single vibronic channel n, overlap only."""
    lt = np.log(gated_overlap(n, MU["T"], A, w))
    lh = np.log(gated_overlap(n, MU["H"], A, w)) - lt
    ld = np.log(gated_overlap(n, MU["D"], A, w)) - lt
    return lh - GSC * ld, lh, ld


# ------------------------------------------------------- full vibronic sum
def summed_offset(A, w, nmax=6, hw_H=3.0, lam=20.0, dG=-5.0, T=298.15,
                  include_gap=True):
    """Offset of the rate summed over channels 0..nmax.

    hw_H is the acceptor quantum in kcal/mol for protium; omega_L ~ mu_L^{-1/2}
    so hw_L = hw_H / sqrt(mu_L/mu_H).  Marcus factor
    exp[-(dG + lam + n hw_L)^2 / (4 lam kB T)].
    """
    tot = {}
    for L in "HDT":
        hw = hw_H / MU[L]                      # MU is sqrt(mu_L/mu_H)
        s = 0.0
        for n in range(nmax + 1):
            fc = gated_overlap(n, MU[L], A, w)
            gap = n * hw if include_gap else 0.0
            marcus = np.exp(-((dG + lam + gap) ** 2) / (4 * lam * KB_KCAL * T))
            s += fc * marcus
        tot[L] = s
    lh = np.log(tot["H"] / tot["T"])
    ld = np.log(tot["D"] / tot["T"])
    return lh - GSC * ld


def dispersion_penalty(A, w, nmax=6, hw_H=3.0, lam=20.0, dG=-5.0, T=298.15):
    """Delta = ln E[Y^gamma] - gamma ln E[Y] over the tritium channel weights."""
    wt, Y = [], []
    for n in range(nmax + 1):
        fcT = gated_overlap(n, MU["T"], A, w)
        fcD = gated_overlap(n, MU["D"], A, w)
        hwT, hwD = hw_H / MU["T"], hw_H / MU["D"]
        mT = np.exp(-((dG + lam + n * hwT) ** 2) / (4 * lam * KB_KCAL * T))
        mD = np.exp(-((dG + lam + n * hwD) ** 2) / (4 * lam * KB_KCAL * T))
        wt.append(fcT * mT)
        Y.append((fcD * mD) / (fcT * mT))
    wt = np.array(wt); wt /= wt.sum(); Y = np.array(Y)
    return np.log(np.sum(wt * Y ** GSC)) - GSC * np.log(np.sum(wt * Y)), wt, Y


def main():
    print("=" * 78)
    print("PART 1  individual vibronic channels, overlap only")
    print("=" * 78)
    print(f"ground-channel bound F0 = {F0:+.6f}")
    print(f"{'A':>7}{'w':>8}" + "".join(f"{'n=%d'%n:>11}" for n in range(4)))
    for A, w in ((1.0, 0.2), (5.0, 1.0), (20.0, 2.0)):
        row = f"{A:7.1f}{w:8.2f}"
        for n in range(4):
            F, _, _ = channel_offset(n, A, w)
            row += f"{F:11.5f}"
        print(row)
    print("\n  n = 0 respects the bound; excited channels do NOT.  Each added")
    print("  quantum shifts the offset upward by about 2|F0|, because the")
    print("  Franck-Condon progression contributes n ln(mu_L/mu_T) to ln K_LT.")

    print("\n" + "=" * 78)
    print("PART 2  the summed rate, with vibronic energy gaps included")
    print("=" * 78)
    print(f"{'A':>6}{'w':>7}{'hw_H':>7}{'F(sum)':>11}{'F(0 only)':>12}"
          f"{'Delta':>10}{'exceeds F0?':>13}")
    for A, w in ((1.0, 0.2), (5.0, 1.0), (20.0, 2.0)):
        for hw in (2.0, 3.0, 5.0):
            Fs = summed_offset(A, w, hw_H=hw)
            F0c, _, _ = channel_offset(0, A, w)
            D, wt, Y = dispersion_penalty(A, w, hw_H=hw)
            print(f"{A:6.1f}{w:7.2f}{hw:7.1f}{Fs:11.5f}{F0c:12.5f}"
                  f"{D:10.5f}{'YES' if Fs > F0 else 'no':>13}")

    print("\n" + "=" * 78)
    print("PART 3  the dispersion criterion")
    print("=" * 78)
    print(f"  Delta ~ (1/2) gamma (gamma-1) Var(ln Y),  gamma = {GSC:.4f}")
    print(f"  coefficient = {0.5*GSC*(GSC-1):.4f}")
    print(f"  bound stays informative while sd(ln Y) << "
          f"{np.sqrt(abs(F0)/(0.5*GSC*(GSC-1))):.4f}")
    for A, w, hw in ((5.0, 1.0, 3.0), (5.0, 1.0, 8.0), (5.0, 1.0, 15.0)):
        D, wt, Y = dispersion_penalty(A, w, hw_H=hw)
        lnY = np.log(Y)
        sd = np.sqrt(np.sum(wt * (lnY - np.sum(wt * lnY)) ** 2))
        print(f"  hw_H={hw:4.1f}: excited weight = {1-wt[0]:.4f}, "
              f"sd(ln Y) = {sd:.4f}, Delta = {D:.5f}, "
              f"approx = {0.5*GSC*(GSC-1)*sd**2:.5f}")


if __name__ == "__main__":
    main()
