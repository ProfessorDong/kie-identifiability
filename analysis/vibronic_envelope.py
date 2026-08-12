"""The envelope of attainable offsets for the full vibronic sum.

The main text proves Theorem 1 for the ground channel and shows by one worked
parameter set that summation over excited channels can carry the offset above
zero.  What it does not do is map where.  This module computes

    B_vib(nu, lambda, dG, T) = sup over (A, w) of F_summed ,

restricted to intrinsic H/T effects above sqrt(mu_T/mu_H) = 1.2689, which is the
regime in which the ground-channel envelope is flat at F0 and therefore the
regime the empirical comparison uses.  Where B_vib < 0 a separation between the
gated and semiclassical families survives summation; where B_vib >= 0 it does
not, and no offset-based exclusion is possible even in principle.

EXACT GATED FRANCK-CONDON FACTORS.  The gating average of the 0 -> n factor is
done in closed form rather than by quadrature, which is what makes a phase
diagram affordable.  With S = mu R^2 and R ~ N(R0, sigma^2),

    <|S^{0n}|^2> = (mu^n / n!) <R^{2n} e^{-mu R^2}>

and completing the square in the exponent gives R ~ N(m, s^2) under the tilted
measure, with

    1/s^2 = 2 mu + 1/sigma^2,   m = R0 s^2/sigma^2,
    <R^{2n} e^{-mu R^2}> = (s/sigma) exp[-mu R0^2/(1 + 2 mu sigma^2)] E[R^{2n}],

where the moments of the tilted normal follow the recursion
E[X^k] = m E[X^{k-1}] + (k-1) s^2 E[X^{k-2}].  This is exact at every n and
costs no integration.  It is checked against the quadrature routine of
vibronic.py below.
"""
from __future__ import annotations

import math

import numpy as np

import masses as M
import vibronic as V

GSC = M.gamma_sc("C")
F0 = M.offset_F0("C")
MU = dict(zip("HDT", M.mu_ratios("C")))
KB_KCAL = 1.98720425e-3
KTHRESH = math.sqrt(M.mu_ratios("C")[2])          # 1.2689
CM_TO_KCAL = 2.8591e-3                            # kcal/mol per cm^-1


def fc_gated(nmax, mu, A, w):
    """Exact <|S^{0n}|^2> over Gaussian gating, n = 0..nmax, vectorized in n."""
    R0, sig2 = math.sqrt(A), w / 2.0
    if sig2 <= 0:
        S = mu * A
        return np.array([math.exp(n * math.log(max(S, 1e-300)) - S
                                  - math.lgamma(n + 1)) for n in range(nmax + 1)])
    s2 = 1.0 / (2.0 * mu + 1.0 / sig2)
    m = R0 * s2 / sig2
    pref = math.sqrt(s2 / sig2) * math.exp(-mu * R0 ** 2 / (1.0 + 2.0 * mu * sig2))
    mom = np.empty(2 * nmax + 1)                  # E[X^k] under N(m, s2)
    mom[0] = 1.0
    if len(mom) > 1:
        mom[1] = m
    for k in range(2, 2 * nmax + 1):
        mom[k] = m * mom[k - 1] + (k - 1) * s2 * mom[k - 2]
    out = np.empty(nmax + 1)
    for n in range(nmax + 1):
        out[n] = pref * math.exp(n * math.log(mu) - math.lgamma(n + 1)) * mom[2 * n]
    return np.maximum(out, 0.0)


def summed_X(A, w, nu_cm, lam, dG, T, nmax=30):
    """(ln K_HT, ln K_DT) for the channel-summed rate."""
    hw_H = nu_cm * CM_TO_KCAL
    tot = {}
    for L in "HDT":
        hw = hw_H / MU[L]                          # omega ~ mu^{-1/2}
        fc = fc_gated(nmax, MU[L], A, w)
        n = np.arange(nmax + 1)
        marcus = np.exp(-((dG + lam + n * hw) ** 2) / (4.0 * lam * KB_KCAL * T))
        tot[L] = float(np.sum(fc * marcus))
    if tot["T"] <= 0:
        return np.nan, np.nan
    return math.log(tot["H"] / tot["T"]), math.log(tot["D"] / tot["T"])


def B_vib(nu_cm, lam, dG, T, nmax=30, nA=26, nw=26, kmin=KTHRESH):
    """sup of F over (A, w), restricted to K_HT^int >= kmin.

    Returns (B, A*, w*, K*).  nan if no grid point meets the restriction.
    """
    best, arg = -np.inf, (np.nan, np.nan, np.nan)
    for A in np.logspace(-1.3, 1.8, nA):
        for w in np.logspace(-2.6, 1.4, nw):
            xh, xd = summed_X(A, w, nu_cm, lam, dG, T, nmax)
            if not (np.isfinite(xh) and np.isfinite(xd)):
                continue
            if math.exp(xh) < kmin:
                continue
            f = xh - GSC * xd
            if f > best:
                best, arg = f, (A, w, math.exp(xh))
    if not np.isfinite(best):
        return np.nan, np.nan, np.nan, np.nan
    return best, arg[0], arg[1], arg[2]


# =============================================================== verification
def _verify():
    print("=" * 78)
    print("1. EXACT GATED FRANCK-CONDON FACTORS vs QUADRATURE")
    print("=" * 78)
    worst = 0.0
    print(f"{'n':>3}{'mu':>8}{'A':>7}{'w':>7}{'exact':>15}{'quad':>15}{'rel':>11}")
    for A, w in ((2.0, 0.35), (5.0, 1.0), (0.5, 0.05)):
        for L in "HDT":
            ex = fc_gated(4, MU[L], A, w)
            for n in (0, 2, 4):
                q = V.gated_overlap(n, MU[L], A, w)
                rel = abs(ex[n] - q) / max(q, 1e-300)
                worst = max(worst, rel)
                if n == 2:
                    print(f"{n:>3}{MU[L]:8.4f}{A:7.2f}{w:7.2f}"
                          f"{ex[n]:15.6e}{q:15.6e}{rel:11.2e}")
    print(f"  worst relative difference {worst:.2e}   "
          f"{'PASS' if worst < 1e-6 else 'FAIL'}")

    print()
    print("=" * 78)
    print("2. REPRODUCE THE WORKED POINT OF THE MAIN TEXT")
    print("=" * 78)
    xh, xd = summed_X(2.011, 0.346, 3626.0, 11.55, -11.42, 284.8, nmax=30)
    print(f"  A=2.011 w=0.346 nu=3626 lam=11.55 dG=-11.42 T=284.8")
    print(f"  K_HT = {math.exp(xh):.5f}  K_DT = {math.exp(xd):.5f}  "
          f"F = {xh - GSC*xd:+.5f}")
    print(f"  main text reports K_HT = 1.51949, K_DT = 1.13128, F = +0.00530")

    print()
    print("=" * 78)
    print("3. GROUND-CHANNEL LIMIT RECOVERS F0")
    print("=" * 78)
    print("  Suppressing excited channels (nmax = 0) must return the")
    print(f"  ground-channel envelope, flat at F0 = {F0:+.6f} above K = {KTHRESH:.4f}:")
    b, A, w, K = B_vib(3000.0, 20.0, -5.0, 298.15, nmax=0)
    print(f"  B_vib(nmax=0) = {b:+.6f}   at A={A:.3g}, w={w:.3g}, K={K:.3f}")
    print(f"  difference from F0: {b - F0:+.2e}")


def _phase():
    print()
    print("=" * 78)
    print("4. PHASE DIAGRAM: WHERE DOES THE SEPARATION SURVIVE?")
    print("=" * 78)
    nus = (1500.0, 2000.0, 2500.0, 2800.0, 3000.0, 3300.0, 3600.0)
    lams = (5.0, 10.0, 20.0, 40.0)
    for dG, T in ((-5.0, 298.15), (-15.0, 298.15)):
        print(f"\n  dG = {dG} kcal/mol, T = {T} K"
              f"      (entries are B_vib; negative = separation survives)")
        print("  " + f"{'nu (cm^-1)':>12}" +
              "".join(f"{'lam='+str(l):>12}" for l in lams))
        for nu in nus:
            row = ""
            for lam in lams:
                b, *_ = B_vib(nu, lam, dG, T)
                row += f"{b:12.4f}" if np.isfinite(b) else f"{'--':>12}"
            print(f"  {nu:12.0f}" + row)
    print()
    print("  Reading: B_vib < F0 would leave the main-text bound intact;")
    print("  F0 <= B_vib < 0 leaves a separation but a smaller one; B_vib >= 0")
    print("  means the gated family reaches the semiclassical locus and no")
    print("  offset-based exclusion is possible at those parameters.")


if __name__ == "__main__":
    _verify()
    _phase()
