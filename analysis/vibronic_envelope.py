"""The envelope of attainable offsets for the full vibronic sum.

The sampled (A,w) grids give lower bounds on a supremum for the chosen
channel model. A positive grid value is a constructive example. A negative
grid maximum cannot establish that all allowable offsets are negative.
We restrict intrinsic H/T to at least 1.2689, the flat part of the
ground-channel envelope used in the manuscript.

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
from scipy.special import gammaln, logsumexp

import masses as M
import vibronic as V

GSC = M.gamma_sc("C")
F0 = M.offset_F0("C")
MU = dict(zip("HDT", M.mu_ratios("C")))
KB_KCAL = 1.98720425e-3
KTHRESH = math.sqrt(M.mu_ratios("C")[2])          # 1.2689
CM_TO_KCAL = 2.859144e-3                          # kcal/mol per cm^-1 (hc N_A)


def log_fc_gated(nmax, mu, A, w):
    """Logarithms of the Gaussian-averaged channel factors, including at w=0.

    Axis zero enumerates quantum numbers. A and w can be scalars or arrays
    broadcastable to the shape of a parameter grid.
    """
    A, w = np.broadcast_arrays(np.asarray(A, dtype=float), np.asarray(w, dtype=float))
    if np.any(A < 0) or np.any(w < 0):
        raise ValueError("Action and width must be nonnegative.")
    denom = 1 + mu * w
    log_prefactor = -0.5 * np.log1p(mu * w) - mu * A / denom
    with np.errstate(divide="ignore"):
        log_mean = 0.5 * np.log(A) - np.log(denom)
        log_variance = np.log(w / 2) - np.log(denom)
    moments = [np.zeros_like(A), log_mean]
    for k in range(2, 2 * nmax + 1):
        moments.append(np.logaddexp(log_mean + moments[-1],
                                    np.log(k - 1) + log_variance + moments[-2]))
    return np.stack([log_prefactor + n * np.log(mu) - gammaln(n + 1) + moments[2*n]
                     for n in range(nmax + 1)])


def fc_gated(nmax, mu, A, w):
    """Gaussian-averaged factors on the ordinary scale, for quadrature comparisons."""
    return np.exp(log_fc_gated(nmax, mu, A, w))


def summed_X(A, w, nu_cm, lam, dG, T, nmax=30):
    """Return log H/T and D/T without ever forming a ratio of absolute rates."""
    log_rates = {}
    for isotope, mu in MU.items():
        factors = log_fc_gated(nmax, mu, A, w)
        channels = np.arange(nmax + 1).reshape((-1,) + (1,) * (factors.ndim - 1))
        gap = dG + lam + channels * nu_cm * CM_TO_KCAL / mu
        log_rates[isotope] = logsumexp(factors - gap**2 / (4 * lam * KB_KCAL * T), axis=0)
    return log_rates["H"] - log_rates["T"], log_rates["D"] - log_rates["T"]


def B_vib(nu_cm, lam, dG, T, nmax=30, nA=26, nw=26, kmin=KTHRESH):
    """Return the sampled maximum and its (action, width, intrinsic H/T)."""
    A, w = np.meshgrid(np.logspace(-1.3, 1.8, nA), np.logspace(-2.6, 1.4, nw),
                       indexing="ij")
    xh, xd = summed_X(A, w, nu_cm, lam, dG, T, nmax)
    offsets = np.where(xh >= np.log(kmin), xh - GSC * xd, -np.inf)
    i = np.unravel_index(np.argmax(offsets), offsets.shape)
    if not np.isfinite(offsets[i]):
        return np.nan, np.nan, np.nan, np.nan
    return offsets[i], A[i], w[i], np.exp(xh[i])


BVIB_CELLS = ((1500, 10, -5), (3000, 20, -5), (2500, 40, -5), (3600, 10, -15))


def bvib_grid_max(nu_cm, lam, dG, T=298.15, A_hi=63.1, w_hi=25.1, nA=60, nw=60):
    """Largest F on a log grid in (A, w), restricted to K_HT^int >= 1.2689.

    A from 0.050 and w from 10^-2.6 upward.  A grid maximum is a LOWER bound on
    the supremum and nothing more.  The defaults are the original 60 x 60 box;
    bvib_extended() widens it.  Returns (F, A*, w*).
    """
    A, w = np.meshgrid(np.logspace(math.log10(0.05), math.log10(A_hi), nA),
                       np.logspace(-2.6, math.log10(w_hi), nw), indexing="ij")
    xh, xd = summed_X(A, w, nu_cm, lam, dG, T, 30)
    offsets = np.where(xh >= np.log(KTHRESH), xh - GSC * xd, -np.inf)
    i = np.unravel_index(np.argmax(offsets), offsets.shape)
    if not np.isfinite(offsets[i]):
        return -np.inf, np.nan, np.nan
    return offsets[i], A[i], w[i]


def bvib_extended():
    """Supplementary Note 14 table: the original box against A <= 1e5, w <= 1e3.

    The extended values still move when the grid is refined or widened, so they
    are reported as grid maxima, lower bounds on a supremum that is not located.
    """
    return {cell: (bvib_grid_max(*cell),
                   bvib_grid_max(*cell, A_hi=1e5, w_hi=1e3, nA=160, nw=100))
            for cell in BVIB_CELLS}


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
    assert np.isclose(xh - GSC * xd, 0.005118213495, rtol=0, atol=2e-10)
    # Absolute rates at this grid point are below double-precision range.
    for cutoff in (30, 60):
        h, d = summed_X(1e5, 26.07123571336003, 1500, 10, -5, 298.15, cutoff)
        assert np.isclose(h - GSC * d, 0.5250265, rtol=0, atol=2e-6)

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
    print("4. SAMPLED MAXIMUM OFFSETS ON THE FINITE GRID")
    print("=" * 78)
    nus = (1500.0, 2000.0, 2500.0, 2800.0, 3000.0, 3300.0, 3600.0)
    lams = (5.0, 10.0, 20.0, 40.0)
    for dG, T in ((-5.0, 298.15), (-15.0, 298.15)):
        print(f"\n  dG = {dG} kcal/mol, T = {T} K"
              f"      (each entry is a finite-grid maximum)")
        print("  " + f"{'nu (cm^-1)':>12}" +
              "".join(f"{'lam='+str(l):>12}" for l in lams))
        for nu in nus:
            row = ""
            for lam in lams:
                b, *_ = B_vib(nu, lam, dG, T)
                row += f"{b:12.4f}" if np.isfinite(b) else f"{'--':>12}"
            print(f"  {nu:12.0f}" + row)
    print()
    print("  Positive values witness positive offsets. Negative grid maxima")
    print("  do not certify separation throughout the parameter domain.")


if __name__ == "__main__":
    _verify()
    _phase()
    print("\nSI Note 14. Sampled offsets on the small and extended grids:")
    print("frequency  reorganization  driving force  small max   extended max    action    width")
    for cell, (small, large) in bvib_extended().items():
        print(*cell, f"{small[0]:+.6f} {large[0]:+.6f} {large[1]:.0f} {large[2]:.6f}")
