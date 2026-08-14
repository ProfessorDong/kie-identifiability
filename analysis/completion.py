"""Completing the identified set, and what the completed offset then measures.

The isotope data alone determine a half-line.  One additional quantity closes
it.  Since the observation map is K = x(1+c)/(x+c), knowing the commitment c
gives the intrinsic effect in closed form,

    x = K c / (1 + c - K),                                             (1)

finite and above unity only for c > K - 1, which is the singular limit where the
observed effect saturates.  The offset F = ln x_H - gamma ln x_D then follows as
a point rather than a bound, with sensitivity

    d ln x / d ln c = (1 - K)/(1 + c - K),
    dF/d ln c = (1-K_H)/(1+c-K_H) - gamma (1-K_D)/(1+c-K_D).           (2)

Two routes to c were compared.  Varying solvent viscosity (Stokes-Einstein,
k_off ~ 1/eta) makes the pair of observations over-determine x and c, but the
inversion is ill conditioned: it needs the isotope effects to 0.006-0.3%, an
order of magnitude beyond current practice, and it is worst exactly where
masking is worst.  Measuring c directly, by isotope trapping or the standard
commitment methods, needs c only to a few tens of percent once c exceeds about
twice the singular value, and is therefore the practical route.

What the completed offset buys is a probe of the gating coordinate.  In the
gated model the offset depends on the dimensionless gating width
w = 2 kappa sqrt(mu_H) sigma^2(T), and for a harmonic promoting mode the exact
thermal variance is w(T) = w0 coth(theta/2T) with theta = hbar omega_g / k_B.
Zero-point motion makes w saturate below theta, so F(T) carries the mode
frequency.  The sensitivity is largest where hbar omega_g ~ 2 k_B T, i.e. near
400 cm^-1 at physiological temperature, which is the promoting-mode range.
"""
from __future__ import annotations

import numpy as np

import masses as M
from qtunnel import gating_w, log_kie_intrinsic

GSC = M.gamma_sc("C")
F0 = M.offset_F0("C")
CM_TO_K = 1.4387769          # hbar omega / k_B, per cm^-1


# --------------------------------------------------- closing the identified set
def intrinsic_from_commitment(K, c):
    """Intrinsic effect from an observed effect and a known commitment, Eq. (1)."""
    K = np.asarray(K, float); c = np.asarray(c, float)
    return np.where(c > K - 1.0, K * c / (1.0 + c - K), np.nan)


def dF_dlnc(KH, KD, c, gamma=GSC):
    """Sensitivity of the completed offset to the commitment, Eq. (2)."""
    return (1.0 - KH)/(1.0 + c - KH) - gamma*(1.0 - KD)/(1.0 + c - KD)


def commitment_precision_needed(KH, KD, c, target=abs(F0)/2, gamma=GSC):
    """Relative precision on c that fixes F to +/- `target`."""
    d = abs(dF_dlnc(KH, KD, c, gamma))
    return np.inf if d == 0 else target/d


# ------------------------------------------------------- the two-viscosity route
def intrinsic_from_two_viscosities(K1, K2, r):
    """Intrinsic effect from observations at viscosity ratio r (c2 = c1/r)."""
    return (r*K1*K2 - K1*K2 - r*K1 + K2)/(r*K2 - r - K1 + 1.0)


# ------------------------------------------------- what the completed offset sees
def offset_gated(A, w):
    lnKH, lnKD = log_kie_intrinsic(A, w)
    return lnKH - GSC*lnKD


def offset_vs_T(A, w0, nu_cm, temps):
    """F(T) for a harmonic promoting mode of wavenumber nu_cm."""
    theta = nu_cm*CM_TO_K
    return np.array([offset_gated(A, gating_w(T, w0, theta)) for T in temps])


def sensitivity_peak(A, w0, T_lo=250.0, T_hi=350.0, nus=None):
    """Wavenumber at which F(T) responds most strongly, and the span in |F0|."""
    nus = np.arange(30.0, 1200.0, 10.0) if nus is None else np.asarray(nus, float)
    span = np.array([abs(offset_vs_T(A, w0, nu, [T_hi])[0]
                         - offset_vs_T(A, w0, nu, [T_lo])[0]) for nu in nus])
    i = int(np.argmax(span))
    return nus[i], span[i]/abs(F0)


def run():
    print("Closing the half-line with one commitment measurement\n")
    print(f"{'system':22s} {'K_HT':>7s} {'K_DT':>7s} {'c/(K_H-1)':>10s} "
          f"{'F_int':>9s} {'gain vs bound':>14s} {'c needed':>9s}")
    for lab, KH, KD in [("YADH", 7.13, 1.73), ("BSAO", 35.2, 3.07),
                        ("LADH F93W", 7.755, 1.858)]:
        Fobs = np.log(KH) - GSC*np.log(KD)
        for mult in (2.0, 3.0, 5.0):
            c = (KH - 1.0)*mult
            xh = intrinsic_from_commitment(KH, c); xd = intrinsic_from_commitment(KD, c)
            Fi = np.log(xh) - GSC*np.log(xd)
            need = commitment_precision_needed(KH, KD, c)
            print(f"{lab:22s} {KH:7.2f} {KD:7.3f} {mult:10.1f} {Fi:+9.4f} "
                  f"{(Fi-Fobs)/abs(F0):11.1f}|F0| {need*100:8.1f}%")

    print("\nThe viscosity route, for comparison: precision needed on each effect")
    for c in (1.0, 3.0, 10.0):
        xH, xD = 2.0**GSC, 2.0
        K = lambda x, cc: x*(1 + cc)/(x + cc)
        r = 4.0
        # finite-difference conditioning of the closed-form inversion
        eps = 1e-6
        base = intrinsic_from_two_viscosities(K(xH, c), K(xH, c/r), r)
        pert = intrinsic_from_two_viscosities(K(xH, c)*(1 + eps), K(xH, c/r), r)
        amp = abs((pert - base)/base)/eps
        print(f"   c={c:5.1f}: error amplification {amp:8.1f}x"
              f"  -> {abs(F0)/2/amp*100:.3f}% needed")

    print("\nWhat the completed offset then measures: the promoting mode")
    print(f"{'A':>6s} {'w0':>6s} {'peak nu_g':>11s} {'span over 250-350 K':>21s}")
    for A in (1.0, 5.0, 20.0):
        for w0 in (0.1, 0.5):
            nu, sp = sensitivity_peak(A, w0)
            print(f"{A:6.1f} {w0:6.2f} {nu:9.0f} cm-1 {sp:16.2f}|F0|")
    print(f"\n   quantum-classical crossover at 300 K: nu = {2*300/CM_TO_K:.0f} cm-1")


if __name__ == "__main__":
    run()
