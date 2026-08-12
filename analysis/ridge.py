"""Corrected theory: the ridge limit, the offset invariant, the finite-scale bound.

Replaces the superseded "supremum = boundary" reading.

Geometry.  Write the intrinsic log-KIE pair as (X_D, X_H) = (ln K_DT, ln K_HT).

  semiclassical   X_H = gamma_SC X_D                     (a line through 0)
  gated tunneling X_H -> P_H + s C_H,  X_D -> P_D + s C_D  (the ridge, s = A/w^2)

with P_L = (1/2) ln(mu_T/mu_L) and C_L = 1/mu_L - 1/mu_T.  Because C_L is
exactly the zero-point-energy mass factor mu_L^(-1/2) - mu_T^(-1/2), we have
C_H/C_D = gamma_SC identically, so the ridge is a line PARALLEL to the
semiclassical locus, displaced by the fixed offset

  F0 = P_H - gamma_SC P_D.

Physical content.  Gating inverts the isotope mass power of the overlap from
sqrt(mu) (rigid) to 1/sqrt(mu) (strongly gated), and 1/sqrt(mu) is the
zero-point-energy scaling.  The exponents therefore coincide asymptotically for
a mechanistic reason.

Why the Swain-Schaad exponent fails.  gamma = X_H/X_D is the angle from the
origin.  Parallel lines share a direction, so the angle converges as the KIE
grows: the exponent is structurally blind to a fixed affine offset.  The
invariant that sees it is

  F = X_H - gamma_SC X_D,

which is 0 on the semiclassical locus and F0 on the ridge, at ANY KIE magnitude.
"""
from __future__ import annotations

import numpy as np

import masses as M


def ridge_coeffs(donor: str = "C"):
    PH, PD = M.ridge_P(donor)
    CH, CD = M.ridge_C(donor)
    return PH, PD, CH, CD


def ridge_logkies(s, donor: str = "C"):
    """Intrinsic (ln K_HT, ln K_DT) on the ridge at scale parameter s = A/w^2."""
    PH, PD, CH, CD = ridge_coeffs(donor)
    s = np.asarray(s, dtype=float)
    return PH + s * CH, PD + s * CD


def offset(lnKHT, lnKDT, donor: str = "C"):
    """The discriminating invariant F = ln K_HT - gamma_SC ln K_DT."""
    return np.asarray(lnKHT) - M.gamma_sc(donor) * np.asarray(lnKDT)


def gamma_infty(lnKHT, donor: str = "C"):
    """Largest exponent reachable at a given intrinsic H/T magnitude.

    Solving X_H = P_H + s C_H for s and substituting gives the finite-scale
    bound; it is strictly below gamma_SC and tends to it only as X_H -> inf.
    """
    PH, PD, CH, CD = ridge_coeffs(donor)
    X = np.asarray(lnKHT, dtype=float)
    s = (X - PH) / CH
    return np.where(s > 0, X / (PD + s * CD), np.nan)


def observed_from_intrinsic(x, c):
    """Competitive observation map, shared by both model classes."""
    return x * (1.0 + c) / (x + c)


def observed_offset_semiclassical(u, c, donor: str = "C"):
    """F measured on OBSERVED effects when the intrinsic pair is semiclassical.

    This is not zero.  Commitment masking drives the observed offset away from
    the intrinsic one, which is why an observed F cannot be read directly as
    evidence for or against either model class.
    """
    g = M.gamma_sc(donor)
    u = np.asarray(u, dtype=float)
    return ((1 - g) * np.log1p(c) - np.log(u ** g + c) + g * np.log(u + c))


if __name__ == "__main__":
    for donor in ("C", "N", "O"):
        PH, PD, CH, CD = ridge_coeffs(donor)
        print(f"donor {donor}:  gamma_SC={M.gamma_sc(donor):.5f}  "
              f"C_H/C_D={CH/CD:.5f}  F0={M.offset_F0(donor):+.6f}")
    print()
    print("finite-scale bound on the exponent:")
    for k in (2, 3, 5, 11, 100, 1e4):
        print(f"  K_HT^int={k:8g}  gamma_inf={gamma_infty(np.log(k)):.5f}")
    print(f"  gamma_SC = {M.gamma_sc():.5f} (reached only as K_HT^int -> inf)")
    print()
    print("commitment masking of the offset, semiclassical intrinsic pair:")
    print(f"{'u':>6}{'c=0.1':>10}{'c=1':>10}{'c=10':>10}{'c=1e3':>10}")
    for u in (1.5, 2.0, 3.0):
        row = "".join(f"{observed_offset_semiclassical(u, c):10.4f}"
                      for c in (0.1, 1, 10, 1e3))
        print(f"{u:6.1f}{row}")
    print(f"\nF0 for comparison: {M.offset_F0():+.6f}")
