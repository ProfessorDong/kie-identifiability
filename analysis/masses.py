"""The isotope mass convention, defined once and imported everywhere.

Rationale.  The semiclassical Swain-Schaad exponent comes from the zero-point
energy of the transferred-hydrogen stretch, hbar*omega/2 with
omega = sqrt(k_f/mu).  The relevant mass is therefore the REDUCED mass of the
X-H oscillator, not the bare isotope mass.  The same reduced mass sets the
width of the transferred-particle wavefunction in the tunneling overlap,
sqrt(hbar/(2 mu omega)), because mu*omega = sqrt(k_f*mu).

That matters here: the supremum identity sup gamma_TUN = gamma_SC is an
identity in whatever masses are used, so it is exact only when BOTH models are
evaluated with the same masses.  Mixing a bare-mass tunneling model with a
reduced-mass semiclassical exponent, as an earlier draft did, manufactures a
spurious gap between them.

Atomic masses are CODATA/AME2020 values in unified atomic mass units.
"""
from __future__ import annotations

import numpy as np

M_H_ATOMIC = 1.00782503207
M_D_ATOMIC = 2.01410177785
M_T_ATOMIC = 3.01604928199

# donor heavy-atom masses, for the X-H oscillator reduced mass
HEAVY = {"C": 12.0, "N": 14.0030740048, "O": 15.9949146196, "S": 31.97207100}


def reduced(m_iso: float, donor: str = "C") -> float:
    """Reduced mass of the X-H oscillator for donor heavy atom X."""
    M = HEAVY[donor]
    return m_iso * M / (m_iso + M)


def triple(donor: str = "C"):
    """(mu_H, mu_D, mu_T) reduced masses for the given donor atom."""
    return tuple(reduced(m, donor) for m in
                 (M_H_ATOMIC, M_D_ATOMIC, M_T_ATOMIC))


def mu_ratios(donor: str = "C"):
    """Square-root mass ratios (1, sqrt(mu_D/mu_H), sqrt(mu_T/mu_H))."""
    a, b, c = triple(donor)
    return 1.0, np.sqrt(b / a), np.sqrt(c / a)


def gamma_sc(donor: str = "C") -> float:
    """Semiclassical Swain-Schaad exponent from zero-point energy."""
    a, b, c = triple(donor)
    return (a ** -0.5 - c ** -0.5) / (b ** -0.5 - c ** -0.5)


def gamma_rigid(donor: str = "C") -> float:
    """Rigid-barrier tunneling exponent: sqrt(m) differences."""
    muH, muD, muT = mu_ratios(donor)
    return (muH - muT) / (muD - muT)


def gamma_gated(donor: str = "C") -> float:
    """Strong-gating limit at finite barrier width: log-mass ratio."""
    a, b, c = triple(donor)
    return np.log(a / c) / np.log(b / c)


def ridge_P(donor: str = "C"):
    """Ridge intercepts P_L = (1/2) ln(mu_T/mu_L)."""
    muH, muD, muT = mu_ratios(donor)
    return 0.5 * np.log(muT / muH), 0.5 * np.log(muT / muD)


def ridge_C(donor: str = "C"):
    """Ridge slopes C_L = 1/mu_L - 1/mu_T, the zero-point-energy mass factor."""
    muH, muD, muT = mu_ratios(donor)
    return 1 / muH - 1 / muT, 1 / muD - 1 / muT


def offset_F0(donor: str = "C") -> float:
    """Mass-scaling offset of the tunneling ridge from the semiclassical locus.

    F = ln K_HT - gamma_SC ln K_DT equals 0 on the semiclassical locus and
    F0 on the gated-tunneling ridge, independently of the KIE magnitude.
    """
    PH, PD = ridge_P(donor)
    return PH - gamma_sc(donor) * PD


if __name__ == "__main__":
    print(f"{'donor':>6}{'gamma_SC':>11}{'rigid':>10}{'gated':>10}{'F0':>11}")
    for d in HEAVY:
        print(f"{d:>6}{gamma_sc(d):11.5f}{gamma_rigid(d):10.5f}"
              f"{gamma_gated(d):10.5f}{offset_F0(d):11.6f}")
    CH, CD = ridge_C("C")
    print(f"\nC-H donor: C_H/C_D = {CH/CD:.10f}  gamma_SC = {gamma_sc('C'):.10f}")
    print("(equal by construction: C_L is the zero-point-energy mass factor)")
