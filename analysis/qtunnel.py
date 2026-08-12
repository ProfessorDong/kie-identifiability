"""Gated vibronic tunneling: exact ground-channel rate ratios and their geometry.

Model
-----
Vibronically nonadiabatic hydrogen transfer in the Kuznetsov--Ulstrup /
Hammes-Schiffer--Soudackov form.  For isotopologue L in {H, D, T} the rate at a
fixed donor--acceptor distance R factorises into an isotope-independent
Marcus/electronic part and the squared overlap of the transferred-particle
vibrational wavefunctions,

    k_L(R)  =  K_el(T) * |S_L(R)|^2 ,

and the observed rate is the thermal average over a gating (promoting-mode)
coordinate,

    k_L  =  < |S_L(R)|^2 >_{P(R;T)} .

Ground-channel (0-0) overlap.  Two harmonic wells of force constant k_f whose
minima are separated by R, with a particle of mass m_L, give

    |S_L(R)|^2 = exp[ - sqrt(k_f m_L) R^2 / (2 hbar) ]
               = exp[ - kappa sqrt(m_L) R^2 ],      kappa = sqrt(k_f)/(2 hbar).

Note this is Gaussian in R, not exponential.  The familiar exp[-2 alpha_L (R-R0)]
form is its linearisation about R0, with alpha_L = kappa R0 sqrt(m_L); we do NOT
linearise, because the average below is available in closed form exactly.

Gating distribution.  A harmonic promoting mode of frequency Omega and effective
mass M, in thermal equilibrium, gives a Gaussian P(R;T) with mean R0 and the
exact quantum variance

    sigma^2(T) = (hbar / (2 M Omega)) coth( hbar Omega / (2 k_B T) ),

which carries hbar explicitly and reduces to k_B T/(M Omega^2) classically.

Exact averaged overlap.  With a = 1/(2 sigma^2) and b_L = kappa sqrt(m_L),

    < |S_L|^2 >  =  (1 + 2 b_L sigma^2)^{-1/2} exp[ - b_L R0^2 / (1 + 2 b_L sigma^2) ].

Reduced parameters.  Everything observable depends on two dimensionless groups,

    A     = kappa R0^2                     (rigidity / barrier-width group)
    w(T)  = 2 kappa sigma^2(T)             (gating group),

so with mu_L = sqrt(m_L / m_H),

    ln K_LT^int = -1/2 ln[(1 + w mu_L)/(1 + w mu_T)]
                  - A [ mu_L/(1 + w mu_L) - mu_T/(1 + w mu_T) ].

The temperature dependence of the intrinsic KIEs enters through w(T) alone:

    w(T) = w0 coth( theta / (2 T) ),       theta = hbar Omega / k_B.

Caveats that belong in the paper, not in a docstring footnote:
  * ground-channel only; excited vibronic channels open at high T and raise the
    effective exponent.  The predictions below are therefore a lower bound on
    gamma_TUN within this model class.
  * ZPE differences in the donor well are absorbed into the isotope-independent
    activation term; only the overlap carries isotope dependence here.
"""
from __future__ import annotations

import numpy as np

# masses of the transferred particle in units of the protium mass
M_H, M_D, M_T = 1.0, 2.0, 3.0
MU_H, MU_D, MU_T = np.sqrt(M_H), np.sqrt(M_D), np.sqrt(M_T)


def log_overlap_avg(A: float, w: float, mu: float) -> float:
    """ln < |S|^2 > for reduced mass-root mu, rigidity A, gating w. Exact."""
    return -0.5 * np.log1p(w * mu) - A * mu / (1.0 + w * mu)


def log_kie_intrinsic(A, w):
    """(ln K_HT^int, ln K_DT^int) under the gated tunneling model."""
    lt = log_overlap_avg(A, w, MU_T)
    return (log_overlap_avg(A, w, MU_H) - lt,
            log_overlap_avg(A, w, MU_D) - lt)


def gamma_tunneling(A, w):
    """Effective Swain-Schaad exponent ln K_HT^int / ln K_DT^int."""
    lh, ld = log_kie_intrinsic(A, w)
    return lh / ld


def gating_w(T_K, w0, theta_K):
    """w(T) = w0 coth(theta / 2T): exact harmonic-oscillator thermal variance."""
    return w0 / np.tanh(theta_K / (2.0 * np.asarray(T_K, dtype=float)))


# ---------------------------------------------------------------- limits
def gamma_rigid() -> float:
    """w -> 0.  gamma = (mu_H - mu_T)/(mu_D - mu_T), from the overlap exponent."""
    return (MU_H - MU_T) / (MU_D - MU_T)


def gamma_fully_gated() -> float:
    """w -> infinity.  Only the prefactor survives: ln(m_H/m_T)/ln(m_D/m_T)."""
    return np.log(M_H / M_T) / np.log(M_D / M_T)


# ------------------------------------------------- commitment (from v2)
def observed_from_intrinsic(x, c):
    """Competitive observation map K^obs = x(1+c)/(x+c); c = k_off/k_ref."""
    return x * (1.0 + c) / (x + c)


def predict_observed(T_K, A, w0, theta_K, c0, Ec_kcal, R_KCAL=1.98720425e-3):
    """Full forward model: gated tunneling + Arrhenius reciprocal commitment."""
    T_K = np.asarray(T_K, dtype=float)
    w = gating_w(T_K, w0, theta_K)
    lh, ld = log_kie_intrinsic(A, w)
    c = c0 * np.exp(-Ec_kcal / (R_KCAL * T_K))
    return (observed_from_intrinsic(np.exp(lh), c),
            observed_from_intrinsic(np.exp(ld), c))
