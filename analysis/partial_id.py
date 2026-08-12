"""Partial identification of the mass-scaling offset, done exactly.

Replaces three earlier results that were wrong or incomplete:

  * the flat bound F <= F0, which fails for intrinsic effects below
    sqrt(mu_T/mu_H) = 1.2689 and is replaced here by a scale-dependent
    envelope B(K);
  * the numerical minimization of F(c), replaced by a closed form;
  * the single-commitment observation map, extended here to the two-commitment
    Northrop form with an equilibrium isotope effect, which is the form used in
    the enzymology literature.

Conventions follow masses.py: reduced masses of the X-H oscillator with CODATA
atomic masses, giving gamma_SC = 3.34887 for a carbon donor.  That is the value
obtained by Kohen and Jensen from the reduced mass of 12C and the hydrogen
isotopes; the alternative bare-mass value 3.2628 is also in use and every result
here is reported for a stated convention.
"""
from __future__ import annotations

import numpy as np

import masses as M

GSC = M.gamma_sc("C")
F0 = M.offset_F0("C")
MUH, MUD, MUT = M.mu_ratios("C")          # 1, d, t   (square-root mass ratios)


# ===================================================================== theory
def envelope_B(K):
    """Largest offset the ground-channel gated model can reach at fixed K_HT.

    The supremum over (A, w) at unconstrained scale is 0, attained as the
    isotope effects collapse to unity, so no flat negative bound exists.  At
    fixed intrinsic H/T magnitude the envelope is

        B(K) = ln K - (gamma/2) ln[ K^2 (t-1) / (t - d + K^2 (d-1)) ]   K < sqrt(t)
             = F0                                                       K >= sqrt(t)

    with t = mu_T/mu_H and d = mu_D/mu_H in square-root mass ratios.
    """
    K = np.asarray(K, dtype=float)
    t, d = MUT, MUD
    inner = K ** 2 * (t - 1.0) / (t - d + K ** 2 * (d - 1.0))
    below = np.log(K) - 0.5 * GSC * np.log(inner)
    return np.where(K < np.sqrt(t), below, F0)


def scale_threshold():
    """K_HT above which the envelope is flat at F0."""
    return np.sqrt(MUT)


# ------------------------------------------- one-commitment identified set
def F_min_exact(kht, kdt, gamma=GSC):
    """Closed-form lower endpoint of the identified set for the shared map.

    With a = K_HT - 1, b = K_DT - 1 and c the reciprocal commitment,

        dF/dc = [ c (gamma b - a) - a b (gamma - 1) ] / [ c (c - a) (c - b) ].

    If a >= gamma b the derivative is negative throughout and the infimum is
    F_obs, approached only as c -> infinity, so the identified set is the OPEN
    half-line (F_obs, inf).  Otherwise there is an interior stationary point at
    c* = a b (gamma - 1) / (gamma b - a) and the set is closed at F_min.

    Returns (F_min, c_star, interior, closed).
    """
    a, b = kht - 1.0, kdt - 1.0
    F_obs = np.log(kht) - gamma * np.log(kdt)
    if a >= gamma * b:                       # L_H >= gamma
        return F_obs, np.inf, False, False
    c_star = a * b * (gamma - 1.0) / (gamma * b - a)
    xh = kht * c_star / (1.0 + c_star - kht)
    xd = kdt * c_star / (1.0 + c_star - kdt)
    return np.log(xh) - gamma * np.log(xd), c_star, True, True


def F_min_vec(kht, kdt, gamma=GSC):
    """Vectorized F_min_exact.

    Both branches of the scalar routine are closed form, so the whole
    calculation is a pair of algebraic expressions selected by the L_H
    criterion.  Inadmissible pairs (those violating x_H > x_D > 1, which the
    inversion requires) return -inf rather than being dropped: a replicate that
    carries no information about the offset belongs in the lower tail of the
    sampling distribution, and discarding it would bias a lower confidence
    bound upward.

    On the interior branch the endpoint is evaluated through the closed forms

        x_H* = K_HT b (gamma-1) / (a-b),
        x_D* = K_DT a (gamma-1) / [gamma (a-b)],

    which follow from substituting c* into x(c) and require a > b, i.e. the
    normal ordering K_HT > K_DT.  That is exactly the admissibility mask.
    """
    kht = np.asarray(kht, dtype=float)
    kdt = np.asarray(kdt, dtype=float)
    a, b = kht - 1.0, kdt - 1.0
    ok = (kht > kdt) & (kdt > 1.0)
    with np.errstate(divide="ignore", invalid="ignore"):
        F_obs = np.log(kht) - gamma * np.log(kdt)
        xh = kht * b * (gamma - 1.0) / (a - b)
        xd = kdt * a * (gamma - 1.0) / (gamma * (a - b))
        F_int = np.log(xh) - gamma * np.log(xd)
    out = np.where(ok & (a < gamma * b), F_int, F_obs)
    return np.where(ok, out, -np.inf)


# ------------------------------------------- two-commitment identified set
def observed_two(x, cf, cr, eie):
    """Northrop form: K_obs = (x + Cf + Cr*EIE) / (1 + Cf + Cr)."""
    return (x + cf + cr * eie) / (1.0 + cf + cr)


def F_min_two_commitment(kht, kdt, gamma=GSC, eie_dt_max=1.4, n=240):
    """Lower endpoint of the identified set under the two-commitment map.

    Nuisance parameters are the forward commitment Cf >= 0, the reverse
    commitment Cr >= 0, and the D/T equilibrium isotope effect, with the H/T
    equilibrium effect tied by mass scaling.  Admissibility requires
    x_H > x_D > 1.
    """
    best = np.inf
    grid_c = np.concatenate([[0.0], np.logspace(-3, 3, n)])
    grid_e = np.linspace(1.0, eie_dt_max, 24)
    for cf in grid_c:
        for cr in grid_c[::4]:
            for edt in grid_e:
                eht = edt ** gamma
                S = 1.0 + cf + cr
                xh = kht * S - cf - cr * eht
                xd = kdt * S - cf - cr * edt
                if not (xh > xd > 1.0):
                    continue
                F = np.log(xh) - gamma * np.log(xd)
                if F < best:
                    best = F
    return best


# =============================================================== verification
def _verify():
    print("=" * 74)
    print("1. ENVELOPE B(K): closed form vs direct maximization over (A, w)")
    print("=" * 74)

    def sup_numeric(K, n=4000):
        L = np.log(K)
        best = -np.inf
        for w in np.logspace(-8, 8, n):
            pre = lambda mu: -0.5 * np.log1p(w * mu)
            acoef = -(MUH / (1 + w * MUH) - MUT / (1 + w * MUT))
            const = pre(MUH) - pre(MUT)
            if abs(acoef) < 1e-300:
                continue
            A = (L - const) / acoef
            if A <= 0:
                continue
            f = lambda mu: -0.5 * np.log1p(w * mu) - A * mu / (1.0 + w * mu)
            lt = f(MUT)
            best = max(best, (f(MUH) - lt) - GSC * (f(MUD) - lt))
        return best

    print(f"threshold sqrt(mu_T/mu_H) = {scale_threshold():.6f}")
    print(f"{'K_HT':>9}{'closed form':>14}{'numeric':>13}{'|diff|':>11}")
    worst = 0.0
    for K in (1.01, 1.05, 1.10, 1.20, 1.2689, 1.5, 3.0, 11.0, 100.0):
        cf, nm = float(envelope_B(K)), sup_numeric(K)
        worst = max(worst, abs(cf - nm))
        print(f"{K:9.4f}{cf:14.6f}{nm:13.6f}{abs(cf-nm):11.2e}")
    print(f"  worst discrepancy {worst:.2e}   "
          f"{'PASS' if worst < 1e-5 else 'FAIL'}")

    print()
    print("=" * 74)
    print("2. CLOSED-FORM F_min vs grid search, one-commitment map")
    print("=" * 74)
    rng = np.random.default_rng(3)
    worst, n_int = 0.0, 0
    for _ in range(600):
        kdt = 1 + 10 ** rng.uniform(-2.5, 0.9)
        kht = kdt * (1 + 10 ** rng.uniform(-2.5, 1.2))
        Fc, cs, interior, closed = F_min_exact(kht, kdt)
        n_int += interior
        c = (kht - 1) * (1 + np.logspace(-9, 9, 40000))
        x = lambda K: K * c / (1.0 + c - K)
        Fg = (np.log(x(kht)) - GSC * np.log(x(kdt))).min()
        worst = max(worst, abs(Fc - Fg))
    print(f"  600 random admissible pairs, {n_int} with an interior minimum")
    print(f"  worst |closed form - grid| = {worst:.2e}   "
          f"{'PASS' if worst < 1e-6 else 'FAIL'}")

    print()
    print("=" * 74)
    print("3. TWO-COMMITMENT MAP: does the reverse commitment lower the bound?")
    print("=" * 74)
    print(f"{'K_HT':>7}{'K_DT':>7}{'1-commit':>11}{'2-commit':>11}{'shift':>10}")
    for kht, kdt in ((5.04, 1.65), (2.362, 1.49), (6.44, 1.92), (1.93, 1.29)):
        one = F_min_exact(kht, kdt)[0]
        two = F_min_two_commitment(kht, kdt)
        print(f"{kht:7.2f}{kdt:7.2f}{one:11.4f}{two:11.4f}{two-one:10.4f}")
    print("  The 2-commitment column is a GRID EDGE, not a minimum. Under the")
    print("  additive parameterization x_L = K_LT(1+Cf+Cr) - Cf - Cr*EIE, both")
    print("  intrinsic effects diverge together as Cf grows, so")
    print("  F ~ (1-gamma) ln Cf -> -infinity and the set is unbounded below.")
    print("  That divergence is driven by Cf, NOT by the reverse commitment Cr:")
    print("  setting Cr = 0 reproduces it. It is the same parameterization")
    print("  dependence discussed in the Supplemental Material, and it is why")
    print("  the paper uses the shared-c map of Eq. (1), which follows from the")
    print("  kinetic scheme, rather than a shared additive commitment.")


if __name__ == "__main__":
    _verify()
