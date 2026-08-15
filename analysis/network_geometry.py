"""Kinetic-network topology determines the geometry of isotope identification.

For a steady-state mechanism with a single isotope-sensitive step, every
King-Altman term is a product of distinct rate constants, so the isotope-
sensitive constant k appears to at most first power in the numerator and in the
denominator of V/K:

    (V/K)(k) = (alpha k + beta) / (gamma k + delta),   alpha..delta >= 0.

Referencing to tritium, x = k/k_T, gives a Moebius observation map fixing x = 1,

    K(x) = (A x + B)/(C x + D),
    A = alpha k_T (delta + gamma k_T),   B = beta (delta + gamma k_T),
    C = gamma k_T (alpha k_T + beta),    D = delta (alpha k_T + beta),

and in log-rate coordinates h(t) = ln K(e^t) has

    h''(t) = -(AD - BC)(AC e^{2t} - BD) e^t / [(A e^t + B)^2 (C e^t + D)^2].

With AD > BC (masking monotone), the sign of h'' is the sign of BD - AC e^{2t}:

    B = 0            -> concave everywhere      -> identified set opens ABOVE
    C = 0            -> convex everywhere       -> identified set opens BELOW
    B, C > 0         -> convex below x*, concave above it,
                        x* = sqrt(BD/(AC)) = sqrt(beta delta/(alpha gamma)) / k_T.

The two coefficients have direct topological meaning.  beta is the part of the
V/K numerator not carrying k: it is non-zero exactly when some route to product
bypasses the isotope-sensitive step.  gamma is the coefficient of k in the
denominator: it is non-zero exactly when the isotope-sensitive step partitions
against another step, which is what a commitment is.

Hence: a committed step with no isotope-blind bypass gives an upper half-line;
an isotope-blind bypass with no competition gives a lower half-line; and with
both, the direction is decided by whether the intrinsic effect exceeds x*.
"""
from __future__ import annotations

import numpy as np
import sympy as sp

import masses as M

GSC = M.gamma_sc("C")


# ---------------------------------------------------------------- coefficients
def moebius_coeffs(vk_expr, k, k_T_val=None):
    """(alpha, beta, gamma, delta) for V/K written as (alpha k + beta)/(gamma k + delta)."""
    num, den = sp.fraction(sp.cancel(sp.together(vk_expr)))
    pn, pd = sp.Poly(sp.expand(num), k), sp.Poly(sp.expand(den), k)
    if pn.degree() > 1 or pd.degree() > 1:
        raise ValueError("V/K is not first order in the isotope-sensitive constant")
    return (sp.simplify(pn.coeff_monomial(k)), sp.simplify(pn.coeff_monomial(1)),
            sp.simplify(pd.coeff_monomial(k)), sp.simplify(pd.coeff_monomial(1)))


def classify(alpha, beta, gamma, delta, k_T):
    """Curvature class and the concave/convex threshold x*."""
    beta = float(beta); gamma = float(gamma)
    alpha = float(alpha); delta = float(delta); k_T = float(k_T)
    if beta == 0.0 and gamma > 0.0:
        return "concave", 0.0
    if gamma == 0.0 and beta > 0.0:
        return "convex", np.inf
    if beta == 0.0 and gamma == 0.0:
        return "linear", np.nan
    return "mixed", np.sqrt(beta * delta / (alpha * gamma)) / k_T


def h_of_t(t, alpha, beta, gamma, delta, k_T):
    """log-rate masking map for the given V/K coefficients."""
    x = np.exp(t)
    vk = lambda k: (alpha * k + beta) / (gamma * k + delta)
    return np.log(vk(x * k_T) / vk(k_T))


def curvature(t, *p, eps=1e-5):
    return (h_of_t(t + eps, *p) - 2 * h_of_t(t, *p) + h_of_t(t - eps, *p)) / eps**2


# ------------------------------------------------------------------ mechanisms
def mechanisms():
    """Concrete steady-state mechanisms, as (name, V/K expression, symbols)."""
    k1, k2, kc, kb, kr = sp.symbols("k1 k2 k_c k_b k_r", positive=True)
    return [
        ("committed step, no bypass", k1 * kc / (k2 + kc), kc,
         {k1: 10.0, k2: 5.0}),
        ("bypass + competition", k1 * (kc + kb) / (k2 + kc + kb), kc,
         {k1: 10.0, k2: 5.0, kb: 0.3}),
        ("strong bypass + competition", k1 * (kc + kb) / (k2 + kc + kb), kc,
         {k1: 10.0, k2: 5.0, kb: 8.0}),
        ("bypass, no competition", (kc + kb) / kr, kc,
         {kb: 0.5, kr: 2.0}),
    ]


def run(k_T=1.0):
    print("Curvature class from network topology, and the identified-set direction\n")
    print(f"{'mechanism':30s} {'beta':>8s} {'gamma':>7s} {'class':>8s} {'x*':>8s}  h'' check")
    for name, vk, k, subs in mechanisms():
        a, b, g, d = (float(z.subs(subs)) for z in moebius_coeffs(vk, k))
        cls, xstar = classify(a, b, g, d, k_T)
        # numerical curvature at a point on each side of x*
        probes = [1.5 * max(xstar, 1e-3), 0.5 * xstar] if cls == "mixed" else [2.0]
        got = []
        for xp in probes:
            if not np.isfinite(xp) or xp <= 0:
                continue
            got.append(np.sign(curvature(np.log(xp), a, b, g, d, k_T)))
        agree = {
            "concave": all(s < 0 for s in got),
            "convex": all(s > 0 for s in got),
            "mixed": len(got) == 2 and got[0] < 0 < got[1],
            "linear": True,
        }[cls]
        xs = "--" if cls in ("concave", "convex", "linear") else f"{xstar:.3f}"
        print(f"{name:30s} {b:8.3f} {g:7.3f} {cls:>8s} {xs:>8s}  {'OK' if agree else 'FAIL'}")

    print("\nConsequence for the offset (F_int = 0 on the reference ray):")
    for name, vk, k, subs in mechanisms():
        a, b, g, d = (float(z.subs(subs)) for z in moebius_coeffs(vk, k))
        tD = 0.6
        F = h_of_t(GSC * tD, a, b, g, d, k_T) - GSC * h_of_t(tD, a, b, g, d, k_T)
        cls, _ = classify(a, b, g, d, k_T)
        print(f"  {name:30s} class {cls:8s} F_obs = {F:+.5f}"
              f"   {'below' if F < 0 else 'above'} the ray")



# ----------------------------------------- reduction of the general map to series
def to_series(x, phi):
    """The contracted intrinsic effect Y = (x + phi)/(1 + phi).

    The general one-isotope-sensitive-step map factors exactly:

        K(x; phi, q) = K_series(Y; c),   Y = (x+phi)/(1+phi),  c = q/(1+phi),

    verified symbolically in verify_derivation.py.  A bypass therefore does not
    change the FORM of the masking, it changes what is being masked: the
    intrinsic effect is contracted affinely toward unity before the ordinary
    series map acts.  Since Y - 1 = (x-1)/(1+phi), the ratio (Y_H-1)/(Y_D-1)
    equals (x_H-1)/(x_D-1), so the L_H case split is bypass invariant.
    """
    return (np.asarray(x, float) + phi) / (1.0 + phi)


def from_series(Y, phi):
    """Inverse of to_series: the intrinsic effect behind a contracted one."""
    return np.asarray(Y, float) * (1.0 + phi) - phi


def endpoint(KH, KD, phi, r=1.0, n=200000):
    """Endpoint of the identified set for F at bypass fraction phi.

    The two tritium references differ at the non-transferred position, so their
    chemical rates differ by r = c_D/c_H (Proposition S4).  An isotope-blind
    bypass k_b and competing step k_2 are then referenced to DIFFERENT tritium
    rates in the two experiments, giving

        phi_D = r phi_H,    q_D = r q_H,

    so the bypass and the reference asymmetry cannot be treated separately.
    Profiling the commitment exactly and taking c -> infinity gives

        E_r(phi) = ln[K_HT + (K_HT-1)phi] - gamma ln[K_DT + (K_DT-1) r phi],

    which reduces to the shared-reference expression at r = 1.
    """
    cmin = max(KH, KD) - 1.0
    cs = cmin * (1.0 + np.geomspace(1e-10, 1e12, n))
    YH, YD = KH * cs / (1 + cs - KH), KD * cs / (1 + cs - KD)
    xh, xd = YH * (1 + phi) - phi, YD * (1 + r * phi) - r * phi
    ok = (xh > xd) & (xd > 1.0)
    F = np.where(ok, np.log(np.where(ok, xh, 2.0)) - GSC * np.log(np.where(ok, xd, 2.0)),
                 np.inf)
    return float(min(F.min(), endpoint_closed(KH, KD, phi, r)))


def endpoint_closed(KH, KD, phi, r=1.0):
    """Closed-form endpoint E_r(phi), exact where the infimum is at c -> oo."""
    return float(np.log(KH + (KH - 1.0) * phi)
                 - GSC * np.log(KD + (KD - 1.0) * r * phi))


def bypass_to_destroy(KH, KD, target, r=1.0):
    """Smallest bypass fraction phi_H dragging the endpoint down to `target`.

    Returns 0.0 when the endpoint already lies below the target, i.e. when the
    system excludes nothing and there is no inference for a bypass to overturn.
    """
    if endpoint(KH, KD, 0.0, r) < target:
        return 0.0
    lo, hi = 0.0, 80.0
    for _ in range(60):
        mid = 0.5 * (lo + hi)
        lo, hi = (mid, hi) if endpoint(KH, KD, mid, r) > target else (lo, mid)
    return 0.5 * (lo + hi)


# (label, K_HT, K_DT, r).  r is the reference asymmetry of Table S4 where a
# source reports secondary effects, and 1.0 where it is unmeasured; the ecDHFR
# light-enzyme entry is the 25 C record of Wang et al., an earlier version of
# this list carried (3.10, 1.75), which is not any record in the benchmark.
BENCH = [("YADH (Cha 1989)", 7.13, 1.73, 1.31),
         ("BSAO (Grant 1989)", 35.2, 3.07, 1.14),
         ("MAO-B pH 6.1, 10 C", 28.89, 2.896, 1.14),
         ("LADH F93W", 7.755, 1.858, 1.27),
         ("ecDHFR light, 25 C", 4.85, 1.66, 1.0)]


def report_tolerance():
    F0 = M.offset_F0("C")
    print("\nHow far an isotope-blind bypass moves the identified-set endpoint\n")
    print(f"{'system':22s} {'L_H':>6} {'endpoint':>9} {'phi to reach F0':>16}")
    for lab, KH, KD, r in BENCH:
        e0 = endpoint(KH, KD, 0.0, r)
        p = bypass_to_destroy(KH, KD, F0, r)
        tag = f"{100*p:14.1f}%" if p > 0 else "   (no exclusion)"
        print(f"{lab:22s} {(KH-1)/(KD-1):6.2f} {e0:+9.4f} {tag}   r={r}")
    print("\n  The endpoint falls monotonically with phi and has no threshold: a")
    print("  bypass does not have to reverse the curvature to matter.  An earlier")
    print("  version of this module asserted a safe tolerance from the curvature")
    print("  switch x* = sqrt(phi(q+phi)) evaluated at the SMALLEST admissible")
    print("  competing rate q = (M-1)(1+phi).  That is the infimum of x* over")
    print("  admissible q, not a bound above it, so it proved nothing.")


if __name__ == "__main__":
    run()
    report_tolerance()
