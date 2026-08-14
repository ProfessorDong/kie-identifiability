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



# ------------------------------------------------- how much bypass is tolerable
def bypass_tolerance(M, xD):
    """Largest isotope-blind bypass k_b/k_T leaving the map concave over x >= xD.

    A bypass competing alongside the isotope-sensitive step saturates the
    observable at sup_x K = M, which forces k2 = (M-1)(k_T+k_b); the threshold is
    then x* = sqrt(phi[(M-1)(1+phi)+phi]) with phi = k_b/k_T.  Concavity over the
    admissible range needs x* <= xD.
    """
    f = lambda phi: np.sqrt(phi * ((M - 1.0) * (1.0 + phi) + phi))
    lo, hi = 0.0, 10.0
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        lo, hi = (mid, hi) if f(mid) <= xD else (lo, mid)
    return lo


BENCH = [("YADH (Cha 1989)", 7.13, 1.73), ("BSAO (Grant 1989)", 35.2, 3.07),
         ("MAO-B pH 6.1, 10 C", 28.89, 2.896), ("LADH F93W", 7.755, 1.858),
         ("ecDHFR light, 25 C", 3.10, 1.75)]


def report_tolerance():
    print("\nIsotope-blind bypass tolerated before the map stops being concave:")
    for lab, M, xD in BENCH:
        print(f"  {lab:22s} K_HT={M:6.2f}  K_DT={xD:5.3f}  ->  "
              f"k_b/k_T up to {bypass_tolerance(M, xD)*100:.1f}%")


if __name__ == "__main__":
    run()
    report_tolerance()
