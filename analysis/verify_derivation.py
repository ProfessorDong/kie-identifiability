"""Numerical verification of every closed form in qtunnel.py.

Nothing in the manuscript may rest on an algebraic step that is not checked here.
"""
from __future__ import annotations

import numpy as np
from scipy import integrate

import qtunnel as Q

FAIL = 0


def check(name, ok, detail=""):
    global FAIL
    ok = bool(np.all(ok))
    print(f"  {'ok  ' if ok else 'FAIL'}  {name:56s} {detail}")
    if not ok:
        FAIL += 1


def close(a, b, tol=1e-8):
    return np.all(np.abs(np.asarray(a) - np.asarray(b)) < tol)


print("1. exact Gaussian average vs numerical quadrature")
for A, w, mu in [(1.0, 0.3, Q.MU_H), (3.5, 0.05, Q.MU_T), (0.7, 2.0, Q.MU_D),
                 (12.0, 0.01, Q.MU_D), (2.0, 8.0, Q.MU_H), (40.0, 0.02, Q.MU_H)]:
    kappa = 1.0
    R0 = np.sqrt(A / kappa)
    sig2 = w / (2.0 * kappa)
    sig = np.sqrt(sig2)
    b = kappa * mu
    # the integrand is a product of two Gaussians; its peak sits at
    # R* = R0/(1+2 b sigma^2).  Integrate on a wide window around both peaks,
    # and tell quad where the structure is, or it misses narrow ones entirely.
    Rstar = R0 / (1.0 + 2 * b * sig2)
    lo = min(Rstar, R0) - 40 * sig - 1.0
    hi = max(Rstar, R0) + 40 * sig + 1.0

    # Factor out exp(-b R0^2) before integrating: the raw integrand underflows
    # (values ~1e-17) and quad then reports garbage in the last digits.
    def integrand(R):
        return (np.exp(-(R - R0) ** 2 / (2 * sig2)) / np.sqrt(2 * np.pi * sig2)
                * np.exp(-b * (R ** 2 - R0 ** 2)))

    num, _ = integrate.quad(integrand, lo, hi, points=[Rstar, R0], limit=400)
    got = np.log(num) - b * R0 ** 2
    want = Q.log_overlap_avg(A, w, mu)
    check(f"A={A} w={w} mu={mu:.4f}", close(got, want, 1e-7),
          f"{got:.10f} vs {want:.10f}")

print("\n2. gamma_TUN limits in w at fixed A")
print(f"     rigid  gamma_TUN(w->0)   = {Q.gamma_rigid():.6f}")
print(f"     gated  gamma_TUN(w->inf) = {Q.gamma_fully_gated():.6f}   (A fixed, finite)")
check("w->0 limit", close(Q.gamma_tunneling(2.0, 1e-9), Q.gamma_rigid(), 1e-6))
check("w->inf limit", close(Q.gamma_tunneling(2.0, 1e9), Q.gamma_fully_gated(), 1e-5))

print("\n3. the large-A limit is a Mobius factor times the rigid exponent")
# ln K ~ -A[f(mu_L)-f(mu_T)], f(mu)=mu/(1+w mu), and
# f(mu_H)-f(mu_T) = (mu_H-mu_T)/((1+w mu_H)(1+w mu_T)), so
# gamma(A->inf, w) = gamma_rigid * (1+w mu_D)/(1+w mu_H).
for w in [0.0, 0.1, 1.0, 10.0, 1e3]:
    pred = Q.gamma_rigid() * (1 + w * Q.MU_D) / (1 + w * Q.MU_H)
    # the asymptotic needs A >> w, since the A-term scales as A/w^2
    got = Q.gamma_tunneling(1e6 * max(1.0, w) ** 2, w)
    check(f"w={w:<8g} large-A formula", close(got, pred, 1e-4),
          f"{got:.6f} vs {pred:.6f}")

print("\n4. supremum of gamma_TUN over the whole parameter space")
sup = Q.gamma_rigid() * Q.MU_D / Q.MU_H          # w->inf of the large-A formula
gamma_sc_atomic = ((Q.M_H ** -0.5 - Q.M_T ** -0.5) /
                   (Q.M_D ** -0.5 - Q.M_T ** -0.5))
print(f"     sup gamma_TUN                    = {sup:.10f}")
print(f"     semiclassical (atomic mass) gamma = {gamma_sc_atomic:.10f}")
check("sup gamma_TUN EQUALS the semiclassical atomic-mass exponent",
      close(sup, gamma_sc_atomic, 1e-12), f"difference {sup-gamma_sc_atomic:.2e}")

print("\n4b. is that identity special to masses 1:2:3, or general? (symbolic)")
import sympy as sp

a, b_, c = sp.symbols("a b c", positive=True)          # sqrt-masses
g_sc = (1 / a - 1 / c) / (1 / b_ - 1 / c)              # ZPE / semiclassical
g_sup = ((a - c) / (b_ - c)) * (b_ / a)                # sup of gated tunneling
check("sup(gamma_TUN) - gamma_SC simplifies to exactly 0",
      sp.simplify(g_sc - g_sup) == 0,
      f"simplify -> {sp.simplify(g_sc - g_sup)}")
print(f"     both equal {sp.simplify(g_sc)}  in sqrt-mass variables:")
print("     an exact identity for ANY isotope triple, not a coincidence of 1:2:3.")

print("\n5. physically admissible region: gating amplitude below the barrier width")
# w/2 = kappa sigma^2 and A = kappa R0^2, so sigma/R0 = sqrt(w/(2A)).
for frac in [0.1, 0.2, 1/3]:
    best = -np.inf
    for A in np.logspace(-2, 4, 400):
        w = 2 * A * frac ** 2
        best = max(best, Q.gamma_tunneling(A, w))
    print(f"     sigma/R0 <= {frac:.3f}:  max gamma_TUN = {best:.6f}")

print("\n6. monotonicity in w is NOT universal")
for A in [0.2, 1.0, 5.0, 20.0, 200.0]:
    g = np.array([Q.gamma_tunneling(A, w) for w in np.logspace(-4, 4, 600)])
    print(f"     A={A:7.1f}: gamma in [{g.min():.4f}, {g.max():.4f}]  "
          f"monotone={bool(np.all(np.diff(g) > -1e-12))}")

print("\n7. sanity: KIEs ordered and above unity everywhere")
bad = 0
for A in np.logspace(-2, 3, 60):
    for w in np.logspace(-4, 4, 60):
        lh, ld = Q.log_kie_intrinsic(A, w)
        if not (lh > ld > 0):
            bad += 1
check("K_HT^int > K_DT^int > 1 everywhere", bad == 0, f"{bad} violations")

print("\n8. commitment map (inherited from the classical analysis)")
check("c->inf recovers the intrinsic effect",
      close(Q.observed_from_intrinsic(12.0, 1e12), 12.0, 1e-5))
check("c->0 fully masks it", close(Q.observed_from_intrinsic(12.0, 1e-12), 1.0, 1e-6))
check("masking is strict for finite c", 1.0 < Q.observed_from_intrinsic(12.0, 2.5) < 12.0)

print("\n9. network reduction (Proposition S6)")
_x, _phi, _q = sp.symbols("x phi q", positive=True)
_N = (_q + 1 + _phi) / (1 + _phi)
_Kgen = _N * (_x + _phi) / (_x + _phi + _q)
_Y, _c = (_x + _phi) / (1 + _phi), _q / (1 + _phi)
check("K(x;phi,q) = K_series(Y;c) symbolically",
      sp.simplify(sp.expand(_Kgen - _Y * (1 + _c) / (_Y + _c))) == 0)
check("the map fixes x=1", sp.simplify(_Kgen.subs(_x, 1) - 1) == 0)
check("contraction gives Y-1 = (x-1)/(1+phi), so L_H is bypass invariant",
      sp.simplify((_Y - 1) - (_x - 1) / (1 + _phi)) == 0)

print(f"\nfailures: {FAIL}")
raise SystemExit(1 if FAIL else 0)
