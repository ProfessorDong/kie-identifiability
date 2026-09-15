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

print("\n10. unequal-reference endpoint (Corollary S4.2)")
_c, _g, _a, _b, _r = sp.symbols("c gamma a b r", positive=True)
_xH = (_a + 1) * _c / (1 + _c - (_a + 1))
_xD = (_b + 1) * (_r * _c) / (1 + _r * _c - (_b + 1))
_F = sp.log(_xH) - _g * sp.log(_xD)
_stated = (_c * (_g * _b - _a * _r) - _a * _b * (_g - 1)) / (_c * (_c - _a) * (_r * _c - _b))
check("dF_r/dc closed form", sp.simplify(sp.expand(sp.diff(_F, _c) - _stated)) == 0)
_cs = _a * _b * (_g - 1) / (_g * _b - _a * _r)
check("c* is stationary", sp.simplify(sp.diff(_F, _c).subs(_c, _cs)) == 0)
check("x_H* closed form",
      sp.simplify(_xH.subs(_c, _cs) - (_a + 1) * _b * (_g - 1) / (_a * _r - _b)) == 0)
check("x_D* closed form",
      sp.simplify(_xD.subs(_c, _cs) - (_b + 1) * _r * _a * (_g - 1) / (_g * (_a * _r - _b))) == 0)
check("reduces to Proposition S2 at r=1",
      sp.simplify(_stated.subs(_r, 1)
                  - (_c * (_g * _b - _a) - _a * _b * (_g - 1)) / (_c * (_c - _a) * (_c - _b))) == 0)

print("\n11. binding-aware map (Proposition S8)")
_al,_be,_x,_c,_u = sp.symbols("alpha beta x c u", positive=True)
_K = _al*_x*(1+_c)/(_x+_be*_c)
check("reduces to the series map at alpha=beta=1",
      sp.simplify(_K.subs({_al:1,_be:1}) - _x*(1+_c)/(_x+_c)) == 0)
check("K/alpha = series map on u = x/beta",
      sp.simplify(_K/_al - (_u*(1+_c)/(_u+_c)).subs(_u, _x/_be)) == 0)
_aH,_aD,_bH,_bD,_xH,_xD,_gm = sp.symbols(
    "alpha_H alpha_D beta_H beta_D x_H x_D gamma", positive=True)
_Fo = sp.log(_aH*_xH/_bH) - _gm*sp.log(_aD*_xD/_bD)
_Fi = sp.log(_xH) - _gm*sp.log(_xD)
_Fb = sp.log(_aH/_bH) - _gm*sp.log(_aD/_bD)
check("F_obs = F_int + F_bind in the commitment-free limit",
      sp.simplify(_Fo - _Fi - _Fb) == 0)
check("F_bind vanishes when binding effects are mass scaled",
      sp.simplify(_Fb.subs(_aH/_bH, (_aD/_bD)**_gm)) == 0 or
      sp.simplify((sp.log((_aD/_bD)**_gm) - _gm*sp.log(_aD/_bD))) == 0)

print("\n12. envelope convergence rate at fixed K (Supplementary Note 2), 60 digits")
# For K > sqrt(t) the constrained supremum approaches F0 as 1/w with coefficient
# (d-1)/d (ln K - ln t/2); only at K = sqrt(t) does that coefficient vanish and the
# gap fall as 1/w^2.  The supplement said w^-2 for every K until 2026-09-15; it is
# the DERIVATIVE, eq. (sm-dXD), that falls as w^-2.
import mpmath as mp
import masses as _M
mp.mp.dps = 60
_mC = mp.mpf(12)
_red = lambda m: _mC * mp.mpf(m) / (_mC + mp.mpf(m))
_d = mp.sqrt(_red(_M.M_D_ATOMIC) / _red(_M.M_H_ATOMIC))
_t = mp.sqrt(_red(_M.M_T_ATOMIC) / _red(_M.M_H_ATOMIC))
_gg = (1 - 1 / _t) / (1 / _d - 1 / _t)
_F0 = mp.log(_t) / 2 - _gg * mp.log(_t / _d) / 2
_p = lambda w, m: mp.log((1 + w * _t) / (1 + w * m)) / 2
_q = lambda w, m: (_t - m) / ((1 + w * m) * (1 + w * _t))
def _gap(K, w):
    L = mp.log(K)
    A = (L - _p(w, 1)) / _q(w, 1)
    return _F0 - (L - _gg * (_p(w, _d) + A * _q(w, _d)))
for K in (mp.mpf(7), mp.mpf(2)):
    coef = (_d - 1) / _d * (mp.log(K) - mp.log(_t) / 2)
    w = mp.mpf(10) ** 8
    check(f"K={float(K):g}: w (F0-F) -> (d-1)/d (ln K - ln t/2)",
          abs(w * _gap(K, w) / coef - 1) < 1e-6,
          f"{mp.nstr(w * _gap(K, w), 10)} vs {mp.nstr(coef, 10)}")
_ks = mp.sqrt(_t)
_r8, _r6 = [mp.mpf(10) ** e * mp.mpf(10) ** e * _gap(_ks, mp.mpf(10) ** e) for e in (8, 6)]
check("K=sqrt(t): w^2 (F0-F) converges (gap falls as w^-2)",
      abs(_r8 / _r6 - 1) < 1e-3 and _r8 > 0, f"{mp.nstr(_r6, 8)} -> {mp.nstr(_r8, 8)}")

print("\n13. bounded-bypass endpoint: stationarity quadratic and both branches")
_KH, _KD, _ph, _rr, _c, _gm2 = sp.symbols("K_H K_D phi r c gamma", positive=True)
_a2, _b2 = _KH - 1, _KD - 1
_rho = _rr * (1 + _ph) / (1 + _rr * _ph)
_PH, _PD = _KH + _a2 * _ph, _KD + _b2 * _rr * _ph
_xH2 = (_c * _PH + _ph * _a2) / (_c - _a2)
_xD2 = (_rho * _c * _PD + _rr * _ph * _b2) / (_rho * _c - _b2)
# the pair reproduces both observations through the bypass map (q = c(1+phi))
_Kmap = lambda x, ph, q: (x + ph) / (1 + ph) * (q + 1 + ph) / (q + x + ph)
check("intrinsic pair reproduces K_HT and K_DT",
      sp.simplify(_Kmap(_xH2, _ph, _c * (1 + _ph)) - _KH) == 0
      and sp.simplify(_Kmap(_xD2, _rr * _ph, _rr * _c * (1 + _ph)) - _KD) == 0)
_dF = sp.diff(sp.log(_xH2), _c) - _gm2 * sp.diff(sp.log(_xD2), _c)
_stat = (_a2 * _KH * (_PD * _rho * _c + _rr * _ph * _b2) * (_rho * _c - _b2)
         - _gm2 * _rr * _b2 * _KD * (_PH * _c + _ph * _a2) * (_c - _a2))
_den = (_PH * _c + _ph * _a2) * (_c - _a2) * (_PD * _rho * _c + _rr * _ph * _b2) * (_rho * _c - _b2)
check("dF/dc = -(1+phi) * quadratic / positive denominator",
      sp.simplify(_dF + (1 + _ph) * _stat / _den) == 0)
import network_geometry as _NG
_n, _worst, _int, _crit, _n1 = _NG.check_endpoint_exact(n=3000)
check("endpoint_exact = direct profile on random admissible triples",
      _worst < 1e-7, f"{_n} triples, {_int} interior, max diff {_worst:.1e}")
check("c -> infinity criterion agrees with the profile (r = 1)", _crit == 0,
      f"{_crit} disagreements in {_n1} triples")

print("\n14. yeast reference asymmetry: the joint relations reduce to a quadratic")
_a3, _hh, _dd, _SH3, _SD3, _r3 = sp.symbols("a h d S_H S_D r", positive=True)
_c3 = _hh * (1 + _a3) - 1
_s3 = (_SH3 * (1 + _a3) - 1) / _a3
_qq = (_r3 * _c3 - (_dd - 1)) / _dd
_sD3 = (_SD3 * (_qq + 1) - 1) / _qq
# q d (r sigma_D - sigma_H) is exactly the quadratic, so the two share their roots
check("r sigma_D = sigma_H  <=>  S_D c r^2 + (S_D - d - s c) r + s(d-1) = 0",
      sp.simplify((_r3 * _sD3 - _s3) * _qq * _dd
                  - (_SD3 * _c3 * _r3**2 + (_SD3 - _dd - _s3 * _c3) * _r3
                     + _s3 * (_dd - 1))) == 0)
import completion as _C
_nq, _nanq, _dq, _pq = _C.check_joint_r_quadratic()
check("quadratic root = bisection over the Monte Carlo inputs",
      _nanq == 0 and _dq < 1e-10 and _pq < 1, f"{_nq} draws, max diff {_dq:.1e}, root product <= {_pq:.3f}")

print("\n15. mixed-labeling secondary design: forward commitments are ordered")
# c_i = k_off/k_iT with reference molecules HT and DT (transferred isotope first),
# so c_D/c_H = k_HT/k_DT, the primary H/D effect: >= 1 for a normal primary effect.
import mixed_label as _ML
_n, _bad, _worst, _lift = _ML.check_ordered()
check("ordered maps never raise F above F_int (any reference exponent > 1)",
      _bad == 0, f"{_bad} of {_n}, largest excess {_worst:.1e}")
check("anti-ordered maps do lift pairs on the ray", _lift > 0, f"{_lift} of {_n}")
check("4.8 locus: ordered c_H=1, c_D=5 deflates; anti-ordered 5, 1 inflates",
      _ML.observed_exponent(1.10, 4.8, 1.0, 5.0) < 4.8 < _ML.observed_exponent(1.10, 4.8, 5.0, 1.0),
      f"{_ML.observed_exponent(1.10, 4.8, 1.0, 5.0):.2f} and {_ML.observed_exponent(1.10, 4.8, 5.0, 1.0):.2f}")

print("\n16. binding correction: exact shift on the open branch only")
from partial_id import F_min_binding as _Fb, F_min_exact as _Fme
import masses as _M2
_G = _M2.gamma_sc("C")
def _profile_binding(KH, KD, aH, bH, aD, bD, n=2000000):
    kh, kd = KH / aH, KD / aD
    cs = (max(kh, kd) - 1.0) * (1 + np.geomspace(1e-12, 1e12, n))
    xh, xd = bH * kh * cs / (1 + cs - kh), bD * kd * cs / (1 + cs - kd)
    ok = (xh > xd) & (xd > 1)
    return float(np.min(np.log(xh[ok]) - _G * np.log(xd[ok])))
_args = (1.05, 1.05, 0.95, 0.95)          # alpha_H, beta_H, alpha_D, beta_D: F_bind = 0
_closed = _Fb(3.0, 1.8, *_args) - _Fme(3.0, 1.8)[0]
check("closed form = direct profile with binding, closed branch",
      abs(_Fb(3.0, 1.8, *_args) - _profile_binding(3.0, 1.8, *_args)) < 1e-6)
check("closed branch (3, 1.8): F_bind = 0 yet the endpoint moves",
      abs(_closed) > 0.1, f"shift {_closed:+.6f}")
_open = _Fb(7.13, 1.73, *_args) - _Fme(7.13, 1.73)[0]
check("open branch (yeast): F_bind = 0 leaves the endpoint unchanged",
      abs(_open) < 1e-12, f"shift {_open:+.2e}")
_rng = np.random.default_rng(16)
_bad = 0
for _ in range(20000):
    KD = 1 + _rng.uniform(0.05, 2.0); KH = KD + _rng.uniform(0.05, 20.0)
    aH, bH, aD, bD = _rng.uniform(0.9, 1.1, 4)
    kh, kd = KH / aH, KD / aD
    if not (kh > kd > 1) or (kh - 1) / (kd - 1) < _G:
        continue
    Fb = np.log(aH / bH) - _G * np.log(aD / bD)
    _bad += abs(_Fb(KH, KD, aH, bH, aD, bD) - (np.log(KH) - _G * np.log(KD) - Fb)) > 1e-10
check("open branch: endpoint = F_obs - F_bind exactly", _bad == 0, f"{_bad} violations")

print(f"\nfailures: {FAIL}")
raise SystemExit(1 if FAIL else 0)
