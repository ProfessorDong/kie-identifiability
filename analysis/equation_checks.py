"""Every displayed equation of the main text and the supplement, checked.

An equation audit (2026-10-22) went through the 8 displayed equations of the
main text, the 68 of the supplement and the inline formulas the text relies on.
Each identity below is verified symbolically with sympy (exact zero after
simplification), or numerically where it is a limit or an inequality.  The
labels are the LaTeX labels of the supplement (eq:sm-...) or of the main text.
Numerical results quoted beside the equations are checked by the scripts named
in the supplement and by audit_numbers.py.

Run from analysis/:   python equation_checks.py      (exits non-zero on a failure)
"""
from __future__ import annotations

import sys

import mpmath as mp
import sympy as sp

import masses as M

FAILS = []


def ok(name, cond):
    print(("PASS " if cond else "FAIL ") + name)
    if not cond:
        FAILS.append(name)


def zero(expr):
    return sp.simplify(expr) == 0


def proofs():
    """Supplementary Note 2: overlap, envelope, identified set."""
    R, R0, a, b, sig, w, A, t, d, mu, K, c, x, KH, KD, g, L = sp.symbols(
        "R R0 a b sigma w A t d mu K c x K_H K_D gamma L", positive=True)
    I = sp.sqrt(a / sp.pi) * sp.integrate(sp.exp(-a * (R - R0) ** 2 - b * R ** 2), (R, -sp.oo, sp.oo))
    ok("eq:sm-overlap (Proposition S1)",
       zero(I.subs(a, 1 / (2 * sig ** 2))
            - (1 + 2 * b * sig ** 2) ** sp.Rational(-1, 2) * sp.exp(-b * R0 ** 2 / (1 + 2 * b * sig ** 2))))
    ok("completing the square",
       zero(sp.expand(-a * (R - R0) ** 2 - b * R ** 2
                      + (a + b) * (R - a * R0 / (a + b)) ** 2 + a * b * R0 ** 2 / (a + b))))
    kap, muH = sp.symbols("kappa mu_H", positive=True)

    def lnS(m):
        bL = kap * sp.sqrt(muH) * m
        return -sp.log(1 + 2 * bL * sig ** 2) / 2 - bL * R0 ** 2 / (1 + 2 * bL * sig ** 2)
    sub = {sig: sp.sqrt(w / (2 * kap * sp.sqrt(muH))), R0: sp.sqrt(A / (kap * sp.sqrt(muH)))}
    XLp = sp.log((1 + w * t) / (1 + w * mu)) / 2 + A * (t - mu) / ((1 + w * mu) * (1 + w * t))
    ok("eq:sm-XL / main eq:XL", zero(sp.expand_log((lnS(mu) - lnS(t)).subs(sub) - XLp, force=True)))
    p = lambda m: sp.log((1 + w * t) / (1 + w * m)) / 2
    q = lambda m: (t - m) / ((1 + w * m) * (1 + w * t))
    G = d * (t - 1) / (t - d)
    ok("eq:sm-dFdA", zero(q(1) - G * q(d) + (t - 1) * (d - 1) / ((1 + w) * (1 + w * d) * (1 + w * t))))
    ok("eq:sm-pq", all(zero(sp.diff(p(m), w) - q(m) / 2) for m in (1, d, t)))
    r = q(d) / q(1)
    ok("eq:sm-r", zero(r - (t - d) * (1 + w) / ((t - 1) * (1 + w * d)))
       and zero(sp.diff(r, w) + (t - d) * (d - 1) / ((t - 1) * (1 + w * d) ** 2)))
    Aw = (L - p(1)) / q(1)
    ok("eq:sm-dXD", zero(sp.diff(p(d) + Aw * q(d), w) - Aw * q(1) * sp.diff(r, w)))
    ws = (K ** 2 - 1) / (t - K ** 2)
    ok("w*, 1 + w* t, 1 + w* d, p_D(w*)",
       zero(sp.exp(2 * p(1)).subs(w, ws) - K ** 2)
       and zero(1 + ws * t - K ** 2 * (t - 1) / (t - K ** 2))
       and zero(1 + ws * d - (t - d + K ** 2 * (d - 1)) / (t - K ** 2))
       and zero(sp.exp(2 * p(d)).subs(w, ws) - K ** 2 * (t - 1) / (t - d + K ** 2 * (d - 1))))
    ok("ridge: zeta_H = gamma zeta_D", zero((1 - 1 / t) - G * (1 / d - 1 / t)))
    s = sp.symbols("s", positive=True)
    XLw = lambda m: sp.log((1 + w * t) / (1 + w * m)) / 2 + s * w ** 2 * (t - m) / ((1 + w * m) * (1 + w * t))
    ok("ridge limit X_L -> P_L + s zeta_L",
       zero(sp.limit(XLw(d) - (sp.log(t / d) / 2 + s * (1 / d - 1 / t)), w, sp.oo)))
    xc = lambda KK: KK * c / (1 + c - KK)
    F = sp.log(xc(KH)) - g * sp.log(xc(KD))
    aa, bb = KH - 1, KD - 1
    ok("eq:sm-dFdc", zero(sp.diff(F, c) - (c * (g * bb - aa) - aa * bb * (g - 1)) / (c * (c - aa) * (c - bb))))
    cs = aa * bb * (g - 1) / (g * bb - aa)
    ok("eq:sm-xstar / main eq:xstar",
       zero(xc(KH).subs(c, cs) - KH * bb * (g - 1) / (aa - bb))
       and zero(xc(KD).subs(c, cs) - KD * aa * (g - 1) / (g * (aa - bb))))
    ok("eq:sm-xfromc inverts main eq:map", zero((x * (1 + c) / (x + c)).subs(x, xc(K)) - K))
    tv = sp.symbols("tv")
    h = sp.log(sp.exp(tv) * (1 + c) / (sp.exp(tv) + c))
    ok("eq:sm-hprops", zero(sp.diff(h, tv) - c / (sp.exp(tv) + c))
       and zero(sp.diff(h, tv, 2) + c * sp.exp(tv) / (sp.exp(tv) + c) ** 2))
    ok("eq:sm-dhdc", zero(sp.diff(h, c) - (sp.exp(tv) - 1) / ((1 + c) * (sp.exp(tv) + c))))
    Cf = sp.symbols("C_f", positive=True)
    ok("additive = reciprocal at C_L = x/c (eq:sm-chacd)", zero((x + x / c) / (1 + x / c) - x * (1 + c) / (x + c)))
    XH, XD = sp.symbols("X_H X_D", positive=True)
    ok("main: exponent = gamma + F/X_D", zero(XH / XD - (g + (XH - g * XD) / XD)))
    # numerical constants with the masses of masses.py
    mH, mD, mT, mC = M.M_H_ATOMIC, M.M_D_ATOMIC, M.M_T_ATOMIC, 12.0
    mu_ = lambda m: m * mC / (m + mC)
    tt = mp.sqrt(mu_(mT) / mu_(mH))
    dd = mp.sqrt(mu_(mD) / mu_(mH))
    GG = dd * (tt - 1) / (tt - dd)
    ok("gamma_SC = d(t-1)/(t-d)", abs(GG - M.gamma_sc("C")) < 1e-9)
    F0 = (mp.log(tt) - GG * mp.log(tt / dd)) / 2
    ok("F0 = P_H - gamma P_D", abs(F0 - M.offset_F0("C")) < 1e-9)
    env = lambda KK: mp.log(KK) - GG / 2 * mp.log(KK ** 2 * (tt - 1) / (tt - dd + KK ** 2 * (dd - 1)))
    ok("eq:sm-env continuous at sqrt t", abs(env(mp.sqrt(tt)) - F0) < 1e-12)


def networks():
    """Supplementary Note 3."""
    x, c, phi, q, t = sp.symbols("x c phi q t", positive=True)
    n1, n0, d1, d0, kT, kon, koff, k = sp.symbols("n1 n0 d1 d0 k_T k_on k_off k", positive=True)
    A, B, C, D = sp.symbols("A B C D", positive=True)
    ok("eq:sm-series", zero(1 / (kon * k / (k + koff)) - (1 / kon + koff / (kon * k))))
    VK = lambda kk: (n1 * kk + n0) / (d1 * kk + d0)
    Aa = n1 * kT * (d0 + d1 * kT); Bb = n0 * (d0 + d1 * kT)
    Cc = d1 * kT * (n1 * kT + n0); Dd = d0 * (n1 * kT + n0)
    ok("eq:sm-moebius", zero(VK(x * kT) / VK(kT) - (Aa * x + Bb) / (Cc * x + Dd)) and zero(Aa + Bb - Cc - Dd))
    hh = sp.log((A * sp.exp(t) + B) / (C * sp.exp(t) + D))
    ok("eq:sm-hcurv", zero(sp.diff(hh, t, 2) + (A * D - B * C) * (A * C * sp.exp(2 * t) - B * D) * sp.exp(t)
                           / ((A * sp.exp(t) + B) ** 2 * (C * sp.exp(t) + D) ** 2)))
    ok("x* = sqrt(BD/AC) = sqrt(n0 d0/(n1 d1))/k_T", zero(sp.sqrt(Bb * Dd / (Aa * Cc)) - sp.sqrt(n0 * d0 / (n1 * d1)) / kT))
    ser = lambda Y, cc: Y * (1 + cc) / (Y + cc)
    ok("Proposition S5: B = 0 gives main eq:map", zero((C + D) * x / (C * x + D) - ser(x, D / C)))
    K = sp.symbols("K", positive=True)
    ok("Proposition S5: C = 0 inversion", zero(((x + phi) / (1 + phi)).subs(x, K + phi * (K - 1)) - K))
    Dm = A + B - C
    ok("Proposition S5: mixed = series on the contracted effect",
       zero((A * x + B) / (C * x + Dm) - ser((x + B / A) / (1 + B / A), A / C - 1)))
    N = (q + 1 + phi) / (1 + phi)
    ok("eq:sm-reduction (Proposition S6)", zero((N * x + N * phi) / (x + phi + q) - ser((x + phi) / (1 + phi), q / (1 + phi))))
    KH, KD, g, r = sp.symbols("K_H K_D gamma r", positive=True)
    a, b = KH - 1, KD - 1
    Y = lambda KK: KK * c / (1 + c - KK)
    xH = Y(KH) * (1 + phi) - phi
    xD = Y(KD) * (1 + phi) - phi
    ok("bypass x_L(c)", zero(xH - ((KH + a * phi) * c + phi * a) / (c - a)))
    F = sp.log(xH) - g * sp.log(xD)
    PH, PD = KH + a * phi, KD + b * phi
    stat = g * KD * b * (PH * c + phi * a) * (c - a) - KH * a * (PD * c + phi * b) * (c - b)
    ratio = sp.simplify(sp.numer(sp.together(sp.diff(F, c))) / stat)
    ok("eq:sm-bypassstat (numerator = (1 + phi) x quadratic)", zero(ratio - (1 + phi)))
    ok("eq:sm-Ephi", zero(sp.limit(F, c, sp.oo) - (sp.log(KH * (1 + phi) - phi) - g * sp.log(KD * (1 + phi) - phi))))
    lead = sp.limit(sp.diff(F, c) * c ** 2, c, sp.oo)
    ok("eq:sm-limcrit (leading coefficient = -(1 + phi) x criterion)",
       zero(lead + (1 + phi) * (KH * a / (KH + a * phi) - g * KD * b / (KD + b * phi))))
    xc = lambda KK, cc: KK * cc / (1 + cc - KK)
    Fr = sp.log(xc(KH, c)) - g * sp.log(xc(KD, r * c))
    ok("eq:sm-dFdc-r", zero(sp.diff(Fr, c) - (c * (g * b - a * r) - a * b * (g - 1)) / (c * (c - a) * (r * c - b))))
    cs = a * b * (g - 1) / (g * b - a * r)
    ok("eq:sm-xstar-r", zero(xc(KH, cs) - KH * b * (g - 1) / (a * r - b))
       and zero(xc(KD, r * cs) - KD * r * a * (g - 1) / (g * (a * r - b))))


def binding_and_protocol():
    """Supplementary Notes 4, 5 and 6."""
    x, c, al, be, u, r, KH, KD, g = sp.symbols("x c alpha beta u r K_H K_D gamma", positive=True)
    ok("eq:sm-bindreduce (Proposition S7)", zero((al * x * (1 + c) / (x + be * c) / al).subs(x, be * u) - u * (1 + c) / (u + c)))
    kon, koff, kb, kT, k = sp.symbols("k_on k_off k_b k_T k", positive=True)
    aH, bH, q, ph = sp.symbols("a_H b_H q phi", positive=True)
    vk = lambda kc, ko, kf: kf * (kc + kb) / (kc + kb + ko)
    Ki = vk(x * kT, bH * koff, aH * kon) / vk(kT, koff, kon)
    joint = aH * (x + ph) / (1 + ph) * (q + 1 + ph) / (bH * q + x + ph)
    ok("eq:sm-joint (first passage)", zero(Ki.subs({kb: ph * kT, koff: q * kT}) - joint))
    pH, pD, sH, sD, sT, d = sp.symbols("p_H p_D s_H s_D s_T d", positive=True)
    v = lambda p, s: p * s / (d + p * s)
    Robs = (v(pH, sH) / v(pH, sT)) / (v(pD, sD) / v(pD, sT))
    ok("eq:sm-Robs", zero(Robs - (sH / sD) * (d + pH * sT) / (d + pH * sH) * (d + pD * sD) / (d + pD * sT)))
    xH = KH * c / (1 + c - KH)
    Mp = lambda xx, cc: xx * (1 + cc) / (xx + cc)
    ok("eq:sm-chi", zero(Mp(xH, r * c) - KH * (1 + r * c) / (r * (1 + c) - (r - 1) * KH)))


def completion_and_rest():
    """Supplementary Notes 8, 9, 12 and 13."""
    c, r, KH, KD, g, a, x = sp.symbols("c r K_H K_D gamma a x", positive=True)
    F = sp.log(KH * c / (1 + c - KH)) - g * sp.log(KD * r * c / (1 + r * c - KD))
    ok("eq:sm-dFdlnc", zero(c * sp.diff(F, c) - ((1 - KH) / (1 + c - KH) - g * (1 - KD) / (1 + r * c - KD))))
    ok("eq:sm-cfroma", zero(sp.solve(sp.Eq(c, a * KH * c / (c - (KH - 1))), c)[0] - ((KH - 1) + KH * a))
       if len(sp.solve(sp.Eq(c, a * KH * c / (c - (KH - 1))), c)) == 1 else
       any(zero(s_ - ((KH - 1) + KH * a)) for s_ in sp.solve(sp.Eq(c, a * KH * c / (c - (KH - 1))), c)))
    K1x = x * (1 + c) / (x + c)
    K2x = x * (1 + c / r) / (x + c / r)
    ok("eq:sm-twovisc (c_2 = c_1/r)", zero((r * K1x * K2x - K1x * K2x - r * K1x + K2x) / (r * K2x - r - K1x + 1) - x))
    ok("eq:sm-pool", zero(1 - (x + c) / (x * (1 + c)) - (1 - 1 / x) * c / (1 + c)))
    gam = M.gamma_sc("C")
    ok("eq:sm-delta-approx coefficient 3.933", abs(gam * (gam - 1) / 2 - 3.933) < 5e-4)
    k1, k2, k3, k5, Keq, E = sp.symbols("k1 k2 k3 k5 K_eq E", positive=True)
    VK = lambda a3, a4: k1 * a3 * k5 / (k2 * a4 + k2 * k5 + a3 * k5)
    Kobs = VK(x * k3, x * k3 / (Keq * E)) / VK(k3, k3 / Keq)
    Cf, Cr = k3 / k2, (k3 / Keq) / k5
    ok("eq:sm-rev", zero(Kobs - x * (1 + Cf + Cr) / (1 + Cf * x + Cr * x / E)))
    cf, cr, K = sp.symbols("C_f C_r K", positive=True)
    rev = x * (1 + cf + cr) / (1 + cf * x + cr * x / E)
    ok("reversible map = Cleland form at C_f^i = C_f x, C_r^i = C_r x/E",
       zero(rev - (x + cf * x + cr * x) / (1 + cf * x + cr * x / E)))
    Dfac = 1 + cf * (1 - K) + cr * (1 - K / E)
    ok("eq:sm-Frev inversion x = K/D", zero(rev.subs(x, K / Dfac) - K))
    EH, ED = sp.symbols("E_H E_D", positive=True)
    aa, bb = KH - 1, KD - 1
    u, v = KH / EH - 1, KD / ED - 1
    cfs, crs = (v - u) / (aa * v - bb * u), (aa - bb) / (aa * v - bb * u)
    ok("C_f*, C_r* make D_H = D_D = 0",
       zero(1 + cfs * (1 - KH) + crs * (1 - KH / EH)) and zero(1 + cfs * (1 - KD) + crs * (1 - KD / ED)))
    ok("delta = u(a - b) at the lower window edge", zero((aa * v - bb * u).subs(ED, KD / (u + 1)) - u * (aa - bb)))


def main():
    proofs()
    networks()
    binding_and_protocol()
    completion_and_rest()
    print(f"\n{len(FAILS)} failures")
    return 1 if FAILS else 0


if __name__ == "__main__":
    sys.exit(main())
