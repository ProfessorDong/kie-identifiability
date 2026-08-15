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


def dF_dlnc(KH, KD, c, r=1.0, gamma=GSC):
    """Sensitivity of the completed offset to the commitment.

    c is c_H; the D/T comparison carries r*c, so d ln x_D/d ln c_H picks up the
    denominator 1 + r c - K_DT rather than 1 + c - K_DT.  Setting r = 1 recovers
    the shared-reference expression.
    """
    return (1.0 - KH)/(1.0 + c - KH) - gamma*(1.0 - KD)/(1.0 + r*c - KD)


def commitment_precision_needed(KH, KD, c, r=1.0, target=abs(F0)/2, gamma=GSC):
    """Relative precision on c_H that fixes F to +/- `target`."""
    d = abs(dF_dlnc(KH, KD, c, r, gamma))
    return np.inf if d == 0 else target/d


# ------------------------------------- systems completed from published kinetics
def yadh_completion(a, KH=7.13, KD=1.73, r=1.31):
    """Yeast ADH closed by Klinman's dissociation ratio a = k_-1/k_cat.

    a is referenced to the protiated substrate, our commitment to tritium, so
    c = k_off/k_T = a x_H; combining with x_H = K_HT c/(c - (K_HT - 1)) gives
    c = (K_HT - 1) + K_HT a with no iteration.  Isotope effects are Cha, Murray
    & Klinman (1989); a = 1.3 to 7.3 is Klinman (1976), Table IV.
    """
    c = (KH - 1.0) + KH*a
    # the two tritium references differ at the non-transferred position, so the
    # D/T comparison carries commitment r*c, not c (Proposition S4).  Using one
    # shared c here would contradict the result the paper proves.
    xh, xd = intrinsic_from_commitment(KH, c), intrinsic_from_commitment(KD, r*c)
    return dict(c=c, cD=r*c, xH=xh, xD=xd, gamma=np.log(xh)/np.log(xd),
                F=np.log(xh) - GSC*np.log(xd))


def bsao_completion(mask, KH=35.2, KD=3.07, r=1.14):
    """Bovine serum amine oxidase closed by a measured masking factor.

    Grant & Klinman (1989) Table IV reports the steady-state and pre-steady-state
    H/D effects side by side, and their ratio is the masking factor directly:
    with K = x(1+c)/(x+c) for both isotopes, (K_HT/K_DT)/(x_H/x_D) collapses to
    (x_D + c)/(x_H + c), a function of c alone.  Inverting it is the whole step.
    """
    g = lambda c: ((1 + c)*(intrinsic_from_commitment(KD, r*c) + r*c)
                   / ((intrinsic_from_commitment(KH, c) + c)*(1 + r*c)) - mask)
    lo, hi = KH + 1e-9, 1e12          # g is increasing in c, -> 1 as c -> inf
    for _ in range(400):
        mid = 0.5*(lo + hi)
        lo, hi = (mid, hi) if g(mid) < 0 else (lo, mid)
    c = 0.5*(lo + hi)
    xh, xd = intrinsic_from_commitment(KH, c), intrinsic_from_commitment(KD, r*c)
    return dict(c=c, cD=r*c, xH=xh, xD=xd, gamma=np.log(xh)/np.log(xd),
                F=np.log(xh) - GSC*np.log(xd))


def report_completions():
    """Reproduces the two completion tables of the supplement."""
    print("\nSystems completed from published kinetics, no new experiment\n")
    print(f"{'system':26s} {'c':>8s} {'x_H':>7s} {'x_D':>7s} "
          f"{'gamma_int':>10s} {'F_int':>8s}")
    Fobs_y = np.log(7.13) - GSC*np.log(1.73)
    for lab, a in (("YADH  a = 1.3", 1.3), ("YADH  a = 7.3", 7.3)):
        r = yadh_completion(a)
        print(f"{lab:26s} {r['c']:8.1f} {r['xH']:7.2f} {r['xD']:7.3f} "
              f"{r['gamma']:10.2f} {r['F']:+8.3f}")
    print(f"   half-line endpoint was F > {Fobs_y:+.3f}; the completion moves it up by "
          f"{(yadh_completion(7.3)['F']-Fobs_y)/abs(F0):.1f} to "
          f"{(yadh_completion(1.3)['F']-Fobs_y)/abs(F0):.1f} |F0|")
    for lab, m in (("BSAO  -1 sigma", 0.948 - 0.044), ("BSAO  central", 0.948),
                   ("BSAO  +1 sigma", 0.948 + 0.044)):
        r = bsao_completion(m)
        print(f"{lab:26s} {r['c']:8.0f} {r['xH']:7.2f} {r['xD']:7.3f} "
              f"{r['gamma']:10.3f} {r['F']:+8.3f}")
    print(f"   both intrinsic exponents lie below gamma_SC = {GSC:.3f}, and the "
          f"offsets below F0 = {F0:+.4f}")
    print("   note: c is poorly determined near the no-masking limit (a 5% shift in\n"
          "   the masking factor moves c by a factor 6) while F moves by 0.04 |F0|.\n"
          "   The offset is the well-conditioned coordinate, not the commitment.")



# ------------------------------------- Grant & Klinman 1989 Table IV, verbatim
# (k_H/k_D) at substrate saturation, stopped-flow against steady state.  The
# stopped-flow column used dideuterated benzylamine, so both columns carry the
# primary and secondary effect together; the steady-state column is derived from
# the competitive tritium effects of Tables I and II.
BSAO_TABLE_IV = np.array([
    #  T(C)  stopped-flow  +-     steady-state  +-
    [  0.0,   14.8, 0.7,   19.0, 0.8],
    [  5.0,   15.7, 1.2,   18.7, 0.2],
    [ 15.0,   16.6, 1.4,   16.5, 0.4],
    [ 25.0,   16.1, 1.8,   13.5, 0.4],
    [ 35.0,   13.8, 1.1,   13.0, 0.4],
    [ 45.0,   10.4, 0.9,   10.4, 0.2],
])


def bsao_masking_by_temperature():
    """Per-temperature masking factor m = steady/stopped, with propagated error."""
    T, sf, sfe, ss, sse = BSAO_TABLE_IV.T
    m = ss / sf
    me = m * np.sqrt((sse / ss) ** 2 + (sfe / sf) ** 2)
    return T, m, me


def bsao_homogeneity(rows=None):
    """Inverse-variance mean and Cochran Q for a subset of the temperatures."""
    T, m, me = bsao_masking_by_temperature()
    sel = np.ones_like(T, bool) if rows is None else rows
    w = 1.0 / me[sel] ** 2
    mu = float((w * m[sel]).sum() / w.sum())
    Q = float((w * (m[sel] - mu) ** 2).sum())
    return mu, float(1 / np.sqrt(w.sum())), Q, int(sel.sum()) - 1


def _mc_bsao(n, seed, m_val, m_err, KH=(35.2, 0.8), KD=(3.07, 0.07), r=(1.1370, 0.0193)):
    rng = np.random.default_rng(seed)
    kh, kd, rr = (rng.normal(*KH, n), rng.normal(*KD, n), rng.normal(*r, n))
    mm = rng.normal(m_val, m_err, n)
    Fobs = np.log(kh) - GSC * np.log(kd)
    lo = np.maximum(kh, kd) - 1 + 1e-9
    hi = np.full(n, 1e12)
    mc = np.clip(mm, 1e-6, 1 - 1e-9)
    for _ in range(200):
        mid = np.sqrt(lo * hi)
        f = ((1 + mid) * (intrinsic_from_commitment(kd, rr * mid) + rr * mid)
             / ((intrinsic_from_commitment(kh, mid) + mid) * (1 + rr * mid)) - mc)
        lo, hi = np.where(f < 0, mid, lo), np.where(f < 0, hi, mid)
    c = np.sqrt(lo * hi)
    xh, xd = intrinsic_from_commitment(kh, c), intrinsic_from_commitment(kd, rr * c)
    F = np.where(mm >= 1.0, Fobs, np.log(xh) - GSC * np.log(xd))
    g = np.log(xh) / np.log(xd)
    ok = (kh > kd) & (kd > 1) & np.isfinite(F)
    return F[ok], g[ok]


def _mc_yadh(n, seed, a, KH=(7.13, 0.07), KD=(1.73, 0.02), r=(1.3107, 0.0164)):
    rng = np.random.default_rng(seed)
    kh, kd, rr = (rng.normal(*KH, n), rng.normal(*KD, n), rng.normal(*r, n))
    aa = rng.uniform(*a, n) if isinstance(a, tuple) else np.full(n, a)
    c = (kh - 1) + kh * aa
    xh, xd = intrinsic_from_commitment(kh, c), intrinsic_from_commitment(kd, rr * c)
    F = np.log(xh) - GSC * np.log(xd)
    g = np.log(xh) / np.log(xd)
    ok = (kh > kd) & (kd > 1) & np.isfinite(F)
    return F[ok], g[ok]


def report_uncertainty(n=400000, seed=20260815):
    """Propagate every measured input through both completions."""
    T, m, me = bsao_masking_by_temperature()
    print("\nBSAO masking factor, temperature by temperature (Grant & Klinman Table IV)")
    print(f"  {'T(C)':>5} {'m':>8} {'+-':>7} {'(m-1)/sd':>9}")
    for i in range(len(T)):
        print(f"  {T[i]:5.0f} {m[i]:8.3f} {me[i]:7.3f} {(m[i]-1)/me[i]:+9.2f}")
    mu, se, Q, dof = bsao_homogeneity()
    print(f"  all six pooled : {mu:.4f} +- {se:.4f}   Q = {Q:.2f} on {dof} dof")
    mu2, se2, Q2, d2 = bsao_homogeneity(T >= 15)
    print(f"  15-45 C pooled : {mu2:.4f} +- {se2:.4f}   Q = {Q2:.2f} on {d2} dof")
    print("  masking cannot exceed unity; the 0 C value does so by 3.5 standard")
    print("  deviations, so a single masking factor does not describe these data.")

    print("\nCompletions with every measured input propagated")
    F, g = _mc_bsao(n, seed, m[3], me[3])
    q = np.percentile(F, [2.5, 50, 97.5]); qg = np.percentile(g, [2.5, 50, 97.5])
    print(f"  BSAO at 25 C, m = {m[3]:.3f} +- {me[3]:.3f} (matched condition)")
    print(f"    F_int {q[1]:+.3f}  95% [{q[0]:+.3f},{q[2]:+.3f}]   "
          f"P(F<F0) = {np.mean(F < F0):.2f}")
    print(f"    gamma_int {qg[1]:.2f} 95% [{qg[0]:.2f},{qg[2]:.2f}]   "
          f"P(gamma<gamma_SC) = {np.mean(g < GSC):.2f}")
    print("    -> the interval straddles both thresholds: not informative.")
    for a in (1.3, 2.3, 7.3):
        F, g = _mc_yadh(n, seed, a)
        q = np.percentile(F, [2.5, 50, 97.5]); qg = np.percentile(g, [2.5, 50, 97.5])
        print(f"  YADH a = {a:.1f}: F_int {q[1]:+.3f} 95% [{q[0]:+.3f},{q[2]:+.3f}]"
              f"  gamma_int {qg[1]:.2f} [{qg[0]:.2f},{qg[2]:.2f}]"
              f"  P(F>0) = {np.mean(F > 0):.4f}")
    F, _ = _mc_yadh(n, seed, (1.3, 7.3))
    q = np.percentile(F, [2.5, 50, 97.5])
    print(f"  YADH, a uniform on [1.3,7.3]: F_int {q[1]:+.3f} "
          f"95% [{q[0]:+.3f},{q[2]:+.3f}]  P(F>0) = {np.mean(F > 0):.4f}")


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
    for lab, KH, KD, r in [("YADH", 7.13, 1.73, 1.3107),
                           ("BSAO", 35.2, 3.07, 1.1370),
                           ("LADH F93W", 7.755, 1.858, 1.27)]:
        Fobs = np.log(KH) - GSC*np.log(KD)
        for mult in (2.0, 3.0, 5.0):
            c = (KH - 1.0)*mult
            xh = intrinsic_from_commitment(KH, c); xd = intrinsic_from_commitment(KD, r*c)
            Fi = np.log(xh) - GSC*np.log(xd)
            need = commitment_precision_needed(KH, KD, c, r)
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
    report_completions()
    report_uncertainty()
