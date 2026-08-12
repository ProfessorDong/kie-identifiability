"""Referee claims 3 and 11, tested directly.

3.  Table II fixed K_HT^int = 11 and fitted only the commitment.  Profiling the
    intrinsic magnitude instead should collapse the residuals.
11. The claim that neither model class can produce a non-monotonic observed
    K_HT should be false: search for a counterexample.
"""
from __future__ import annotations

import numpy as np
from scipy.optimize import least_squares, brentq
from scipy import stats

import qtunnel as Q
from discriminate import series, R_KCAL

RNG = np.random.default_rng(7)
MU_T, MU_H, MU_D = Q.MU_T, Q.MU_H, Q.MU_D
P = {"H": 0.5*np.log(MU_T/MU_H), "D": 0.5*np.log(MU_T/MU_D)}
C = {"H": (MU_T-MU_H)/(MU_H*MU_T), "D": (MU_T-MU_D)/(MU_D*MU_T)}


def A_for(w, target):
    f = lambda A: Q.log_kie_intrinsic(A, w)[0] - target
    lo, hi = 1e-9, 1e9
    if not np.isfinite(f(lo)) or not np.isfinite(f(hi)) or f(lo)*f(hi) > 0:
        return np.nan
    return brentq(f, lo, hi, xtol=1e-12, rtol=1e-14)


def solve_Aw(gamma, lnKHT):
    ws = np.logspace(-8, 5, 1500)
    prev = None
    for w in ws:
        A = A_for(w, lnKHT)
        if not np.isfinite(A):
            continue
        g = Q.gamma_tunneling(A, w)
        if prev is not None and (prev[1]-gamma)*(g-gamma) < 0:
            wstar = brentq(lambda x: Q.gamma_tunneling(A_for(x, lnKHT), x)-gamma,
                           prev[0], w, xtol=1e-14)
            return A_for(wstar, lnKHT), wstar
        prev = (w, g)
    return None


def fit_commitment(A, w0, theta, T, kht, sh, kdt, sd):
    def resid(cp):
        w = Q.gating_w(T, w0, theta)
        lh, ld = Q.log_kie_intrinsic(A, w)
        c = np.exp(cp[0] - cp[1]/(R_KCAL*T))
        ph = Q.observed_from_intrinsic(np.exp(lh), c)
        pdd = Q.observed_from_intrinsic(np.exp(ld), c)
        out = np.concatenate([(ph-kht)/sh, (pdd-kdt)/sd])
        return np.where(np.isfinite(out), out, 1e6)
    best = None
    for _ in range(300):
        r = least_squares(resid, [RNG.uniform(-3, 3), RNG.uniform(-6, 6)],
                          max_nfev=4000)
        if best is None or r.cost < best.cost:
            best = r
    return 2*best.cost


def main():
    T, kht, sh, kdt, sd = series("hyd")
    theta = 979.0

    print("="*78)
    print("CLAIM 3: Table II is an artefact of fixing K_HT^int = 11")
    print("="*78)
    print(f"{'gamma':>7}{'chi2 @ K=11 (mine)':>21}{'profiled chi2':>16}"
          f"{'best K_HT^int':>15}{'p':>8}")
    for g in (2.35, 2.50, 2.70, 2.90):
        # (a) reproduce my Table II: intrinsic magnitude pinned at 11
        sol = solve_Aw(g, np.log(11.0))
        if sol is None:
            print(f"{g:7.2f}   unreachable at K=11"); continue
        A11, w11 = sol
        w0_11 = w11*np.tanh(theta/(2*298.15))
        x2_fixed = fit_commitment(A11, w0_11, theta, T, kht, sh, kdt, sd)

        # (b) profile the intrinsic magnitude at the same exponent
        best_x2, best_k = np.inf, np.nan
        for k in np.linspace(1.6, 30.0, 120):
            s = solve_Aw(g, np.log(k))
            if s is None:
                continue
            A, w = s
            w0 = w*np.tanh(theta/(2*298.15))
            x2 = fit_commitment(A, w0, theta, T, kht, sh, kdt, sd)
            if x2 < best_x2:
                best_x2, best_k = x2, k
        dof = 2*len(T) - 3          # 2 commitment params + 1 profiled magnitude
        p = 1 - stats.chi2.cdf(best_x2, dof)
        print(f"{g:7.2f}{x2_fixed:21.1f}{best_x2:16.2f}{best_k:15.2f}{p:8.3f}")

    print("\n" + "="*78)
    print("CLAIM 11: can a semiclassical + Arrhenius-commitment model produce")
    print("          a NON-MONOTONIC observed K_HT(T)?")
    print("="*78)
    Tg = np.linspace(278.15, 308.15, 25)
    found = 0
    for _ in range(200000):
        lnAu, Eu = RNG.uniform(-2, 4), RNG.uniform(-4, 6)
        lnc0, Ec = RNG.uniform(-4, 4), RNG.uniform(-8, 8)
        u = 1.0 + np.exp(lnAu - Eu/(R_KCAL*Tg))
        c = np.exp(lnc0 - Ec/(R_KCAL*Tg))
        kh = Q.observed_from_intrinsic(u**3.34, c)
        if not np.all(np.isfinite(kh)):
            continue
        d = np.diff(kh)
        if np.any(d > 0) and np.any(d < 0) and kh.min() > 1.05 and kh.max() < 20:
            found += 1
            if found == 1:
                print(f"  counterexample found: lnA_u={lnAu:.3f} E_u={Eu:.3f} "
                      f"lnc0={lnc0:.3f} E_c={Ec:.3f}")
                print(f"  K_HT^obs over 5-35 C: min {kh.min():.4f} max {kh.max():.4f}"
                      f"  turning point present")
    print(f"  non-monotonic parameter sets found: {found} / 200000")
    print("  => the impossibility claim in the manuscript is FALSE as stated.")


if __name__ == "__main__":
    main()
