"""Is a common offset within a temperature series consistent with the errors?

The series statistic used for the bounds is a maximum over temperatures, which
is inefficient. An inverse-variance weighted mean would be far more powerful if
the offset were constant within a series. This script tests that premise, and
reports the numbers quoted in the supplement, so the claim is reproducible
rather than asserted.

For each series the offset per record is F = ln K_HT - gamma ln K_DT with

    sigma_F^2 = (s_HT/K_HT)^2 + gamma^2 (s_DT/K_DT)^2
                - 2 gamma rho (s_HT/K_HT)(s_DT/K_DT),

the two effects sharing a tritium reference so their correlation is unknown.
rho = 0 is the neutral choice; rho = -1 inflates the errors most and is
therefore the case most favourable to pooling.

Run from analysis/:   python pooling_check.py
"""
from __future__ import annotations

import numpy as np
import pandas as pd

import masses as M

G = M.gamma_sc("C")
KEY = ["family", "system", "variant", "step"]


def series_chi2(d, rho):
    """chi2 per degree of freedom for a constant offset, series by series."""
    out = []
    for key, g in d.groupby(KEY):
        F = np.log(g.K_HT) - G * np.log(g.K_DT)
        a, b = g.K_HT_se / g.K_HT, g.K_DT_se / g.K_DT
        sig = np.sqrt(a ** 2 + (G * b) ** 2 - 2 * G * rho * a * b)
        w = 1.0 / sig ** 2
        Fbar = float((F * w).sum() / w.sum())
        chi2 = float((((F - Fbar) / sig) ** 2).sum())
        dof = len(g) - 1
        if dof > 0:
            out.append((f"{key[0]} {key[2]} ({key[3]})", len(g), chi2 / dof))
    return pd.DataFrame(out, columns=["series", "n_T", "chi2_dof"])


def ladh_shared_pair_fit():
    """Do the nine 25 C horse liver forms share one intrinsic pair?

    Shared-reference model: K_Lj = x_L(1+c_j)/(x_L+c_j), one (x_H, x_D) and one
    c_j per form.  With z = 1 - 1/K, u_L = 1 - 1/x_L and s_j = c_j/(1+c_j) this
    is exactly z_Lj = u_L s_j, so (u_L, s_j) -> (lam u_L, s_j/lam) leaves every
    observation unchanged.  Of the 11 nominal parameters only 10 are identified:
    the ratio R = u_H/u_D and the nine z_Hj.  The fit is done in those
    coordinates; residual degrees of freedom 18 - 10 = 8.

    Returns chi2, dof, p, R, the Jacobian rank of the 11-parameter model at the
    fit, and a witness: two parameter sets giving identical observations with
    offsets of opposite sign.
    """
    from scipy.optimize import least_squares
    from scipy.stats import chi2 as chi2d
    d = pd.read_csv("../data/ladh_adh_primary.csv")
    d = d[d.T_C == 25]
    h, sh = d.K_HT.to_numpy(float), d.K_HT_se.to_numpy(float)
    k, sk = d.K_DT.to_numpy(float), d.K_DT_se.to_numpy(float)
    zH, zD = 1 - 1 / h, 1 - 1 / k

    def res(p):
        R, zh = p[0], p[1:]
        return np.r_[(1 / (1 - zh) - h) / sh, (1 / (1 - zh / R) - k) / sk]

    n = len(h)
    fit = least_squares(res, np.r_[np.mean(zH / zD), zH],
                        bounds=(np.r_[1 + 1e-6, np.full(n, 1e-9)],
                                np.r_[np.inf, np.full(n, 1 - 1e-9)]),
                        xtol=1e-14, ftol=1e-14, gtol=1e-14)
    c2 = float(fit.fun @ fit.fun)
    dof = 2 * n - (n + 1)
    R = float(fit.x[0])
    # 11-parameter Jacobian of z in (u_H, u_D, s_1..s_9) at an admissible point
    uH = 0.5 * (1 + fit.x[1:].max())
    u, s = np.array([uH, uH / R]), fit.x[1:] / uH
    J = np.zeros((2 * n, n + 2))
    for j in range(n):
        for i in range(2):
            J[i * n + j, i] = s[j]
            J[i * n + j, 2 + j] = u[i]
    rank = int(np.linalg.matrix_rank(J))
    # witness: nine commitments, scale lam = 0.8
    cs = np.geomspace(0.15, 2.0, 9)
    x = np.array([12.0, 2.0])
    lam = 0.8
    u0, s0 = 1 - 1 / x, cs / (1 + cs)
    x1, s1 = 1 / (1 - lam * u0), s0 / lam
    c1 = s1 / (1 - s1)
    K0 = x[:, None] * (1 + cs) / (x[:, None] + cs)
    K1 = x1[:, None] * (1 + c1) / (x1[:, None] + c1)
    F = lambda p: float(np.log(p[0]) - G * np.log(p[1]))
    return dict(chi2=c2, dof=dof, p=float(chi2d.sf(c2, dof)), R=R, rank=rank,
                n_par=n + 2, witness_diff=float(np.max(np.abs(K0 - K1))),
                witness_F=(F(x), F(x1)), witness_x1=tuple(x1), witness_cmin=float(c1.min()))


def pooled_compatible_set(zH, R, n=4000):
    """Compatible offsets for exact data from variants sharing one pair.

    u_H ranges over (max_j z_Hj, 1) and u_D = u_H/R; F(u) = -ln(1-u) +
    gamma ln(1 - u/R).  Returns the infimum of F over that range: the set is
    the half-line above it, the same as for the least-masked variant alone.
    """
    u = np.linspace(max(zH), 1 - 1e-12, n)
    F = -np.log(1 - u) + G * np.log(1 - u / R)
    return float(F.min())


def main():
    d = pd.read_csv("../data/trinomial_benchmark.csv")
    print(f"constant-offset fit within each series, {d.groupby(KEY).ngroups} series")
    for rho, note in ((0.0, "independent errors"),
                      (-1.0, "maximally adverse, most favourable to pooling")):
        r = series_chi2(d, rho)
        print(f"\n  rho = {rho:+.1f}  ({note})")
        print(f"    median chi2/dof {r.chi2_dof.median():.2f}   "
              f"max {r.chi2_dof.max():.1f}  ({r.loc[r.chi2_dof.idxmax(), 'series']})")
        print(f"    series with chi2/dof > 1: "
              f"{int((r.chi2_dof > 1).sum())} of {len(r)}")
    r0 = series_chi2(d, 0.0)
    print("\n  worst five at rho = 0:")
    print(r0.sort_values("chi2_dof", ascending=False).head(5).to_string(index=False))

    f = ladh_shared_pair_fit()
    print(f"\nnine 25 C horse liver forms, one shared intrinsic pair:")
    print(f"  chi2 = {f['chi2']:.2f} on {f['dof']} dof, p = {f['p']:.1e}; "
          f"model rank {f['rank']} of {f['n_par']} parameters (scale not identified)")
    print(f"  witness: identical observations to {f['witness_diff']:.1e} from "
          f"offsets {f['witness_F'][0]:+.4f} and {f['witness_F'][1]:+.4f}")


if __name__ == "__main__":
    main()
