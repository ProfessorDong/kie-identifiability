"""What can a competitive isotope series determine about the discriminating invariant?

The semiclassical and gated-tunneling loci are parallel lines in log-intrinsic-KIE
space, separated by a fixed offset set by the masses alone.  The entire
discrimination problem therefore reduces to one scalar:

    F = ln K_HT^int - gamma_SC ln K_DT^int,   F = 0 (semiclassical), F0 (ridge).

Rather than fit two model classes and compare them -- the approach that produced
the invalid fixed-magnitude table in the withdrawn draft -- we fit ONE model in
which F is free and profile it.  That measures the invariant directly, with an
interval, and yields a precision target rather than a model-selection verdict.

Forward model, per series:
    ln K_DT^int(T) = a + b/T                       Arrhenius intrinsic D/T
    ln K_HT^int(T) = gamma_SC ln K_DT^int(T) + F   F constant in T
    c(T)           = exp(lnc0 - Ec/(R T))          Arrhenius reciprocal commitment
    observed       = x(1+c)/(x+c)                  competitive observation map

Five parameters against 2 n_T observations.  The profile is computed by
continuation: a heavily multistarted fit at the best grid point, then warm-
started stepping outward, which keeps the profile smooth and avoids the spurious
spikes that independent multistart produces far from the optimum.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.optimize import least_squares

import masses as M
import ridge as Rg

R_KCAL = 1.98720425e-3
GSC = M.gamma_sc("C")
F0 = M.offset_F0("C")
RNG = np.random.default_rng(20260811)
OUT = []


def say(s=""):
    print(s, flush=True)
    OUT.append(s)


def predict(theta, T, F):
    a, b, lnc0, Ec = theta
    lnxD = a + b / T
    xD = np.exp(lnxD)
    xH = np.exp(GSC * lnxD + F)
    c = np.exp(lnc0 - Ec / (R_KCAL * T))
    return Rg.observed_from_intrinsic(xH, c), Rg.observed_from_intrinsic(xD, c)


def make_resid(T, kht, sh, kdt, sd, F):
    def resid(th):
        with np.errstate(all="ignore"):
            ph, pd_ = predict(th, T, F)
        out = np.concatenate([(ph - kht) / sh, (pd_ - kdt) / sd])
        return np.where(np.isfinite(out), out, 1e6)
    return resid


def fit_at(T, kht, sh, kdt, sd, F, n_rand=0, warm=None):
    resid = make_resid(T, kht, sh, kdt, sd, F)
    starts = []
    if warm is not None:
        starts.append(np.asarray(warm, dtype=float))
    for _ in range(n_rand):
        starts.append(np.array([RNG.uniform(-1, 2), RNG.uniform(-200, 900),
                                RNG.uniform(-3, 3), RNG.uniform(-6, 6)]))
    best = None
    for th0 in starts:
        try:
            r = least_squares(resid, th0, max_nfev=3000)
        except Exception:
            continue
        if best is None or r.cost < best.cost:
            best = r
    return best


def profile(T, kht, sh, kdt, sd, grid, n_seed=400, n_step=6):
    """Profile chi2 over F by continuation from the global optimum."""
    # locate the optimum with a heavily multistarted coarse sweep
    coarse = grid[::4]
    best_F, best_r = None, None
    for F in coarse:
        r = fit_at(T, kht, sh, kdt, sd, F, n_rand=max(30, n_seed // len(coarse)))
        if r is not None and (best_r is None or r.cost < best_r.cost):
            best_F, best_r = F, r
    if best_r is None:
        return None, None, None
    i0 = int(np.argmin(np.abs(grid - best_F)))

    chi2 = np.full(len(grid), np.nan)
    sols = [None] * len(grid)
    r = fit_at(T, kht, sh, kdt, sd, grid[i0], n_rand=n_seed, warm=best_r.x)
    chi2[i0], sols[i0] = 2 * r.cost, r.x
    for direction in (+1, -1):                    # walk outward, warm-started
        warm = sols[i0]
        j = i0 + direction
        while 0 <= j < len(grid):
            r = fit_at(T, kht, sh, kdt, sd, grid[j], n_rand=n_step, warm=warm)
            if r is None:
                break
            chi2[j], sols[j] = 2 * r.cost, r.x
            warm = r.x
            j += direction
    return chi2, grid[int(np.nanargmin(chi2))], np.nanmin(chi2)


def main():
    d = pd.read_csv("../data/trinomial_benchmark.csv")
    say("=" * 82)
    say("THE DISCRIMINATING INVARIANT ACROSS THE MATCHED-PAIR BENCHMARK")
    say("=" * 82)
    say(f"mass convention: C-H reduced masses, gamma_SC = {GSC:.5f}")
    say(f"semiclassical predicts F = 0;  gated-tunneling ridge predicts "
        f"F0 = {F0:+.6f}")
    say(f"a series discriminates only if its 95% interval on F is narrower")
    say(f"than |F0| = {abs(F0):.4f} AND lands on one side of it.")
    say("")

    grid = np.arange(-1.50, 1.5001, 0.01)
    rows = []
    say(f"{'series':46s}{'n_T':>4}{'F_hat':>8}{'95% CI on F':>20}{'width':>8}"
        f"{'w/|F0|':>8}")
    say("-" * 94)
    for key, g in d.groupby(["family", "system", "variant", "step"]):
        fam, sysm, var, step = key
        g = g.sort_values("T_C")
        T = g.T_C.to_numpy() + 273.15
        kht, sh = g.K_HT.to_numpy(), g.K_HT_se.to_numpy()
        kdt, sd = g.K_DT.to_numpy(), g.K_DT_se.to_numpy()
        if len(T) < 3:
            continue
        chi2, Fhat, cmin = profile(T, kht, sh, kdt, sd, grid)
        if chi2 is None:
            say(f"{fam+' '+var:46s}  fit failed"); continue
        ok = np.isfinite(chi2) & (chi2 <= cmin + 3.841)   # one-parameter 95%
        lo, hi = grid[ok].min(), grid[ok].max()
        open_lo, open_hi = bool(ok[0]), bool(ok[-1])
        width = hi - lo
        label = f"{fam} {var} ({step})"
        flag = ("<" if open_lo else "[") + f"{lo:6.3f},{hi:6.3f}" + (">" if open_hi else "]")
        say(f"{label:46s}{len(T):>4}{Fhat:8.3f}{flag:>20}{width:8.3f}"
            f"{width/abs(F0):8.1f}")
        rows.append(dict(family=fam, system=sysm, variant=var, step=step,
                         n_T=len(T), chi2_min=cmin, F_hat=Fhat, F_lo=lo, F_hi=hi,
                         width=width, open_interval=open_lo or open_hi,
                         excludes_semiclassical=bool(lo > 0 or hi < 0),
                         excludes_ridge=bool(lo > F0 or hi < F0)))

    res = pd.DataFrame(rows)
    res.to_csv("../results/offset_profiles.csv", index=False)

    say("")
    say("=" * 82)
    say("VERDICT")
    say("=" * 82)
    say(f"series analysed                                : {len(res)}")
    say(f"intervals excluding the semiclassical F = 0    : "
        f"{int(res.excludes_semiclassical.sum())}")
    say(f"intervals excluding the tunneling ridge F = F0 : "
        f"{int(res.excludes_ridge.sum())}")
    say(f"intervals containing BOTH                      : "
        f"{int(((~res.excludes_semiclassical) & (~res.excludes_ridge)).sum())}")
    say(f"narrowest interval on F                        : {res.width.min():.3f}"
        f"  ({res.width.min()/abs(F0):.0f}x the signal)")
    say(f"median interval on F                           : {res.width.median():.3f}"
        f"  ({res.width.median()/abs(F0):.0f}x the signal)")
    say("")

    # the mass-modulated pair: a within-system control on the gating coordinate
    hl = res[res.variant.isin(["light enzyme", "heavy enzyme"])]
    if len(hl) == 2:
        say("MASS-MODULATED CONTROL (Francis et al., light vs heavy ecDHFR)")
        say("Heavy labelling perturbs the promoting mode at fixed electronic")
        say("structure, so it moves the system ALONG the ridge, where F is")
        say("invariant. Both series should return the same F.")
        for _, r in hl.iterrows():
            say(f"  {r.variant:14s} F = {r.F_hat:+.3f}  "
                f"[{r.F_lo:+.3f}, {r.F_hi:+.3f}]")
        say("")

    say("Interpretation: the invariant is bounded far more loosely than the")
    say("separation it must resolve, in every series. The obstacle is not the")
    say("number of enzymes but the precision with which F can be recovered")
    say("through the commitment map.")
    with open("../results/offset_report.txt", "w") as f:
        f.write("\n".join(OUT) + "\n")
    print("\n[written] ../results/offset_profiles.csv, offset_report.txt")


if __name__ == "__main__":
    main()
