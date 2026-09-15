"""Uncertainty-aware one-sided bounds on the mass-scaling offset.

The deterministic bounds are transformations of point estimates.  They are
reported here with sampling uncertainty, which changes the empirical
conclusion: the apparent near-miss of the closest series is well inside one
standard error and is not an inferential statement.

Three features of the design matter.

  * The H/T and D/T effects share a tritium reference, so their errors are
    correlated and the correlation is not reported in the sources; we bracket
    it over rho in {-1, 0, 0.5, 0.9}; rho = -1 is the maximally adverse case.

  * The target of a series is psi = max_T F_min(K(T)), the sharpest exclusion
    its temperatures afford.  Its lower bound is max_T L_T, where L_T is the
    one-sided lower bound for temperature T at level alpha/n_T.  That is valid:
    max_T L_T > psi requires L_T > F_min(K(T)) for some T, so by the union bound
    P(max_T L_T > psi) <= n_T (alpha/n_T) = alpha.

    Until 2026-09-14 the bound was instead the alpha-quantile of the maximum of
    the RESAMPLED endpoints.  That is not a confidence bound for a maximum: the
    resampled maximum is biased upward, and simulated coverage for the ecDHFR
    W133F series as recorded gave P(bound > psi) = 0.20 at a nominal 0.05, and
    0.80 for five temperatures with equal endpoints.  The Bonferroni maximum gave
    0.02-0.04 in the same simulations.  The old bounds were therefore too high;
    no series crossed F0 under either, so no verdict changes.

  * A replicate whose draw violates admissibility (x_H > x_D > 1) carries no
    information about the offset and is retained at -inf rather than dropped.
    Dropping such replicates removes mass from the lower tail and biases a
    lower confidence bound upward, which is the anti-conservative direction.

Quantiles use the order statistic itself (method="lower") rather than linear
interpolation, so the presence of -inf replicates cannot propagate a NaN.  The
Monte Carlo uncertainty of every reported bound is estimated directly, by
running independent replications, so the quoted precision is justified rather
than assumed.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

import masses as M
from partial_id import F_min_vec, F_min_exact

GSC = M.gamma_sc("C")
F0 = M.offset_F0("C")
N_DRAW = 1_000_000        # replicates per replication (tail quantiles at alpha/n_T)
N_REP = 5                # independent replications, for the Monte Carlo s.d.
ALPHA = 0.05
BASE_SEED = 20260812
OUT = []


def say(s=""):
    print(s)
    OUT.append(s)


def draw(kht, sh, kdt, sd, rho, n, rng):
    """Correlated draws of the observed pair, lognormal on the log scale."""
    m = np.array([np.log(kht), np.log(kdt)])
    s = np.array([sh / kht, sd / kdt])          # delta-method log-scale s.d.
    cov = np.array([[s[0] ** 2, rho * s[0] * s[1]],
                    [rho * s[0] * s[1], s[1] ** 2]])
    z = rng.multivariate_normal(m, cov, size=n)
    return np.exp(z[:, 0]), np.exp(z[:, 1])


def one_replication(g, rho, n, rng):
    """One simultaneous-in-T replication.

    Returns (lower confidence bound, sampling s.d. of the statistic,
    fraction of replicates carrying no information).
    """
    reps = []
    for _, r in g.iterrows():
        a, b = draw(r.K_HT, r.K_HT_se, r.K_DT, r.K_DT_se, rho, n, rng)
        reps.append(F_min_vec(a, b))
    reps = np.vstack(reps)                        # temperatures x replicates
    stat = np.max(reps, axis=0)                   # resampled series maximum
    dead = float(np.mean(~np.isfinite(stat)))
    # Bonferroni over temperatures: per-temperature bounds at alpha/n_T, then
    # their maximum (see module docstring for why not a quantile of stat)
    lcb = float(np.max(np.quantile(reps, ALPHA / reps.shape[0], axis=1,
                                   method="lower")))
    above = float(np.mean(stat > F0))            # resampling frequency of clearing F0
    finite = stat[np.isfinite(stat)]
    return lcb, float(finite.std()), dead, above


def series_bound(g, rho, n=N_DRAW, reps=N_REP, seed=BASE_SEED):
    """Bound averaged over independent replications, with its Monte Carlo s.d."""
    vals, sds, deads, abv = [], [], [], []
    for k in range(reps):
        rng = np.random.default_rng(seed + 1009 * k)
        l, sd, dd, ab = one_replication(g, rho, n, rng)
        vals.append(l); sds.append(sd); deads.append(dd); abv.append(ab)
    return (float(np.mean(vals)), float(np.std(vals, ddof=1)),
            float(np.mean(sds)), float(np.mean(deads)), float(np.mean(abv)))


def assignment_sensitivity(rho=-1.0):
    """How far the unpaired 10 C ecDHFR record can move the two series bounds.

    Wang et al. distinguish light and heavy rows only by colour, so the record
    is assigned by interpolation.  Returns (largest move of either bound over
    the three assignments -- light, heavy, dropped -- and the smallest distance
    of any of those bounds below F0).
    """
    d = pd.read_csv("../data/trinomial_benchmark.csv")
    eco = d.system.str.contains("Escherichia", na=False) & \
        d.variant.astype(str).isin(["light enzyme", "heavy enzyme"])
    r10 = eco & d.T_C.eq(10.0)
    out = {}
    for tag in ("light enzyme", "heavy enzyme", None):
        dd = d[eco].copy()
        if tag is None:
            dd = dd[~r10[eco]]
        else:
            dd.loc[r10[eco], "variant"] = tag
        out[tag] = {v: series_bound(g, rho)[0] for v, g in dd.groupby("variant")}
    moves = [abs(out[a][v] - out[b][v]) for v in ("light enzyme", "heavy enzyme")
             for a in out for b in out]
    gap = min(F0 - x for o in out.values() for x in o.values())
    return max(moves), gap


def coverage_check(family="ecDHFR", variant="W133F", rho=-1.0, reps=20000,
                   n=4000, seed=20260914):
    """Simulated probability that each procedure's bound exceeds the true maximum.

    Truth is the recorded series; observed data are drawn from it, and both the
    quantile-of-the-resampled-maximum bound (used until 2026-09-14) and the
    Bonferroni maximum are computed from each simulated data set.  A valid 95%
    lower bound exceeds the true maximum with probability at most 0.05.
    """
    d = pd.read_csv("../data/trinomial_benchmark.csv")
    g = d[(d.family == family) & (d.variant == variant)]
    kh, sh, kd, sd = (g[c].to_numpy(float) for c in ("K_HT", "K_HT_se", "K_DT", "K_DT_se"))
    psi = float(np.max(F_min_vec(kh, kd)))
    rng = np.random.default_rng(seed)
    nT = len(kh)
    fail_old = fail_new = 0
    for _ in range(reps):
        obs = [draw(kh[j], sh[j], kd[j], sd[j], rho, 1, rng) for j in range(nT)]
        oh = np.array([o[0][0] for o in obs]); od = np.array([o[1][0] for o in obs])
        E = np.vstack([F_min_vec(*draw(oh[j], sh[j], od[j], sd[j], rho, n, rng))
                       for j in range(nT)])
        fail_old += np.quantile(E.max(axis=0), ALPHA, method="lower") > psi
        fail_new += np.max(np.quantile(E, ALPHA / nT, axis=1, method="lower")) > psi
    return fail_old / reps, fail_new / reps


def plugin_coverage(kht, sh, kdt, sd, alpha, n=1_000_000, seed=20260915, plugin=True):
    """Repeated-sampling failure rate of the single-record rho = -1 lower bound.

    Truth: log K ~ N(log K_true, s/K_true), perfectly anti-correlated (the model
    this module draws from).  The bound as implemented re-expresses the reported
    absolute errors at the OBSERVED effects (s/K_obs) and takes the alpha quantile
    of the endpoint over draws centred there.  At rho = -1 the draws lie on a line
    along which the endpoint is monotone, so that quantile is the endpoint at the
    z_alpha point and no inner Monte Carlo is needed.

    Returns P(bound > true endpoint); a valid bound has at most alpha.  With
    plugin=False the log-scale standard deviations are held at their true values,
    and the rate is alpha exactly.  The plug-in bound is only approximate: until
    2026-09-15 the Methods called it exact at rho = -1, which holds only for known
    log-scale standard deviations (4.8% for the yeast errors, 7.6% for 2.0 +- 0.2,
    1.1 +- 0.001, at a nominal 5%).
    """
    from scipy.stats import norm
    rng = np.random.default_rng(seed)
    truth = F_min_exact(kht, kdt)[0]
    s_h, s_d = sh / kht, sd / kdt
    z = rng.standard_normal(n)
    lh, ld = np.log(kht) + s_h * z, np.log(kdt) - s_d * z
    if plugin:
        u_h, u_d = sh / np.exp(lh), sd / np.exp(ld)
    else:
        u_h, u_d = s_h, s_d
    za = norm.ppf(alpha)
    bh, bd = np.exp(lh + u_h * za), np.exp(ld - u_d * za)
    ok = (bh > bd) & (bd > 1)
    lcb = np.where(ok, F_min_vec(np.where(ok, bh, 2.0), np.where(ok, bd, 1.5)), -np.inf)
    return float(np.mean(lcb > truth))


def main():
    d = pd.read_csv("../data/trinomial_benchmark.csv")
    say("=" * 96)
    say("UNCERTAINTY-AWARE ONE-SIDED BOUNDS ON THE OFFSET")
    say("=" * 96)
    say(f"gamma_SC = {GSC:.5f} (C-H reduced masses),  F0 = {F0:+.6f}")
    say(f"{N_DRAW} replicates x {N_REP} independent replications per cell;")
    say(f"one-sided {100*(1-ALPHA):.0f}% bound on the series maximum: per-temperature")
    say(f"bounds at {ALPHA}/n_T (Bonferroni), then their maximum.")
    say("")
    say("A series refutes the ground-channel gated model only if its lower")
    say("confidence bound exceeds the mechanism envelope at the relevant scale.")
    say(f"Every intrinsic H/T effect here exceeds "
        f"{np.sqrt(M.mu_ratios('C')[2]):.4f}, so that envelope is the flat "
        f"value F0 = {F0:+.6f}.")
    say("")
    say("The H/T and D/T effects share a tritium reference. The correlation is")
    say("not reported, so results are bracketed over rho in {-1, 0, 0.5, 0.9}.")
    say("")

    rows = []
    for rho in (-1.0, 0.0, 0.5, 0.9):
        say(f"--- rho = {rho} " + "-" * 76)
        say(f"{'series':42s}{'point':>9}{'95% LCB':>10}{'MC sd':>8}"
            f"{'sd(stat)':>10}{'P(>F0)':>9}{'refutes?':>10}")
        for key, g in d.groupby(["family", "variant", "step"]):
            fam, var, step = key
            pt = max(F_min_exact(r.K_HT, r.K_DT)[0] for _, r in g.iterrows())
            lcb, mcsd, sdstat, dead, above = series_bound(g, rho)
            say(f"{fam+' '+var+' ('+step+')':42s}{pt:9.4f}{lcb:10.4f}"
                f"{mcsd:8.4f}{sdstat:10.4f}{above:9.3f}"
                f"{'YES' if lcb > F0 else 'no':>10}")
            rows.append(dict(rho=rho, family=fam, variant=var, step=step,
                             point=pt, lcb=lcb, lcb_mc_sd=mcsd,
                             sd_stat=sdstat, p_above_F0=above,
                             frac_uninformative=dead,
                             refutes=bool(lcb > F0)))
        say("")

    res = pd.DataFrame(rows)
    res.to_csv("../results/bounds_uncertainty.csv", index=False)

    say("=" * 96)
    say("SUMMARY")
    say("=" * 96)
    for rho, sub in res.groupby("rho"):
        best = sub.loc[sub.lcb.idxmax()]
        say(f"rho = {rho}: series refuting the model = {int(sub.refutes.sum())}"
            f" of {len(sub)};  best lower bound {best.lcb:+.4f}"
            f" ({best.family} {best.variant})")
    say("")
    say(f"Largest Monte Carlo s.d. over all cells: {res.lcb_mc_sd.max():.5f}")
    say(f"Largest uninformative-replicate fraction: "
        f"{res.frac_uninformative.max():.2e}")
    mc = res.lcb_mc_sd.max()
    secure = 3 if mc < 5e-4 else 2
    say(f"Bounds are therefore secure to {secure} decimals everywhere "
        f"(Monte Carlo s.d. {mc:.1e}).")
    say("")
    w = res[res.rho == 0.0]
    b = w.loc[w.point.idxmax()]
    say(f"Closest series on the point estimate: {b.family} {b.variant}.")
    say(f"  point bound      {b.point:+.4f}  (short of F0 by {F0-b.point:.4f})")
    say(f"  95% lower bound  {b.lcb:+.4f}  (short of F0 by {F0-b.lcb:.4f})")
    say(f"  s.d. of the series statistic {b.sd_stat:.4f}")
    say(f"  resampling frequency of clearing F0: {b.p_above_F0:.3f}")
    say("")
    say("The resampled series maximum is right skewed, which is why its quantile")
    say("is not used as the bound; the bound is the Bonferroni maximum above.")
    say("")
    bl = w.loc[w.lcb.idxmax()]
    say(f"Ranking by confidence bound instead promotes {bl.family} {bl.variant}"
        f" (LCB {bl.lcb:+.4f}),")
    say("which has more temperatures and tighter errors.")

    old, new = coverage_check()
    say("")
    say("Coverage, ecDHFR W133F as recorded taken as truth, rho = -1:")
    say(f"  P(bound > true maximum): quantile of resampled maximum {old:.3f}, "
        f"Bonferroni maximum {new:.3f} (nominal 0.05)")
    mv, gap = assignment_sensitivity()
    say("")
    say("The unpaired 10 C ecDHFR record (assigned by interpolation): over light,")
    say(f"heavy and dropped, a bound moves by at most {mv:.4f} at rho = -1 and every")
    say(f"bound stays at least {gap:.4f} below F0.")

    with open("../results/bounds_uncertainty.txt", "w") as f:
        f.write("\n".join(OUT) + "\n")
    print("\n[written] ../results/bounds_uncertainty.{csv,txt}")


if __name__ == "__main__":
    main()
