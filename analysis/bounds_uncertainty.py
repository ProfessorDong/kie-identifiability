"""Uncertainty-aware one-sided bounds on the mass-scaling offset.

The earlier bounds were deterministic transformations of point estimates.  They
are reported here with sampling uncertainty, which changes the empirical
conclusion: the apparent near-miss of the closest series is well inside one
standard error and is not an inferential statement.

Two features of the design matter.  The H/T and D/T effects share a tritium
reference, so their errors are correlated and the correlation is not reported
in the sources; we therefore bracket it.  And the per-series bound is a maximum
over temperatures, which is a selection, so a simultaneous procedure is used
rather than a pointwise one.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

import masses as M
from partial_id import F_min_exact

GSC = M.gamma_sc("C")
F0 = M.offset_F0("C")
RNG = np.random.default_rng(20260812)
OUT = []


def say(s=""):
    print(s)
    OUT.append(s)


def draw(kht, sh, kdt, sd, rho, n):
    """Correlated draws of the observed pair on the log scale."""
    m = np.array([np.log(kht), np.log(kdt)])
    s = np.array([sh / kht, sd / kdt])
    cov = np.array([[s[0] ** 2, rho * s[0] * s[1]],
                    [rho * s[0] * s[1], s[1] ** 2]])
    z = RNG.multivariate_normal(m, cov, size=n)
    return np.exp(z[:, 0]), np.exp(z[:, 1])


def series_bound(g, rho, n=20000, alpha=0.05):
    """One-sided lower confidence bound on max_T F_min(T), simultaneous in T."""
    reps = []
    for _, r in g.iterrows():
        a, b = draw(r.K_HT, r.K_HT_se, r.K_DT, r.K_DT_se, rho, n)
        f = np.empty(n)
        for i in range(n):
            if a[i] <= b[i] or b[i] <= 1.0:
                f[i] = -np.inf
            else:
                f[i] = F_min_exact(a[i], b[i])[0]
        reps.append(f)
    # the series statistic is the maximum over temperatures; take its
    # alpha-quantile across the joint replicates (temperatures independent)
    stat = np.max(np.vstack(reps), axis=0)
    stat = stat[np.isfinite(stat)]
    return np.quantile(stat, alpha), stat.mean()


def main():
    d = pd.read_csv("../data/trinomial_benchmark.csv")
    say("=" * 86)
    say("UNCERTAINTY-AWARE ONE-SIDED BOUNDS ON THE OFFSET")
    say("=" * 86)
    say(f"gamma_SC = {GSC:.5f} (C-H reduced masses)")
    say(f"a series refutes the ground-channel gated model only if its lower")
    say(f"confidence bound exceeds the envelope at the relevant scale; for")
    say(f"intrinsic effects above {np.sqrt(M.mu_ratios('C')[2]):.4f} that envelope is "
        f"F0 = {F0:+.6f}.")
    say("")
    say("The H/T and D/T effects share a tritium reference. The correlation is")
    say("not reported, so results are bracketed over rho in {0, 0.5, 0.9}.")
    say("")

    rows = []
    for rho in (0.0, 0.5, 0.9):
        say(f"--- rho = {rho} " + "-" * 62)
        say(f"{'series':44s}{'point':>9}{'95% LCB':>10}{'SE':>8}{'refutes?':>10}")
        for key, g in d.groupby(["family", "variant", "step"]):
            fam, var, step = key
            pt = max(F_min_exact(r.K_HT, r.K_DT)[0] for _, r in g.iterrows())
            lcb, mean = series_bound(g, rho, n=6000)
            se = (pt - lcb) / 1.645
            say(f"{fam+' '+var+' ('+step+')':44s}{pt:9.4f}{lcb:10.4f}{se:8.4f}"
                f"{'YES' if lcb > F0 else 'no':>10}")
            rows.append(dict(rho=rho, family=fam, variant=var, step=step,
                             point=pt, lcb=lcb, se=se, refutes=bool(lcb > F0)))
        say("")

    res = pd.DataFrame(rows)
    res.to_csv("../results/bounds_uncertainty.csv", index=False)

    say("=" * 86)
    say("SUMMARY")
    say("=" * 86)
    for rho, sub in res.groupby("rho"):
        best = sub.loc[sub.lcb.idxmax()]
        say(f"rho = {rho}: series refuting the model = {int(sub.refutes.sum())}"
            f" of {len(sub)};  best lower bound {best.lcb:+.4f}"
            f" ({best.family} {best.variant})")
    say("")
    w = res[res.rho == 0.0]
    b = w.loc[w.point.idxmax()]
    say(f"The closest series on the point estimate is {b.family} {b.variant},")
    say(f"point bound {b.point:+.4f}, short of F0 by {F0-b.point:.4f}.")
    say(f"Its 95% lower confidence bound is {b.lcb:+.4f}, short by "
        f"{F0-b.lcb:.4f}, i.e. {abs(F0-b.point)/b.se:.2f} standard errors.")
    say("The apparent near-miss is not an inferential statement.")

    with open("../results/bounds_uncertainty.txt", "w") as f:
        f.write("\n".join(OUT) + "\n")
    print("\n[written] ../results/bounds_uncertainty.{csv,txt}")


if __name__ == "__main__":
    main()
