"""Out-of-sample stability of the per-system verdicts.

The blind test one would most like to run is: classify each system as
informative or not from independently published kinetics alone, freeze the
prediction, then reveal the isotope data.  The curated corpus does not carry a
per-system kinetic characterization, so that test is not available here and is
not attempted; what it would require is stated in the supplement.

What the data do support is an out-of-sample check on the verdicts themselves.
A temperature series contributes one identified set per temperature, and the
row plotted is the largest endpoint over the series.  A verdict that depends on
a single temperature is an artifact of which temperatures happened to be
measured.  Leaving each temperature out in turn and recomputing answers that.

Run from analysis/:   python holdout.py
"""
from __future__ import annotations

import numpy as np
import pandas as pd

import masses as M
from partial_id import F_min_exact

G = M.gamma_sc("C")
F0 = M.offset_F0("C")


def series_endpoint(rows):
    """Largest endpoint over a set of (K_HT, K_DT) records, as the paper plots."""
    vals = [F_min_exact(kh, kd)[0] for kh, kd in rows if kh > kd > 1]
    return max(vals) if vals else np.nan


def leave_one_out():
    """Per series: the full verdict, and whether dropping one temperature flips it."""
    tb = pd.read_csv("../data/trinomial_benchmark.csv")
    out = []
    for (sysname, variant, step), g in tb.groupby(["system", "variant", "step"]):
        rows = [(r.K_HT, r.K_DT) for _, r in g.iterrows()]
        rows = [r for r in rows if r[0] > r[1] > 1]
        if len(rows) < 2:
            continue
        full = series_endpoint(rows)
        verdict = full > F0
        flips, held = 0, []
        for i in range(len(rows)):
            sub = rows[:i] + rows[i + 1:]
            e = series_endpoint(sub)
            held.append(e)
            if (e > F0) != verdict:
                flips += 1
        out.append(dict(system=f"{sysname} {variant} ({step})", n=len(rows),
                        full=full, verdict=verdict, flips=flips,
                        spread=float(np.nanmax(held) - np.nanmin(held))))
    return pd.DataFrame(out)


def confidence_deletion(rho=-1.0, only=None):
    """Deletion sensitivity of the CONFIDENCE decision, not just the point.

    For each series whose full-data bound clears F0, drop each temperature in
    turn and recompute the one-sided series bound with the paper's procedure
    (bounds_uncertainty.series_bound: per-temperature bounds at 0.05/n_T over the
    remaining temperatures, maximum taken).  Until 2026-10-15 the supplement
    said the two F93W methanol series "exclude under every leave-one-out", but
    leave_one_out() above compares point endpoints only; an external audit
    showed that both lose the confidence-level exclusion without their 3 C
    record.  Returns one row per (series, omitted temperature).
    """
    import bounds_uncertainty as BU
    tb = pd.read_csv("../data/trinomial_benchmark.csv")
    out = []
    for (sysname, variant, step), g in tb.groupby(["system", "variant", "step"]):
        if only is not None and variant not in only:
            continue
        full = BU.series_bound(g, rho)[0]
        if not full > F0:
            continue
        for i, (_, r) in enumerate(g.iterrows()):
            sub = g.drop(index=r.name)
            pts = [F_min_exact(a, b)[0] for a, b in zip(sub.K_HT, sub.K_DT)
                   if a > b > 1]
            b = BU.series_bound(sub, rho)[0]
            out.append(dict(variant=variant, omitted_T=float(r.T_C),
                            full_bound=full, point=max(pts), bound=b,
                            excludes=bool(b > F0)))
    return pd.DataFrame(out)


def structural_prediction_check():
    """Consistency check: a positive observed offset implies an open set.

    Theorem 1 implies that any record whose observed offset is positive has
    L_H >= gamma, so its identified set is the open half-line.  This cannot fail
    on any data: ln K_HT > gamma ln K_DT gives K_HT - 1 > K_DT^gamma - 1 >=
    gamma (K_DT - 1) by Bernoulli's inequality.  It checks the implementation,
    not the theory.  (Described until 2026-09-14 as an out-of-sample test.)
    """
    import corpus
    frames = [corpus.series(), corpus.single_condition()]
    n = viol = 0
    for f in frames:
        for _, r in f.iterrows():
            kh, kd = r.K_HT, r.K_DT
            if not (kh > kd > 1):
                continue
            n += 1
            F_obs = np.log(kh) - G * np.log(kd)
            _, _, _, closed = F_min_exact(kh, kd)
            if F_obs > 0 and closed:
                viol += 1
    return n, viol


def main():
    print("A. Out-of-sample stability of the per-series verdicts")
    df = leave_one_out()
    print(f"   temperature series tested: {len(df)}")
    print(f"   verdicts that flip when any single temperature is dropped: "
          f"{int((df.flips > 0).sum())}")
    print(f"   informative under the full series: {int(df.verdict.sum())}")
    print(f"\n   {'series':52s} {'n':>2s} {'endpoint':>9s} {'flips':>5s} {'spread':>7s}")
    for _, r in df.sort_values("full", ascending=False).iterrows():
        print(f"   {r.system[:52]:52s} {r.n:2d} {r.full:+9.4f} {r.flips:5d} {r.spread:7.4f}")

    print("\nA2. Deletion sensitivity of the confidence-level exclusions (rho = -1)")
    cd = confidence_deletion()
    for _, r in cd.iterrows():
        print(f"   {r.variant:22s} omit {r.omitted_T:+6.0f} C: point {r.point:+.4f}, "
              f"bound {r.bound:+.4f} ({'excludes' if r.excludes else 'does not exclude'}; "
              f"full bound {r.full_bound:+.4f})")
    cd.to_csv("../results/deletion_sensitivity.csv", index=False)

    n, viol = structural_prediction_check()
    print(f"\nB. Structural prediction tested out of sample")
    print(f"   Theorem 1 forecasts: observed offset > 0  =>  identified set open")
    print(f"   records checked {n}, violations {viol}")
    return df


if __name__ == "__main__":
    main()
