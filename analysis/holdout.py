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


def structural_prediction_check():
    """A prediction made before the data: exclusion at stake implies an open set.

    Theorem 1 is proved from the map alone.  It forecasts that any record whose
    observed offset is positive has L_H >= gamma, so its identified set is the
    open half-line.  Checking that against every record is a genuine
    out-of-sample test of a structural claim, since no record entered the proof.
    """
    frames = [pd.read_csv(f) for f in ("../data/trinomial_benchmark.csv",
                                       "../data/ladh_adh_primary.csv",
                                       "../data/bsao_grant1989.csv")]
    ya = pd.read_csv("../data/cha1989_yadh.csv")
    frames.append(ya[ya.note.str.contains("average")])
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

    n, viol = structural_prediction_check()
    print(f"\nB. Structural prediction tested out of sample")
    print(f"   Theorem 1 forecasts: observed offset > 0  =>  identified set open")
    print(f"   records checked {n}, violations {viol}")
    return df


if __name__ == "__main__":
    main()
