"""Synthetic validation: the whole procedure run forward and back on known truths.

An external audit (2026-10-21) asked for one end-to-end demonstration under a
specified kinetic model: start from a known intrinsic pair and commitment,
simulate the assay, compute the compatible set and its bound exactly as the
paper does, and show what an independent commitment measurement adds.  Nothing
here is experimental confirmation; every number is a property of the procedure.

Truths.  A shared-reference scheme, Eq. (1), with intrinsic D/T effect 1.8 and
  * "above":  intrinsic offset +0.25, a pair the ground-channel model cannot
    produce (it lies above both F0 and the semiclassical reference);
  * "model":  intrinsic offset F0 - 0.005, a pair inside the model's range
    (at x_H > 1.2689 the model fills offsets up to, not including, F0).
Each is observed at reciprocal commitment c = 2 (strong masking), 8 and 32
(weak), with independent lognormal errors of 1% (close to the published median,
0.84% and 1.05%) or 0.3% on each effect.

The analyst sees only the two effects and their errors and reports the paper's
per-unit bound: the exact 95% quantile of the endpoint at the adverse
correlation rho = -1, with plug-in errors (bounds_uncertainty.exact_bound).
Reported per truth and commitment: the endpoint of the noise-free observation;
the fraction of simulated experiments whose bound exceeds F0 (a false exclusion
for the "model" truth, power for the "above" truth); and the fraction whose
bound lies at or below the true intrinsic offset (coverage).

Calibration.  The analyst then also measures c, with 5% relative error, and
inverts Eq. (1) at that value: x = K c/(1 + c - K).  Reported: the 5th, 50th
and 95th percentiles of the recovered offset over the simulated experiments.

Run from analysis/:   python synthetic_validation.py
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from scipy import stats

import masses as M
from partial_id import F_min_vec, F_min_exact
import bounds_uncertainty as BU

G = M.gamma_sc("C")
F0 = M.offset_F0("C")
XD = 1.8
TRUTHS = (("above", 0.25), ("model", F0 - 0.005))
C_GRID = (2.0, 8.0, 32.0)
PREC = (0.01, 0.003)
S_C = 0.05
N = 20_000
SEED = 20261021
ALPHA = 0.05


def observe(x, c):
    """Eq. (1): the observed effect of intrinsic effect x at reciprocal commitment c."""
    return x * (1.0 + c) / (x + c)


def bound_vec(kh, sh, kd, sd, alpha=ALPHA):
    """bounds_uncertainty.exact_bound on arrays (no transform): the exact
    alpha-quantile of the endpoint at rho = -1 under plug-in log-scale errors."""
    a, b = sh / kh, sd / kd
    u2 = np.log(kd) / b
    u1 = np.log(kd / kh) / (a + b)
    upper = stats.norm.sf(u2)
    ok = (u1 < u2) & (alpha > stats.norm.cdf(u1) + upper)
    u = stats.norm.ppf(np.clip(alpha - upper, 1e-300, 1.0))
    h, d = kh * np.exp(a * u), kd * np.exp(-b * u)
    return np.where(ok, F_min_vec(np.where(ok, h, 2.0), np.where(ok, d, 1.5)), -np.inf)


def truth_pair(F_int):
    return float(np.exp(F_int + G * np.log(XD))), XD


def screen(F_int, c, s, n=N, seed=SEED):
    """Simulate n experiments; return (noise-free endpoint, P(bound > F0),
    P(bound > 0), P(bound <= F_int))."""
    xh, xd = truth_pair(F_int)
    KH, KD = observe(xh, c), observe(xd, c)
    rng = np.random.default_rng(seed)
    kh = KH * np.exp(s * rng.standard_normal(n))
    kd = KD * np.exp(s * rng.standard_normal(n))
    lb = bound_vec(kh, s * kh, kd, s * kd)
    return (float(F_min_exact(KH, KD)[0]), float(np.mean(lb > F0)),
            float(np.mean(lb > 0.0)), float(np.mean(lb <= F_int)))


def calibrate(F_int, c, s, s_c=S_C, n=N, seed=SEED):
    """Recovered offset at a measured commitment: its 5th, 50th, 95th percentiles,
    and the fraction of experiments with no admissible inversion."""
    xh, xd = truth_pair(F_int)
    KH, KD = observe(xh, c), observe(xd, c)
    rng = np.random.default_rng(seed + 1)
    kh = KH * np.exp(s * rng.standard_normal(n))
    kd = KD * np.exp(s * rng.standard_normal(n))
    ch = c * np.exp(s_c * rng.standard_normal(n))
    ok = ch > kh - 1.0
    with np.errstate(divide="ignore", invalid="ignore"):
        xh_hat = kh * ch / (1.0 + ch - kh)
        xd_hat = kd * ch / (1.0 + ch - kd)
        F = np.where(ok, np.log(xh_hat) - G * np.log(xd_hat), np.inf)
    q = np.quantile(F, [0.05, 0.5, 0.95])
    return float(q[0]), float(q[1]), float(q[2]), float(np.mean(~ok))


def check_vectorized(n=200, seed=7):
    """bound_vec against the scalar exact_bound; returns the largest difference."""
    rng = np.random.default_rng(seed)
    kd = 1.0 + rng.uniform(0.05, 1.5, n)
    kh = kd * (1.0 + rng.uniform(0.05, 6.0, n))
    s = rng.uniform(0.002, 0.05, n)
    v = bound_vec(kh, s * kh, kd, s * kd)
    w = np.array([BU.exact_bound(a, s_ * a, b, s_ * b) for a, b, s_ in zip(kh, kd, s)])
    fin = np.isfinite(v) & np.isfinite(w)
    assert (np.isfinite(v) == np.isfinite(w)).all()
    return float(np.max(np.abs(v[fin] - w[fin])))


def table():
    rows = []
    for name, F_int in TRUTHS:
        for c in C_GRID:
            row = dict(truth=name, F_int=F_int, c=c)
            for s in PREC:
                e, p0, p_ref, cov = screen(F_int, c, s)
                tag = f"{100 * s:g}".replace(".", "")
                row.update({"endpoint": e, f"excl_{tag}": p0, f"excl_ref_{tag}": p_ref,
                            f"cover_{tag}": cov})
            lo, med, hi, bad = calibrate(F_int, c, PREC[1])
            row.update(cal_lo=lo, cal_med=med, cal_hi=hi, cal_bad=bad)
            rows.append(row)
    return pd.DataFrame(rows)


def tex(d):
    """Body of the supplement's synthetic-validation table (tab:synth)."""
    lines = ["% Generated by analysis/synthetic_validation.py -- do not edit; edit the generator."]
    lab = {"above": "$+0.250$", "model": f"${F0 - 0.005:+.3f}$"}
    for r in d.itertuples():
        lines.append(f"{lab[r.truth]} & ${r.c:g}$ & ${r.endpoint:+.3f}$ & "
                     f"${r.excl_1:.3f}$ & ${r.excl_03:.3f}$ & "
                     f"${min(r.cover_1, r.cover_03):.3f}$ & "
                     f"$[{r.cal_lo:+.3f},\\,{r.cal_hi:+.3f}]$\\\\")
    return "\n".join(lines + ["\\bottomrule"]) + "\n"


def main():
    print(f"vectorized bound vs scalar exact_bound: max |diff| {check_vectorized():.1e}")
    d = table()
    pd.set_option("display.width", 200)
    print(d.to_string(index=False, float_format=lambda v: f"{v:+.4f}"))
    d.to_csv("../results/synthetic_validation.csv", index=False)
    with open("../results/sm_table_synth.tex", "w") as fh:
        fh.write(tex(d))
    print("[written] ../results/synthetic_validation.csv, sm_table_synth.tex")


if __name__ == "__main__":
    main()
