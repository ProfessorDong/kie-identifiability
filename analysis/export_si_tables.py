"""Emit two tables into the Supplementary Information.

The correlation bracket of the six highest series bounds, and the observed
secondary records of the alcohol dehydrogenases.  Generated from the results
and data files, never retyped.
"""
from __future__ import annotations

import pathlib

import numpy as np
import pandas as pd

import masses as M
from export_benchmark_table import _pair

F0 = M.offset_F0("C")
# The manuscript tree is not part of the public package; fall back to results/
# so the tables are still reproducible from the shipped code and data.
_MS = pathlib.Path("../manuscript/si")
OUT = str(_MS / "si_extra_tables.tex") if _MS.is_dir() else "../results/si_extra_tables.tex"
L = []

b = pd.read_csv("../results/bounds_uncertainty.csv")
import bounds_uncertainty as _bu
_ONES = ("zero one two three four five six seven eight nine ten eleven twelve "
         "thirteen fourteen fifteen sixteen seventeen eighteen nineteen").split()
_TENS = "_ _ twenty thirty forty fifty sixty seventy eighty ninety".split()


def _words(n: int) -> str:
    """An integer below 1000 in words, as the captions print counts."""
    if n < 20:
        return _ONES[n]
    if n < 100:
        t, o = divmod(n, 10)
        return _TENS[t] + (f"-{_ONES[o]}" if o else "")
    h, r = divmod(n, 100)
    return f"{_ONES[h]} hundred" + (f" and {_words(r)}" if r else "")


class _Words(dict):
    def __missing__(self, n):
        return _words(n)


_WORD = _Words()


def _sci(x, nd=0):
    e = int(np.floor(np.log10(abs(x))))
    m = x / 10 ** e
    return (f"$10^{{{e}}}$" if nd == 0 and abs(m - 1) < 1e-9
            else f"${m:.{nd}f}\\times10^{{{e}}}$")
_top = b[b.rho == 0.0].sort_values("point", ascending=False).head(6)
top = list(zip(_top.family, _top.variant))
_NT = 6
# how many of the listed series clear the threshold at the adverse rho = -1;
# the caption states it rather than asserting a fixed answer
_clear = [(f, v) for f, v in top
          if b[(b.family == f) & (b.variant == v) & (b.rho == -1.0)].lcb.iloc[0] > F0]
_WORD.update({0: "no entry reaches it", 1: "one entry clears it at every $\\rho$",
              2: "two entries clear it at every $\\rho$"})
_CLEARWORD = (_WORD[len(_clear)] if len(_clear) <= 2
              else f"{_words(len(_clear))} entries clear it at every $\\rho$")
L += [r"\begin{table}[!htbp]",
      r"\caption{Nominal one-sided $95\%$ lower confidence bounds on the offset for the "
            r"six series with the highest point estimates, bracketed over the unknown "
      r"correlation $\rho$ between the errors of the H/T and D/T effects, which no "
      r"source reports. The threshold for excluding the ground-channel model "
            r"is $\Fz=-0.0421$; " + _CLEARWORD + r". Each bound is on the series "
      r"maximum of the endpoint: the one-sided bound for each temperature at "
      r"level $0.05/n_T$ ($n_T$ counts every recorded temperature; draws, and "
      r"records, with no admissible intrinsic pair enter at $-\infty$), then the "
      r"largest of these, which by the union bound has "
      r"coverage at least $95\%$ whenever each per-temperature bound has its "
      r"nominal coverage. That coverage is approximate: the errors are re-expressed "
      r"on the log scale at the observed effects, so the bound would be exact only "
      r"if the log-scale standard deviations were known and the quantile stayed on "
      r"the admissible branch. Exclusion by such a bound has type-I error at most "
      r"$5\%$ only when that coverage holds. Bounds average " + _WORD[_bu.N_REP] +
      r" replications of " + _sci(_bu.N_DRAW) + r" draws; the largest Monte Carlo "
      r"standard deviation over all " + _WORD[4 * len(b) // 4] + r" cells is " +
      _sci(b.lcb_mc_sd.max(), 1) + r". The bracket includes $\rho=-1$, the "
      r"maximally adverse case.}",
      r"\label{tab:unc-si}", r"\begin{center}\small",
      r"\begin{tabular}{lccccc}", r"\toprule",
            r"Series & point & $\rho=-1$ & $\rho=0$ & $\rho=0.5$ & $\rho=0.9$\\",
      r"\midrule"]
for fam, v in top:
    s = b[(b.family == fam) & (b.variant == v)]
    lab = f"{fam} {v}".replace("%", r"\%")
    L.append(f"{lab} & ${s[s.rho==0].point.iloc[0]:.4f}$ & "
             f"${s[s.rho==-1.0].lcb.iloc[0]:.3f}$ & "
             f"${s[s.rho==0].lcb.iloc[0]:.3f}$ & ${s[s.rho==0.5].lcb.iloc[0]:.3f}$ & "
             f"${s[s.rho==0.9].lcb.iloc[0]:.3f}$\\\\")
L += [r"\bottomrule", r"\end{tabular}", r"\end{center}", r"\end{table}", ""]

# The summed-vibronic table (tab:bvib-si) was withdrawn together with its data.
# Its grid maxima sat at the upper edge of the A range in 27 of 32 cells and
# reverse sign when that range is extended, so the entries were wrong in sign
# rather than merely uncertain.  See the SI section "The summed vibronic
# envelope is permissive".

s = pd.read_csv("../data/secondary_adh.csv")
s["gamma_obs"] = np.log(s.K_HT) / np.log(s.K_DT)
s = s.sort_values("gamma_obs", ascending=False)
L += [r"\begin{table}[!htbp]",
      r"\caption{Observed secondary competitive effects for alcohol "
      r"dehydrogenase, tabulated for reference only. "
      r"$\gamma^{\mathrm{obs}}=\ln\KHT/\ln\KDT$ reproduces the published "
      r"exponent to $0.05$ in every case. These measurements come from the "
      r"mixed-labeling design, in which the H/T effect accompanies C--H cleavage "
      r"and the D/T effect accompanies C--D cleavage, so the two ratios are "
      r"referred to different molecules. Their forward commitments are ordered "
      r"by the primary effect, so forward masking can only deflate these "
      r"exponents (Supplementary Note~\ref{sm:mixed}); reverse commitments and "
      r"secondary equilibrium effects are not bounded, and the ground-channel "
      r"envelope has no prediction for a secondary exponent, so no envelope "
      r"inference is drawn from these values.}",
      r"\label{tab:sec-si}", r"\begin{center}\small",
      r"\begin{tabular}{llcccc}", r"\toprule",
      r"Enzyme & Form & $\KHT$ & $\KDT$ & $\gamma^{\mathrm{obs}}$ & published\\",
      r"\midrule"]
for _, r in s.iterrows():
    L.append(f"{r.family} & {r.variant} & "
             f"${r.K_HT:.3f}\\pm{r.K_HT_se:.3f}$ & ${r.K_DT:.3f}\\pm{r.K_DT_se:.3f}$ & "
             f"${r.gamma_obs:.2f}$ & {_pair(r.exp_pub, r.exp_pub_se)}\\\\")
L += [r"\bottomrule", r"\end{tabular}", r"\end{center}", r"\end{table}", ""]

open(OUT, "w").write("% Generated by export_si_tables.py -- do not edit.\n"
                     + "\n".join(L) + "\n")
print(f"wrote {OUT}: 2 tables, {len(top)} bound rows, {len(s)} secondary rows")
