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
from corpus import display_family

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
# a series is (family, variant, step): Y209W thymidylate synthase has a hydride
# and a proton series, and selecting by family and variant alone mixed the two
top = list(zip(_top.family, _top.variant, _top.step))
_NT = 6
# Which of the six clear the upper limit of the model at every correlation
# examined, and which only at some.  Until 2026-11-20 the caption said that
# three series clear it "at each of the four assumed correlations"; at
# rho = 0.9 a fourth does.
_RHOS = (-1.0, 0.0, 0.5, 0.9)
_all, _some = [], []
for f, v, st in top:
    sub = b[(b.family == f) & (b.variant == v) & (b.step == st)].set_index("rho").lcb
    above = [rho for rho in _RHOS if sub[rho] > F0]
    if len(above) == len(_RHOS):
        _all.append((f, v, st))
    elif above:
        _some.append((f, v, st, above, sub))
_only = ""
for f, v, st, above, sub in _some:
    lab = f"{display_family(f, v)} {v} ({st})".replace("%", r"\%")
    at = " and ".join(f"$\\rho={rho:g}$" for rho in above)
    vals = ", ".join(f"${sub[rho]:.3f}$" for rho in above)
    _only += f" {lab} clears it only at {at} ({vals})."
_sim = b[b.rho != -1.0]
_adv = b[b.rho == -1.0]
L += [r"\begin{table}[!htbp]",
      r"\caption{Nominal one-sided $95\%$ lower bounds on the largest "
      r"per-temperature endpoint for the six leading series, at the four "
      r"assumed error correlations. " + _words(len(_all)).capitalize() +
      r" series clear the upper limit of the ground-channel model, "
      r"$\Fz=-0.0421$, at all four correlations." + _only +
      r" Each series uses the largest temperature-specific lower bound at "
      r"$0.05/n_T$, counting all measured temperatures. Inadmissible draws enter "
      r"at $-\infty$. This Bonferroni construction requires the individual "
      r"bounds to have their nominal coverage (Supplementary Note~\ref{sm:bounds}). "
      r"Entries at $\rho=-1$ are closed-form quantiles. The others average "
      + _words(_bu.N_REP) + r" replications of " + _sci(_bu.N_DRAW) +
      r" draws; the largest Monte Carlo standard deviation over the "
      + str(len(_sim)) + r" simulated cells of the " + str(b.groupby(["family", "variant", "step"]).ngroups) +
      r" series was " + _sci(_sim.lcb_mc_sd.max(), 1) +
      r", and at $\rho=-1$ the simulated bounds agree with the closed form within "
      + f"${np.ceil((_adv.lcb_mc - _adv.lcb).abs().max() * 1e4) / 1e4:.4f}$" + r".}",
      r"\label{tab:unc-si}", r"\begin{center}\small",
      r"\begin{tabular}{lccccc}", r"\toprule",
            r"Series & point & $\rho=-1$ & $\rho=0$ & $\rho=0.5$ & $\rho=0.9$\\",
      r"\midrule"]
for fam, v, st in top:
    s = b[(b.family == fam) & (b.variant == v) & (b.step == st)]
    lab = f"{display_family(fam, v)} {v}".replace("%", r"\%")
    L.append(f"{lab} & ${s[s.rho==0].point.iloc[0]:.4f}$ & "
             f"${s[s.rho==-1.0].lcb.iloc[0]:.3f}$ & "
             f"${s[s.rho==0].lcb.iloc[0]:.3f}$ & ${s[s.rho==0.5].lcb.iloc[0]:.3f}$ & "
             f"${s[s.rho==0.9].lcb.iloc[0]:.3f}$\\\\")
L += [r"\bottomrule", r"\end{tabular}", r"\end{center}", r"\end{table}", ""]

# The summed-vibronic table (tab:bvib-si) was removed together with its data.
# Its grid maxima sat at the upper edge of the A range in 27 of 32 cells and
# reverse sign when that range is extended, so the entries were wrong in sign
# rather than merely uncertain.  See the SI section "The summed vibronic
# envelope is permissive".

_split = len(L)
s = pd.read_csv("../data/secondary_adh.csv")
s["gamma_obs"] = np.log(s.K_HT) / np.log(s.K_DT)
s = s.sort_values("gamma_obs", ascending=False)
L += [r"\begin{table}[!htbp]",
      r"\caption{Secondary competitive effects for alcohol dehydrogenases with "
      r"separately published Swain--Schaad exponents. The log ratio computed "
      r"from the two printed effects differs from the published exponent by at most "
      + f"${np.abs(s.gamma_obs - s.exp_pub).max():.2f}$" +
      r" and lies within the quoted uncertainty in all " + _words(len(s)) +
      r" rows. The twelve secondary exponents that Tsai and Klinman print for "
      r"their $3\,^{\circ}$C and cryosolvent conditions~\cite{TsaiKlinman2001} "
      r"are not included. H/T accompanies C--H transfer and D/T accompanies C--D transfer. "
      r"Interpretation therefore depends on the tracer backgrounds, reverse "
      r"commitments and secondary equilibrium isotope effects "
      r"(Supplementary Note~\ref{sm:mixed}).}",
      r"\label{tab:sec-si}", r"\begin{center}\small",
      r"\begin{tabular}{llcccc}", r"\toprule",
      r"Enzyme & Form & $2^{\circ}$ H/T & $2^{\circ}$ D/T & $\gamma^{\mathrm{obs}}$ & published\\",
      r"\midrule"]
for _, r in s.iterrows():
    variant = r.variant
    if r.source == "Rucker1992" and "chlorobenzyl" in variant:
        variant = f"p-Cl-benzyl alcohol, ${int(r.T_C)}\\,^{{\\circ}}$C"
    L.append(f"{r.family} & {variant} & "
             f"${r.K_HT:.3f}\\pm{r.K_HT_se:.3f}$ & ${r.K_DT:.3f}\\pm{r.K_DT_se:.3f}$ & "
             f"${r.gamma_obs:.2f}$ & {_pair(r.exp_pub, r.exp_pub_se)}\\\\")
L += [r"\bottomrule", r"\end{tabular}", r"\end{center}", r"\end{table}", ""]

# Float order is a layout choice: the two tables share a note with no text of
# its own.  The secondary-effects table went first while it fitted below the
# preceding note (2026-10-22, morning); after the corpus grew the bounds table,
# the shorter, has to lead or the note's heading is left alone on a half-empty
# page (2026-10-22, later).  Check the page after any change upstream.
open(OUT, "w").write("% Generated by export_si_tables.py -- do not edit.\n"
                     + "\n".join(L) + "\n")
print(f"wrote {OUT}: 2 tables, {len(top)} bound rows, {len(s)} secondary rows")
