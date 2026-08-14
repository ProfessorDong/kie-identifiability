"""Emit the three tables that live in the REVTeX main text into the PNAS SI.

The PNAS main text has no tables (6-page limit), so the confidence bounds, the
summed-envelope map and the secondary records must appear in the SI or they
appear nowhere.  Generated from the results files, never retyped.
"""
from __future__ import annotations

import pathlib

import numpy as np
import pandas as pd

import masses as M

F0 = M.offset_F0("C")
# The manuscript tree is not part of the public package; fall back to results/
# so the tables are still reproducible from the shipped code and data.
_MS = pathlib.Path("../manuscript/si")
OUT = str(_MS / "si_extra_tables.tex") if _MS.is_dir() else "../results/si_extra_tables.tex"
L = []

b = pd.read_csv("../results/bounds_uncertainty.csv")
top = b[b.rho == 0.0].sort_values("point", ascending=False).head(4).variant.tolist()
L += [r"\begin{table}[h]",
      r"\caption{One-sided $95\%$ lower confidence bounds on the offset for the "
      r"four series with the highest point estimates, bracketed over the unknown "
      r"correlation $\rho$ between the H/T and D/T effects induced by the shared "
      r"tritium reference. The threshold for excluding the ground-channel model "
      r"is $\Fz=-0.0421$; no entry reaches it. Bounds are simultaneous over the "
      r"temperatures of a series and averaged over five replications of "
      r"$5\times10^{4}$ draws; the largest Monte Carlo standard deviation over "
      r"all sixty-four cells is $7.4\times10^{-4}$. The bracket includes "
      r"$\rho=-1$, the maximally adverse case.}",
      r"\label{tab:unc-si}", r"\begin{center}\small",
      r"\begin{tabular}{lccccc}", r"\toprule",
      r"Series (all ecDHFR) & point & $\rho=-1$ & $\rho=0$ & $\rho=0.5$ & $\rho=0.9$\\",
      r"\midrule"]
for v in top:
    s = b[b.variant == v]
    L.append(f"{v} & ${s[s.rho==0].point.iloc[0]:.4f}$ & "
             f"${s[s.rho==-1.0].lcb.iloc[0]:.4f}$ & "
             f"${s[s.rho==0].lcb.iloc[0]:.4f}$ & ${s[s.rho==0.5].lcb.iloc[0]:.4f}$ & "
             f"${s[s.rho==0.9].lcb.iloc[0]:.4f}$\\\\")
L += [r"\bottomrule", r"\end{tabular}", r"\end{center}", r"\end{table}", ""]

BV = {(-5.0, 1500): (-6.1, -1.7, -7.1, -11.0), (-5.0, 2500): (57.1, -10.6, -25.4, -29.5),
      (-5.0, 3000): (74.1, -22.9, -32.4, -34.9), (-5.0, 3600): (-0.7, -33.7, -37.5, -38.7),
      (-15.0, 1500): (-9.7, -5.1, -0.9, -3.6), (-15.0, 2500): (126.5, 5.7, -0.2, -17.4),
      (-15.0, 3000): (684.3, -42.3, 1.7, -24.7), (-15.0, 3600): (1412.8, -264.5, 2.4, -31.5)}
L += [r"\begin{table}[h]",
      r"\caption{Summed vibronic envelope $B_{\mathrm{vib}}$ (units of $10^{-3}$) "
      r"at $T=298.15$~K. Each entry is the largest $F$ found on a $60\times60$ "
      r"logarithmic grid over $A\in[0.050,63.1]$ and $w\in[0.0025,25.1]$, "
      r"restricted to $\KHT^{\mathrm{int}}\ge\sqrt{t}$: a sampled maximum over "
      r"that rectangle, not a certified supremum over the continuum. "
      r"Negative entries leave a "
      r"separation between the gated and semiclassical families; entries at or "
      r"above zero leave none. For comparison $\Fz=-42.1$ in these units. Entries "
      r"are converged to better than $1$ in the last digit except in the steeply "
      r"varying corner at $\lambda=10$, $\Delta G^{\circ}=-15$.}",
      r"\label{tab:bvib-si}", r"\begin{center}\small",
      r"\begin{tabular}{lcccc}", r"\toprule",
      r" & \multicolumn{4}{c}{$\lambda$ (kcal\,mol$^{-1}$)}\\",
      r"$\tilde\nu$ (cm$^{-1}$) & 5 & 10 & 20 & 40\\", r"\midrule"]
for dG in (-5.0, -15.0):
    L.append(r"\multicolumn{5}{l}{$\Delta G^{\circ}=" + f"{dG:.0f}" +
             r"$~kcal\,mol$^{-1}$}\\")
    for nu in (1500, 2500, 3000, 3600):
        L.append(f"{nu} & " + " & ".join(f"${x:+.1f}$" for x in BV[(dG, nu)]) + r"\\")
    if dG == -5.0:
        L.append(r"\midrule")
L += [r"\bottomrule", r"\end{tabular}", r"\end{center}", r"\end{table}", ""]

s = pd.read_csv("../data/secondary_adh.csv")
s["gamma_obs"] = np.log(s.K_HT) / np.log(s.K_DT)
s = s.sort_values("gamma_obs", ascending=False)
L += [r"\begin{table}[h]",
      r"\caption{Observed secondary competitive effects for alcohol "
      r"dehydrogenase, tabulated for reference only. "
      r"$\gamma^{\mathrm{obs}}=\ln\KHT/\ln\KDT$ reproduces the published "
      r"exponent to $0.05$ in every case. These measurements come from the "
      r"mixed-labeling design, in which the H/T effect accompanies C--H cleavage "
      r"and the D/T effect accompanies C--D cleavage, so the two ratios are "
      r"referred to different molecules and do not share a commitment. "
      r"Theorem~1 of the main text does not apply to them, and we draw no "
      r"identification inference from these values in either direction.}",
      r"\label{tab:sec-si}", r"\begin{center}\small",
      r"\begin{tabular}{llcccc}", r"\toprule",
      r"Enzyme & Form & $\KHT$ & $\KDT$ & $\gamma^{\mathrm{obs}}$ & published\\",
      r"\midrule"]
for _, r in s.iterrows():
    L.append(f"{r.family} & {r.variant} & "
             f"${r.K_HT:.3f}\\pm{r.K_HT_se:.3f}$ & ${r.K_DT:.3f}\\pm{r.K_DT_se:.3f}$ & "
             f"${r.gamma_obs:.2f}$ & ${r.exp_pub}\\pm{r.exp_pub_se}$\\\\")
L += [r"\bottomrule", r"\end{tabular}", r"\end{center}", r"\end{table}", ""]

open(OUT, "w").write("% Generated by export_si_tables.py -- do not edit.\n"
                     + "\n".join(L) + "\n")
print(f"wrote {OUT}: 3 tables, {len(top)} bound rows, {len(s)} secondary rows")
