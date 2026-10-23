"""Data for the asymmetry figure: primary identified sets against the two
thresholds.  The secondary panel was removed: the mixed-labeling design does not
satisfy the shared-commitment condition, so no inference is drawn from it."""
from __future__ import annotations

import numpy as np
import pandas as pd

import masses as M
from partial_id import F_min_exact
import corpus

GSC = M.gamma_sc("C"); F0 = M.offset_F0("C")
OUT = "../figures/tikz/data/"

# ---- panel (a): primary identified sets, one row per series / enzyme form
rows = []
RHO = -1.0            # maximally adverse; every bound drawn is the worst case


def _lcb(kht, kht_se, kdt, kdt_se, rho=RHO):
    """One-sided 95% lower bound on a single record's endpoint at `rho`.

    The quantile of the endpoint itself (bounds_uncertainty.single_bound), not
    F_min - 1.645 sd(F_obs): the shortcut is the quantile only on the open
    branch, and it drew the two closed-branch Hong 2007 whiskers 0.06 and 0.10
    too short until 2026-10-21.
    """
    import bounds_uncertainty as BU
    return BU.single_bound(kht, kht_se, kdt, kdt_se, rho)

_TB = pd.read_csv("../data/trinomial_benchmark.csv")
_FAM = {"FDH": "formate dehydrogenase", "TSase": "thymidylate synthase",
        "ecDHFR": "Escherichia coli dihydrofolate reductase",
        "hsDHFR": "Homo sapiens dihydrofolate reductase",
        "R67DHFR": "R67 dihydrofolate reductase",
        "MAOB": "monoamine oxidase",
        "LADH": "horse liver alcohol dehydrogenase",
        "YADH": "Saccharomyces cerevisiae alcohol dehydrogenase"}


def _series_closed(r):
    """Is the plotted endpoint of this series attained at finite commitment?

    A series row plots the largest endpoint over its temperatures, so whether
    that endpoint is attained is a property of the single record achieving it.
    The identified set is closed exactly when L_H < gamma; otherwise the
    endpoint is approached only as the commitment diverges and is excluded.
    """
    sub = _TB[_TB.system.str.contains(_FAM.get(r.family, r.family), case=False,
                                      na=False)
              & (_TB.variant.astype(str) == str(r.variant))]
    if r.family == "TSase" and "step" in _TB.columns:
        sub = sub[sub.step.astype(str) == str(r.step)]
    vals = [F_min_exact(x.K_HT, x.K_DT) for _, x in sub.iterrows()
            if x.K_HT > x.K_DT > 1]
    if not vals:
        raise ValueError(f"no benchmark record for {r.family} {r.variant}")
    best = max(vals, key=lambda v: v[0])
    if abs(best[0] - r.point) > 1e-6:
        raise ValueError(f"endpoint mismatch for {r.family} {r.variant}: "
                         f"{r.point} vs {best[0]}")
    return bool(best[3])


b = pd.read_csv("../results/bounds_uncertainty.csv")
# Conservative throughout: every bound is the maximally adverse rho = -1 case,
# so any system drawn as clearing F0 clears it under any correlation.
b = b[b.rho == RHO]
# TSase contributes both a hydride- and a proton-transfer series; without the
# step in the label the two WT rows are indistinguishable in the figure.
multi = {f for f, g in b.groupby("family") if g.step.nunique() > 1}
# Thymidylate synthase appears as two species, and since the wild-type E. coli
# record was admitted (2026-09-14) "TSase WT" no longer names one enzyme.  The
# species is read from the benchmark's own system column, not assumed.
_SPECIES = {"Escherichia coli": "ec", "Homo sapiens": "hs"}


def _tsase_prefix(variant):
    sysname = _TB[(_TB.family == "TSase")
                  & (_TB.variant.astype(str) == str(variant))].system.iloc[0]
    return next(v for k, v in _SPECIES.items() if sysname.startswith(k))


for _, r in b.iterrows():
    if r.family == "TSase":
        lab = f"{_tsase_prefix(r.variant)}TSase {r.variant} ({r.step})"
    else:
        lab = f"{r.family} {r.variant}"
        if r.family in multi:
            lab += f" ({r.step})"
    g = _TB[(_TB.family == r.family) & (_TB.variant.astype(str) == str(r.variant))
            & (_TB.step.astype(str) == str(r.step))]
    rows.append((lab, r.point, r.lcb, "series", _series_closed(r), g))

single = corpus.single_condition()
# A form measured by two studies at different temperatures is two analysis
# units, as the two amine oxidase pH values already are; the temperature
# disambiguates the label only where a form appears more than once.
_dup = {(fam, v) for (fam, v), g in single.groupby(["family", "variant"]) if len(g) > 1}
_SHORT = {"wild type": "WT"}
for _, r in single.iterrows():
    f, _, _, cl = F_min_exact(r.K_HT, r.K_DT)
    if r.grp == "ectsase":
        lab = f"ecTSase {_SHORT.get(r.variant, r.variant)} ({r.step})"
    else:
        lab = f"{r.family} {r.variant}"
        if (r.family, r.variant) in _dup:
            lab += f" {int(r.T_C)}$^\\circ$C"
    rows.append((lab, f, _lcb(r.K_HT, r.K_HT_se, r.K_DT, r.K_DT_se), r.grp, cl,
                 pd.DataFrame([r])))

df = pd.DataFrame(rows, columns=["label", "point", "lcb", "grp", "closed", "recs"])
# Two marks, both computed (2026-10-15, after an external audit):
#   dagger   - a series with a recorded pair outside the normal-order domain;
#              its row summarizes the admissible temperatures only
#   asterisk - a unit whose rho = -1 bound still clears F0 after adjustment
#              across all analysis units (bounds_uncertainty.multiplicity_bound)
import bounds_uncertainty as BU
_U = corpus.counts()["units"]
df["dagger"] = [bool(((g.K_HT <= g.K_DT) | (g.K_DT <= 1)).any()) for g in df.recs]
df["star"] = [bool(lcb > F0 and BU.multiplicity_bound(g, _U, RHO) > F0)
              for lcb, g in zip(df.lcb, df.recs)]
df = df.sort_values("point").reset_index(drop=True)
df["y"] = np.arange(len(df)) + 1
inf_mask = df.lcb > F0

# ---- single-column panel.  The full endpoint range spans 1.0 log units but the
# comparison is decided within 0.05 of F0, so the axis is cut at XMIN and the two
# systems below it are marked off-scale with their endpoints printed in the row
# label.  At this scale |F0| spans ~3.5 mm, so no magnified inset is needed.
XMIN, XMAX = -0.50, 0.20
RAY = 0.050                      # short ray: the set is open above
off = df[df.point < XMIN]

# Row labels: the figure is two panels side by side, so labels are kept short.
# The caption and Supplementary Table 7 give the full names.
_ABBR = (("wild type (p-chlorobenzyl alcohol)", "WT, p-Cl-benzyl alcohol"),
         ("wild type (p-methoxybenzyl alcohol)", "WT, p-MeO-benzyl alcohol"),
         ("BSAO p-methoxybenzylamine", "BSAO WT, p-MeO-benzylamine"),
         ("BSAO wild type", "BSAO WT, benzylamine"),
         ("wild type", "WT"), ("Mg2+", "Mg$^{2+}$"),
         ("(dihydrobiopterin)", "(dihydrobiopterin)"),
         ("heavy enzyme", "heavy"), ("light enzyme", "light"))


def _short(lab):
    for a, b in _ABBR:
        lab = lab.replace(a, b)
    return lab


# A form measured both as a series and as a single record (R67 DHFR: Yahashiri
# et al. 2008, Sen et al. 2011) gets the same short label twice; name the source
# on both rows, as the caption promises.  Series sources come from the labeling
# table, single records carry their own.
import re as _re
_TD = pd.read_csv("../data/tracer_design.csv")


def _source(r):
    recs = r.recs
    if "source" in recs.columns:
        key = str(recs.source.iloc[0])
    else:
        x = recs.iloc[0]
        key = _TD[(_TD.family == x.family) & (_TD.variant == x.variant)].source.iloc[0]
    return _re.sub(r"([A-Za-z])(\d{4})", r"\1 \2", key.split(";")[0].strip())


_shortlab = df.label.map(lambda s: _short(s))
for i in df.index[_shortlab.duplicated(keep=False)]:
    df.at[i, "label"] = f"{df.at[i, 'label']} ({_source(df.loc[i])})"
df["label"] = df.label.str.replace("R67DHFR", "R67 DHFR", regex=False)


def _write_panel(fh, sub):
    """Rays, whiskers, endpoints and labels for one panel, rows y = 1, 2, ..."""
    for y, (_, r) in enumerate(sub.iterrows(), start=1):
        col  = "cPro" if r.lcb > F0 else "cGrey"
        wcol = "cPro!45" if r.lcb > F0 else "cGrey!40"
        lw = "0.9pt" if r.lcb > F0 else "0.6pt"
        start = max(r.point, XMIN)
        right = r.point > XMAX            # endpoint beyond the right edge
        if not right:
            fh.write(f"\\draw[{col}, line width={lw}, -{{Straight Barb[length=1.9pt,"
                     f"width=2.8pt]}}] (axis cs:{start:.5f},{y}) -- "
                     f"(axis cs:{min(start+RAY, XMAX):.5f},{y});\n")
        if right:
            # off the right edge: whisker to the edge and a chevron there; the
            # endpoint is printed in the row label, as for the left edge
            fh.write(f"\\draw[{wcol}, line width={lw}] "
                     f"(axis cs:{max(r.lcb,XMIN):.5f},{y}) -- "
                     f"(axis cs:{XMAX:.5f},{y});\n")
            fh.write(f"\\node[font=\\fontsize{{7.6}}{{9}}\\selectfont, color={col}, "
                     f"anchor=east, inner sep=0.6pt] at (axis cs:{XMAX:.5f},{y}) "
                     f"{{$\\gg$}};\n")
        elif r.point >= XMIN:
            fh.write(f"\\draw[{wcol}, line width={lw}] "
                     f"(axis cs:{max(r.lcb,XMIN):.5f},{y}) -- "
                     f"(axis cs:{r.point:.5f},{y});\n")
            # filled when the endpoint is attained at finite commitment,
            # open when it is only approached, which the caption relies on
            if r.closed:
                fh.write(f"\\fill[{col}] (axis cs:{r.point:.5f},{y}) "
                         f"circle (0.85pt);\n")
            else:
                fh.write(f"\\draw[{col}, line width=0.45pt, fill=white] "
                         f"(axis cs:{r.point:.5f},{y}) circle (0.95pt);\n")
        else:
            # endpoint lies off the left edge: chevron instead of a dot
            fh.write(f"\\node[font=\\fontsize{{7.6}}{{9}}\\selectfont, color={col}, "
                     f"anchor=east, inner sep=0.6pt] at (axis cs:{XMIN:.5f},{y}) "
                     f"{{$\\ll$}};\n")
        lab = _short(r.label)
        lab = (f"{lab} (${r.point:+.2f}$)" if r.point > XMAX else
               lab if r.point >= XMIN else f"{lab} (${r.point:.2f}$)")
        lab = lab.replace("%", "\\%")   # solvent rows carry a percent sign
        lab += "$^{\\dagger}$" if r.dagger else ""
        lab += "$^{\\ast}$" if r.star else ""
        # off-scale rows carry a chevron at XMIN, so their labels sit further
        # left to keep clear of it
        dx = 0.032 if r.point >= XMIN else 0.060
        fh.write(f"\\node[rowlab] at (axis cs:{XMIN-dx:.5f},{y}) {{{lab}}};\n")


# Panel A holds the lower half of the endpoints, panel B the upper half (every
# unit whose bound or endpoint reaches the model envelope).  Until 2026-10-20 the
# figure was one panel of 35 rows; the census doubled the corpus.
NA = len(df) // 2
for name, sub in (("fa_rays_a.tex", df.iloc[:NA]), ("fa_rays_b.tex", df.iloc[NA:])):
    with open(OUT + name, "w") as fh:
        fh.write("% generated by export_fig_asym.py -- do not edit\n")
        _write_panel(fh, sub)
with open(OUT + "fa_rows.tex", "w") as fh:
    fh.write("% generated by export_fig_asym.py -- do not edit\n")
    fh.write(f"\\def\\faRowsA{{{NA}}}\\def\\faRowsB{{{len(df) - NA}}}\n")
print(f"panel: {len(df)} rows, {len(off)} off-scale below {XMIN}")
for _, r in off.iterrows():
    print(f"   off-scale: {r.label} at {r.point:+.4f}")

print(f"panel: {int(inf_mask.sum())} informative of {len(df)}")
print(f"panel: dagger on {list(df.label[df.dagger])}; asterisk on {list(df.label[df.star])}")
print(f"panel: {int(df.closed.sum())} closed endpoints, "
      f"{int((~df.closed).sum())} open")
