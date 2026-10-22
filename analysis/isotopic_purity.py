"""Protium in the deuterated comparator, across the whole corpus.

A competitive D/T experiment compares a tritiated tracer with a deuterated,
14C-labeled comparator, and the conversion is read from the 14C.  If part of the
comparator carries protium, that part reacts faster, and the experiment reports
a D/T effect larger than the one for pure deuterated material.  Several sources
correct for this (Cha 1989, Bahnson 1993, Grant & Klinman 1989, Rucker 1992,
Jonsson 1994); the others report uncorrected values, some with a stated purity
and some with none.  data/deuterium_purity.csv records what each source says,
with where it says it, and whether its comparator is deuterated at one position
or two.

Measurement model.  Rates are V/K relative to the D/T experiment's tritiated
tracer.  With component fractions p_i reacting at relative rates k_i, at
fractional conversion f of the comparator and tritium "time" tau,

    1 - f = sum_i p_i exp(-k_i tau),

and the fractional-conversion formula every source uses returns
-ln(1 - f)/tau.  As f -> 0 this is sum_i p_i k_i, the largest overstatement; it
shrinks as conversion proceeds, because the faster molecules are consumed first.

Which rates enter depends on the labeling design (an external audit,
2026-10-21, found the earlier version of this module applied one design to all):

  * Monodeuterated comparator, one tritiated tracer for both experiments (the
    formate dehydrogenase, thymidylate synthase and most dihydrofolate
    reductase sources).  The protium contaminant IS the H/T comparator,
    measured against the same tracer, so it reacts at K_HT exactly:
        components  (1 - p, K_DT), (p, K_HT).
  * Dideuterated comparator with a tracer deuterated at the non-transferred
    position (the alcohol dehydrogenases, amine oxidases and three
    dihydrofolate reductase sources).  Protium can sit at the transferred
    position (an HD molecule) or the other one (DH).  Against the TD tracer the
    HD molecule reacts at chi = M(x_H, r c), not at the measured K_HT = M(x_H, c)
    of the HH/TH experiment; eliminating x_H,
        chi = K_HT (1 + r c) / [r (1 + c) - (r - 1) K_HT]  >=  K_HT  for r >= 1,
    and the DH molecule at K_DT times a secondary H/D effect >= 1.  The
    two-component model with K_HT is therefore the SMALLEST correction a given
    protium level implies; the true correction is larger, and every offset it
    produces is a lower limit.

Correcting lowers K_DT, and the endpoint decreases in K_DT, so correcting an
uncorrected source can only raise its offset: uncorrected protium biases a unit
toward the uninformative side.  The converse holds for a source that HAS
corrected: had its material carried less protium than it corrected for, its
reported D/T effect is too low and its offset too high.  That is the one purity
bias toward exclusion.  Its limit is the uncorrected value, recovered here by
inverting the source's OWN correction with the constants it printed (Cha et al.
1989, note 28: 98.7% DD, 0.65% H at each position, relative rates 1, 3.9 and
1.2), not by a model of ours.  A correction with fixed constants multiplies the
D/T effect by a known factor, so undoing it rescales the value and its error
together, exactly.

Uncertainty.  Every bound is the exact rho = -1 quantile of the endpoint
(bounds_uncertainty.exact_bound) with the correction applied to each draw of the
original pair, so the H/T error enters the corrected D/T effect as it should.
Earlier versions refitted a lognormal to the corrected value and dropped that
term.

Run from analysis/:   python isotopic_purity.py
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.optimize import brentq

import corpus
import masses as M
import bounds_uncertainty as BU
from partial_id import F_min_exact

F0 = M.offset_F0("C")
G = M.gamma_sc("C")
PURITY = pd.read_csv(corpus.DATA / "deuterium_purity.csv")
CORRECTED = {"corrected", "corrected_by_cited_method"}
F_GRID = (0.0, 0.5)                 # the f -> 0 limit and a mid-range conversion
UNITS = corpus.counts()["units"]
STORED = pd.read_csv(corpus.DATA.parent / "results" / "bounds_uncertainty.csv").query("rho == -1.0")
OUT = []

# The corrections the excluding units carry, as their sources state them.  B is
# the primary H/D rate of the protium-at-the-transferred-position molecule and
# C the secondary H/D rate of the other, both relative to the dideuterated
# material.  Cha et al. print every constant (note 28).  Rucker et al. used the
# method of Grant & Klinman 1989 at 98.9% D (1.1% H, split as Grant & Klinman
# split it) and do not print B or C; B is taken as each record's own observed
# H/D ratio and C as Cha's 1.2, whose term changes the factor by 0.1%.
PUBLISHED = {
    "Cha1989": dict(pP=0.0065, pS=0.0065, B=3.9, C=1.2,
                    where="note 28 (Science 243, 1329-1330)"),
    "Rucker1992": dict(pP=0.0055, pS=0.0055, B=None, C=1.2,
                       where="Experimental Procedures; method of Grant & Klinman 1989"),
}
# A correction of Bahnson et al.'s size, for sources that cite their method
# without stating whether, or by how much, they corrected (Chin & Klinman 2000).
CITED_SCENARIO = dict(pP=0.0065, pS=0.0065, B=None, C=1.2)


def say(s=""):
    print(s)
    OUT.append(s)


# ------------------------------------------------------------------ the model
def mixture_rate(comps, f=0.0):
    """Rate a fractional-conversion experiment reports for a comparator made of
    components [(fraction, relative rate), ...], at conversion f."""
    if f == 0.0:
        return sum(p * k for p, k in comps)
    g = lambda tau: sum(p * np.exp(-k * tau) for p, k in comps) - (1 - f)
    kmin = min(k for _, k in comps)
    tau = brentq(g, 0.0, 60.0 / kmin)
    return -np.log(1.0 - f) / tau


def measured_dt(kht, kdt, p, f=0.0):
    """D/T effect reported at conversion f by material with protium fraction p at
    the transferred position, two-component model."""
    if p == 0.0:
        return kdt
    return mixture_rate([(1 - p, kdt), (p, kht)], f)


def true_dt(kht, kdt_meas, p, f=0.0):
    """D/T effect of pure deuterated material, given the reported value (the
    inverse of measured_dt; works on arrays at f = 0)."""
    if p == 0.0:
        return kdt_meas
    if f == 0.0:
        return (kdt_meas - p * kht) / (1.0 - p)
    # A draw so far in the tail that even infinitely slow deuterated material
    # would report more than kdt_meas has no pure value; 0 marks it inadmissible.
    lo, hi = 1e-9, 2.0 * kdt_meas / (1.0 - p) + 1.0
    h = lambda k: measured_dt(kht, k, p, f) - kdt_meas
    if h(lo) >= 0:
        return 0.0
    return brentq(h, lo, hi, xtol=1e-12)


def correction_factor(pP, pS, B, C, f=0.0):
    """Factor by which a correction with fixed constants divided the reported D/T
    effect: rate of the mixture over that of pure dideuterated material."""
    return mixture_rate([(1 - pP - pS, 1.0), (pP, B), (pS, C)], f)


def chi_contaminant(kht, c, r):
    """V/K of the HD contaminant against the TD tracer, in the separable model."""
    return kht * (1 + r * c) / (r * (1 + c) - (r - 1) * kht)


# ------------------------------------------------------------------ the units
def units():
    """Every analysis unit as (label, source row, frame of its records).

    Series are matched to the purity table by DOI and single-condition records
    by source name.  A source whose series and single records used different
    deuterated material (Rucker 1992) has one row for each, marked in the
    `applies` column.
    """
    by_doi = {r.doi: r for r in PURITY.itertuples() if r.applies in ("all", "series")}
    by_name = {r.source: r for r in PURITY.itertuples() if r.applies in ("all", "single")}
    out = []
    for (fam, var, step), g in corpus.series().groupby(["family", "variant", "step"]):
        out.append((f"{corpus.display_family(fam, var)} {var} ({step})",
                    by_doi[g.source_DOI.iloc[0]], g.copy()))
    for _, r in corpus.single_condition().iterrows():
        out.append((f"{r.family} {r.variant}, {r.T_C:g} C ({r.source})",
                    by_name[r.source], pd.DataFrame([r])))
    return out


def evaluate(g, p=0.0, f=0.0, units=1):
    """(point endpoint, one-sided bound at rho = -1) after removing protium
    fraction p from every record, the correction applied draw by draw."""
    pts = [F_min_exact(r.K_HT, true_dt(r.K_HT, r.K_DT, p, f))[0] for r in g.itertuples()]
    if p == 0.0:
        tr = None
    elif f == 0.0:
        tr = lambda h, d: (d - p * h) / (1.0 - p)
    else:
        tr = lambda h, d: true_dt(h, d, p, f)
    return max(pts), BU.exact_series_bound(g, units, transform=tr)


def klass(pt, lb):
    if lb > F0:
        return "excludes"
    return "precision" if pt > F0 else "uninformative"


def threshold(g, which, f=0.0, hi=0.05):
    """Smallest p (to 1e-5) at which the point endpoint (which='point') or the
    bound (which='bound') rises above F0; None if not reached by hi."""
    idx = 0 if which == "point" else 1
    if evaluate(g, 0.0, f)[idx] > F0:
        return 0.0
    if evaluate(g, hi, f)[idx] <= F0:
        return None
    lo, up = 0.0, hi
    while up - lo > 1e-5:
        mid = 0.5 * (lo + up)
        lo, up = (lo, mid) if evaluate(g, mid, f)[idx] > F0 else (mid, up)
    return up


def scenario(p_of, f=0.0):
    """Class counts and changes when each unit is corrected by p_of(source row)."""
    rows = []
    for lab, src, g in units():
        base = klass(*evaluate(g))
        p = p_of(src)
        pt, lb = evaluate(g, p, f)
        rows.append((lab, src.source, src.comparator, p, base, klass(pt, lb), pt, lb))
    return pd.DataFrame(rows, columns=["unit", "source", "comparator", "p", "base",
                                       "now", "point", "bound"])


def as_stated(src):
    """Each uncorrected source at its own measured contamination or purity limit."""
    if src.status in CORRECTED or np.isnan(src.p_transferred):
        return 0.0
    return float(src.p_transferred)


def uniform(p):
    return lambda src: 0.0 if src.status in CORRECTED else p


def undo(g, consts, f=0.0, units=1):
    """(point, bound) for corrected records with the correction removed.

    The correction is taken with its stated constants, so the reported D/T
    effect and its error are multiplied by the same known factor; B = None
    means each record's own observed H/D ratio.
    """
    g = g.copy()
    fac = [correction_factor(consts["pP"], consts["pS"],
                             consts["B"] if consts["B"] else r.K_HT / r.K_DT,
                             consts["C"], f) for r in g.itertuples()]
    g["K_DT"] = g.K_DT.values * fac
    g["K_DT_se"] = g.K_DT_se.values * fac
    pt = max(F_min_exact(r.K_HT, r.K_DT)[0] for r in g.itertuples())
    return pt, BU.exact_series_bound(g, units)


def exclusions():
    """The units that exclude as published, with their purity status and the
    bounds per unit and across all units, with and without their correction.

    Each row: unit, source, status, then bounds (per unit, across UNITS) as
    published and, where a correction was or may have been applied, with it
    removed at f -> 0 and at f = 0.5.
    """
    rows = []
    for lab, src, g in units():
        pt, lb = evaluate(g)
        if lb <= F0:
            continue
        if len(g) > 1:
            # a series: the Monte Carlo bound every other table and figure prints
            # (the exact quantile agrees with it within 0.0004)
            st = STORED[(STORED.family == g.family.iloc[0]) & (STORED.variant == g.variant.iloc[0])
                        & (STORED.step == g.step.iloc[0])]
            lb = float(st.lcb.iloc[0])
        row = dict(unit=lab, source=src.source, status=src.status, point=pt,
                   bound=lb, bound_all=evaluate(g, units=UNITS)[1])
        consts = PUBLISHED.get(src.source) if src.status == "corrected" else (
            CITED_SCENARIO if src.status == "corrected_by_cited_method" else None)
        if src.status == "corrected" and consts is None:
            raise ValueError(f"no stated constants for the correction of {src.source}")
        for f in F_GRID:
            k = "0" if f == 0.0 else "5"
            if consts is None:
                row[f"undo{k}"] = row[f"undo{k}_all"] = np.nan
            else:
                row[f"undo{k}"] = undo(g, consts, f)[1]
                row[f"undo{k}_all"] = undo(g, consts, f, UNITS)[1]
        rows.append(row)
    return pd.DataFrame(rows)


def protium_needed(lab, f=0.0, units=None):
    """Smallest protium content (total, both positions) of a corrected source's
    material for which its bound, adjusted across `units` units at rho = -1,
    still clears F0; None if it clears even with the correction removed.

    The source corrected for its measured level p; had the true level been
    lam p, the pure D/T effect is the reported one times g(p) / g(lam p).
    """
    units = units or UNITS
    src, g = {l: (s, gg) for l, s, gg in units_()}[lab]
    c = PUBLISHED[src.source]

    def bound(lam):
        h = g.copy()
        fac = [correction_factor(c["pP"], c["pS"], c["B"] or r.K_HT / r.K_DT, c["C"], f)
               / correction_factor(lam * c["pP"], lam * c["pS"], c["B"] or r.K_HT / r.K_DT,
                                   c["C"], f) for r in g.itertuples()]
        h["K_DT"] = g.K_DT.values * fac
        h["K_DT_se"] = g.K_DT_se.values * fac
        return BU.exact_series_bound(h, units)
    if bound(0.0) > F0:
        return None
    lam = brentq(lambda x: bound(x) - F0, 0.0, 1.0, xtol=1e-6)
    return lam * (c["pP"] + c["pS"])


units_ = units


def _tex(label):
    """A unit label as the supplement prints it."""
    import re
    label = label.replace("%", "\\%").replace("Mg2+", "Mg$^{2+}$")
    label = re.sub(r", (-?\d+) C \((\w+)\)", r", \1\\,$^{\\circ}$C (\2)", label)
    label = re.sub(r"\(([A-Za-z]+?)(Klinman)?(\d{4})\)", lambda m: f"({m.group(1)} {m.group(3)})", label)
    return label.replace(" (hydride)", "").replace("wild type", "WT")


def purity_table(allrows):
    """Body of the supplement's purity-scenario table (tab:purity)."""
    lines = ["% Generated by analysis/isotopic_purity.py -- do not edit; edit the generator."]
    short = {"uninformative": "unin.", "precision": "prec.", "excludes": "excl."}
    for (f, name), s in allrows.groupby(["f", "scenario"], sort=False):
        c = s.now.value_counts()
        up = s[(s.base != s.now)]
        moved = "; ".join(f"{_tex(u.unit)}, {short[u.base]}~$\\to$ {short[u.now]}"
                          for u in up.itertuples())
        lab = "$f\\to0$" if f == 0 else f"$f={f:g}$"
        lines.append(f"{name.replace('%', '\\%')} & {lab} & {c.get('uninformative', 0)} & "
                     f"{c.get('precision', 0)} & {c.get('excludes', 0)} & {moved or 'none'}\\\\")
    return "\n".join(lines + ["\\bottomrule"]) + "\n"


_EXCL_NAMES = {
    "YADH wild type, 25 C (Cha1989)": "YADH, benzyl alcohol",
    "YADH wild type (p-chlorobenzyl alcohol) (hydride)": "YADH, $p$-Cl-benzyl alcohol",
    "LADH F93W (40% MeOH) (hydride)": "LADH F93W, $40\\%$ MeOH",
    "LADH F93W (60% MeOH) (hydride)": "LADH F93W, $60\\%$ MeOH",
    "LADH F93W-I224G, 25 C (ChinKlinman2000)": "LADH F93W:I224G",
    "LADH F93W-I224V, 25 C (ChinKlinman2000)": "LADH F93W:I224V",
}
_STATUS = {"corrected": "corrected", "corrected_by_cited_method": "cites a correction",
           "uncorrected_measured": "uncorrected"}


def exclusion_table(ex):
    """Body of the table of the nominal exclusions (tab:excl): per-unit and
    all-unit bounds at rho = -1, with each correction as published and removed."""
    def cell(v):
        return "n/a" if np.isnan(v) else f"${v:+.3f}$"
    lines = ["% Generated by analysis/isotopic_purity.py -- do not edit; edit the generator."]
    order = list(_EXCL_NAMES)
    if set(ex.unit) != set(order):
        raise ValueError(f"the nominal exclusions changed: {sorted(ex.unit)}")
    ex = ex.set_index("unit").loc[order].reset_index()
    for r in ex.itertuples():
        lines.append(f"{_EXCL_NAMES[r.unit]} & {_STATUS[r.status]} & ${r.bound:+.3f}$ & "
                     f"${r.bound_all:+.3f}$ & {cell(r.undo0)} & {cell(r.undo0_all)} & "
                     f"{cell(r.undo5_all)}\\\\")
    return "\n".join(lines + ["\\bottomrule"]) + "\n"


def main():
    say("=" * 96)
    say("PROTIUM IN THE DEUTERATED COMPARATOR")
    say("=" * 96)
    say(f"F0 = {F0:+.4f}; bounds one-sided 95% at rho = -1, exact quantiles of the")
    say("endpoint with every correction applied draw by draw (series: Bonferroni over T).")
    say("")
    say("Sources (data/deuterium_purity.csv):")
    for r in PURITY.itertuples():
        pv = "" if np.isnan(r.p_transferred) else f"p = {r.p_transferred:.4f}"
        say(f"  {r.source:20s} {r.comparator:15s} {r.status:26s} {pv}")
    say("")
    kh, kd = 7.0, 1.8
    say("Size of the bias, K_HT = 7.0, K_DT = 1.8, p = 0.01:")
    for f in (0.0, 0.2, 0.5, 0.8):
        say(f"  f = {f:.1f}: reported D/T {measured_dt(kh, kd, 0.01, f):.4f}")
    say("")
    say("Dideuterated designs: the HD contaminant against the TD tracer, chi = M(x_H, r c).")
    for c, r in ((22.529, 1.456134), (1e6, 1.456134), (22.529, 1.0)):
        say(f"  K_HT 7.13, c {c:g}, r {r:.4f}: chi = {chi_contaminant(7.13, c, r):.4f}")
    say("  chi >= K_HT whenever r >= 1, so the two-component model is the smallest correction.")
    say("")
    say("Cha et al.'s own correction (note 28), as a factor on the D/T effect:")
    for f in (0.0, 0.2, 0.5, 0.8):
        say(f"  f = {f:.1f}: {correction_factor(0.0065, 0.0065, 3.9, 1.2, f):.5f}")
    say("")
    rows = []
    for f in F_GRID:
        for name, pof in (("as stated", as_stated), ("0.5% where uncorrected", uniform(0.005)),
                          ("1% where uncorrected", uniform(0.01))):
            s = scenario(pof, f)
            c = s.now.value_counts()
            say(f"f = {f:.1f}, {name}: uninformative {c.get('uninformative', 0)}, "
                f"precision-limited {c.get('precision', 0)}, excluding {c.get('excludes', 0)}")
            for r in s[s.base != s.now].itertuples():
                say(f"    {r.unit:52s} {r.comparator:15s} {r.base:>13s} -> {r.now:13s} "
                    f"point {r.point:+.4f} bound {r.bound:+.4f}")
            rows.append(s.assign(f=f, scenario=name))
    allrows = pd.concat(rows)
    allrows.to_csv("../results/isotopic_purity.csv", index=False)
    with open("../results/sm_table_purity.tex", "w") as fh:
        fh.write(purity_table(allrows))
    say("")
    say("Smallest protium fraction at the transferred position that moves each")
    say("uncorrected unit (f -> 0 | f = 0.5); for dideuterated comparators an upper limit:")
    say(f"{'unit':52s}{'point > F0':>18s}{'bound > F0':>18s}")
    for lab, src, g in units():
        if src.status in CORRECTED:
            continue
        cells = []
        for which in ("point", "bound"):
            v = [threshold(g, which, f) for f in F_GRID]
            cells.append(" | ".join("  -  " if x is None else f"{100 * x:.2f}%" for x in v))
        say(f"{lab:52s}{cells[0]:>18s}{cells[1]:>18s}")
    say("")
    ex = exclusions()
    ex.to_csv("../results/purity_exclusions.csv", index=False)
    with open("../results/sm_table_exclusions.tex", "w") as fh:
        fh.write(exclusion_table(ex))
    say(f"The nominal exclusions: bounds per unit and across {UNITS} units (rho = -1), as")
    say("published and with any correction removed (published constants; for a source")
    say("that only cites a correction, a correction of Bahnson et al.'s size):")
    say(f"{'unit':52s}{'status':>27s}{'unit':>8s}{'all':>8s}"
        f"{'undo0':>8s}{'all':>8s}{'undo.5':>8s}{'all':>8s}")
    for r in ex.itertuples():
        say(f"{r.unit:52s}{r.status:>27s}{r.bound:+8.3f}{r.bound_all:+8.3f}"
            f"{r.undo0:+8.3f}{r.undo0_all:+8.3f}{r.undo5:+8.3f}{r.undo5_all:+8.3f}")
    say("")
    say(f"Protium the corrected material must have carried for the bound across {UNITS} units")
    say("at rho = -1 to clear F0 (the source measured 1.3% for benzyl alcohol, 1.1% for")
    say("p-chlorobenzyl alcohol):")
    for lab in ("YADH wild type, 25 C (Cha1989)",
                "YADH wild type (p-chlorobenzyl alcohol) (hydride)"):
        for f in F_GRID:
            v = protium_needed(lab, f)
            say(f"  {lab:52s} f = {f:.1f}: "
                + ("clears with no protium" if v is None else f"at least {100 * v:.2f}%"))
    with open("../results/isotopic_purity.txt", "w") as fh:
        fh.write("\n".join(OUT) + "\n")
    say("[written] ../results/isotopic_purity.txt, isotopic_purity.csv, purity_exclusions.csv,")
    say("          sm_table_purity.tex, sm_table_exclusions.tex")


if __name__ == "__main__":
    main()
