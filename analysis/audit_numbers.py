"""Sweep every derived number in the manuscript against the code that produces it.

Motivation.  A number written into prose while composing a sentence is
indistinguishable, to the author, from a number recalled from a computation.
Whitelist verification -- checking the numbers one thinks to check -- cannot
catch a fabricated value, because a fabricated value is exactly the one that
does not come to mind as needing a check.  This module inverts that: it
recomputes every derived quantity from the code and asserts that the manuscript
contains it, and flags any place where the manuscript instead contains a
*different* value for the same quantity (the stale-value failure).

Run from analysis/:   python audit_numbers.py
Exit status is non-zero if any derived quantity is missing or contradicted.
"""
from __future__ import annotations

import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

import masses as M
from completion import (bsao_homogeneity, bsao_masking_by_temperature,
                        yadh_completion)
from network_geometry import BENCH, bypass_to_destroy, endpoint, endpoint_closed
from partial_id import F_min_exact

GSC = M.gamma_sc("C")
F0 = M.offset_F0("C")
ROOT = Path("/home/dong/Workspace/WritePaper/MDPI_Entropy_QuantumB/v3_quantum/manuscript")
DOCS = {"main": ROOT / "pnas_manuscript.tex", "si": ROOT / "si/si_body.tex"}


def derived():
    """(label, value, decimals, which documents must contain it)."""
    out = [
        ("gamma_SC", GSC, 5, ("main", "si")),
        ("F0", F0, 6, ("main", "si")),
    ]
    for a, lab in ((1.3, "hi"), (7.3, "lo")):
        r = yadh_completion(a)
        out.append((f"YADH F_int {lab}", r["F"], 3, ("main", "si")))
        out.append((f"YADH gamma_int {lab}", r["gamma"], 2, ("main", "si")))
    # The BSAO completion detail moved to the supplement when the main text was
    # trimmed; these are SI-only by editorial decision, not by oversight.
    T, m, me = bsao_masking_by_temperature()
    out.append(("BSAO m(25C)", m[3], 3, ("si",)))
    out.append(("BSAO sd m(25C)", me[3], 3, ("si",)))
    mu, se, Q, dof = bsao_homogeneity()
    out.append(("BSAO Q all six", Q, 1, ("si",)))
    out.append(("bypass phi* r=1.31 vs F0",
                bypass_to_destroy(7.13, 1.73, F0, 1.31), 3, ("si",)))
    out.append(("bypass phi* r=1.31 vs 0",
                bypass_to_destroy(7.13, 1.73, 0.0, 1.31), 3, ("si",)))
    out.append(("E_r(0.15) r=1.31", endpoint_closed(7.13, 1.73, 0.15, 1.31), 4, ("si",)))
    # The envelope threshold: mu_ratios() returns SQUARE-ROOT mass ratios, so the
    # code's t is sqrt(mu_T/mu_H) and the threshold sqrt(t) is (mu_T/mu_H)^(1/4).
    # Reading mu_ratios as a plain mass ratio once led to "fixing" the exponent
    # from 1/4 to 1/2, which is wrong; this pins the printed value.
    # binding tolerance (Proposition S8): F_bind below which the exclusion survives
    _Fo = 0.128726
    out.append(("binding tolerance vs F0", _Fo - F0, 3, ("main", "si")))
    out.append(("binding tolerance vs zero", _Fo, 3, ("main", "si")))
    # joint binding-and-bypass tolerance: the marginal numbers are correct but
    # only their conjunction fails, so the joint boundary is what must be pinned
    from joint_nuisance import endpoint as _je, marginal_tolerances as _mt
    _fb, _ph = _mt()
    out.append(("joint: marginal binding tolerance", _fb, 3, ("main", "si")))
    out.append(("joint: marginal bypass tolerance", _ph, 3, ("si",)))
    out.append(("joint: endpoint at (0.090, 0.080)", _je(F_bind=0.090, phi_H=0.080),
                4, ("si",)))
    # the placement effect: quoted both at the reference point and as the
    # regional maximum, because the point value alone understates it
    _pb = _je(F_bind=0.090, phi_H=0.080, put_binding_in="beta")
    _pa = _je(F_bind=0.090, phi_H=0.080, put_binding_in="alpha")
    out.append(("joint: placement effect at the point", abs(_pa - _pb), 3, ("si",)))
    # the profiled worst-case boundary, which is the one that may be quoted
    from joint_nuisance import endpoint_profiled as _jp
    out.append(("joint: profiled endpoint at (0.085,0.080)", _jp(0.085, 0.080, n=5),
                3, ("si",)))
    mC = M.HEAVY["C"]
    _red = lambda m: mC * m / (mC + m)
    _ratio = _red(M.M_T_ATOMIC) / _red(M.M_H_ATOMIC)
    out.append(("envelope threshold (mu_T/mu_H)^(1/4)", _ratio ** 0.25, 4, ("main", "si")))
    # sensitivity illustration, now carrying r
    from completion import commitment_precision_needed, intrinsic_from_commitment
    for m, dp in ((3.0, 3), (5.0, 3)):
        c = m * (7.13 - 1)
        xh = intrinsic_from_commitment(7.13, c)
        xd = intrinsic_from_commitment(1.73, 1.3107 * c)
        Fo = float(np.log(7.13) - GSC * np.log(1.73))
        gain = (float(np.log(xh) - GSC * np.log(xd)) - Fo) / abs(F0)
        out.append((f"YADH gain at m={m:.0f}", gain, 1, ("si",)))
        out.append((f"YADH c precision m={m:.0f}",
                    100 * commitment_precision_needed(7.13, 1.73, c, 1.3107), 1,
                    ("main", "si")))
    out.append(("E(0.35) r=1", endpoint_closed(7.13, 1.73, 0.35, 1.0), 4, ("si",)))
    # the two aqueous horse liver records added from Tsai and Klinman
    from partial_id import F_min_exact as _fme
    out.append(("LADH WT 3C endpoint", _fme(9.0, 2.2)[0], 3, ("si",)))
    out.append(("LADH F93W 3C endpoint", _fme(10.3, 2.0)[0], 3, ("main", "si")))
    for lab, KH, KD, r in BENCH:
        out.append((f"endpoint {lab}", endpoint(KH, KD, 0.0, r), 4, ("si",)))
    ya = pd.read_csv("../data/cha1989_yadh.csv")
    a = ya[ya.note.str.contains("average")].iloc[0]
    out.append(("YADH F_obs", F_min_exact(a.K_HT, a.K_DT)[0], 3, ("main", "si")))
    return out




# The supplement cannot \\ref into the main document, so its references to
# main-text figures are hard numbers.  Inserting a figure renumbers everything
# after it and silently repoints them, which is exactly what happened when the
# atlas figure was added.  This pins each one to a keyword of the caption it is
# meant to cite.
# The supplement also cites main-text equations by number, ten times for Eq. (2).
# Numbered equations are counted in source order whether or not they carry a
# label, so inserting one earlier silently repoints every such reference.
# The general-network census is expensive to recompute on every audit run, so
# the counts are pinned here and checked against the supplement.  Regenerate
# with atlas_general_census.py if the enumeration or its well-posedness
# criteria ever change.
# The corpus size appears in prose in both documents and is not otherwise a
# derived quantity, so it grew stale when records were added.  Pinned here.
CORPUS = {"records": 96, "units": 31}


GENERAL_CENSUS = {"networks": 764, "edges": 5538, "mixed": 3415,
                  "inverse": 1026, "blind": 626, "concave": 471}


SI_EQREFS = {
    2: "eq:map",                          # the observation map
}


SI_FIGREFS = {
    3: "Identified sets for the primary",   # endpoints, the benchmark figure
}


def check_si_figrefs():
    """Every 'Fig.~N of the main text' must point at the intended caption."""
    main = DOCS["main"].read_text()
    caps = re.findall(r"\\caption\{(.{0,120})", main, re.S)
    si = DOCS["si"].read_text()
    bad = 0
    for n in sorted(set(int(m) for m in
                        re.findall(r"Fig\.~(\d+) of the main text", si))):
        want = SI_FIGREFS.get(n)
        got = caps[n - 1].replace("\n", " ") if n <= len(caps) else "<out of range>"
        ok = want is not None and want in got
        print(f"  SI cites main-text Fig. {n}: {got[:52].strip()!r} "
              f"{'OK' if ok else 'MISMATCH'}")
        if not ok:
            bad += 1
    return bad



def check_si_eqrefs():
    """Main-text equation numbers cited by the supplement must still resolve."""
    main = DOCS["main"].read_text()
    order = []
    for m in re.finditer(r"\\begin\{equation\}(.*?)\\end\{equation\}", main, re.S):
        lab = re.search(r"\\label\{([^}]+)\}", m.group(1))
        order.append(lab.group(1) if lab else None)
    si = DOCS["si"].read_text()
    cited = set(int(x) for x in re.findall(r"main-text Eq\.~\((\d+)\)", si))
    bad = 0
    for n in sorted(cited):
        got = order[n - 1] if n <= len(order) else "<out of range>"
        want = SI_EQREFS.get(n)
        ok = want is not None and got == want
        print(f"  SI cites main-text Eq. ({n}) -> {got} {'OK' if ok else 'MISMATCH'}")
        if not ok:
            bad += 1
    return bad


def _si_section_titles(text):
    """Section titles with proper brace matching, so multi-line titles survive."""
    out = []
    for m in re.finditer(r"\\(?:sub)*section\*?\{", text):
        i, d = m.end(), 0
        for j in range(i, len(text)):
            if text[j] == "{":
                d += 1
            elif text[j] == "}":
                if d == 0:
                    break
                d -= 1
        out.append(" ".join(text[i:j].split()))
    return out


def check_si_pointers():
    """Every 'SI Appendix, Some Section' pointer must name a real SI section.

    A pointer may shorten a long title, so a prefix counts as a match; anything
    that is not even a prefix is a pointer to nothing.
    """
    titles = _si_section_titles(DOCS["si"].read_text())
    main = DOCS["main"].read_text()
    ptrs = re.findall(r"\\textit\{SI Appendix\},\s*\\textit\{([^}]*)\}", main)
    bad = 0
    for p in sorted(set(ptrs)):
        q = " ".join(p.split())
        if not any(t.lower().startswith(q.lower()) for t in titles):
            print(f"  SI pointer names no section: {q!r}")
            bad += 1
    print(f"  SI Appendix pointers: {len(set(ptrs))} distinct, "
          f"{len(set(ptrs)) - bad} resolve")
    return bad


def check_general_census():
    """The general-network counts quoted in the supplement must be the ones
    atlas_general_census.py produces."""
    txt = DOCS["si"].read_text()
    a = txt.index(r"\emph{Beyond chains.}")
    sec = txt[a:a + 3000]
    bad = 0
    for lab, v in GENERAL_CENSUS.items():
        if f"${v}$" not in sec:
            print(f"  FAIL general census: {lab} = {v} not in the atlas section")
            bad += 1
    print(f"  general census: {len(GENERAL_CENSUS) - bad}/{len(GENERAL_CENSUS)} "
          f"counts present")
    return bad


# Author's standing instruction: no reference to AI, language models or
# machine-generated writing appears anywhere in the article, the supplement or
# the repository metadata.  Reviewers have repeatedly proposed adding one; this
# check fails the audit if any such text is ever introduced.
AI_TERMS = ("chatgpt", "openai", "anthropic", "copilot", "gemini",
            "large language model", "language model", "generative ai",
            "artificial intelligence", "ai-assisted", "ai assisted",
            "machine-generated", "llm")


def check_cited_scripts():
    """Every script named in the documents must exist in analysis/.

    The SI names a script wherever the computation does not live in
    verify_derivation.py.  Filenames are the least stable thing the text can
    cite, and nothing else here would notice a rename, so they are resolved
    against the directory.
    """
    import re
    from pathlib import Path
    here = Path(__file__).resolve().parent
    bad, seen = 0, set()
    for tag, path in (("main", DOCS["main"]), ("si", DOCS["si"])):
        text = path.read_text()
        for raw in re.findall(r"\\texttt\{([A-Za-z0-9\\_]+\.py)\}", text):
            name = raw.replace("\\_", "_")
            seen.add(name)
            if not (here / name).is_file():
                print(f"  FAIL scripts: {tag} cites {name}, which does not exist")
                bad += 1
    if not seen:
        print("  FAIL scripts: no cited scripts found; the pattern may have "
              "stopped matching")
        return 1
    if not bad:
        print(f"  cited scripts: {len(seen)} named in the text, all present")
    return bad


def check_no_ai_mentions():
    """No AI/model reference may appear in the manuscript or the supplement."""
    from pathlib import Path
    bad = 0
    targets = [("main", DOCS["main"]), ("si", DOCS["si"])]
    root = Path(__file__).resolve().parent.parent
    for extra in ("README.md", "CITATION.cff"):
        f = root / extra
        if f.exists():
            targets.append((extra, f))
    for tag, path in targets:
        t = path.read_text().lower()
        hits = [w for w in AI_TERMS if w in t]
        if hits:
            print(f"  FAIL ai-mention: {tag} contains {hits}")
            bad += 1
    if not bad:
        print(f"  no AI/model reference in {len(targets)} documents, as required")
    return bad


def check_methods_sources():
    """Every data file's own bib key must be cited in Materials and Methods.

    Adding the two aqueous 3 C records brought TsaiKlinman2001 into the corpus
    but not into the Methods source list, so the provenance of two records went
    uncited.  The data files that name their source are checked against the
    citation here.
    """
    import pandas as pd, glob, re
    main = DOCS["main"].read_text()
    m = re.search(r"Isotope effects were taken from the published tables of\s*"
                  r"refs\.?\\?\s*\\cite\{([^}]*)\}", main)
    if not m:
        print("  FAIL methods: could not locate the source citation")
        return 1
    cited = {k.strip() for k in m.group(1).split(",")}
    needed = set()
    for f in glob.glob("../data/*.csv"):
        try:
            d = pd.read_csv(f)
        except Exception:
            continue
        if "source" in d.columns:
            needed |= {str(v) for v in d["source"].dropna().unique()}
    # only files that actually feed the benchmark matter
    needed &= {"Grant1989", "Cha1989", "Bahnson1993", "Bahnson1997",
               "TsaiKlinman2001"}
    missing = sorted(needed - cited)
    if missing:
        print(f"  FAIL methods: data sources not cited in Methods: {missing}")
        return 1
    print(f"  Methods source list: all {len(needed)} named data sources cited")
    return 0


def check_ladh_zero_claim():
    """How many horse liver sets contain zero must match what the text claims.

    Adding the two aqueous 3 C records put one LADH endpoint above zero, which
    silently falsified the standing claim that every horse liver set contains
    zero.  Nothing caught it, so the count is checked against the data here.
    """
    import pandas as pd
    from partial_id import F_min_exact
    NUM = {10: "ten", 11: "eleven", 9: "nine", 8: "eight"}
    d = pd.read_csv("../data/ladh_adh_primary.csv")
    n_tot = len(d)
    n_zero = 0
    for _, r in d.iterrows():
        f, _, _, closed = F_min_exact(r.K_HT, r.K_DT)
        if (f <= 0) if closed else (f < 0):
            n_zero += 1
    bad = 0
    if n_zero == n_tot:
        want = "every"
    else:
        want = f"{NUM.get(n_zero, n_zero)} of the {NUM.get(n_tot, n_tot)}"
    for tag, path in (("main", DOCS["main"]), ("si", DOCS["si"])):
        t = " ".join(path.read_text().split()).lower()
        if want not in t:
            print(f"  FAIL LADH: {tag} does not say '{want}' "
                  f"({n_zero}/{n_tot} sets contain zero)")
            bad += 1
    if not bad:
        print(f"  horse liver: {n_zero} of {n_tot} sets contain zero, "
              f"stated as '{want}' in both documents")
    return bad


def check_ecdhfr_assignment():
    """The disputed 10 C ecDHFR row must not be able to move either series point.

    Ref. Wang2014 distinguishes the light and heavy enzyme rows only by color,
    which the text layer loses, so the unpaired 10 C row is assigned by
    interpolation.  The series statistic is a maximum over temperatures and
    that record is interior to both series, so the assignment cannot select it.
    This test fails if new data ever makes it the maximum.
    """
    import pandas as pd
    from partial_id import F_min_exact
    d = pd.read_csv("../data/trinomial_benchmark.csv")
    m = (d.system.str.contains("Escherichia", na=False)
         & d.variant.astype(str).str.contains("heavy|light", case=False, na=False))
    r10 = m & d.T_C.eq(10.0)
    if r10.sum() != 1:
        print(f"  FAIL ecDHFR: expected one 10 C row, found {int(r10.sum())}")
        return 1
    pts = []
    for tag in ("light enzyme", "heavy enzyme", None):
        dd = d.copy()
        if tag is None:
            dd = dd[~r10]
        else:
            dd.loc[r10, "variant"] = tag
        mm = (dd.system.str.contains("Escherichia", na=False)
              & dd.variant.astype(str).str.contains("heavy|light", case=False,
                                                    na=False))
        pts.append({v: round(max(F_min_exact(r.K_HT, r.K_DT)[0]
                                 for _, r in g.iterrows()), 10)
                    for v, g in dd[mm].groupby("variant")})
    if pts[0] != pts[1] or pts[0] != pts[2]:
        print(f"  FAIL ecDHFR: the 10 C assignment moves a series point: {pts}")
        return 1
    print("  ecDHFR 10 C row: point endpoints invariant under all three "
          "assignments")
    return 0


def check_envelope():
    """The quoted joint envelope must survive the profiled worst case."""
    import joint_nuisance as J
    return J.check_envelope_safe()


def check_joint_counterexamples():
    """The conservative boundary must reject both constructed counterexamples."""
    from joint_nuisance import check_counterexamples
    bad = check_counterexamples(verbose=False)
    print(f"  joint counterexamples rejected by the bound: "
          f"{2 - bad}/2")
    return bad


def check_corpus():
    """The record and unit counts must match the data files, not just each other."""
    import pandas as pd
    n = (len(pd.read_csv("../data/trinomial_benchmark.csv"))
         + len(pd.read_csv("../data/ladh_adh_primary.csv"))
         + len(pd.read_csv("../data/bsao_grant1989.csv")) + 1)
    bad = 0
    if n != CORPUS["records"]:
        print(f"  FAIL corpus: data files hold {n} records, prose says "
              f"{CORPUS['records']}")
        bad += 1
    for d in ("main", "si"):
        t = DOCS[d].read_text()
        if str(CORPUS["records"]) not in t:
            print(f"  FAIL corpus: {CORPUS['records']} absent from {d}")
            bad += 1
    # the README drifted once, holding 94/29/28 after the corpus grew
    readme = (ROOT.parent / "github-package" / "README.md")
    if not readme.exists():
        readme = Path(__file__).resolve().parent.parent / "README.md"
    if readme.exists():
        rt = readme.read_text()
        for q in (str(CORPUS["records"]), str(CORPUS["units"])):
            if q not in rt:
                print(f"  FAIL corpus: README.md does not carry {q}")
                bad += 1
        print(f"  README corpus figures: {CORPUS['records']} records, "
              f"{CORPUS['units']} units")
    print(f"  corpus: {n} records in the data, quoted consistently")
    return bad


def check_repo_metadata():
    """CITATION.cff and README must carry the article's title.

    CITATION.cff is what Zenodo reads when a release is cut, so a stale title
    here becomes a stale archive record.  It drifted once already, through two
    retitlings that were not propagated.
    """
    from pathlib import Path
    main = DOCS["main"].read_text()
    i = main.index("\\title{") + 7
    d = 0
    for j in range(i, len(main)):
        if main[j] == "{":
            d += 1
        elif main[j] == "}":
            if d == 0:
                break
            d -= 1
    title = " ".join(main[i:j].split())
    root = Path(__file__).resolve().parent.parent
    bad = 0
    for name in ("CITATION.cff", "README.md"):
        f = root / name
        if not f.exists():
            continue
        txt = " ".join(f.read_text().split())
        if title not in txt:
            print(f"  FAIL {name} does not carry the article title")
            bad += 1
    print(f"  repo metadata: {2 - bad}/2 files carry the article title")
    return bad


def check_titles():
    """The supplement must carry the article's title.

    PNAS builds the two documents from separate sources, so a retitled article
    leaves the supplement on the old title with nothing to flag it.  That is
    what happened when the paper was reframed around network structure.
    """
    def _title(p):
        t = p.read_text()
        i = t.index(r"\title{")
        depth, out = 0, []
        for ch in t[i + 7:]:
            if ch == "{":
                depth += 1
            elif ch == "}":
                if depth == 0:
                    break
                depth -= 1
            out.append(ch)
        return " ".join("".join(out).split())

    main = _title(DOCS["main"])
    si = _title(ROOT / "si/pnas_si.tex")
    ok = main == si
    print(f"  article title : {main[:58]}...")
    print(f"  supplement    : {'matches' if ok else 'DIFFERS: ' + si[:58]}")
    return 0 if ok else 1


def _texts():
    return {k: p.read_text() for k, p in DOCS.items()}


def run():
    texts = _texts()
    bad = 0
    print(f"{'quantity':32s} {'value':>11s}  documents")
    for lab, val, dp, docs in derived():
        s = f"{abs(val):.{dp}f}"
        for d in docs:
            t = texts[d]
            if s in t:
                verdict = "found"
            else:
                # a near-miss at the same precision is a STALE value, not a gap
                near = re.findall(rf"{re.escape(s[:dp - 1])}\d", t)
                verdict = f"MISSING (near: {sorted(set(near))[:3]})" if near else "MISSING"
                bad += 1
            print(f"  {lab:30s} {val:+11.6f}  {d}: {verdict}")
    bad += check_si_figrefs()
    bad += check_si_eqrefs()
    bad += check_si_pointers()
    bad += check_general_census()
    bad += check_corpus()
    bad += check_cited_scripts()
    bad += check_no_ai_mentions()
    bad += check_methods_sources()
    bad += check_ladh_zero_claim()
    bad += check_ecdhfr_assignment()
    bad += check_envelope()
    bad += check_joint_counterexamples()
    bad += check_repo_metadata()
    bad += check_titles()
    print(f"\n{'FAIL' if bad else 'PASS'}: {bad} derived quantities missing or contradicted")
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(run())
