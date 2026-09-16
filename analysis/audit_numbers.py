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

import os
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

import masses as M
from completion import (bsao_homogeneity, bsao_masking_by_temperature,
                        yadh_completion, yadh_completion_joint, bypass_tolerance,
                        mc_yadh_tolerance, _mc_yadh)
from network_geometry import BENCH, bypass_to_destroy, endpoint, endpoint_closed
from partial_id import F_min_exact

GSC = M.gamma_sc("C")
F0 = M.offset_F0("C")
# Which manuscript to audit: see manuscript_root.py.  KIE_MANUSCRIPT, or the
# untracked manuscript_root.local, names the directory; the main text is its one
# *_manuscript.tex.  A target named but not readable is an error, not a skip --
# see run().
import manuscript_root as MR
ROOT = MR.root()
DOCS = {"main": MR.one(ROOT, "*_manuscript.tex"), "si": ROOT / "si/si_body.tex"}


def _titled(text):
    """Content of \\title{...}, tolerating the optional \\title[short]{...} form.

    The document class takes an optional running-head argument, and a bare
    index("\\title{") raises on it.
    """
    m = re.search(r"\\title(?![a-zA-Z])", text)   # not \titleformat, \titlespacing
    if m is None:
        raise ValueError("no \\title in document")
    j = m.end()
    while j < len(text) and text[j] in " \n\t":
        j += 1
    if j < len(text) and text[j] == "[":
        d = 0
        while j < len(text):
            if text[j] == "[":
                d += 1
            elif text[j] == "]":
                d -= 1
                if d == 0:
                    j += 1
                    break
            j += 1
        while j < len(text) and text[j] in " \n\t":
            j += 1
    if j >= len(text) or text[j] != "{":
        raise ValueError("no brace group after \\title")
    d, out = 0, []
    for ch in text[j + 1:]:
        if ch == "{":
            d += 1
        elif ch == "}":
            if d == 0:
                break
            d -= 1
        out.append(ch)
    return " ".join("".join(out).split())


def derived():
    """(label, value, decimals, which documents must contain it)."""
    out = [
        ("gamma_SC", GSC, 5, ("main", "si")),
        ("F0", F0, 6, ("main", "si")),
    ]
    # Reference exponents quoted in Supplementary Note 1 for comparison with the
    # sources.  These were hand-typed until 2026-09-14, when a reference-review
    # found the bare-mass value given as 3.26278 in one place and 3.26281 in
    # another; masses.py constants give 3.262807.  Each is now recomputed here
    # from the same atomic masses, so a stale or transposed value fails.
    def _exp(mH, mD, mT, heavy=None):
        red = (lambda m: m) if heavy is None else (lambda m: m * heavy / (m + heavy))
        f = lambda m: red(m) ** -0.5
        return (f(mH) - f(mT)) / (f(mD) - f(mT))
    g_bare = _exp(M.M_H_ATOMIC, M.M_D_ATOMIC, M.M_T_ATOMIC)
    out.append(("gamma bare atomic masses", g_bare, 5, ("si",)))
    out.append(("gamma integer 1:2:3", _exp(1.0, 2.0, 3.0), 5, ("si",)))
    # Saunders 1985 quotes 3.34 for reduced masses; integer isotope masses
    # reduced against carbon-12 reproduce it, as the note states.
    out.append(("gamma reduced, integer vs C12", _exp(1.0, 2.0, 3.0, 12.0), 5, ("si",)))
    # Swain et al. 1958 wrote the relation as the 1.442 power of the H/D effect.
    out.append(("gamma/(gamma-1) bare (Swain 1.442)", g_bare / (g_bare - 1), 5, ("si",)))
    for a, lab in ((1.3, "hi"), (7.3, "lo")):
        # r solved jointly with the commitment (completion.yadh_joint_r); the
        # observed secondary ratio 1.31 is only a lower bound on it
        r = yadh_completion_joint(a)
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
    # yeast reference asymmetry and bypass tolerance at the jointly solved r
    # (until 2026-09-14 quoted at r = 1.31, which overstates the tolerance)
    for a in (1.3, 7.3):
        j = yadh_completion_joint(a)
        out.append((f"YADH joint r a={a}", j["r"], 3, ("si",)))
        out.append((f"YADH tolerance vs F0 a={a} (%)",
                    100 * float(bypass_tolerance(7.13, 1.73, j["r"], F0)), 1, ("si",)))
        out.append((f"YADH tolerance vs 0 a={a} (%)",
                    100 * float(bypass_tolerance(7.13, 1.73, j["r"], 0.0)), 1, ("si",)))
    c23 = yadh_completion_joint(2.3)
    out.append(("YADH F_int a=2.3", c23["F"], 3, ("main", "si")))
    out.append(("YADH gamma_int a=2.3", c23["gamma"], 2, ("si",)))
    t = mc_yadh_tolerance()
    for lab in ("F0", "zero"):
        for k in ("median", "lo", "hi"):
            out.append((f"YADH tolerance MC {lab} {k} (%)", 100 * t[lab][k], 1, ("si",)))
    out.append(("YADH tolerance MC F0 P(>10%)", t["F0"]["p10"], 2, ("si",)))
    Fu, _ = _mc_yadh(400000, 20260815, (1.3, 7.3))
    for q in (2.5, 50, 97.5):
        out.append((f"YADH completion MC uniform a, q{q}", float(np.percentile(Fu, q)), 3, ("si",)))
    ya = pd.read_csv("../data/cha1989_yadh.csv")
    a = ya[ya.note.str.contains("average")].iloc[0]
    out.append(("YADH F_obs", F_min_exact(a.K_HT, a.K_DT)[0], 3, ("main", "si")))
    # Cha's errors are labelled standard deviations; until 2026-09-14 both
    # documents claimed that using them as standard errors was conservative by
    # sqrt(5).  The errors on the average in fact match the standard error of the
    # three determinations, which is now what the text says, with these numbers.
    det = ya[ya.note.str.contains("determination")]
    out.append(("YADH SE of 3 determinations, H/T",
                det.K_HT.std(ddof=1) / np.sqrt(len(det)), 3, ("main", "si")))
    out.append(("YADH SE of 3 determinations, D/T",
                det.K_DT.std(ddof=1) / np.sqrt(len(det)), 3, ("main", "si")))
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
CORPUS = {"records": 97, "units": 32}


GENERAL_CENSUS = {"networks": 764, "edges": 5538, "mixed": 3415,
                  "inverse": 1026, "blind": 626, "concave": 471}


SI_EQREFS = {
    2: "eq:map",                          # the observation map
}


SI_FIGREFS = {
    3: "Identified sets for the primary",   # endpoints, the benchmark figure
}


def check_si_figrefs():
    """Every 'Fig.~N of the main text' must point at the intended caption.

    Captions are collected from figure environments only.  Counting every
    \\caption in the file breaks the moment the manuscript gains a table: the
    table's caption shifts the index and every later figure reference is then
    compared against the wrong float.  That is what happened when Table 1 was
    added on 2026-09-11 -- the guard reported Fig. 3 as the atlas figure.
    """
    main = DOCS["main"].read_text()
    caps = []
    for fig in re.findall(r"\\begin\{figure\*?\}(.*?)\\end\{figure\*?\}",
                          main, re.S):
        m = re.search(r"\\caption\{(.{0,120})", fig, re.S)
        if m:
            caps.append(m.group(1))
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


def _si_numbered_titles(text):
    """(supplementary-note number, title) for every SI heading.

    The note number is the index of the enclosing \\section, so a pointer at a
    subsection resolves to the note that contains it.  Brace matching so that
    multi-line titles survive.
    """
    out, note = [], 0
    for m in re.finditer(r"\\(sub)?section\*?\{", text):
        i, d = m.end(), 0
        for j in range(i, len(text)):
            if text[j] == "{":
                d += 1
            elif text[j] == "}":
                if d == 0:
                    break
                d -= 1
        if not m.group(1):
            note += 1
        out.append((note, " ".join(text[i:j].split())))
    return out


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
    """Every pointer into the supplement must name a real heading.

    The form used is "Supplementary Note N, \\textit{Title}",
    which carries both a number and a title, so both are checked: the title must
    name a heading and N must be the note that heading actually sits in.  A
    mismatched pair is the failure this guards -- renumbering the supplement
    while leaving the main text pointing at the old numbers.

    A document carrying no pointers at all fails: a renamed supplement once
    turned this guard into a silent 0-of-0 pass.
    """
    numbered = _si_numbered_titles(DOCS["si"].read_text())
    main = DOCS["main"].read_text()

    def resolve(q):
        return [n for n, t in numbered if t.lower().startswith(q.lower())]

    bad = 0
    dotted = re.findall(
        r"Supplementary Note~?(\d+),\s*\\textit\{([^}]*)\}", main)
    legacy = re.findall(
        r"Supplementary Information,\s*\\textit\{([^}]*)\}", main)

    if not dotted and not legacy:
        print("  FAIL SI pointers: main text names no supplement heading at all")
        return 1

    seen = set()
    for num, title in dotted:
        q = " ".join(title.split())
        seen.add((num, q))
    for num, q in sorted(seen):
        hits = resolve(q)
        if not hits:
            print(f"  FAIL SI pointer names no heading: {q!r}")
            bad += 1
        elif len(set(hits)) > 1:
            print(f"  FAIL SI pointer is ambiguous: {q!r} -> notes {sorted(set(hits))}")
            bad += 1
        elif hits[0] != int(num):
            print(f"  FAIL SI pointer {q!r} says Note {num}, heading is in Note {hits[0]}")
            bad += 1

    seenl = {" ".join(x.split()) for x in legacy}
    for q in sorted(seenl):
        if not resolve(q):
            print(f"  FAIL SI pointer names no heading: {q!r}")
            bad += 1

    total = len(seen) + len(seenl)
    print(f"  Supplementary pointers: {total} distinct, {total - bad} resolve "
          f"(number and title checked on {len(seen)})")
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
# The hard rule this enforces is that no AI or model reference appears in the
# manuscript, the supplement, or the repository metadata.  The list below is
# deliberately over-inclusive: a false positive costs one look, a false
# negative defeats the guard.  It omitted "claude" and "co-authored-by" until
# 2026-09-09, so a perturbation test injecting "Claude" passed silently -- the
# exact failure this suite exists to prevent.
#
# Two deliberate exclusions, both to avoid false positives that would train the
# reader to ignore the guard: bare "gpt" (GPT is glutamic-pyruvic transaminase,
# an enzyme name that could legitimately appear here, so the hyphenated
# "gpt-" is used instead) and "assistant" (a research assistant may properly be
# thanked in an acknowledgement).
AI_TERMS = ("chatgpt", "openai", "anthropic", "claude", "copilot", "gemini",
            "gpt-", "large language model", "language model", "generative ai",
            "artificial intelligence", "ai-assisted", "ai assisted",
            "ai-generated", "ai generated", "co-authored-by", "chatbot",
            "machine-generated", "llm")


def _corpus_records():
    import corpus
    return corpus.counts()["records"]


def check_provenance_split():
    """The SI's provenance counts must match the benchmark's own source DOIs.

    Five sources are parsed from a publisher PDF, one from a published
    spreadsheet, and one was transcribed by hand from scanned pages.  Only the
    first two classes are evidence of transcription; quoting a larger parsed
    count than the data supports would overstate the audit.
    """
    import pandas as pd
    PDF = {"10.1021/ja411998h", "10.1021/acs.biochem.1c00558",
           "10.1021/acscatal.9b03345", "10.1021/ja501936d",
           "10.3390/ijms16047304"}
    XLSX = {"10.1371/journal.pone.0196506"}
    HAND = {"10.1021/bi00253a026"}
    d = pd.read_csv("../data/trinomial_benchmark.csv")
    n_pdf = int(d.source_DOI.isin(PDF).sum())
    n_xls = int(d.source_DOI.isin(XLSX).sum())
    n_hand = int(d.source_DOI.isin(HAND).sum())
    parsed = n_pdf + n_xls
    bad = 0
    if parsed + n_hand != len(d):
        print(f"  FAIL provenance: {parsed}+{n_hand} != {len(d)}; "
              "a source DOI is unclassified")
        bad += 1
    si = " ".join(DOCS["si"].read_text().split())
    main = " ".join(DOCS["main"].read_text().split())
    for label, need, doc, text in (
            ("SI parsed total", f"${parsed}$ of the ${len(d)}$", "si", si),
            ("SI value count", f"${parsed * 4}$ values", "si", si),
            ("SI pdf values", f"{n_pdf * 4} values", "si", si),
            ("main parsed total", f"{parsed} of the {len(d)} temperature-resolved records",
             "main", main),
            # every record not parsed was transcribed by hand, including all the
            # single-condition records; the count follows from the corpus
            ("main hand total", f"The remaining ${_corpus_records() - parsed}$ were "
             f"transcribed by hand", "main", main)):
        if need not in text:
            print(f"  FAIL provenance: {doc} does not state {need!r} ({label})")
            bad += 1
    if not bad:
        print(f"  provenance: {n_pdf} PDF ({n_pdf*4} values) + {n_xls} spreadsheet "
              f"({n_xls*4}) = {parsed} parsed, {n_hand} hand-transcribed")
    return bad


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


def check_precision_limited_count():
    """How many units clear F0 on the point estimate but not on the bound.

    This is the group the SI describes when the two aqueous records are folded
    in.  Adding them took the group from two to three, and the SI carried the
    old count for a while, so it is checked against the data here.
    """
    import pandas as pd, numpy as np
    from partial_id import F_min_exact
    NUM = {2: "two", 3: "three", 4: "four"}
    G = M.gamma_sc("C")
    F0 = M.offset_F0("C")

    def lcb(pt, kh, sh, kd, sd, rho=-1.0):
        a, b = sh / kh, sd / kd
        return pt - 1.645 * np.sqrt(a ** 2 + (G * b) ** 2 - 2 * G * rho * a * b)

    import corpus
    pts = []
    b = pd.read_csv("../results/bounds_uncertainty.csv")
    b = b[b.rho == -1.0]
    pts += [(r.point, r.lcb) for _, r in b.iterrows()]
    # every single-condition record, read from the one list that defines them;
    # until 2026-09-14 this function kept its own list of files
    for _, r in corpus.single_condition().iterrows():
        p = F_min_exact(r.K_HT, r.K_DT)[0]
        pts.append((p, lcb(p, r.K_HT, r.K_HT_se, r.K_DT, r.K_DT_se)))
    n = sum(1 for pt, l in pts if pt > F0 and l <= F0)
    word = NUM.get(n, str(n))
    ORD = {3: "third", 4: "fourth", 5: "fifth"}
    need = {"main": (f"while {word},", f"remaining {word} the endpoint",
                     f"and {word} more limited by precision alone",
                     f"and {word} more are limited by precision alone"),
            "si": (f"a {ORD.get(n, n)} unit whose endpoint exceeds",)}
    bad = 0
    for tag, phrases in need.items():
        t = " ".join(DOCS[tag].read_text().split()).lower()
        for ph in phrases:
            if ph not in t:
                print(f"  FAIL precision-limited: {tag} lacks '{ph}' "
                      f"({n} units clear on the point estimate but not the bound)")
                bad += 1
    if not bad:
        print(f"  precision-limited units: {n}, stated as '{word}' in both documents")
    return bad


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


from functools import lru_cache


@lru_cache(maxsize=1)
def _envelope_r_result():
    import joint_nuisance as J
    return J.check_envelope_r(verbose=False)


def check_envelope():
    """Both quoted joint lines must survive the profiled worst case.

    The r = 1.31 line at that r; the r-range line at every stipulated r from
    1.31 to the largest jointly solved value (added 2026-09-15).
    """
    import joint_nuisance as J
    bad = J.check_envelope_safe()
    fails, tight, where = _envelope_r_result()
    print(f"  r-range joint line: {fails} failures, tightest slack {tight:.2e} at phi={where:.3f}")
    return bad + fails


def check_joint_counterexamples():
    """The conservative boundary must reject both constructed counterexamples."""
    from joint_nuisance import check_counterexamples
    bad = check_counterexamples(verbose=False)
    print(f"  joint counterexamples rejected by the bound: "
          f"{2 - bad}/2")
    return bad


def check_corpus():
    """The record and unit counts must match the data files, not just each other."""
    import corpus
    c = corpus.counts()
    n = c["records"]
    bad = 0
    if n != CORPUS["records"]:
        print(f"  FAIL corpus: data files hold {n} records, prose says "
              f"{CORPUS['records']}")
        bad += 1
    if c["units"] != CORPUS["units"]:
        print(f"  FAIL corpus: data files give {c['units']} analysis units, prose "
              f"says {CORPUS['units']}")
        bad += 1
    # Phrases, not bare numbers: "97" occurs inside "1997", so the old bare test
    # passed on a document that still said 96 (found 2026-09-14).
    R, U = CORPUS["records"], CORPUS["units"]
    for d, phrases in (("main", (f"benchmark of {R} matched", f"{R} matched records in {U} analysis units",
                                 f"We assembled {R} matched primary records", f"giving {U} analysis units",
                                 f"${R}$ records in ${U}$ analysis units")),
                       ("si", (f"all ${R}$ records",))):
        t = " ".join(DOCS[d].read_text().split())
        for ph in phrases:
            if ph not in t:
                print(f"  FAIL corpus: {d} lacks {ph!r}")
                bad += 1
    # the README drifted once, holding 94/29/28 after the corpus grew
    readme = (ROOT.parent / "github-package" / "README.md")
    if not readme.exists():
        readme = Path(__file__).resolve().parent.parent / "README.md"
    if readme.exists():
        rt = " ".join(readme.read_text().split())
        # The whole phrase, not the bare numbers: "97" and "32" occur inside
        # page numbers and years, so a bare-number test passed on a README that
        # still said 96 and 31 (found 2026-09-14).
        q = f"{CORPUS['records']} records over {CORPUS['units']} analysis units"
        if q not in rt:
            print(f"  FAIL corpus: README.md does not carry {q!r}")
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
    title = _titled(DOCS["main"].read_text())
    root = Path(__file__).resolve().parent.parent
    bad = 0
    for name in ("CITATION.cff", "README.md"):
        f = root / name
        if not f.exists():
            continue
        txt = " ".join(f.read_text().split())
        # case-insensitive: a journal may require Title Case, the deposit keeps sentence
        # case, and what must not drift is the wording
        if title.lower() not in txt.lower():
            print(f"  FAIL {name} does not carry the article title")
            bad += 1
    # a version-specific archive citation must name the release CITATION.cff describes
    cff = (root / "CITATION.cff").read_text() if (root / "CITATION.cff").exists() else ""
    mv = re.search(r"^version:\s*v?([\d.]+)", cff, re.M)
    cited = set(re.findall(r"version analyzed here is v([\d.]+)|version v([\d.]+), doi", " ".join(DOCS["main"].read_text().split())))
    cited = {a or b for a, b in cited}
    if mv and cited and cited != {mv.group(1)}:
        print(f"  FAIL manuscript cites deposit version {sorted(cited)}, CITATION.cff is v{mv.group(1)}")
        bad += 1
    print(f"  repo metadata: {2 - bad}/2 files carry the article title")
    return bad


def check_reference_asymmetry():
    """Every r stated in the text must equal r computed from the secondary effects.

    Until 2026-09-14 each r was hand-typed, and monoamine oxidase B at pH 6.1 was
    tabulated as 1.14 -- the value of bovine serum amine oxidase and of MAO-B at
    pH 7.5 -- when its own data give 1.13 to 1.21.  This checks, against
    reference_asymmetry.py: every "inferred" cell of the protocol table (a single
    value, or a range where secondary effects were measured at several
    temperatures), the r column of the bypass table, the r of each
    network_geometry.BENCH entry, and the overall range quoted in both documents.
    """
    import reference_asymmetry as RA
    import network_geometry as NG
    d = RA.table()
    bad = 0
    fam = {"YADH": "YADH", "LADH": "LADH", "BSAO": "BSAO", "MAOB": "MAOB"}

    def expect(family, variant):
        lo, hi = RA.r_range(family, variant)
        a, b = f"{lo:.2f}", f"{hi:.2f}"
        return a if a == b else f"{a}--{b}"

    si = DOCS["si"].read_text()
    rows = re.findall(r"^(YADH|LADH|BSAO|MAOB) (.+?) & [^&\n]+ & yes & ([\d.\-]+) & inferred",
                      si, re.M)
    if len(rows) != d.groupby(["family", "variant"]).ngroups:
        print(f"  FAIL reference asymmetry: protocol table has {len(rows)} inferred "
              f"rows, data hold {d.groupby(['family', 'variant']).ngroups} forms")
        bad += 1
    for f, v, cell in rows:
        want = expect(fam[f], v)
        if cell != want:
            print(f"  FAIL reference asymmetry: protocol table {f} {v} says {cell}, "
                  f"data give {want}")
            bad += 1

    bypass = {"Yeast alcohol dehydrogenase": ("YADH", "wild type", 25.0),
              "Horse liver F93W": ("LADH", "F93W", 25.0),
              "Bovine serum amine oxidase": ("BSAO", "wild type", 25.0),
              "Monoamine oxidase~B, pH~6.1, $10\\,^{\\circ}$C": ("MAOB", "pH 6.1", 10.0)}
    for label, (f, v, T) in bypass.items():
        m = re.search(re.escape(label) + r" & \$[\d.]+\$ & \$([\d.]+)\$", si)
        want = f"{RA.r_at(f, v, T):.2f}"
        if f == "YADH":
            # the yeast row carries r solved jointly with the published commitment,
            # tabulated as the range over Klinman's a; its lower end is checked here
            from completion import yadh_joint_r
            want = f"{float(yadh_joint_r(7.3)):.2f}"
        if m is None or m.group(1) != want:
            print(f"  FAIL reference asymmetry: bypass table {label!r} r = "
                  f"{m.group(1) if m else 'absent'}, data give {want}")
            bad += 1

    # a shared-tracer row carries r = 1 exactly, from tracer_design.csv
    m = re.search(r"DHFR light enzyme, \$25\\,\^\{\\circ\}\$C & \$[\d.]+\$ & \$?([\d.\-]+)\$?", si)
    want = f"{RA.r_design('ecDHFR', 'light enzyme'):g}"
    if m is None or m.group(1) != want:
        print(f"  FAIL reference asymmetry: bypass table DHFR light r = "
              f"{m.group(1) if m else 'absent'}, tracer design gives {want}")
        bad += 1
    bench = {"YADH (Cha 1989)": ("YADH", "wild type", 25.0),
             "BSAO (Grant 1989)": ("BSAO", "wild type", 25.0),
             "MAO-B pH 6.1, 10 C": ("MAOB", "pH 6.1", 10.0),
             "LADH F93W": ("LADH", "F93W", 25.0)}
    for lab, _, _, r in NG.BENCH:
        if lab == "ecDHFR light, 25 C" and r != RA.r_design("ecDHFR", "light enzyme"):
            print(f"  FAIL reference asymmetry: BENCH {lab} r = {r}, tracer design gives 1")
            bad += 1
        if lab in bench:
            want = RA.r_at(*bench[lab])
            if abs(r - want) > 0.005:
                print(f"  FAIL reference asymmetry: network_geometry BENCH {lab} r = "
                      f"{r:.4f}, data give {want:.4f}")
                bad += 1

    lo, hi = f"{d.r.min():.2f}", f"{d.r.max():.2f}"
    for doc in ("main", "si"):
        txt = " ".join(DOCS[doc].read_text().split())
        # Two phrasings.  The second ("$r$ running from $a$, for ..., to $b$")
        # escaped the first pattern, and a stale "1.14 ... to 1.31" survived in
        # it after the MAO-B correction (found 2026-09-14).
        # Since 2026-09-15 the observed ratio is written R (a lower bound on r),
        # so both letters are matched, and a document quoting no range fails
        # rather than passing on nothing.
        found = re.findall(r"\$[rR]=(\d\.\d\d)\$ to \$(\d\.\d\d)\$", txt)
        found += re.findall(r"\$[rR]\$ (?:runs|running) from \$(\d\.\d\d)\$[^$]*?"
                            r"(?:\$[^$]*\$[^$]*?)*? to \$(\d\.\d\d)\$", txt)
        if not found:
            print(f"  FAIL reference asymmetry: {doc} quotes no range for the observed ratio")
            bad += 1
        # a quoted range must be the overall one or that of a single family
        # (the monoamine oxidase series is quoted on its own)
        fams = {(f"{g.r.min():.2f}", f"{g.r.max():.2f}") for _, g in d.groupby("family")}
        for a, b in found:
            if (a, b) != (lo, hi) and (a, b) not in fams:
                print(f"  FAIL reference asymmetry: {doc} quotes r = {a} to {b}, "
                      f"data give {lo} to {hi} overall and no family has that range")
                bad += 1
    if not bad:
        print(f"  reference asymmetry: {len(rows)} protocol rows, 4 bypass rows, "
              f"BENCH and the quoted range {lo} to {hi} agree with the data")
    return bad


def check_promoting_mode():
    r"""The promoting-mode numbers must be those completion.py computes.

    The temperature spans of F, the crossover wavenumber and the range over
    which the sensitivity peak moves were hand-typed.  On 2026-09-14 the
    reference review found the text arguing from the 50 cm^-1 span, so each is
    now recomputed: every "$s|\Fz|$ at $nu$~cm$^{-1}$" pair in either document
    must match offset_vs_T at A = 5, w0 = 0.5 over 250-350 K, all four quoted
    wavenumbers must appear in both, and the crossover and peak range must agree
    with sensitivity_peak over the sweep A in [1, 50], w0 in [0.1, 2].
    """
    from completion import CM_TO_K, F0 as _F0, offset_vs_T, sensitivity_peak
    span = {nu: abs(np.diff(offset_vs_T(5.0, 0.5, nu, [250.0, 350.0]))[0]) / abs(_F0)
            for nu in (50, 200, 400, 800)}
    peaks = [sensitivity_peak(A, w0)[0] for A in (1.0, 2.0, 5.0, 10.0, 20.0, 50.0)
             for w0 in (0.1, 0.25, 0.5, 1.0, 2.0)]
    lo, hi, cross = f"{min(peaks):.0f}", f"{max(peaks):.0f}", f"{2 * 300 / CM_TO_K:.0f}"
    bad = 0
    pair = re.compile(r"\$(\d+\.\d\d)\|\\Fz\|\$\s+(?:for an?|at)\s+\$(\d+)\$~cm")
    for doc in ("main", "si"):
        t = DOCS[doc].read_text()
        seen = set()
        for val, nu in pair.findall(t):
            nu = int(nu)
            seen.add(nu)
            if nu not in span or val != f"{span[nu]:.2f}":
                want = f"{span[nu]:.2f}" if nu in span else "not computed"
                print(f"  FAIL promoting mode: {doc} gives {val}|F0| at {nu} cm-1, "
                      f"completion.py gives {want}")
                bad += 1
        for nu in span:
            if nu not in seen:
                print(f"  FAIL promoting mode: {doc} does not quote the span at {nu} cm-1")
                bad += 1
        if not re.search(rf"is\s+\${cross}\$~cm", t):
            print(f"  FAIL promoting mode: {doc} does not give the crossover as {cross} cm-1")
            bad += 1
        if not re.search(rf"between\s+\${lo}\$\s+and\s+\${hi}\$~cm", t):
            print(f"  FAIL promoting mode: {doc} does not give the peak range {lo}-{hi} cm-1")
            bad += 1
    if not bad:
        print("  promoting mode: spans " + ", ".join(f"{span[n]:.2f} at {n}" for n in span)
              + f"; crossover {cross}; peak {lo}-{hi} cm-1, in both documents")
    return bad


def check_corpus_consequences():
    r"""Statements that follow from the size and content of the corpus.

    Admitting the Agrawal 2004 wild-type E. coli thymidylate synthase record
    (2026-09-14) changed the number of units, and so the multiplicity correction,
    the number of forms in the protocol table, and every "all N records" claim.
    None of those was computed by any deposited code: the Bonferroni bounds were
    hand-typed.  Each is now derived here and matched to the sentence that
    states it, not merely to a number appearing somewhere.
    """
    import corpus
    from scipy import stats
    from holdout import structural_prediction_check
    _W = {0: "zero", 13: "thirteen", 16: "sixteen", 17: "seventeen", 29: "twenty-nine",
          30: "thirty", 31: "thirty-one", 32: "thirty-two", 33: "thirty-three"}

    class _Words(dict):
        # a count with no word is written as digits, so an unexpected count
        # fails the phrase test instead of crashing the whole suite
        def __missing__(self, k):
            return str(k)
    W = _Words(_W)
    main = " ".join(DOCS["main"].read_text().split())
    si = " ".join(DOCS["si"].read_text().split())
    bad = 0

    def need(doc, text, what):
        nonlocal bad
        if text not in (main if doc == "main" else si):
            print(f"  FAIL corpus consequences: {doc} lacks {what}: {text!r}")
            bad += 1

    c = corpus.counts()
    u, n = c["units"], c["records"]
    # multiplicity, for the yeast record, at rho = 0 and rho = -1
    ya = corpus.single_condition()
    a = ya[ya.grp == "yadh"].iloc[0]
    z = stats.norm.ppf(1 - 0.05 / u)
    F = np.log(a.K_HT) - GSC * np.log(a.K_DT)
    x, y = a.K_HT_se / a.K_HT, GSC * a.K_DT_se / a.K_DT
    b0 = F - z * np.hypot(x, y)
    b1 = F - z * (x + y)
    need("si", f"correction across the {W[u]} gives ${b0:+.3f}$ under independence "
               f"and ${b1:+.3f}$ at $\\rho=-1$", "the Bonferroni bounds")
    need("si", f"the margin over $\\Fz$ is ${b1 - F0:.3f}$", "the simultaneous margin")
    need("si", f"Bonferroni-corrected bound at that correlation by ${b1 - F0:.3f}$",
         "the simultaneous margin (joint-nuisance note)")
    need("main", f"singled out from {W[u]} examined", "the multiplicity denominator")
    need("si", f"singled out from {W[u]} examined", "the multiplicity denominator")
    need("main", f"applied uniformly to all {W[u]} systems", "the unit count")
    # the monotonicity condition of E(phi), and the Theorem 1 forecast
    rows = pd.concat([corpus.series(), corpus.single_condition()])
    holds = sum(1 for _, r in rows.iterrows()
                if (r.K_HT - 1) - GSC * (r.K_DT - 1) < (r.K_HT - 1) * (r.K_DT - 1) * (GSC - 1))
    if holds == n:
        need("si", f"holds for all ${n}$ records here", "the E(phi) monotonicity count")
    else:
        print(f"  FAIL corpus consequences: E(phi) monotone for {holds} of {n}")
        bad += 1
    chk, viol = structural_prediction_check()
    if viol or chk != n:
        print(f"  FAIL corpus consequences: Theorem 1 forecast {chk} checked, {viol} violations")
        bad += 1
    need("si", f"Checked against all ${chk}$ records, there are no violations",
         "the Theorem 1 forecast count")
    # the protocol table: one row per distinct series or form
    tab = DOCS["si"].read_text()
    tab = tab[tab.index(r"\label{tab:protocol}"):]
    tab = tab[:tab.index(r"\end{tabular}")]
    lines = [ln.rstrip() for ln in tab.splitlines()]
    inferred = sum(1 for ln in lines if ln.endswith(r"inferred$^{a}$\\"))
    shared = sum(1 for ln in lines if ln.endswith(r"& 1 & shared$^{b}$\\"))
    forms = inferred + shared
    # the table's statuses must be those deposited in tracer_design.csv, which
    # records the tracer each source used in the H/T and the D/T experiment
    import reference_asymmetry as RA
    td = RA.tracer_design()
    # row by row, so a swap of two statuses cannot pass on the counts alone
    cells = {}
    for ln in lines:
        if ln.endswith(r"inferred$^{a}$\\") or ln.endswith(r"& 1 & shared$^{b}$\\"):
            cells[ln.split(" & ")[0]] = "labeled" if "inferred" in ln else "shared"
    from corpus import display_family
    for _, x in td.iterrows():
        if x.design == "labeled":
            stem = f"{x.family} {x.variant}"
        else:
            fam = display_family(x.family, x.variant) if x.family == "TSase" else x.family
            var = "WT" if x.family == "ecTSase" else x.variant
            stem = f"{fam} {var} ({x.step}"
        # labeled rows are named exactly; shared rows carry "(step" and may add a
        # temperature, so they match on that prefix (which includes the step)
        hit = [k for k in cells if (k == stem if x.design == "labeled" else k.startswith(stem))]
        if len(hit) != 1 or cells[hit[0]] != x.design:
            print(f"  FAIL corpus consequences: protocol table row for {stem!r} is "
                  f"{[cells[h] for h in hit] or 'absent'}, tracer_design.csv says {x.design}")
            bad += 1
    if (inferred, shared) != (int((td.design == "labeled").sum()), int((td.design == "shared").sum())):
        print(f"  FAIL corpus consequences: protocol table {inferred} inferred / {shared} shared, "
              f"tracer_design.csv {(td.design == 'labeled').sum()} / {(td.design == 'shared').sum()}")
        bad += 1
    if forms != u - 2:
        print(f"  FAIL corpus consequences: protocol table has {forms} forms, "
              f"units minus the two repeat 3 C records give {u - 2}")
        bad += 1
    need("si", f"outcome for the {W[forms]} distinct series and enzyme forms", "the form count")
    need("si", f"{W[inferred].capitalize()} of the {W[forms]} report secondary", "the inferred count")
    need("si", f"the remaining {W[shared]} need no inference", "the shared-tracer count")
    need("si", f"For the remaining {W[shared]} the question does not arise", "the shared-tracer count")
    need("main", f"across {W[inferred]} of the {W[forms]} forms", "the form count")
    need("main", f"For the remaining {W[shared]} the question does not arise", "the shared-tracer count")
    # the confidence procedure for temperature series (corrected 2026-09-14):
    # its simulated failure rate, Monte Carlo precision and assignment check
    import bounds_uncertainty as BU
    bu = pd.read_csv("../results/bounds_uncertainty.csv")
    import inspect
    reps = inspect.signature(BU.coverage_check).parameters["reps"].default
    old_rate, new_rate = BU.coverage_check()
    se = 100 * np.sqrt(old_rate * (1 - old_rate) / reps)
    need("main", f"true maximum in ${100 * old_rate:.0f}\\%$ of $2\\times10^{{4}}$ simulated "
                 f"trials (binomial standard error ${se:.1f}\\%$)",
         "the simulated failure rate of the old procedure, with its trial count")
    if reps != 20000:
        print(f"  FAIL corpus consequences: coverage_check runs {reps} trials, text says 2x10^4")
        bad += 1
    need("main", f"and the corrected bound in ${100 * new_rate:.1f}\\%$",
         "the simulated failure rate of the corrected procedure")
    mc = bu.lcb_mc_sd.max()
    e = int(np.floor(np.log10(mc)))
    need("main", f"largest Monte Carlo standard deviation ${mc / 10 ** e:.1f}\\times10^{{{e}}}$",
         "the Monte Carlo precision of the bounds")
    mv, gap = BU.assignment_sensitivity()
    need("si", f"moves a bound by at most ${mv:.4f}$, and every bound stays at least ${gap:.4f}$",
         "the ecDHFR assignment sensitivity")
    # the yeast bypass tolerance quoted in the main text, at the jointly solved r
    from completion import yadh_completion_joint, bypass_tolerance
    t13 = 100 * float(bypass_tolerance(7.13, 1.73, yadh_completion_joint(1.3)["r"], F0))
    t73 = 100 * float(bypass_tolerance(7.13, 1.73, yadh_completion_joint(7.3)["r"], F0))
    need("main", f"less than ${t13:.0f}\\%$ to ${t73:.0f}\\%$ of the", "the yeast bypass tolerance")
    # the Discussion now gives both thresholds; this one is against F0
    need("main", f"${t13:.0f}\\%$ to ${t73:.0f}\\%$ to reach the gated envelope",
         "the yeast bypass tolerance against F0 (Discussion)")
    # the profile-likelihood stratification must say what its own table shows
    # (until 2026-09-14 it said no adequately fitting series discriminates)
    from scipy import stats as _st
    pr = pd.read_csv("../results/offset_profiles.csv")
    pr["p_fit"] = _st.chi2.sf(pr.chi2_min, (2 * pr.n_T - 5).clip(lower=1))
    excl = pr[(pr.F_lo > 0) | (pr.F_hi < 0)]
    adq = excl[excl.p_fit > 0.05]
    above = adq[adq.F_lo > F0]
    WS = {3: "three", 5: "five", 7: "seven", 12: "twelve", 18: "eighteen"}
    need("si", f"profiles exclude $F=0$ in {WS[len(excl)]} of {WS[len(pr)]} series",
         "the profile exclusion count")
    need("si", f"{WS[len(excl) - len(adq)]} of those {WS[len(excl)]} belong to series the model",
         "the misfit exclusion count")
    need("si", f"The other {WS[len(adq)]} exclusions are in adequately fitting series",
         "the adequate exclusion count")
    need("si", f"{WS[len(above)]} ecDHFR and hsDHFR series", "the count wholly above F0")
    # the admitted E. coli thymidylate synthase record
    e = ya[ya.grp == "ectsase"].iloc[0]
    fe = F_min_exact(e.K_HT, e.K_DT)[0]
    se = np.sqrt((e.K_HT_se / e.K_HT) ** 2 + (GSC * e.K_DT_se / e.K_DT) ** 2
                 + 2 * GSC * (e.K_HT_se / e.K_HT) * (e.K_DT_se / e.K_DT))
    lh = (e.K_HT - 1) / (e.K_DT - 1)
    need("main", f"endpoint of ${fe:+.3f}$", "the ecTSase endpoint")
    need("main", f"bound of ${fe - 1.645 * se:.3f}$", "the ecTSase bound")
    need("si", f"$L_H={lh:.2f}>\\gSC$", "the ecTSase L_H")
    need("si", f"$F_{{\\mathrm{{obs}}}}={fe:+.3f}$", "the ecTSase endpoint")
    need("si", f"at $\\rho=-1$ is ${fe - 1.645 * se:.3f}$", "the ecTSase bound")
    # the yeast completion table, row by row: the prose beside it quotes F and
    # gamma, so a stale table cell would otherwise pass the value-somewhere test
    from completion import yadh_completion_joint
    for a in (1.3, 7.3):
        j = yadh_completion_joint(a)
        xh = f"{float(j['xH']):.2f}"
        xh = xh if len(xh) == 5 else f"\\phantom{{0}}{xh}"
        need("si", f"${a}$ & ${j['r']:.3f}$ & ${j['c']:.1f}$ & ${j['cD']:.1f}$ & ${xh}$ & "
                   f"${float(j['xD']):.3f}$ & ${j['gamma']:.2f}$ & ${j['F']:+.3f}$\\\\",
             f"the yeast completion table row a={a}")
    if not bad:
        print(f"  corpus consequences: Bonferroni over {u} units ({b0:+.3f}, {b1:+.3f}, "
              f"margin {b1 - F0:.3f}); {forms} protocol forms; all {n} records "
              f"monotone and forecast; ecTSase endpoint {fe:+.3f}, bound {fe - 1.645 * se:.3f}")
    return bad


def check_scope_corrections():
    r"""Numbers added when the scope of four results was corrected (2026-09-15).

    The joint binding-and-bypass line had been quoted without its r (it holds
    only at r = 1.31); the bounded-bypass endpoint E(phi) had been given for the
    interior branch, where it is only an upper bound; the envelope's convergence
    rate was stated as w^-2 for every K; and the reversible extension was said
    to open above for any equilibrium isotope effect.  Each replacement statement
    carries numbers, and each is recomputed here and matched to its sentence.
    """
    import joint_nuisance as J
    import network_geometry as NG
    from completion import (yadh_completion_joint, bypass_tolerance, yadh_joint_r,
                            check_joint_r_quadratic, mc_yadh_tolerance,
                            bsao_completion, bsao_masking_by_temperature)
    main = " ".join(DOCS["main"].read_text().split())
    si = " ".join(DOCS["si"].read_text().split())
    bad = 0

    def need(doc, text, what):
        nonlocal bad
        if text not in (main if doc == "main" else si):
            print(f"  FAIL scope corrections: {doc} lacks {what}: {text!r}")
            bad += 1

    # the joint line over the stipulated range of r
    ic, sl = J.envelope_r(0.0), J.envelope_r(0.0) - J.envelope_r(1.0)
    rhi = float(yadh_joint_r(1.3))
    fails, tight, where = _envelope_r_result()   # failures are counted by check_envelope
    need("main", f"$r$ stipulated anywhere from ${J.R_RANGE[0]:.2f}$ to ${rhi:.2f}$",
         "the stipulated r range")
    need("main", f"$F_{{\\mathrm{{bind}}}}\\le{ic:.3f}-{sl:.3f}\\,\\phi$", "the r-range joint line")
    need("si", f"\\ge\\ {ic:.3f}-{sl:.3f}\\,\\phi_{{\\mathrm H}}", "the r-range joint line")
    need("si", f"nonvacuous for $\\phi_{{\\mathrm H}}\\le{ic / sl:.3f}$", "where the line is nonvacuous")
    e = int(np.floor(np.log10(tight)))
    need("si", f"smallest slack of ${tight / 10 ** e:.1f}\\times10^{{{e}}}$ at "
               f"$\\phi_{{\\mathrm H}}={where:.3f}$", "the tightest slack of the r-range line")
    ce = J._endpoint_ab(1.0, 1.0, 1.0, 1.0, 0.14, rhi)
    need("si", f"the profiled endpoint there is ${ce:.4f}$", "the r = 1.57 counterexample")
    # the yeast tolerance against zero, now distinguished from the one against F0
    z13 = 100 * float(bypass_tolerance(7.13, 1.73, yadh_completion_joint(1.3)["r"], 0.0))
    z73 = 100 * float(bypass_tolerance(7.13, 1.73, yadh_completion_joint(7.3)["r"], 0.0))
    need("main", f"have to carry ${z13:.0f}\\%$ to ${z73:.0f}\\%$ of the isotope-sensitive rate to do so",
         "the yeast bypass tolerance against the semiclassical locus")
    t = mc_yadh_tolerance()
    need("si", f"intact with probability ${t['F0']['p10']:.2f}$, and one of a twentieth with "
               f"probability ${t['F0']['p05']:.3f}$", "the tolerance probabilities")
    # the bounded-bypass endpoint in both branches
    n, worst, interior, crit, n1 = NG.check_endpoint_exact(n=3000)
    if worst > 1e-7 or crit:
        print(f"  FAIL scope corrections: endpoint_exact differs by {worst:.1e}, "
              f"criterion disagreements {crit}")
        bad += 1
    e = int(np.floor(np.log10(worst)))
    need("si", f"On ${n}$ random admissible triples, ${interior}$ of them with an interior "
               f"minimum, it agrees with direct profiling to ${worst / 10 ** e:.0f}\\times10^{{{e}}}$",
         "the endpoint_exact validation")
    need("si", f"on all ${n1}$ of those with a single reference", "the r = 1 triple count")
    need("si", f"$E={NG.endpoint_closed(2.0, 1.5, 0.1):.3f}$ while the profiled endpoint is "
               f"${NG.endpoint_exact(2.0, 1.5, 0.1):.3f}$", "the interior-branch example")
    # the yeast r quadratic
    _, nan_mismatch, diff, prod = check_joint_r_quadratic()
    if nan_mismatch or diff > 1e-10:
        print("  FAIL scope corrections: the r quadratic disagrees with bisection")
        bad += 1
    need("si", f"is below ${np.ceil(prod * 1000) / 1000:.3f}$ over $2\\times10^{{5}}$",
         "the largest product of the r roots")
    # BSAO: the observed ratio is only a lower bound on r
    mm = bsao_masking_by_temperature()[1][3]
    dF = bsao_completion(mm, r=2.0)["F"] - bsao_completion(mm, r=1.1370)["F"]
    need("si", f"raising it to $2.0$ moves the point value at $m={mm:.3f}$ by ${dF:+.3f}$",
         "the BSAO sensitivity to r")
    # the envelope's convergence coefficient at K = 7
    d, tt = M.mu_ratios("C")[1:]        # already square-root mass ratios
    coef = (d - 1) / d * (np.log(7.0) - 0.5 * np.log(tt))
    need("si", f"with coefficient ${coef:.4f}$ at $K=7$", "the 1/w convergence coefficient")
    # the reversible corner path inside the Proposition S7 window
    H, K, E = 5.04, 1.65, 1.64
    EH = E ** GSC
    a, b, u, v = H - 1, K - 1, H / EH - 1, K / E - 1
    cf, cr = (v - u) / (a * v - b * u), (a - b) / (a * v - b * u)
    eps = 0.01
    CF, CR = (1 - eps) * cf, (1 - eps) * cr
    DH, DD = 1 + CF * (1 - H) + CR * (1 - H / EH), 1 + CF * (1 - K) + CR * (1 - K / E)
    Fr = np.log(H / DH) - GSC * np.log(K / DD)
    need("si", f"which is ${Fr:.2f}$ at $\\varepsilon=0.01$", "the reversible corner path")
    if not bad:
        print(f"  scope corrections: r-range line {ic:.3f} - {sl:.3f} phi (slack {tight:.1e}), "
              f"endpoint_exact on {n} triples, reversible corner {Fr:.2f}")
    return bad


def check_source_fidelity():
    r"""Statements corrected against the primary sources on 2026-09-15.

    A detailed audit read Cha et al. (1989), Klinman (1976), Kohen and Jensen
    (2002), Northrop and Duggleby (1990) and Kohen et al. (1999) against the text
    and found: the mixed-labeling counterexample used an anti-ordered pair of
    commitments that the design cannot have; Cha's Table 2 has a fourth, D/T-only
    determination the error discussion omitted; the commitment-precision
    percentages were a commitment-only budget presented as the whole; the
    combined multiplicity-and-correlation margin was not stated; an Arrhenius
    prefactor combination was called the offset; and a Monte Carlo fraction was
    printed as probability one.  Each replacement carries numbers, recomputed
    here and matched to its sentence, and the deposit README must carry the
    current corpus counts rather than an earlier release's.
    """
    from scipy.stats import norm
    import corpus
    import mixed_label as ML
    from completion import isotope_error_budget, _mc_yadh
    from pathlib import Path
    main = " ".join(DOCS["main"].read_text().split())
    si = " ".join(DOCS["si"].read_text().split())
    bad = 0

    def need(doc, text, what):
        nonlocal bad
        if text not in (main if doc == "main" else si):
            print(f"  FAIL source fidelity: {doc} lacks {what}: {text!r}")
            bad += 1

    # combined error-inflation margin: Bonferroni over the units, at rho = -1
    ya = pd.read_csv("../data/cha1989_yadh.csv")
    avg = ya[ya.note.str.contains("average")].iloc[0]
    u = corpus.counts()["units"]
    Fh = np.log(avg.K_HT) - GSC * np.log(avg.K_DT)
    s1 = avg.K_HT_se / avg.K_HT + GSC * avg.K_DT_se / avg.K_DT
    k = (Fh - F0) / (norm.ppf(1 - 0.05 / u) * s1)
    need("main", f"or by ${100 * (k - 1):.0f}\\%$ once the multiplicity adjustment", "the combined margin")
    need("si", f"an understatement of only ${100 * (k - 1):.1f}\\%$ brings the bound down to $\\Fz$",
         "the combined margin")
    # Cha Table 2: four D/T determinations, three paired
    det = ya[ya.note.str.contains("determination")]
    d4 = pd.read_csv("../data/cha1989_yadh_dt_only.csv")
    dt = np.r_[det.K_DT.to_numpy(), d4.K_DT.to_numpy()]
    need("si", f"the fourth, ${d4.K_DT.iloc[0]:.2f}\\pm{d4.K_DT_se.iloc[0]:.2f}$, has no H/T partner",
         "the fourth D/T determination")
    need("si", f"the three H/T determinations, ${det.K_HT.std(ddof=1) / np.sqrt(len(det)):.3f}$",
         "the H/T standard error")
    need("si", f"the three paired D/T determinations, ${det.K_DT.std(ddof=1) / np.sqrt(len(det)):.3f}$",
         "the paired D/T standard error")
    need("si", f"exceeds that of all four, ${dt.std(ddof=1) / np.sqrt(len(dt)):.3f}$",
         "the four-determination D/T standard error")
    # the isotope-error budget of the commitment design
    b3, b5 = (isotope_error_budget(avg.K_HT, avg.K_DT, avg.K_HT_se, avg.K_DT_se,
                                   m * (avg.K_HT - 1), 1.3107) for m in (3, 5))
    for doc in ("main", "si"):
        need(doc, f"${b3[0]:.3f}$ at $\\rho=0$ and ${b3[1]:.3f}$ at $\\rho=-1$", "the isotope-error budget at m = 3")
    need("si", f"${b5[0]:.3f}$ and ${b5[1]:.3f}$ at $m=5$", "the isotope-error budget at m = 5")
    need("si", f"cut to ${100 * abs(F0) / 2 / b3[0]:.0f}\\%$ ($\\rho=0$) or ${100 * abs(F0) / 2 / b3[1]:.0f}\\%$ ($\\rho=-1$)",
         "the isotope precision required")
    need("main", f"already above $|\\Fz|/2={abs(F0) / 2:.3f}$", "the design target")
    # mixed labeling: ordered maps deflate, anti-ordered inflate
    n, nbad, worst, lifted = ML.check_ordered()
    if nbad:
        print(f"  FAIL source fidelity: {nbad} ordered maps raise F_obs above F_int")
        bad += 1
    need("si", f"$c_{{\\mathrm D}}=5$ return an observed exponent of ${ML.observed_exponent(1.10, 4.8, 1.0, 5.0):.2f}$",
         "the ordered mixed-label example")
    need("si", f"$c_{{\\mathrm D}}=1$ return ${ML.observed_exponent(1.10, 4.8, 5.0, 1.0):.2f}$",
         "the anti-ordered mixed-label example")
    need("si", f"Over $4\\times10^{{5}}$ random ordered maps", "the ordered-map count")
    if n != 400000:
        print(f"  FAIL source fidelity: check_ordered draws {n}, text says 4x10^5")
        bad += 1
    # the Arrhenius prefactor combination, which is not the offset
    FA = np.log(4.3) - GSC * np.log(1.73)
    sFA = np.hypot(0.6 / 4.3, GSC * 0.26 / 1.73)
    need("si", f"is ${FA:.2f}\\pm{sFA:.2f}$, but it is not the offset", "the prefactor combination")
    # zero nonpositive completed offsets, and the 95% bound that implies
    zero = all(int(np.sum(_mc_yadh(400000, 20260815, a)[0] <= 0)) == 0
               for a in (1.3, 2.3, 7.3, (1.3, 7.3)))
    if not zero:
        print("  FAIL source fidelity: a completed yeast offset at or below zero was drawn")
        bad += 1
    e = 3 / 400000
    ex = int(np.floor(np.log10(e)))
    need("si", f"probability below ${e / 10 ** ex:.1f}\\times10^{{{ex}}}$ at $95\\%$ confidence",
         "the zero-count bound")
    # the deposit README carries the current corpus, not an earlier release's
    c = corpus.counts()
    readme = " ".join((Path(__file__).resolve().parent.parent / "README.md").read_text().split())
    for stale in ("94 matched", "29 independent systems", "28 systems fall",
                  "Singled out from 29", "Bonferroni across 29", "completes bovine serum"):
        if stale in readme:
            print(f"  FAIL source fidelity: README still says {stale!r}")
            bad += 1
    for want in (f"{c['records']} matched", f"{c['units']} analysis units"):
        if want not in readme:
            print(f"  FAIL source fidelity: README lacks {want!r}")
            bad += 1
    sources = (Path(__file__).resolve().parent.parent / "external_data" / "SOURCES.md").read_text()
    if "completes bovine serum amine oxidase" in sources:
        print("  FAIL source fidelity: SOURCES.md still lists the withdrawn BSAO completion")
        bad += 1
    if not bad:
        print(f"  source fidelity: combined margin {100 * (k - 1):.1f}%, D/T SE of four "
              f"{dt.std(ddof=1) / 2:.3f}, isotope budget {b3[0]:.3f}/{b3[1]:.3f}, ordered maps "
              f"{nbad} of {n}, F_A {FA:.2f}")
    return bad


def check_readiness_scope():
    r"""Scope corrections from a readiness review (2026-09-15).

    The binding shift was stated as exact on any branch (it holds on the open
    branch only); the single-record bound was called exact at rho = -1 (it is a
    plug-in bound, exact only for known log-scale standard deviations); the
    Conclusions still called a model-based commitment estimate a measurement;
    and the Smedarchina-Siebrand comparison credited their calculation with the
    value of F0.  The numbers that replace those statements are recomputed here,
    and the retired phrasings must not reappear.
    """
    from scipy.stats import t as tdist
    import corpus
    import bounds_uncertainty as BU
    from partial_id import F_min_binding, F_min_exact as _fme
    main = " ".join(DOCS["main"].read_text().split())
    si = " ".join(DOCS["si"].read_text().split())
    bad = 0

    def need(doc, text, what):
        nonlocal bad
        if text not in (main if doc == "main" else si):
            print(f"  FAIL readiness scope: {doc} lacks {what}: {text!r}")
            bad += 1

    def forbid(doc, text, what):
        nonlocal bad
        # case-insensitive: a retired claim must not return at the start of a sentence
        if text.lower() in (main if doc == "main" else si).lower():
            print(f"  FAIL readiness scope: {doc} still says {what}: {text!r}")
            bad += 1

    # binding on the closed branch
    args = (1.05, 1.05, 0.95, 0.95)
    e0, e1 = _fme(3.0, 1.8)[0], F_min_binding(3.0, 1.8, *args)
    need("main", f"still lower the endpoint by ${e0 - e1:.3f}$", "the closed-branch binding shift")
    need("si", f"the endpoint falls from ${e0:.3f}$ to ${e1:.3f}$", "the closed-branch binding example")
    if abs(F_min_binding(7.13, 1.73, *args) - _fme(7.13, 1.73)[0]) > 1e-12:
        print("  FAIL readiness scope: yeast endpoint moves under an F_bind = 0 perturbation")
        bad += 1
    # plug-in coverage of the single-record bound
    y = BU.plugin_coverage(7.13, 0.07, 1.73, 0.02, 0.05)
    h = BU.plugin_coverage(2.0, 0.2, 1.1, 0.001, 0.05)
    need("main", f"fails in ${100 * y:.1f}\\%$ of trials for the yeast errors", "the yeast coverage")
    need("main", f"but in ${100 * h:.1f}\\%$ for a hypothetical pair", "the counterexample coverage")
    # t-tail sensitivity of the multiplicity-adjusted yeast bound
    u = corpus.counts()["units"]
    F = np.log(7.13) - GSC * np.log(1.73)
    s1 = 0.07 / 7.13 + GSC * 0.02 / 1.73
    tb = F - tdist.ppf(1 - 0.05 / u, 10) * s1
    need("main", f"bound at $\\rho=-1$ falls to ${tb:.3f}$, below $\\Fz$", "the t-tail sensitivity")
    # house style: no em dashes in the manuscript or the supplement
    for doc in ("main", "si"):
        raw = DOCS[doc].read_text()
        n_em = raw.count("---") + raw.count("\\textemdash") + raw.count("\u2014")
        if n_em:
            print(f"  FAIL readiness scope: {doc} contains {n_em} em dash(es)")
            bad += n_em
    # masking direction: an intrinsic pair ABOVE the reference can be carried below it
    xD, Fi, c = 2.0, 0.10, 1.0
    xH = xD ** GSC * np.exp(Fi)
    K = lambda x: x * (1 + c) / (x + c)
    Fo = float(np.log(K(xH)) - GSC * np.log(K(xD)))
    need("main", f"an intrinsic offset of ${Fi:+.3f}$ at", "the masking-direction example")
    need("main", f"is observed as ${Fo:.3f}$ at a commitment of $c={c:.0f}$",
         "the masking-direction example")
    # retired phrasings
    forbid("main", "toward the semiclassical reference, never past it", "masking never passing the reference")
    forbid("main", "past the semiclassical reference\nat all", "masking never passing the reference")
    forbid("main", "monotone observation map with one isotope-sensitive step and a commitment",
           "the half-line under monotonicity alone")
    forbid("main", "endpoint moves by exactly", "the unconditional binding shift")
    forbid("main", "however large", "binding effects harmless however large")
    forbid("si", "moves by exactly $F_{\\mathrm{bind}}$", "the unconditional binding shift")
    forbid("main", "commitment measurements already published", "a commitment estimate as a measurement")
    forbid("main", "has exact coverage under that sampling model", "exact coverage of the plug-in bound")
    forbid("si", "property of a published calculation", "F0 credited to the instanton calculation")
    forbid("main", "settled by the mechanism assumed, not by the precision achieved",
           "discrimination settled by topology alone")
    if not bad:
        print(f"  readiness scope: closed-branch binding shift {e0 - e1:.3f}, plug-in coverage "
              f"{100 * y:.1f}%/{100 * h:.1f}%, t10 Bonferroni bound {tb:.3f}")
    return bad


def check_generated_tables():
    """Generated table inputs beside the manuscript must match the deposit's outputs.

    Scripts write these files to results/ and they are copied into the manuscript
    tree by its build; nothing re-derives them there.  A regenerated caption or
    row therefore reaches the PDF only if the copy is refreshed, and until
    2026-09-16 two corrected captions in si_extra_tables.tex sat in results/ while
    the supplement still printed the superseded ones.  Clean builds and the
    numerical guards both miss that: the audit reads si_body.tex, not the files it
    \\input.
    """
    from pathlib import Path as _P
    res = _P(__file__).resolve().parent.parent / "results"
    bad = 0
    names = [f.name for f in sorted(res.glob("sm_table_*.tex"))] + ["si_extra_tables.tex"]
    for name in names:
        src = res / name
        if not src.exists():
            continue
        for dst in (ROOT / name, ROOT / "si" / name):
            if not dst.exists():
                continue
            if dst.read_text() != src.read_text():
                print(f"  FAIL generated tables: {dst} differs from results/{name}")
                bad += 1
    if not bad:
        print(f"  generated tables: {len(names)} inputs match the deposited outputs")
    return bad


def check_titles():
    """The supplement must carry the article's title.

    The article and its supplement are built from separate sources, so a
    retitled article leaves the supplement on the old title with nothing to
    flag it.  That has happened before.
    """
    def _title(p):
        return _titled(p.read_text())

    main = _title(DOCS["main"])
    _si_doc = MR.one(ROOT / "si", "*_si.tex")
    si = _title(_si_doc)
    ok = main.lower() == si.lower()          # wording must match; case may follow house style
    print(f"  article title : {main[:58]}...")
    print(f"  supplement    : {'matches' if ok else 'DIFFERS: ' + si[:58]}")
    return 0 if ok else 1


HAVE_DOCS = all(p.exists() for p in DOCS.values())


def _texts():
    """Manuscript sources, or empty when they are not present.

    The public deposit carries the analysis, the data and the results, but not
    the manuscript: the article is not ours to redistribute before publication.
    Most guards here are cross-checks between a computed quantity and the text
    that quotes it, so they cannot run from the archive alone.  Rather than
    fail with a traceback -- which is what a reviewer downloading the archive
    used to get -- the suite reports which checks it can and cannot perform.
    """
    if not HAVE_DOCS:
        return {k: "" for k in DOCS}
    return {k: p.read_text() for k, p in DOCS.items()}


# Guards that compare a computed quantity against the manuscript text.  Without
# the manuscript they have nothing to compare against and are skipped, loudly.
_NEEDS_DOCS = {
    "check_si_figrefs", "check_si_eqrefs", "check_si_pointers",
    "check_general_census", "check_corpus", "check_provenance_split",
    "check_cited_scripts", "check_no_ai_mentions", "check_methods_sources",
    "check_precision_limited_count", "check_ladh_zero_claim",
    "check_repo_metadata", "check_titles", "check_reference_asymmetry",
    "check_promoting_mode", "check_corpus_consequences", "check_scope_corrections",
    "check_source_fidelity", "check_readiness_scope", "check_generated_tables",
}


def run():
    # A target named by the caller must be readable.  Without this check,
    # KIE_MANUSCRIPT pointed at a tree whose files are named differently falls
    # through to the deposit path below: thirteen guards print SKIPPED, nothing
    # fails, and the suite exits 0 -- which reads as a pass on documents that
    # were never opened.  The graceful path exists for the public archive,
    # where no manuscript is expected; it must not cover a requested target.
    _requested = MR.requested()
    if _requested and not HAVE_DOCS:
        missing = "\n  ".join(str(p) for p in DOCS.values() if not p.exists())
        sys.stderr.write(
            f"manuscript directory {_requested}\nbut these documents are missing:\n"
            f"  {missing}\nRefusing to report a pass on documents not read.\n")
        return 2
    texts = _texts()
    bad = 0
    if HAVE_DOCS:
        print(f"auditing: {DOCS['main']}")
        print(f"          {DOCS['si']}\n")
    else:
        print("manuscript sources not found; auditing the deposit alone.")
        print(f"  looked in: {ROOT}")
        print("  Set KIE_MANUSCRIPT to a manuscript directory to cross-check the")
        print("  text as well.  Everything below is computed from the deposited")
        print("  data and code, and is verified in full.\n")
    print(f"{'quantity':32s} {'value':>11s}  documents")
    for lab, val, dp, docs in derived():
        s = f"{abs(val):.{dp}f}"
        if not HAVE_DOCS:
            # nothing to compare against; report the computed value and move on
            print(f"  {lab:30s} {val:+11.6f}  (computed)")
            continue
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
            # A CONFLICTING COPY is invisible to the containment test above:
            # if the correct value appears once, a wrong copy of the same
            # quantity elsewhere passes.  That is exactly how 3.26278 sat beside
            # 3.26281 in the supplement.  Look for numbers written to the same
            # precision that agree with the correct value to all but the last
            # two decimals and yet differ from it.  Only at >= 5 decimals, where
            # agreement to three places is meaningful and distinct quantities
            # (3.34887 vs 3.34278) do not collide.
            if dp >= 5:
                head = s[: s.index(".") + 1 + (dp - 2)]
                pat = rf"(?<![\d.]){re.escape(head)}\d{{2}}(?![\d])"
                clash = sorted({m for m in re.findall(pat, t) if m != s})
                if clash:
                    print(f"  {'':30s} {'':11s}  {d}: CONFLICTING COPY {clash}")
                    bad += 1
    for _fn in (check_si_figrefs, check_si_eqrefs, check_si_pointers,
                check_general_census, check_corpus, check_provenance_split,
                check_cited_scripts, check_no_ai_mentions,
                check_methods_sources, check_precision_limited_count,
                check_ladh_zero_claim, check_ecdhfr_assignment,
                check_envelope, check_joint_counterexamples,
                check_repo_metadata, check_titles, check_reference_asymmetry,
                check_promoting_mode, check_corpus_consequences,
                check_scope_corrections, check_source_fidelity,
                check_readiness_scope, check_generated_tables):
        if not HAVE_DOCS and _fn.__name__ in _NEEDS_DOCS:
            print(f"  SKIPPED (needs the manuscript): {_fn.__name__}")
            continue
        bad += _fn()
    if HAVE_DOCS:
        print(f"\n{'FAIL' if bad else 'PASS'}: {bad} derived quantities "
              f"missing or contradicted")
    else:
        print(f"\n{'FAIL' if bad else 'PASS'}: {bad} failures in the checks that "
              f"run without the manuscript.")
        print("  The text cross-checks above were skipped, not failed.")
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(run())
