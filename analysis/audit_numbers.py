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
    # The text gives this one as a percentage ("to 13.2%").  Checked as 0.132
    # until 2026-10-15, it passed only by matching the unrelated -0.1324 of the
    # horse liver table; when that table moved to a generated file it failed.
    out.append(("bypass phi* r=1.31 vs 0 (%)",
                100 * bypass_to_destroy(7.13, 1.73, 0.0, 1.31), 1, ("si",)))
    out.append(("E_r(0.15) r=1.31", endpoint_closed(7.13, 1.73, 0.15, 1.31), 4, ("si",)))
    # The envelope threshold: mu_ratios() returns SQUARE-ROOT mass ratios, so the
    # code's t is sqrt(mu_T/mu_H) and the threshold sqrt(t) is (mu_T/mu_H)^(1/4).
    # Reading mu_ratios as a plain mass ratio once led to "fixing" the exponent
    # from 1/4 to 1/2, which is wrong; this pins the printed value.
    # binding tolerance (Proposition S7): F_bind below which the exclusion survives
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
CORPUS = {"records": 231, "units": 73}


GENERAL_CENSUS = {"networks": 764, "edges": 5538, "mixed": 3415,
                  "inverse": 1026, "blind": 626, "concave": 471}


SI_EQREFS = {
    2: "eq:map",                          # the observation map
}


SI_FIGREFS = {
    # accepted phrasings for the endpoints figure; a pointer check that admits
    # only one wording silently stops checking when the caption is reworded.
    # Renumbered 2026-10-14 when the overview and envelope figures were added
    # ahead of these two (atlas 2 -> 3, compatible sets 3 -> 5).
    3: ("Topology classifies what can be identified",
        "Reaction pathways constrain the possible bounds"),
    5: ("Identified sets for the primary",
        "Compatible sets for the primary"),
    # the joint binding-bypass contour table takes its cells from Fig. 6 data
    6: ("A bound becomes a measured offset",),
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
        m = re.search(r"\\caption\{(.{0,220})", fig, re.S)
        if m:
            caps.append(m.group(1))
    si = DOCS["si"].read_text()
    bad = 0
    # both phrasings: "Fig.~N of the main text" and "main-text Fig.~N"; the
    # second escaped this guard until 2026-10-14
    for n in sorted(set(int(m) for m in
                        re.findall(r"Fig\.~(\d+) of the main text", si)
                        + re.findall(r"main-text Fig\.~(\d+)", si))):
        want = SI_FIGREFS.get(n)
        got = caps[n - 1].replace("\n", " ") if n <= len(caps) else "<out of range>"
        ok = want is not None and any(w in got for w in want)
        print(f"  SI cites main-text Fig. {n}: {got[:52].strip()!r} "
              f"{'OK' if ok else 'MISMATCH'}")
        if not ok:
            bad += 1
    # The number-to-caption map above cannot tell that a reference meant for one
    # figure now names another valid figure: after the 2026-10-14 renumbering, a
    # stale "Fig.~3" in a sentence about endpoints would still match the atlas
    # caption.  Each reference is therefore also tied, by its own sentence, to
    # the label of the figure it intends, and every reference must be covered.
    order = [m for fig in re.findall(r"\\begin\{figure\*?\}(.*?)\\end\{figure\*?\}",
                                       main, re.S)
             for m in re.findall(r"\\label\{(fig:[^}]*)\}", fig)]
    flat = " ".join(si.split())
    covered = 0
    for pat, label in SI_FIGCONTEXT:
        hits = re.findall(pat, flat)
        want = order.index(label) + 1 if label in order else None
        if len(hits) != 1 or want is None or int(hits[0]) != want:
            print(f"  FAIL SI figure context {pat!r}: found {hits}, "
                  f"{label} is Fig. {want}")
            bad += 1
        covered += len(hits)
    total = len(re.findall(r"Fig\.~\d+ of the main text|main-text Fig\.~\d+", flat))
    if covered != total:
        print(f"  FAIL SI figure context: {total - covered} main-text figure "
              f"reference(s) not tied to a figure label")
        bad += 1
    return bad


# Each SI pointer to a main-text figure, by the sentence that carries it, and
# the label of the figure it intends.  Numbers follow from the label order.
SI_FIGCONTEXT = [
    (r"mechanisms drawn in Fig\.~(\d+) of the main text", "fig:atlas"),
    (r"mechanisms of main-text Fig\.~(\d+)\.", "fig:atlas"),
    (r"Fig\.~(\d+) of the main text report endpoints", "fig:asym"),
    (r"plotted in Fig\.~(\d+) of the main text", "fig:asym"),
    (r"the data of main-text Fig\.~(\d+)\)", "fig:design"),
]


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
    """The provenance counts in both documents must match the benchmark's DOIs.

    Four classes: parsed from a publisher PDF, recovered from a published
    spreadsheet, transcribed from a PDF supplement whose text layer the
    transcription audit checks (the sources admitted by the census of
    2026-10-20), and transcribed by hand from scanned pages or word-processor
    supplements.  Quoting a larger audited count than the data support would
    overstate the audit.  The classes are stated here independently and must
    also equal the builder's own sets.
    """
    import pandas as pd
    import build_trinomial as BT
    PDF = {"10.1021/ja411998h", "10.1021/acs.biochem.1c00558",
           "10.1021/acscatal.9b03345", "10.1021/ja501936d", "10.1021/ja031683w",
           "10.3390/ijms16047304"}
    XLSX = {"10.1371/journal.pone.0196506"}
    AUDITED = {"10.1021/ja103010b", "10.1021/ja307859m", "10.1021/ja400761x",
               "10.1021/ja060196o", "10.1021/ja504341g", "10.1021/acscatal.5b01332",
               "10.1021/acscatal.8b02554", "10.1039/c0cc02988b", "10.1021/acscatal.5b00331",
               "10.1073/pnas.1415940111", "10.1021/acs.biochem.9b00460",
               "10.1021/ja902120t", "10.1021/acs.biochem.6b00181", "10.1021/jacs.7b08359"}
    HAND = {"10.1021/bi00253a026", "10.1021/bi002075l", "10.1021/bi00161a030",
            "10.1021/bi0518242", "10.1098/rstb.2006.1871"}
    bad = 0
    if set(BT.PDF_SOURCES) != PDF or set(BT.TRANSCRIBED) != AUDITED or BT.HAND_TRANSCRIBED != HAND:
        print("  FAIL provenance: the classes here differ from build_trinomial.py's")
        bad += 1
    d = pd.read_csv("../data/trinomial_benchmark.csv")
    n_pdf = int(d.source_DOI.isin(PDF).sum())
    n_xls = int(d.source_DOI.isin(XLSX).sum())
    n_aud = int(d.source_DOI.isin(AUDITED).sum())
    n_hand = int(d.source_DOI.isin(HAND).sum())
    parsed = n_pdf + n_xls
    if parsed + n_aud + n_hand != len(d):
        print(f"  FAIL provenance: {parsed}+{n_aud}+{n_hand} != {len(d)}; "
              "a source DOI is unclassified")
        bad += 1
    by_hand = _corpus_records() - parsed - n_aud
    si = " ".join(DOCS["si"].read_text().split())
    main = " ".join(DOCS["main"].read_text().split())
    for label, need, doc, text in (
            ("SI parsed total", f"${parsed}$ of the ${len(d)}$ benchmark records", "si", si),
            ("SI value count", f"${parsed * 4}$ values", "si", si),
            ("SI pdf values", f"{n_pdf * 4} values", "si", si),
            ("SI audited records", f"${n_aud}$ records, every one of the ${n_aud * 4}$ values", "si", si),
            ("SI audit total", f"the audit covers ${(n_pdf + n_aud) * 4}$ values", "si", si),
            ("SI series by hand", f"{_nw(n_hand)} temperature-resolved records were transcribed by hand", "si", si),
            # since 2026-10-21 the main text keeps the accounting in Methods only
            ("methods parsed total", f"{parsed} records in all", "main", main),
            ("methods audited", f"A further {n_aud} records were transcribed from supplementary PDFs", "main", main),
            ("methods audited values", f"each of their {n_aud * 4} values", "main", main),
            ("methods hand", f"The remaining {by_hand} were transcribed by hand", "main", main)):
        if need not in text:
            print(f"  FAIL provenance: {doc} does not state {need!r} ({label})")
            bad += 1
    if not bad:
        print(f"  provenance: {n_pdf} PDF ({n_pdf*4} values) + {n_xls} spreadsheet "
              f"({n_xls*4}) = {parsed} parsed, {n_aud} audited transcriptions "
              f"({n_aud*4} values), {by_hand} by hand")
    return bad


def _nw(n):
    """A count in words, as the documents write counts below one hundred."""
    ones = ("zero one two three four five six seven eight nine ten eleven twelve thirteen "
            "fourteen fifteen sixteen seventeen eighteen nineteen").split()
    tens = "_ _ twenty thirty forty fifty sixty seventy eighty ninety".split()
    if n < 20:
        return ones[n]
    t, o = divmod(n, 10)
    return tens[t] + (f"-{ones[o]}" if o else "")


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
    """Every source of every corpus record must be cited in the main text.

    Adding the two aqueous 3 C records once brought TsaiKlinman2001 into the
    corpus but not into the Methods source list, so the provenance of two records
    went uncited.  Since the census of 2026-10-20 the main text's data section
    cites every source; this resolves each single-condition record's source key
    and each series' DOI (through the bibliography) and requires a \\cite of it.
    """
    import re
    import corpus
    main = DOCS["main"].read_text()
    cited = set()
    for grp in re.findall(r"\\cite\{([^}]*)\}", main):
        cited |= {k.strip() for k in grp.split(",")}
    bib = (ROOT / "references.bib").read_text()
    doi2key = {}
    # entries close in more than one style ("}" on its own line or "}}"), so
    # split at each entry rather than match its closing brace
    for chunk in re.split(r"\n(?=@)", bib):
        k = re.match(r"\s*@\w+\{([^,]+),", chunk)
        m = re.search(r"\bdoi\s*=\s*\{([^}]*)\}", chunk)
        if k and m:
            doi2key[m.group(1).strip().lower()] = k.group(1).strip()
    needed = set(corpus.single_condition().source)
    missing_doi = []
    for doi in corpus.series().source_DOI.unique():
        k = doi2key.get(doi.lower())
        if k is None:
            missing_doi.append(doi)
        else:
            needed.add(k)
    missing = sorted(needed - cited)
    if missing or missing_doi:
        print(f"  FAIL methods: corpus sources not cited in the main text: {missing}"
              f"{'; series DOIs with no bibliography entry: ' + str(missing_doi) if missing_doi else ''}")
        return 1
    print(f"  data sources: all {len(needed)} sources of the corpus are cited in the main text")
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
    above = sum(1 for pt, l in pts if pt > F0)          # precision-limited + excluding
    word, rest = _nw(n), _nw(above)
    # each entry is a tuple of accepted phrasings: the count must be stated, in
    # one wording or another, in every place the documents state it (the figure
    # caption, the results paragraph, the abstract and the conclusions; the SI's
    # Agrawal paragraph)
    need = {"main": ((f"{word} fail on precision alone",),
                     (f"for {word} of the remaining {rest} the endpoint",),
                     (f"{word} more are limited by precision",),
                     (f"and {word} more are limited by precision",)),
            "si": ((f"one of the {word} units whose endpoint exceeds",),)}
    bad = 0
    for tag, groups in need.items():
        t = " ".join(DOCS[tag].read_text().split()).lower()
        for variants in groups:
            if not any(ph in t for ph in variants):
                print(f"  FAIL precision-limited: {tag} lacks '{variants[0]}' "
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
    # every aqueous horse liver record, as adh_primary.py and the SI table hold them
    from adh_primary import ladh_records
    d = ladh_records()
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
        want = f"{_nw(n_zero)} of the {_nw(n_tot)}"
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
    """The ecDHFR light-enzyme series must come from its original source.

    Ref. Wang2014 reprints the light series beside its heavy one, rounded and
    distinguished only by color, so until 2026-10-15 the unpaired 10 C row was
    assigned to the light series by interpolation and three values carried the
    reprint's rounding.  The light series is now parsed from Sikorski et al.
    2004, SI Table 1, which reports that series alone, and the transcription
    audit finds every value there; this check fails if a light record is
    sourced elsewhere or if the heavy series acquires a 10 C row.
    """
    import pandas as pd
    d = pd.read_csv("../data/trinomial_benchmark.csv")
    light = d[(d.family == "ecDHFR") & (d.variant == "light enzyme")]
    heavy = d[(d.family == "ecDHFR") & (d.variant == "heavy enzyme")]
    bad = 0
    if set(light.source_DOI) != {"10.1021/ja031683w"} or len(light) != 6:
        print(f"  FAIL ecDHFR: light series must be Sikorski 2004's six records, "
              f"found {len(light)} from {sorted(set(light.source_DOI))}")
        bad += 1
    if 10.0 in set(heavy.T_C) or 10.0 not in set(light.T_C):
        print("  FAIL ecDHFR: the 10 C record is the light enzyme's (Sikorski SI Table 1)")
        bad += 1
    if not bad:
        print("  ecDHFR light series: six records from Sikorski 2004 SI Table 1, "
              "10 C included")
    return bad


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
    # The phrasings since the census of 2026-10-20 (abstract, the data section's
    # opening count and unit count, and the SI's data note).
    for d, phrases in (("main", (f"all {R} eligible matched records, in {U} analysis units",
                                 f"The census yields {R} matched records",
                                 f"single-condition record is one of {U} analysis units")),
                       ("si", (f"the ${R}$ matched records form {U} analysis units",))):
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
    # a version number ends at its last digit: "v1.8.6." at the end of a
    # sentence once read as "1.8.6." and failed against a matching CITATION.cff
    cited = set(re.findall(r"version analyzed here is v(\d+(?:\.\d+)+)|version v(\d+(?:\.\d+)+)",
                           " ".join(DOCS["main"].read_text().split())))
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
    fam = {"YADH": "YADH", "LADH": "LADH", "BSAO": "BSAO", "MAOB": "MAOB", "ecDHFR": "ecDHFR"}

    def expect(family, variant):
        lo, hi = RA.r_range(family, variant)
        a, b = f"{lo:.2f}", f"{hi:.2f}"
        return a if a == b else f"{a}--{b}"

    si = DOCS["si"].read_text()
    rows = re.findall(r"^(YADH|LADH|BSAO|MAOB|ecDHFR) (.+?) & [^&\n]+ & yes & ([\d.\-]+) & inferred",
                      si, re.M)
    if len(rows) != d.groupby(["family", "variant"]).ngroups:
        print(f"  FAIL reference asymmetry: protocol table has {len(rows)} inferred "
              f"rows, data hold {d.groupby(['family', 'variant']).ngroups} forms")
        bad += 1
    for f, v, cell in rows:
        want = expect(fam[f], v.replace("\\%", "%"))   # LaTeX-escaped solvent rows
        if cell != want:
            print(f"  FAIL reference asymmetry: protocol table {f} {v} says {cell}, "
                  f"data give {want}")
            bad += 1

    bypass = {"Yeast alcohol dehydrogenase": ("YADH", "wild type", 25.0),
              "Horse liver F93W, $25\\,^{\\circ}$C": ("LADH", "F93W", 25.0),
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

    # the light-enzyme row carries R from Sikorski's secondary effects; it was
    # tabulated as a shared tracer, r = 1, until 2026-10-15
    m = re.search(r"DHFR light enzyme, \$25\\,\^\{\\circ\}\$C & \$[\d.]+\$ & \$?([\d.\-]+)\$?", si)
    want = f"{RA.r_at('ecDHFR', 'light enzyme', 25.0):.2f}"
    if m is None or m.group(1) != want:
        print(f"  FAIL reference asymmetry: bypass table DHFR light r = "
              f"{m.group(1) if m else 'absent'}, secondary effects give {want}")
        bad += 1
    bench = {"YADH (Cha 1989)": ("YADH", "wild type", 25.0),
             "BSAO (Grant 1989)": ("BSAO", "wild type", 25.0),
             "MAO-B pH 6.1, 10 C": ("MAOB", "pH 6.1", 10.0),
             "LADH F93W": ("LADH", "F93W", 25.0),
             "ecDHFR light, 25 C": ("ecDHFR", "light enzyme", 25.0)}
    for lab, _, _, r in NG.BENCH:
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
          30: "thirty", 31: "thirty-one", 32: "thirty-two", 33: "thirty-three",
          34: "thirty-four", 35: "thirty-five"}

    class _Words(dict):
        # a count with no word is written as digits, so an unexpected count
        # fails the phrase test instead of crashing the whole suite
        def __missing__(self, k):
            return str(k)
    W = _Words(_W)
    main = " ".join(DOCS["main"].read_text().split())
    si = " ".join(DOCS["si"].read_text().split())
    bad = 0

    def need(doc, text, what, alternatives=()):
        nonlocal bad
        if not any(candidate in (main if doc == "main" else si)
                   for candidate in (text, *alternatives)):
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
    need("main", f"applied uniformly to all {W[u]} units", "the unit count")
    # the monotonicity condition of E(phi), and the Theorem 1 forecast
    # over the admissible records (K_HT > K_DT > 1); since 2026-10-13 one
    # record (wild type in 50% methanol at -50 C, K_DT = 0.99) is not, and it
    # has no compatible set under the map
    rows = pd.concat([corpus.series(), corpus.single_condition()])
    rows = rows[(rows.K_HT > rows.K_DT) & (rows.K_DT > 1)]
    n_adm = len(rows)
    holds = sum(1 for _, r in rows.iterrows()
                if (r.K_HT - 1) - GSC * (r.K_DT - 1) < (r.K_HT - 1) * (r.K_DT - 1) * (GSC - 1))
    if holds == n_adm:
        need("si", f"holds for all ${n_adm}$ admissible records here",
             "the E(phi) monotonicity count")
    else:
        print(f"  FAIL corpus consequences: E(phi) monotone for {holds} of {n_adm}")
        bad += 1
    chk, viol = structural_prediction_check()
    if viol or chk != n_adm:
        print(f"  FAIL corpus consequences: Theorem 1 forecast {chk} checked, {viol} violations")
        bad += 1
    need("si", f"Checked against all ${chk}$ admissible records, there are no violations",
         "the Theorem 1 forecast count")
    # the protocol table: one row per distinct series or form.  Since the census
    # of 2026-10-20 it is a longtable with three statuses: inferred (labeled
    # design, secondary effects reported), assumed (labeled design, none
    # reported) and shared (one tritiated isotopologue, r = 1 exactly).
    tab = DOCS["si"].read_text()
    tab = tab[tab.index(r"\label{tab:protocol}"):]
    tab = tab[:tab.index(r"\end{longtable}")]
    lines = [ln.rstrip() for ln in tab.splitlines()]
    STATUS = {r"inferred$^{a}$\\": "inferred", r"assumed$^{c}$\\": "assumed",
              r"& 1 & shared$^{b}$\\": "shared"}
    cells = {}
    for ln in lines:
        for end, st in STATUS.items():
            if ln.endswith(end):
                # row labels are LaTeX: undo the escapes of the solvent and Mg rows
                lab = ln.split(" & ")[0].replace("\\%", "%").replace("Mg$^{2+}$", "Mg2+")
                cells[lab] = st
    import reference_asymmetry as RA
    td = RA.tracer_design()
    sec = RA.table()
    with2 = set(map(tuple, sec[["family", "variant"]].drop_duplicates().values))
    from corpus import display_family
    # row by row, so a swap of two statuses cannot pass on the counts alone
    for _, x in td.iterrows():
        if x.design == "labeled":
            stem = f"{x.family} {x.variant}"
            want = "inferred" if (x.family, x.variant) in with2 else "assumed"
        else:
            fam = display_family(x.family, x.variant) if x.family == "TSase" else x.family
            var = "WT" if (x.family == "ecTSase" and x.variant == "wild type") else x.variant
            stem = f"{fam} {var} ({x.step}"
            want = "shared"
        # labeled rows are named exactly; shared rows carry "(step" and may add a
        # temperature, so they match on that prefix (which includes the step)
        hit = [k for k in cells if (k == stem if x.design == "labeled" else k.startswith(stem))]
        if len(hit) != 1 or cells[hit[0]] != want:
            print(f"  FAIL corpus consequences: protocol table row for {stem!r} is "
                  f"{[cells[h] for h in hit] or 'absent'}, the data say {want}")
            bad += 1
    inferred = sum(1 for v in cells.values() if v == "inferred")
    assumed = sum(1 for v in cells.values() if v == "assumed")
    shared = sum(1 for v in cells.values() if v == "shared")
    forms = inferred + assumed + shared
    if (inferred + assumed, shared) != (int((td.design == "labeled").sum()),
                                        int((td.design == "shared").sum())):
        print(f"  FAIL corpus consequences: protocol table {inferred}+{assumed} labeled / "
              f"{shared} shared, tracer_design.csv {(td.design == 'labeled').sum()} / "
              f"{(td.design == 'shared').sum()}")
        bad += 1
    if forms != u - 2:
        print(f"  FAIL corpus consequences: protocol table has {forms} forms, "
              f"units minus the two repeat 3 C records give {u - 2}")
        bad += 1
    need("si", f"outcome for the {_nw(forms)} distinct series and enzyme forms", "the form count")
    need("si", f"{_nw(inferred).capitalize()} of the {_nw(forms)} report secondary", "the inferred count")
    need("si", f"{_nw(shared)} need no inference", "the shared-tracer count")
    need("si", f"For the {_nw(shared)} the question does not arise", "the shared-tracer count")
    need("si", f"and {_nw(assumed)} cannot be checked at all", "the assumed count")
    need("si", f"The {_nw(inferred + assumed)} forms whose protocols place deuterium", "the labeled count")
    need("main", f"Across {_nw(inferred)} of the {_nw(forms)} enzyme forms", "the form count")
    need("main", f"For {_nw(shared)} forms the same tritiated molecule", "the shared-tracer count")
    need("main", f"In the remaining {_nw(assumed)}, dihydrofolate reductase mutants", "the assumed count")
    # the confidence procedure for temperature series (corrected 2026-09-14):
    # its simulated failure rate and Monte Carlo precision
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
    # Isotopic purity (Supplementary Note on the deuterated comparators, added
    # 2026-10-20; rebuilt 2026-10-21 after an external audit found the contaminant
    # rate wrong for dideuterated designs, the H/T error dropped from corrected
    # values, and per-unit robustness quoted as if it held across all units).
    # Every number is recomputed from isotopic_purity.py.
    import isotopic_purity as IP
    U_ = {lab: (src, g) for lab, src, g in IP.units()}
    ex = IP.exclusions().set_index("unit")
    yb, yc = "YADH wild type, 25 C (Cha1989)", "YADH wild type (p-chlorobenzyl alcohol) (hydride)"
    rb, rc = ex.loc[yb], ex.loc[yc]
    # the stated conjunctions must be what the numbers say
    if not (rb.undo0 > F0 and rb.undo5 > F0 and rb.undo0_all < F0 and rb.undo5_all < F0):
        print("  FAIL corpus consequences: the benzyl alcohol purity/multiplicity pattern changed")
        bad += 1
    if not (rc.undo0_all > F0 and rc.undo5_all > F0):
        print("  FAIL corpus consequences: p-chlorobenzyl alcohol no longer survives all three checks")
        bad += 1
    if not all(ex.loc[u].undo0 < F0 for u in ex.index if "ChinKlinman2000" in u):
        print("  FAIL corpus consequences: the double-mutant correction scenario changed sides")
        bad += 1
    need("si", f"falls from $+0.049$ to ${rb.undo0:+.3f}$ at small conversion and ${rb.undo5:+.3f}$ at $f=0.5$",
         "the benzyl alcohol bound without the protium correction")
    need("si", f"it falls to ${rb.undo0_all:+.3f}$ and ${rb.undo5_all:+.3f}$, below $\\Fz$",
         "the benzyl alcohol bound without its correction, across all units",
         (f"to ${rb.undo0_all:+.3f}$ and ${rb.undo5_all:+.3f}$, below $\\Fz$",))
    need("si", f"from $+0.170$ to ${rc.undo0:+.3f}$ and ${rc.undo5:+.3f}$, and across the $73$ units to "
               f"${rc.undo0_all:+.3f}$ and ${rc.undo5_all:+.3f}$",
         "the p-chlorobenzyl alcohol bounds without the correction",
         (f"${rc.undo0_all:+.3f}$ and ${rc.undo5_all:+.3f}$",))
    need("main", f"leaves the benzyl alcohol bound at ${rb.undo0:+.3f}$",
         "the benzyl alcohol bound without its correction (main)")
    need("main", f"well it falls to ${rb.undo0_all:+.3f}$", "the combined benzyl alcohol bound (main)")
    need("main", f"and removal of its correction together (${rc.undo0_all:+.3f}$)",
         "the combined p-chlorobenzyl alcohol bound (main)")
    pn0, pn5 = IP.protium_needed(yb, 0.0), IP.protium_needed(yb, 0.5)
    need("si", f"at least ${100 * pn0:.2f}\\%$ protium (small conversion) or ${100 * pn5:.2f}\\%$",
         "the protium level the combined benzyl alcohol claim needs",
         (f"${100 * pn0:.2f}\\%$ protium (small conversion)",))
    need("main", f"at least ${100 * pn0:.2f}\\%$ protium", "the protium level (main)")
    for u in ex.index:
        if "ChinKlinman2000" in u:
            need("si", f"${ex.loc[u].undo0:+.3f}$", f"the double-mutant scenario bound {u}")
    # the exclusion table the main text inputs must be the generated one
    import os as _os
    tex_tab = open("../results/sm_table_exclusions.tex").read()
    if IP.exclusion_table(IP.exclusions()) != tex_tab:
        print("  FAIL corpus consequences: results/sm_table_exclusions.tex is stale")
        bad += 1
    ta = IP.threshold(U_["ecTSase wild type, 20 C (Agrawal2004)"][1], "bound", 0.0)
    tw = IP.threshold(U_["ecTSase Y209W (hydride)"][1], "bound", 0.0)
    need("si", f"clear $\\Fz$ once protium of ${100 * ta:.2f}\\%$ and ${100 * tw:.2f}\\%$ is",
         "the protium thresholds of the two units that move")
    chi = IP.chi_contaminant(7.13, 22.529, 1.456134)
    need("si", f"$\\chi={chi:.2f}$ against $\\KHT=7.13$", "the dideuterated contaminant rate example")
    fac0 = IP.correction_factor(0.0065, 0.0065, 3.9, 1.2, 0.0) - 1
    fac5 = IP.correction_factor(0.0065, 0.0065, 3.9, 1.2, 0.5) - 1
    need("si", f"here ${100 * fac0:.1f}\\%$ at small conversion and ${100 * fac5:.1f}\\%$ at $f=0.5$",
         "Cha's correction factor", (f"${100 * fac0:.1f}\\%$ at small conversion and ${100 * fac5:.1f}\\%$",))
    pur = pd.read_csv("../results/isotopic_purity.csv")
    st = pur[(pur.f == 0.0) & (pur.scenario == "as stated")]
    n_ex = int((st.now == "excludes").sum())
    n_ex0 = int((st.base == "excludes").sum())
    if not ((st.base != st.now) <= ((st.now == "excludes") | (st.base == "uninformative"))).all() \
            or (st[(st.base == "excludes")].now != "excludes").any():
        print("  FAIL corpus consequences: a purity correction moved a unit toward the uninformative side")
        bad += 1
    need("si", f"the count of exclusions rises from {_nw(n_ex0)} to {_nw(n_ex)} and none falls",
         "the purity scenario count",
         (f"the count of exclusions rises from {_nw(n_ex0)} to {_nw(n_ex)} and none falls",))
    # retired claims: robustness checks passed separately are not a joint pass
    for bad_phrase in ("neither depends on a correction for isotopic impurity",
                       "No exclusion reported here therefore depends on a correction",
                       "neither rests on a correction for isotopic impurity",
                       "The published errors are therefore adequate for the offset"):
        for doc_name, doc_txt in (("main", main), ("si", si)):
            if " ".join(bad_phrase.split()) in doc_txt:
                print(f"  FAIL corpus consequences: retired claim in {doc_name}: {bad_phrase!r}")
                bad += 1
    # the yeast bypass tolerance quoted in the main text, at the jointly solved r
    from completion import yadh_completion_joint, bypass_tolerance
    t13 = 100 * float(bypass_tolerance(7.13, 1.73, yadh_completion_joint(1.3)["r"], F0))
    t73 = 100 * float(bypass_tolerance(7.13, 1.73, yadh_completion_joint(7.3)["r"], F0))
    need("main", f"less than ${t13:.0f}\\%$ to ${t73:.0f}\\%$ of the", "the yeast bypass tolerance")
    # the Discussion now gives both thresholds; this one is against F0
    need("main", f"${t13:.0f}\\%$ to ${t73:.0f}\\%$ to reach the gated envelope",
         "the yeast bypass tolerance against F0 (Discussion)",
         (f"Removing exclusion of the gated envelope requires ${t13:.0f}\\%$ to ${t73:.0f}\\%$",))
    # the profile-likelihood stratification must say what its own table shows
    # (until 2026-09-14 it said no adequately fitting series discriminates)
    from scipy import stats as _st
    pr = pd.read_csv("../results/offset_profiles.csv")
    # Edge profiles use the smallest objective sampled through F=40, with the
    # lower end recomputed against it (legacy column name chi2_inf).
    ed = pd.read_csv("../results/profile_edges.csv")
    key = ["family", "system", "variant", "step"]
    pr = pr.merge(ed[key + ["chi2_inf", "F_lo"]].rename(columns={"F_lo": "F_lo_inf"}),
                  on=key, how="left")
    hit = pr.chi2_inf.notna()
    pr.loc[hit, "chi2_min"] = pr.loc[hit, "chi2_inf"]
    pr.loc[hit, "F_lo"] = pr.loc[hit, "F_lo_inf"]
    pr.loc[hit, "F_hi"] = np.inf
    pr["p_fit"] = _st.chi2.sf(pr.chi2_min, (2 * pr.n_T - 5).clip(lower=1))
    excl = pr[(pr.F_lo > 0) | (pr.F_hi < 0)]
    adq = excl[excl.p_fit > 0.05]
    above = adq[adq.F_lo > F0]
    below = adq[adq.F_hi < F0]
    need("si", f"profiles exclude $F=0$ in {_nw(len(excl))} of {_nw(len(pr))} series",
         "the profile exclusion count")
    need("si", f"{_nw(len(excl) - len(adq))} of those {_nw(len(excl))} belong to series the model",
         "the misfit exclusion count")
    need("si", f"The other {_nw(len(adq))} exclusions are in adequately fitting series",
         "the adequate exclusion count")
    need("si", f"{_nw(len(above)).capitalize()} ecDHFR and hsDHFR series", "the count wholly above F0")
    need("si", f"{_nw(len(below)).capitalize()} give intervals wholly below $\\Fz$", "the count wholly below F0")
    # the admitted E. coli thymidylate synthase record
    e = ya[ya.grp == "ectsase"].iloc[0]
    fe = F_min_exact(e.K_HT, e.K_DT)[0]
    se = np.sqrt((e.K_HT_se / e.K_HT) ** 2 + (GSC * e.K_DT_se / e.K_DT) ** 2
                 + 2 * GSC * (e.K_HT_se / e.K_HT) * (e.K_DT_se / e.K_DT))
    lh = (e.K_HT - 1) / (e.K_DT - 1)
    # the single-record bound is the exact quantile of the endpoint
    # (bounds_uncertainty.exact_bound); on the open branch, as here, it equals
    # the closed form for F_obs, and both are checked
    eb = BU.exact_bound(e.K_HT, e.K_HT_se, e.K_DT, e.K_DT_se)
    if abs(eb - (fe - 1.645 * se)) > 5e-4:
        print("  FAIL corpus consequences: ecTSase exact bound departs from the open-branch form")
        bad += 1
    need("main", f"endpoint $+{fe:.3f}$", "the ecTSase endpoint", (f"endpoint of ${fe:+.3f}$",))
    need("main", f"bound ${eb:.3f}$", "the ecTSase bound", (f"bound of ${eb:.3f}$",))
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
              f"margin {b1 - F0:.3f}); {forms} protocol forms; all {n_adm} admissible records "
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

    def need(doc, text, what, alternatives=()):
        nonlocal bad
        if not any(candidate in (main if doc == "main" else si)
                   for candidate in (text, *alternatives)):
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
    # phi is a ratio to the TRITIUM chemical rate (2026-10-15 audit): the text
    # names that rate rather than "the isotope-sensitive rate"
    need("main", f"a bypass faster than ${z13:.0f}\\%$ to ${z73:.0f}\\%$ of the tritium chemical rate "
                 f"removes the yeast exclusion",
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
    # A derivative sign at infinity is not a global endpoint criterion.
    limiting, actual = NG.asymptotic_criterion_counterexample()
    need("si", f"holds and $E={limiting:+.4f}$, but the endpoint is ${actual:.4f}$",
         "the counterexample to the local bypass criterion")
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
    # the reversible corner path inside the Proposition S9 window
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

    def need(doc, text, what, alternatives=()):
        nonlocal bad
        if not any(candidate in (main if doc == "main" else si)
                   for candidate in (text, *alternatives)):
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

    def need(doc, text, what, alternatives=()):
        nonlocal bad
        if not any(candidate in (main if doc == "main" else si)
                   for candidate in (text, *alternatives)):
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
    import yadh_robustness as YR
    if abs(YR.t_tail_bound(10, u) - tb) > 1e-12:
        print("  FAIL readiness scope: yadh_robustness.t_tail_bound disagrees with the audit")
        bad += 1
    tp = YR.pchloro_series(u)["t10"]
    need("main", f"bounds at $\\rho=-1$ then fall to ${tb:.3f}$ and ${tp:.3f}$,", "the t-tail sensitivity")
    # the main text cites Note 18 for this number; until 2026-10-16 the SI gave none
    need("si", f"at $\\rho=-1$ falls to ${tb:.3f}$, below $\\Fz$", "the t-tail sensitivity (SI)")
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
    # "classifies", never "determines": topology fixes which geometries are
    # available and the rates select among them; a rewrite put it back in a heading
    for doc in ("main", "si"):
        for w in ("pathways determine which", "topology determines", "topology determine ",
                  "network determines which", "structure determines which"):
            forbid(doc, w, "topology determining the geometry")
    forbid("main", "settled by the mechanism assumed, not by the precision achieved",
           "discrimination settled by topology alone")
    if not bad:
        print(f"  readiness scope: closed-branch binding shift {e0 - e1:.3f}, plug-in coverage "
              f"{100 * y:.1f}%/{100 * h:.1f}%, t10 Bonferroni bound {tb:.3f}")
    return bad


def check_external_audit_2026_10_12():
    r"""Numbers corrected or added after the pre-submission audit of 2026-10-12.

    That audit found the promoting-mode spans and the Fisher example computed
    with bare integer masses in qtunnel.py against the reduced-mass gamma_SC and
    F0; the horse liver pooling fit asserted, not deposited, and given one
    degree of freedom too few (an exact scale symmetry removes a parameter);
    the BSAO bound printed at rho = 0 without saying so; a stale vibronic
    example; profile optima pinned at the grid edge reported as optima; and an
    argument for the yeast r root that did not prove uniqueness.  Every
    replacement number is recomputed here and matched to its sentence.
    """
    import completion as C
    import pooling_check as P
    import vibronic_envelope as VE
    import yadh_robustness as YR
    from partial_id import F_min_exact as _fme
    main = " ".join(DOCS["main"].read_text().split())
    si = " ".join(DOCS["si"].read_text().split())
    bad = 0

    def need(doc, text, what, alternatives=()):
        nonlocal bad
        if not any(candidate in (main if doc == "main" else si)
                   for candidate in (text, *alternatives)):
            print(f"  FAIL audit 2026-10-12: {doc} lacks {what}: {text!r}")
            bad += 1

    # one mass convention: the gated model must reproduce masses.py's limits
    import qtunnel as Q
    if abs(Q.gamma_rigid() - M.gamma_rigid("C")) > 1e-12 or \
       abs(Q.gamma_fully_gated() - M.gamma_gated("C")) > 1e-12:
        print("  FAIL audit 2026-10-12: qtunnel.py and masses.py use different masses")
        bad += 1
    sp = C.spans_250_350()
    g = [v[1] for v in sp.values()]
    need("si", f"moves by ${min(g):.3f}$ to ${max(g):.3f}$, against a distance of",
         "the exponent span")
    need("main", f"moves by at most ${max(g):.3f}$ over the same interval", "the exponent span")
    sv, cond = C.fisher_example()
    e1, e2 = (int(np.floor(np.log10(v))) for v in sv[1:])
    ec = int(np.floor(np.log10(cond)))
    need("si", f"${sv[0]:.2f}$, ${sv[1] / 10 ** e1:.1f}\\times10^{{{e1}}}$ and "
               f"${sv[2] / 10 ** e2:.1f}\\times10^{{{e2}}}$", "the Fisher singular values")
    need("si", f"\\simeq{cond / 10 ** ec:.1f}\\times10^{{{ec}}}$", "the Fisher condition number")
    lo, hi, (q1, q3) = C.variance_sensitivity()
    need("si", f"between ${lo:.2f}$ and ${hi:.1f}|\\Fz|$", "the 10% variance sweep")
    need("si", f"interquartile range of ${q1:.2f}$ to ${q3:.2f}|\\Fz|$", "the 10% variance sweep")
    # the yeast r root: negative at r = 1 in every draw, positive leading term
    n, at1, lead = C.joint_r_value_at_one()
    if not (at1 < 0 < lead):
        print("  FAIL audit 2026-10-12: the r quadratic is not negative at r = 1 in every draw")
        bad += 1
    need("si", f"(its largest value there is ${at1:.2f}$)", "the value at r = 1")
    # variant pooling
    f = P.ladh_shared_pair_fit()
    if f["rank"] != f["n_par"] - 1 or f["witness_diff"] > 1e-12:
        print("  FAIL audit 2026-10-12: the pooling symmetry did not reproduce")
        bad += 1
    pe = int(np.floor(np.log10(f["p"])))
    for doc in ("main", "si"):
        need(doc, f"$\\chi^{{2}}={f['chi2']:.1f}$ on ${f['dof']}$ degrees of freedom", "the pooled fit")
    need("si", f"$p={f['p'] / 10 ** pe:.0f}\\times10^{{{pe}}}$", "the pooled fit p value")
    F1, F2 = f["witness_F"]
    need("main", f"offsets ${F1:+.3f}$ and ${F2:+.3f}$", "the pooling witness")
    need("si", f"offsets are ${F1:+.3f}$ and ${F2:+.3f}$", "the pooling witness")
    need("si", f"commitments from ${f['witness_cmin']:.3f}$ to $5.0$", "the pooling witness")
    # within-series pooling of the observed offsets (README once held 7.8 and 88)
    import pandas as _pd
    _bm = _pd.read_csv("../data/trinomial_benchmark.csv")
    r0, r1 = P.series_chi2(_bm, 0.0), P.series_chi2(_bm, -1.0)
    need("si", f"per degree of freedom of ${r0.chi2_dof.median():.1f}$, with values up to "
               f"${r0.chi2_dof.max():.0f}$", "the within-series pooling statistic")
    need("si", f"the median is ${r1.chi2_dof.median():.1f}$ and {int((r1.chi2_dof > 1).sum())} of the {len(r1)}",
         "the within-series pooling statistic at rho = -1")
    # the correlation-free yeast bounds
    import corpus
    u = corpus.counts()["units"]
    b1, bu = YR.correlation_free_bound(), YR.correlation_free_bound(units=u)
    need("main", f"two effects gives ${b1:+.3f}$", "the correlation-free bound")
    need("main", f"no assumption about the correlation at all (${bu:.3f}$)",
         "the correlation-free multiplicity bound")
    need("si", f"For the yeast average this gives ${b1:+.3f}$", "the correlation-free bound")
    need("si", f"it gives ${bu:.3f}$, still above $\\Fz$ by ${bu - F0:.3f}$",
         "the correlation-free multiplicity bound")
    # BSAO at both correlations
    G_ = GSC
    fb = _fme(35.2, 3.07)[0]
    a, b = 0.8 / 35.2, G_ * 0.07 / 3.07
    need("si", f"bound of ${fb - 1.645 * (a + b):.3f}$ at the maximally adverse $\\rho=-1$ used "
               f"throughout (${fb - 1.645 * np.hypot(a, b):.3f}$ under independence)", "the BSAO bound")
    # the positive summed-vibronic example, from the printed parameters
    xh, xd = VE.summed_X(2.011, 0.346, 3626.0, 11.55, -11.42, 284.8, nmax=30)
    need("si", f"$\\KHT^{{\\mathrm{{int}}}}={np.exp(xh):.4f}$, $\\KDT^{{\\mathrm{{int}}}}={np.exp(xd):.4f}$",
         "the positive vibronic example")
    need("si", f"F={xh - GSC * xd:+.4f} .", "the positive vibronic example")
    for cell, (gm, ex) in VE.bvib_extended().items():
        e = int(np.floor(np.log10(ex[1])))
        need("si", f"$({cell[0]},{cell[1]},{cell[2]})$ & ${gm[0]:+.3f}$ & ${ex[0]:+.3f}$ & "
                   f"${ex[1] / 10 ** e:.1f}\\times10^{{{e}}}$", f"the extended B_vib row {cell}")
    # profiles whose grid optimum sat on the edge
    ed = pd.read_csv("../results/profile_edges.csv")
    if ed.finite_optimum.any():
        print("  FAIL audit 2026-10-12: an edge profile has a resolved optimum")
        bad += 1
    d26 = ed[ed.variant == "d26-PPLR-d62-PEKN"].iloc[0]
    need("si", f"falls toward ${d26.chi2_inf:.2f}$ over that range, and its explored confidence region "
               f"extends from ${d26.F_lo:+.2f}$ through $+40$ without resolving an upper limit",
         "the d26 edge profile")
    tp = ed[(ed.family == "TSase") & (ed.step == "proton")].iloc[0]
    need("si", f"($\\chi^{{2}}={tp.chi2_inf:.0f}$ on three degrees of freedom", "the hsTSase proton misfit")
    # the four mechanisms of the atlas figure
    import export_fig_atlas as EA
    for tag, name, sub in EA.CLASSES:
        if tag not in ("mixed", "inverse"):
            continue
        co, sub, (A_, B_, C_, D_) = EA.resolve(name, sub)
        from network_atlas import curvature_point
        xs = curvature_point(co, sub)
        fmt = (lambda v: f"{v:.2f}") if tag == "mixed" else (lambda v: f"{v:.1f}")
        need("si", f"${fmt(A_)}$ & ${fmt(B_)}$ & ${fmt(C_)}$ & ${fmt(D_)}$ & ${xs:.2f}$",
             f"the {tag} row of the atlas-rates table")
        if tag == "mixed":
            need("main", f"$x^{{*}}={xs:.2f}$", "the mixed-case x* in the Fig. 2 caption")
    if not bad:
        print(f"  audit 2026-10-12: spans, Fisher cond {cond:.2e}, pooled chi2 {f['chi2']:.1f}/{f['dof']}, "
              f"correlation-free {b1:+.3f}/{bu:+.3f}, vibronic {xh - GSC * xd:+.4f}, edge profiles, atlas rates")
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


FIGURE_EXPORTERS = ("export_fig_overview", "export_fig_envelope", "export_fig_design")


def check_figure_data():
    """The data behind the three figures added on 2026-10-14 must be current.

    Each exporter is rerun into a scratch directory, which also runs its own
    checks against the values the text quotes, and every file it writes must
    equal the copy the figure sources read.  Nothing else re-derives those
    files, so without this a changed computation would leave a figure drawing
    superseded numbers while every text guard still passed.  Runs without the
    manuscript, so the archive checks it too.
    """
    import contextlib
    import importlib
    import io
    import tempfile
    from pathlib import Path as _P
    data = _P(__file__).resolve().parent.parent / "figures" / "tikz" / "data"
    bad = 0
    nfiles = 0
    for name in FIGURE_EXPORTERS:
        mod = importlib.import_module(name)
        with tempfile.TemporaryDirectory() as d:
            saved, mod.OUT = mod.OUT, d + "/"
            buf = io.StringIO()
            try:
                with contextlib.redirect_stdout(buf):
                    rc = mod.main()
            finally:
                mod.OUT = saved
            if rc:
                print(f"  FAIL figure data: {name} reports {rc} failed checks")
                print("    " + buf.getvalue().strip().replace("\n", "\n    "))
                bad += 1
            for f in sorted(_P(d).iterdir()):
                nfiles += 1
                ship = data / f.name
                if not ship.exists() or ship.read_text() != f.read_text():
                    print(f"  FAIL figure data: figures/tikz/data/{f.name} is not "
                          f"what {name}.py writes")
                    bad += 1
    if not bad:
        print(f"  figure data: {len(FIGURE_EXPORTERS)} exporters pass their checks "
              f"and reproduce {nfiles} shipped files")
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
def check_source_review_2026_10_13():
    r"""Corrections from reading every cited source in full (2026-10-13).

    The review admitted the Tsai & Klinman methanol records (their stated
    exclusion contradicted the source, which ran the competitive experiments at
    saturating NAD+, and the paper's own theory, under which a random mechanism
    with no bypass keeps the base map), and found that the closed-form endpoint
    scored an inadmissible pair (K_DT < 1) as +0.822.  This checks the numbers
    that the admission produced, the admissibility rule, and that corrected
    phrasings do not come back.
    """
    import corpus
    from partial_id import F_min_exact as _fme
    from scipy import stats
    main = " ".join(DOCS["main"].read_text().split())
    si = " ".join(DOCS["si"].read_text().split())
    bad = 0

    def need(doc, text, what):
        nonlocal bad
        if text not in {"main": main, "si": si}[doc]:
            print(f"  FAIL source review: {doc} lacks {what}: {text!r}")
            bad += 1

    # the admissibility rule of the scalar endpoint
    if _fme(2.2, 0.99)[0] != -np.inf:
        print("  FAIL source review: F_min_exact scores an inadmissible pair")
        bad += 1
    # the two methanol series that exclude: endpoints and nominal bounds
    b = pd.read_csv("../results/bounds_uncertainty.csv")
    for v, lo, hi in (("F93W (40% MeOH)", "$+0.183$", "40"), ("F93W (60% MeOH)", "$+0.395$", "60")):
        s = b[(b.family == "LADH") & (b.variant == v)]
        pt = s[s.rho == 0.0].point.iloc[0]
        l1, l0 = s[s.rho == -1.0].lcb.iloc[0], s[s.rho == 0.0].lcb.iloc[0]
        need("main", f"${pt:+.3f}$", f"the {hi}% methanol endpoint")
        need("main", f"${l1:+.3f}$", f"the {hi}% methanol bound at rho = -1")
        need("main", f"${l0:+.3f}$", f"the {hi}% methanol bound under independence")
        if not (l1 > F0 and l1 > 0):
            print(f"  FAIL source review: {v} no longer clears both thresholds")
            bad += 1
    # the units whose rho = -1 bound clears F0, counted from the data (three
    # until the census of 2026-10-20: yeast and the two methanol series)
    _b1 = pd.read_csv("../results/bounds_uncertainty.csv")
    _n = int((_b1[_b1.rho == -1.0].lcb > F0).sum())
    for _, r_ in corpus.single_condition().iterrows():
        _p = _fme(r_.K_HT, r_.K_DT)[0]
        _a, _c = r_.K_HT_se / r_.K_HT, GSC * r_.K_DT_se / r_.K_DT
        _n += int(_p - 1.645 * (_a + _c) > F0)
    need("main", f"{_nw(_n).capitalize()} units have nominal per-unit confidence bounds above $\\Fz$",
         "the exclusion count")
    # multiplicity for the methanol series, across all units.  At alpha/(u n_T)
    # a Monte Carlo quantile rests on a few hundred tail draws, so these are
    # computed exactly, as the yeast multiplicity bounds are.  The endpoint never
    # exceeds F_obs, so each record's closed-form F_obs bound is an upper limit on
    # its own bound; when the record attaining the largest of them stays on the
    # open branch (L_H >= gamma) at its perturbed pair, that value IS the series
    # bound.  Records with no admissible central pair contribute -inf.
    u = corpus.counts()["units"]
    s = corpus.series()
    for v in ("F93W (40% MeOH)", "F93W (60% MeOH)"):
        g = s[(s.family == "LADH") & (s.variant == v)]
        z = stats.norm.ppf(1 - 0.05 / (u * len(g)))
        for rho in (0.0, -1.0):
            best = None
            for _, a in g.iterrows():
                if not (a.K_HT > a.K_DT > 1):
                    continue
                x, y = a.K_HT_se / a.K_HT, a.K_DT_se / a.K_DT
                sd = np.hypot(x, GSC * y) if rho == 0.0 else x + GSC * y
                m = np.log(a.K_HT) - GSC * np.log(a.K_DT) - z * sd
                if best is None or m > best[0]:
                    best = (m, a.K_HT * np.exp(-z * x), a.K_DT * np.exp(z * y))
            m, hl, dh = best
            if not (hl > dh > 1 and (hl - 1) / (dh - 1) >= GSC):
                print(f"  FAIL source review: {v} multiplicity bound leaves the open branch")
                bad += 1
            need("main", f"${m:+.3f}$", f"the {v} multiplicity bound at rho = {rho:+.0f}")
    # credit and convention fixes that must not regress
    for doc, text, what in (
            ("main", "Eq.~(7) of Siebrand and Smedarchina", "credit for the averaged overlap"),
            ("si", "Proposition~S1 is Eq.~(7) of Siebrand and", "credit for Proposition S1"),
            ("main", "\\cite{Sikorski2004}", "credit for the light-enzyme data"),
            ("main", "\\cite{KingAltman1956}", "the King-Altman citation"),
            ("main", "reciprocal of the commitment evaluated with the tritium rate",
             "the commitment convention")):
        need(doc, text, what)
    for doc, text in (("main", "inverse of the commitment to catalysis in the usual convention"),
                      ("main", "Northrop's inversion~\\cite{Northrop1975,Northrop2001"),
                      ("main", "SoudackovHammesSchiffer2000,HammesSchifferSoudackov2008}. The transferred"),
                      ("si", "``reduced relative to''"),
                      ("si", "drift between columns"),
                      ("main", "companion supplemental"), ("main", "this revision"),
                      ("si", "excluded on a ground internal to this paper")):
        if text in {"main": main, "si": si}[doc]:
            print(f"  FAIL source review: {doc} still has the retired phrasing {text!r}")
            bad += 1
    bib = (DOCS["main"].parent / "references.bib").read_text()
    for wrong in ("Kim, Kwangsoo", "Hothi, Poonam"):
        if wrong in bib:
            print(f"  FAIL source review: references.bib still has {wrong!r}")
            bad += 1
    if not bad:
        print("  source review 2026-10-13: methanol series and multiplicity, admissibility "
              "rule, credits and retired phrasings")
    return bad


def check_light_example():
    """The ecDHFR light-enzyme 5 C example, as the geometry caption states it.

    Until 2026-10-15 the caption gave c = 4.3 and x_H = 6.1 against "6.2 +- 0.4",
    computed from the rounded reprint (3.10) and checked by nothing.  Sikorski et
    al. 2004 report 3.06 and, at 5 C, an intrinsic H/T effect of 6.18 +- 0.39
    (SI Table 2, 95% interval); their D/T reference carries deuterium at 4S, so
    the crossing is also stated at the observed lower bound R on r.
    """
    import pandas as pd
    from scipy.optimize import brentq
    import corpus
    import reference_asymmetry as RA
    kh, kd = corpus.series_record("light enzyme", 5.0)
    R = RA.r_at("ecDHFR", "light enzyme", 25.0)
    t = RA.table()
    rse = float(t[(t.family == "ecDHFR") & (t.variant == "light enzyme")].r_se.iloc[0])

    def cross(r):
        F = lambda c: (np.log(kh * c / (1 + c - kh))
                       - GSC * np.log(kd * r * c / (1 + r * c - kd)))
        return brentq(F, (kh - 1) * 1.000001, 1e6)
    c1, cR = cross(1.0), cross(R)
    xh = kh * c1 / (1 + c1 - kh)
    main = " ".join(DOCS["main"].read_text().split())
    si = " ".join(DOCS["si"].read_text().split())
    bad = 0
    for doc, text, need_ in (
            ("main", main, f"crosses the ray at $c={c1:.1f}$ (open circle)"),
            ("main", main, f"recovers an intrinsic H/T effect of ${xh:.1f}$, against the "
                           f"$6.18\\pm0.39$ ($95\\%$ interval) the source reports"),
            ("main", main, f"so $r\\ge{R:.3f}$ and the crossing moves to $c\\ge{cR:.1f}$"),
            ("si", si, f"$R={R:.3f}\\pm{rse:.3f}$ at $25\\,^{{\\circ}}$C")):
        if need_ not in text:
            print(f"  FAIL light example: {doc} does not state {need_!r}")
            bad += 1
    # the series' profile fit, quoted in the SI prose (62.6/7 until 2026-10-15)
    from scipy.stats import chi2 as _chi2
    pr = pd.read_csv("../results/offset_profiles.csv")
    row = pr[(pr.family == "ecDHFR") & (pr.variant == "light enzyme")].iloc[0]
    dof = int(row.n_T) + 1
    pv = float(_chi2.sf(row.chi2_min, dof))
    want = (f"the light series does not ($\\chi^{{2}}={row.chi2_min:.1f}/{dof}$, "
            f"$p={pv:.1e}".replace("e-", "\\times10^{-").replace("p=", "p=") + "}$)")
    if want not in si:
        print(f"  FAIL light example: si does not state {want!r}")
        bad += 1
    if not bad:
        print(f"  ecDHFR 5 C: crossing c = {c1:.3f} (x_H {xh:.3f}), at R = {R:.4f} "
              f"c = {cR:.3f}")
    return bad


def check_audit_2026_10_15():
    """Corrections from the external audit of 2026-10-15, each recomputed.

    Every number here comes from a deposited function, and every retired
    phrasing is forbidden, so the corrections cannot silently regress.
    """
    import pandas as pd
    import qtunnel as Q
    import joint_nuisance as J
    import network_geometry as NG
    import bounds_uncertainty as BU
    import corpus
    main = " ".join(DOCS["main"].read_text().split())
    si = " ".join(DOCS["si"].read_text().split())
    bad = 0

    def need(doc, text, what):
        nonlocal bad
        if text not in (main if doc == "main" else si):
            print(f"  FAIL audit 2026-10-15: {doc} lacks {what}: {text!r}")
            bad += 1

    # C. deletion sensitivity of the confidence decision (holdout.py writes it)
    path = "../results/deletion_sensitivity.csv"
    try:
        cd = pd.read_csv(path)
    except FileNotFoundError:
        print(f"  FAIL audit 2026-10-15: {path} missing; run holdout.py")
        return bad + 1
    for v, pct in (("F93W (40% MeOH)", "40"), ("F93W (60% MeOH)", "60")):
        g = cd[cd.variant == v]
        r3 = g[g.omitted_T == 3.0].iloc[0]
        rest = g[g.omitted_T != 3.0]
        if r3.excludes or not rest.excludes.all() or not r3.point > F0:
            print(f"  FAIL audit 2026-10-15: {v} deletion pattern changed")
            bad += 1
        b3 = f"{r3.bound:+.3f}" if abs(r3.bound) < 1 else f"{r3.bound:+.2f}"
        need("si", f"& $3\\,^{{\\circ}}$C & ${r3.point:+.3f}$ & ${b3}$ & no", f"the {pct}% deletion row")
        need("si", f"& ${rest.point.iloc[0]:+.3f}$ & ${rest.bound.max():+.3f}$ & yes",
             f"the {pct}% retained rows")
    lo40 = cd[(cd.variant == "F93W (40% MeOH)") & (cd.omitted_T == 3.0)].bound.iloc[0]
    lo60 = cd[(cd.variant == "F93W (60% MeOH)") & (cd.omitted_T == 3.0)].bound.iloc[0]
    need("main", f"the bounds fall to ${lo40:+.3f}$ and ${lo60:+.2f}$", "the deletion bounds")
    # D. nominal per-unit wording and the multiplicity survivor
    need("main", "Six units give nominal per-unit exclusions", "the abstract qualifier")
    need("main", "Only the two yeast exclusions survive adjustment across all units",
         "the abstract multiplicity")
    u = corpus.counts()["units"]
    sc = corpus.single_condition()
    # which units survive Bonferroni across all units at rho = -1, computed: it
    # must be exactly the two yeast units the abstract names
    surv = [f"single {r_.family} {r_.variant}" for _, r_ in sc.iterrows()
            if BU.multiplicity_bound(pd.DataFrame([r_]), u, -1.0) > F0]
    for (fam, var, st), g in corpus.series().groupby(["family", "variant", "step"]):
        if BU.multiplicity_bound(g, u, -1.0) > F0:
            surv.append(f"series {fam} {var}")
    if sorted(surv) != ["series YADH wild type (p-chlorobenzyl alcohol)", "single YADH wild type"]:
        print(f"  FAIL audit 2026-10-15: multiplicity survivors at rho = -1 are {surv}")
        bad += 1
    # the figure is two panels since 2026-10-20
    rays = "\n".join(open(f"../figures/tikz/data/fa_rays_{k}.tex").read() for k in "ab")
    stars = [l for l in rays.splitlines() if "^{\\ast}" in l]
    dags = [l for l in rays.splitlines() if "^{\\dagger}" in l]
    if (len(stars) != 2 or not all("YADH" in l for l in stars)
            or len(dags) != 1 or "50" not in dags[0]):
        print("  FAIL audit 2026-10-15: Fig. 5 marks are not on the two yeast units (*) "
              "and WT 50% MeOH (dagger)")
        bad += 1
    # F. a finite member of the model, not a pair on the unattained envelope
    bad += J.check_finite_witness(verbose=False)
    wt = J.finite_witness()
    need("main", f"$\\mathcal{{A}}=3.2\\times10^{{5}}$ and $w={wt['w']:.0f}$", "the witness parameters")
    need("main", f"$F^{{\\mathrm{{int}}}}={wt['F_int']:.4f}<\\Fz$", "the witness offset")
    need("si", f"$x_{{\\mathrm H}}={wt['xH']:.4f}$ and $x_{{\\mathrm D}}={wt['xD']:.4f}$", "the witness pair")
    need("si", f"$\\Ag={wt['A']/1e5:.2f}\\times10^{{5}}$ and $w={wt['w']:.0f}$", "the witness parameters (SI)")
    # I. the large-commitment candidate equals the exact endpoint where used
    worst = max(abs(NG.endpoint_exact(7.13, 1.73, ph, r) - NG.endpoint_closed(7.13, 1.73, ph, r))
                for r in np.linspace(1.31, 1.60, 8) for ph in np.linspace(0, 1, 21))
    if worst > 1e-9:
        print(f"  FAIL audit 2026-10-15: E_r^infinity differs from the exact endpoint by {worst:.2e}")
        bad += 1
    # Note 1: the supremum zero needs w -> 0 as well
    lh, ld = Q.log_kie_intrinsic(1e-12, 0.5)
    need("si", f"${lh - GSC * ld:.3f}$ at $w=0.5$", "the fixed-w prefactor offset")
    # replicate uniqueness (needs the fetched spreadsheet)
    try:
        import build_trinomial as BT
        blocks = BT.replicate_uniqueness()
        if blocks and not all(b[1] == 1 for b in blocks):
            print("  FAIL audit 2026-10-15: a spreadsheet block has more than one matching replicate set")
            bad += 1
    except Exception as e:                        # sources not fetched
        print(f"  replicate uniqueness skipped ({type(e).__name__})")
    # retired phrasings
    for doc, text in (("main", "a bounded $k_b$"), ("main", "of the isotope-sensitive rate"),
                      ("main", "a test at level at most"), ("main", "standard errors at $5"),
                      ("main", "What the data support directly"),
                      ("main", "no general steady-state theory"),
                      ("main", "intrinsic pair on the gated envelope reproduces both"),
                      ("main", "envelope then reproduces both observed effects"),
                      ("si", "under every leave-one-out"), ("si", "with no exceptions; equivalently"),
                      ("si", "makes the extraction self-validating"),
                      ("si", "for which the exponent and contrast are not defined"),
                      ("si", "routinely measured to the required precision"),
                      ("si", "that requires profile likelihoods"), ("si", "comparable to or above"),
                      ("si", "largest offset measured"), ("si", "is the contrapositive"),
                      ("si", "the one system in the record that excludes"),
                      ("si", "being excluded in Supplementary Note"),
                      ("si", "Why solvent viscosity does not work"),
                      ("si", "Out-of-sample checks")):
        if text in (main if doc == "main" else si):
            print(f"  FAIL audit 2026-10-15: retired phrasing in {doc}: {text!r}")
            bad += 1
    if not bad:
        print(f"  audit 2026-10-15: deletion bounds {lo40:+.3f}/{lo60:+.2f}, finite witness "
              f"F = {wt['F_int']:+.4f}, E_r exact, marks on yeast and WT 50% MeOH")
    return bad


def check_tables_2026_10_15():
    """Table-by-table audit of 2026-10-15: every corrected table recomputed.

    Four supplementary tables had no deposited generator (the horse liver
    correlation bracket, the Stojkovic intrinsic pairs, the vibronic penalty) or
    printed values at the wrong precision (the benchmark: "13 +- 0.1" for the
    published 13.0, "0.13" for 0.132).  Their bodies are now written by code and
    recomputed here; the text claims drawn from them are tied to the same code.
    The generated-file checks need no manuscript, so they also run on the
    deposit alone.
    """
    import adh_primary as AP
    import published_intrinsic as PI
    import vibronic as VB
    bad = 0

    def same(path, text, what):
        nonlocal bad
        try:
            shipped = open(path).read()
        except FileNotFoundError:
            shipped = None
        if shipped != text:
            print(f"  FAIL tables 2026-10-15: {path} differs from {what}; rerun it")
            bad += 1

    rb = AP.rho_bracket()
    same("../results/sm_table_ladhprim.tex", AP.ladhprim_table(rb), "adh_primary.py")
    same("../results/sm_table_stoj.tex", PI.table(), "published_intrinsic.py")
    same("../results/sm_table_vibronic.tex", VB.penalty_table(), "vibronic.py")
    lad = pd.read_csv("../data/ladh_adh_primary.csv")
    f93 = lad[(lad.variant == "F93W") & (lad.T_C == 25)].iloc[0]
    cr = AP.clearing_rho(f93.K_HT, f93.K_HT_se, f93.K_DT, f93.K_DT_se)
    # the Monte Carlo bound must agree with the closed-form crossing
    import bounds_uncertainty as BU
    g = lad[(lad.variant == "F93W") & (lad.T_C == 25)]
    if not (BU.series_bound(g, cr - 0.03)[0] < F0 < BU.series_bound(g, cr + 0.03)[0]):
        print(f"  FAIL tables 2026-10-15: Monte Carlo bound does not cross F0 near rho = {cr:.3f}")
        bad += 1
    st = PI.summary()
    if not st["within_rounding"] or round(st["ratio_max"]) != 3 or round(st["ratio_sd"]) != 7:
        print(f"  FAIL tables 2026-10-15: Stojkovic summary changed {st}")
        bad += 1
    nu = VB.crossing_wavenumber()
    data = open("../results/sm_table_data.tex").read()
    for cell in ("$13.0\\pm0.1$", "$2.896\\pm0.132$", "$2.00\\pm0.02$", "$3.06\\pm0.10$",
                 "$1.80\\pm0.003$", "& $-50$ &"):
        if cell not in data:
            print(f"  FAIL tables 2026-10-15: benchmark table lacks the as-published {cell!r}")
            bad += 1
    prof = open("../results/sm_table_profiles.tex").read().splitlines()
    for line in prof:
        if "+1.50]" in line and "$\\ge$" not in line:
            print(f"  FAIL tables 2026-10-15: grid-edge width printed as exact: {line[:50]}")
            bad += 1
    if not HAVE_DOCS:
        return bad

    main = " ".join(DOCS["main"].read_text().split())
    si = " ".join(DOCS["si"].read_text().split())

    def need(doc, text, what):
        nonlocal bad
        if text not in (main if doc == "main" else si):
            print(f"  FAIL tables 2026-10-15: {doc} lacks {what}: {text!r}")
            bad += 1

    for f in ("sm_table_ladhprim.tex", "sm_table_stoj.tex", "sm_table_vibronic.tex"):
        need("si", f"\\input{{{f}}}", f"the generated {f}")
    need("main", f"it clears $\\Fz$ only for $\\rho>{cr:.2f}$", "the F93W clearing correlation")
    need("si", f"for $\\rho>{cr:.2f}$ (\\texttt{{adh\\_primary.clearing\\_rho}})",
         "the F93W clearing correlation (SI)")
    need("si", f"$\\langle F\\rangle={st['mean']:.3f}$ and $\\max|F|={st['max_abs']:.3f}$",
         "the Stojkovic summary")
    need("si", f"Below ${nu:.0f}$~cm$^{{-1}}$ the penalty exceeds", "the vibronic crossing")
    need("si", "Supplementary Table~\\ref{tab:unc-si} gives", "a pointer to the correlation bracket")
    need("si", "ecTSase & \\emph{E.\\ coli} thymidylate synthase", "the species-named family")
    # the main text points at the atlas rate table by number
    aux = MR.one(DOCS["si"].parent, "*_si.tex").with_suffix(".aux")
    if aux.exists():
        m = re.search(r"newlabel\{tab:atlasrates\}\{\{(\d+)\}", aux.read_text())
        if not m or f"Supplementary Table~{m.group(1)}." not in main:
            print("  FAIL tables 2026-10-15: main-text pointer to the atlas rate table is stale")
            bad += 1
    for doc, text in (("main", "$\\rho\\ge0.5$ and not below"), ("si", "for $\\rho\\ge0.5$"),
                      ("si", "Below about $2600$"), ("si", "Tables carried from the main text"),
                      ("si", "the paper reports the maximally adverse"),
                      ("si", "are shown to two significant figures, except"),
                      ("si", "where $X_L\\propto\\sqrt{\\mu_L}$")):
        if text in (main if doc == "main" else si):
            print(f"  FAIL tables 2026-10-15: retired phrasing in {doc}: {text!r}")
            bad += 1
    if re.search(r"(?<![a-z])TSase & \\emph\{", si):
        print("  FAIL tables 2026-10-15: series table names a TSase family without its species")
        bad += 1
    if not bad:
        print(f"  tables 2026-10-15: four generated tables reproduce; F93W clears F0 for "
              f"rho > {cr:.3f}; vibronic crossing {nu:.0f} cm^-1; benchmark as published")
    return bad


def check_consistency_2026_10_16():
    """Main text against supplement, 2026-10-16.

    Propositions were hand-numbered and appeared S1-S6, S8, S9, S7; the main
    text cited a Note 18 number the note did not contain; the methanol
    multiplicity bounds were in the main text only; the joint binding-bypass
    contour was typed from an older coarse profile (its r = 1.31 and r = 1.57
    rows disagreed at phi = 0, where r cannot matter); the topology table still
    listed the bypass-only row without its fast-dissociation qualifier.
    """
    import bounds_uncertainty as BU
    import corpus
    raw = DOCS["si"].read_text()
    main = " ".join(DOCS["main"].read_text().split())
    si = " ".join(raw.split())
    bad = 0

    def fail(msg):
        nonlocal bad
        print(f"  FAIL consistency 2026-10-16: {msg}")
        bad += 1

    nums = [int(n) for n in re.findall(r"\\textbf\{Proposition S(\d+)\.\}", raw)]
    if nums != list(range(1, len(nums) + 1)):
        fail(f"propositions appear out of order: {nums}")
    m = re.search(r"\\textbf\{Proposition S(\d+)\.\} \\emph\{For the reversible scheme", raw)
    if not m or f"(Supplementary Proposition~S{m.group(1)})" not in main:
        fail("main text does not cite the reversible-scheme proposition by its number")
    u = corpus.counts()["units"]
    tb = pd.read_csv("../data/trinomial_benchmark.csv")
    vals = []
    for v in ("F93W (40% MeOH)", "F93W (60% MeOH)"):
        g = tb[tb.variant == v]
        vals.append([BU.multiplicity_bound(g, u, r) for r in (0.0, -1.0)])
    want = (f"gives ${vals[0][0]:+.3f}$ and ${vals[1][0]:+.3f}$ under independence, above "
            f"$\\Fz$, but ${vals[0][1]:+.3f}$ and ${vals[1][1]:+.3f}$ at $\\rho=-1$").replace("$+", "$")
    if want not in si:
        fail(f"SI lacks the methanol multiplicity bounds: {want!r}")
    D = "../figures/tikz/data/"
    for lab, f in (("1.31", "fd_joint_rlo.dat"), ("1.57", "fd_joint_rhi.dat")):
        d = np.loadtxt(D + f, skiprows=1)
        mp = {round(p, 2): v for p, v in d}
        cells = " & ".join(f"${mp[p]:.4f}$" if mp[p] > 0 else "none"
                           for p in np.round(np.arange(0, 0.161, 0.02), 2))
        if f"$F_{{\\mathrm{{bind}}}}^{{\\max}}$, $r={lab}$ & {cells}" not in si:
            fail(f"joint contour row r = {lab} differs from {f}")
    if "lower half-line (fast-dissociation limit only)" not in si:
        fail("topology table lacks the fast-dissociation qualifier")
    if re.search(r"\\toprule\s*\\hline", raw) or re.search(r"\\hline\\hline\s*\\bottomrule", raw):
        fail("a table stacks booktabs and \\hline rules")
    for text in ("the tabulated contour gives", "the contour is $0.0787$",
                 "mean the bisection floor was rejected", "that case is the mixed class"):
        if text in si:
            fail(f"retired phrasing in si: {text!r}")
    # the front matter states how many notes and tables the SI holds; it said
    # Tables 1--15 for a supplement with sixteen until 2026-10-16
    front = MR.one(DOCS["si"].parent, "*_si.tex")
    aux = MR.one(DOCS["si"].parent, "*_si.tex").with_suffix(".aux")
    if front.exists() and aux.exists():
        ft = " ".join(front.read_text().split())
        ax = aux.read_text()
        ntab = max(int(n) for n in re.findall(r"newlabel\{tab:[^}]*\}\{\{(\d+)\}", ax))
        nsec = len(re.findall(r"\\contentsline \{section\}",
                              front.with_suffix(".toc").read_text()))
        if f"Supplementary Notes 1--{nsec}, Supplementary Tables 1--{ntab} " not in ft:
            fail(f"SI front matter does not state Notes 1--{nsec} and Tables 1--{ntab}")
    if "\\cite{DongZenodo}" not in raw + (front.read_text() if front.exists() else ""):
        fail("the SI names deposited scripts but does not cite the deposit")
    if not bad:
        print(f"  consistency 2026-10-16: propositions S1-S{len(nums)} in order; "
              f"methanol multiplicity and joint contour tied to code")
    return bad


def check_bounds_2026_10_21():
    """Deposit-only regressions for the external audit of 2026-10-21.

    * Single-record bounds are the exact rho = -1 quantile of the ENDPOINT.  The
      shortcut F_min - 1.645 sd(F_obs) is right only on the open branch; for the
      two closed-branch Hong 2007 records it drew whiskers 0.06 and 0.10 too
      short.  Both records are regression cases: the exact bound must differ
      from the shortcut there, agree with Monte Carlo draws, and be what the
      figure data draw.
    * endpoint_exact admits a large-commitment limit reached with equality.
    * A bypass can raise the endpoint (the main text says so only of the yeast
      pairs).
    """
    import bounds_uncertainty as BU
    import corpus
    import network_geometry as N
    from pathlib import Path as _P
    bad = 0

    def fail(msg):
        nonlocal bad
        print(f"  FAIL bounds 2026-10-21: {msg}")
        bad += 1
    s = corpus.single_condition()
    hong = s[s.source == "Hong2007"]
    fig = (_P(__file__).resolve().parent.parent / "figures" / "tikz" / "data" / "fa_rays_b.tex")
    figtxt = fig.read_text() if fig.exists() else ""
    figtxt += (fig.parent / "fa_rays_a.tex").read_text() if (fig.parent / "fa_rays_a.tex").exists() else ""
    for r in hong.itertuples():
        fmin, _, _, closed = F_min_exact(r.K_HT, r.K_DT)
        x, y = r.K_HT_se / r.K_HT, GSC * r.K_DT_se / r.K_DT
        short = fmin - 1.645 * (x + y)
        ex = BU.exact_bound(r.K_HT, r.K_HT_se, r.K_DT, r.K_DT_se)
        mc = BU.single_bound(r.K_HT, r.K_HT_se, r.K_DT, r.K_DT_se, rho=-0.9999999)
        if not closed or not ex < short - 0.05:
            fail(f"Hong 2007 {r.variant}: the closed-branch regression case no longer differs "
                 f"from the open-branch shortcut ({ex:+.4f} vs {short:+.4f})")
        if abs(ex - mc) > 2e-3:
            fail(f"Hong 2007 {r.variant}: exact bound {ex:+.4f} disagrees with draws {mc:+.4f}")
        if figtxt and f"(axis cs:{ex:.5f}," not in figtxt:
            fail(f"Hong 2007 {r.variant}: the figure does not draw the exact bound {ex:.5f}")
    # every single record: exact quantile agrees with seeded Monte Carlo draws
    worst = max(abs(BU.exact_bound(r.K_HT, r.K_HT_se, r.K_DT, r.K_DT_se)
                    - BU.single_bound(r.K_HT, r.K_HT_se, r.K_DT, r.K_DT_se, rho=-0.9999999))
                for r in s.itertuples())
    if worst > 2e-3:
        fail(f"exact single-record bounds depart from Monte Carlo by {worst:.4f}")
    e, prof, lim = N.equality_limit_case()
    if not (np.isfinite(e) and abs(e - lim) < 1e-9 and abs(prof - lim) < 1e-4):
        fail(f"equality-limit endpoint {e} (profile {prof}, limit {lim})")
    nm = dict(N.bypass_nonmonotone_case())
    if not (nm[0.1] > nm[0.0] and nm[1.0] > nm[0.1] and nm[10.0] < nm[0.0]):
        fail(f"the non-monotone bypass example changed: {nm}")
    # the yeast endpoint does fall with phi over the range used
    ph = np.linspace(0.0, 1.0, 21)
    for r_ in (1.31, 1.57):
        v = [N.endpoint_exact(7.13, 1.73, p_, r_) for p_ in ph]
        if not all(b_ < a_ for a_, b_ in zip(v, v[1:])):
            fail(f"the yeast bypass endpoint is not decreasing in phi at r = {r_}")
    if not bad:
        print(f"  bounds 2026-10-21: exact single-record quantiles (max |exact - MC| {worst:.1e}); "
              "Hong 2007 regression cases drawn; equality limit and non-monotone bypass reproduced")
    return bad


def check_text_2026_10_21():
    """Text guards for the external audit of 2026-10-21 (needs the documents)."""
    import design_power as DP
    import network_geometry as N
    bad = 0
    main = " ".join(DOCS["main"].read_text().split())
    si = " ".join(DOCS["si"].read_text().split())

    def need(doc, text, what):
        nonlocal bad
        if text not in (main if doc == "main" else si):
            print(f"  FAIL text 2026-10-21: {doc} lacks {what}: {text!r}")
            bad += 1
    st = DP.benchmark_stats()
    need("si", f"${100 * st['med_H']:.2f}\\%$ for $\\KHT$ and ${100 * st['med_D']:.2f}\\%$ for $\\KDT$ "
               f"at the median of the ${st['n']}$", "the median published precision")
    need("si", f"zero for {st['n_zero']} of the {st['n_adm']} admissible", "the zero-gap count")
    need("si", f"${st['gaps'].max():.3f}$, or ${st['gaps'].max() / abs(F0):.1f}|\\Fz|$", "the largest gap")
    need("si", f"Proposition~S2 on all {st['n_adm']} admissible temperature-resolved records",
         "the verification scope")
    for stale in ("72 of the 93", "all 93 admissible", "$0.99\\%$ for $\\KHT$"):
        if stale in si:
            print(f"  FAIL text 2026-10-21: the SI keeps a statistic of the old corpus: {stale!r}")
            bad += 1
    nm = dict(N.bypass_nonmonotone_case())
    need("si", f"it rises from ${nm[0.0]:+.3f}$ at $\\phi=0$ to ${nm[0.1]:+.3f}$ at $\\phi=0.1$ and "
               f"${nm[1.0]:+.3f}$ at $\\phi=1$ before falling, to ${nm[10.0]:+.3f}$ at $\\phi=10$",
         "the non-monotone bypass example")
    if "The endpoint falls smoothly with $\\phi$, with no threshold" in main:
        print("  FAIL text 2026-10-21: main states bypass monotonicity as general")
        bad += 1
    need("main", "That is a property of these pairs, not of bypasses", "the bypass qualification")
    # the synthetic validation
    sv = pd.read_csv("../results/synthetic_validation.csv")
    a32 = sv[(sv.truth == "above") & (sv.c == 32)].iloc[0]
    a2 = sv[(sv.truth == "above") & (sv.c == 2)].iloc[0]
    mod = sv[sv.truth == "model"]
    if not ((mod.excl_1 == 0).all() and (mod.excl_03 == 0).all()
            and (sv.cover_1 == 1).all() and (sv.cover_03 == 1).all()):
        print("  FAIL text 2026-10-21: the synthetic validation no longer has zero false "
              "exclusions and full coverage, which both documents state")
        bad += 1
    for doc in ("main", "si"):
        need(doc, f"${100 * a32.excl_1:.1f}\\%$ of experiments at $1\\%$ errors", "the synthetic power")
        need(doc, f"between ${a2.cal_lo:+.3f}$ and ${a2.cal_hi:+.3f}$ in $90\\%$ of experiments",
             "the calibrated synthetic interval")
    need("main", f"endpoint at ${a2.endpoint:+.3f}$, below $\\Fz$", "the synthetic masked endpoint")
    # SM-only references: every key the SI cites and the main text does not must
    # be listed in the main text's Supplemental Material statement (18 census
    # sources were missing until 2026-10-21)
    raw_m, raw_s = DOCS["main"].read_text(), DOCS["si"].read_text()
    keys = lambda t: {k.strip() for g in re.findall(r"\\cite[pt]?\*?\{([^}]*)\}", t) for k in g.split(",")}
    i = raw_m.find("The Supplemental Material also cites")
    j = raw_m.find("which do not appear elsewhere", i)
    listed = keys(raw_m[i:j]) if i >= 0 and j > i else set()
    elsewhere = keys(raw_m[:i] + raw_m[j:]) if i >= 0 else keys(raw_m)
    missing = keys(raw_s) - elsewhere - listed
    extra = listed & elsewhere
    if missing:
        print(f"  FAIL text 2026-10-21: SI-only references not listed in the main text: {sorted(missing)}")
        bad += 1
    if extra:
        print(f"  FAIL text 2026-10-21: listed as SI-only but cited in the main text: {sorted(extra)}")
        bad += 1
    # the census boundary and the repeat-condition caveat
    need("si", "On 17 October 2026 we ran five queries against Europe PMC", "the census date")
    need("si", "the pairing of their H/T and D/T effects with the other two determinations",
         "the Sen 2011 condition caveat")
    need("main", "We use it as a benchmark because the offsets it can attain are known exactly",
         "the benchmark role of the tunneling model")
    if not bad:
        print("  text 2026-10-21: corpus statistics, bypass scope, synthetic validation, "
              "SM-only references, census boundary")
    return bad


def check_tables_2026_10_22():
    """The bypass table of Supplementary Note 3, row by row (table audit of
    2026-10-22: it still said "three systems" and omitted three of the six
    nominal exclusions after the census).  Every row with an observed secondary
    ratio is recomputed: L_H, R, the endpoint at phi = 0 and phi*."""
    import network_geometry as N
    import reference_asymmetry as RA
    import corpus
    bad = 0
    si = DOCS["si"].read_text()
    a = si.index("\\label{tab:bypass}")
    body = si[a:si.index("\\end{tabular}", a)]
    s = corpus.series()
    sc = corpus.single_condition()

    def rec(fam, var, T):
        g = s[(s.family == fam) & (s.variant == var) & (s.T_C == T)]
        if len(g) == 0:
            g = sc[(sc.family == fam) & (sc.variant == var) & (sc.T_C == T)]
        return float(g.K_HT.iloc[0]), float(g.K_DT.iloc[0])
    rows = [("Yeast ADH, $p$-Cl-benzyl alcohol, $25\\,^{\\circ}$C", "YADH", "wild type (p-chlorobenzyl alcohol)", 25.0),
            ("Horse liver F93W:I224G", "LADH", "F93W-I224G", 25.0),
            ("Horse liver F93W:I224V", "LADH", "F93W-I224V", 25.0),
            ("Horse liver F93W, $40\\%$ MeOH, $3\\,^{\\circ}$C", "LADH", "F93W (40% MeOH)", 3.0),
            ("Horse liver F93W, $60\\%$ MeOH, $3\\,^{\\circ}$C", "LADH", "F93W (60% MeOH)", 3.0)]
    for lab, fam, var, T in rows:
        kh, kd = rec(fam, var, T)
        r = RA.r_at(fam, var, T)
        e = N.endpoint_exact(kh, kd, 0.0, r)
        ph = N.bypass_to_destroy(kh, kd, F0, r)
        want = (f"{lab} & ${(kh - 1) / (kd - 1):.2f}$ & ${r:.2f}$ & ${e:+.4f}$ & ${ph:.3f}$")
        if want not in body:
            print(f"  FAIL tables 2026-10-22: bypass row differs from the code: {want!r}")
            bad += 1
    flat = " ".join(si.split())
    if "Three systems in the record have an exclusion" in flat:
        print("  FAIL tables 2026-10-22: the bypass note still counts three exclusions")
        bad += 1
    if "Six units in the record have a nominal exclusion" not in flat:
        print("  FAIL tables 2026-10-22: the bypass note does not count the six exclusions")
        bad += 1
    if not bad:
        print("  tables 2026-10-22: bypass table rows for the census exclusions match the code")
    return bad


def check_equations_2026_10_22():
    """Equation audit of 2026-10-22: every displayed equation re-derived
    (equation_checks.py), and the numbers that sat beside equations without
    code recomputed from the functions that now produce them."""
    import contextlib
    import io
    import equation_checks as EC
    import curvature as CV
    import completion as CO
    import atlas_general as AG
    import reference_asymmetry as RA
    bad = 0
    EC.FAILS.clear()
    with contextlib.redirect_stdout(io.StringIO()):
        EC.proofs(); EC.networks(); EC.binding_and_protocol(); EC.completion_and_rest()
    if EC.FAILS:
        print(f"  FAIL equations 2026-10-22: equation_checks fails {EC.FAILS}")
        bad += len(EC.FAILS)
    ex, (nv, nk, worst, tot), npos = CV.series_steps()
    (k41, v41), (kc, vc, knc, fnc) = CV.sign_and_binding_draws()
    (e1, q1), (e2, m2), err3 = CV.remaining_draw_checks()
    ntr, sw, big = AG.parallel_case()
    if v41 or vc or npos or e1 or e2 or err3 > 1e-8:
        print("  FAIL equations 2026-10-22: a random-draw check the supplement calls exceptionless has exceptions")
        bad += 1
    if not HAVE_DOCS:
        return bad
    si = " ".join(DOCS["si"].read_text().split())

    def need(text, what):
        nonlocal bad
        if text not in si:
            print(f"  FAIL equations 2026-10-22: si lacks {what}: {text!r}")
            bad += 1
    need(f"the observed offset is ${ex:.3f}$", "the two-step series example")
    need(f"exceeds $F^{{\\mathrm{{int}}}}$ in ${100 * nv / nk:.1f}\\%$, by as much as ${worst:.3f}$ at a total offset of ${tot:+.3f}$",
         "the separate-pairs counterexample")
    need(f"$400$ random draws of such a mechanism, ${sw}$ gave an offset that changed sign along the ray, ${big}$ of",
         "the parallel-step count")
    need(f"Verified on ${k41 / 1e5:.1f}\\times10^{{5}}$ random admissible draws with $r\\in[1,3]$",
         "the Corollary S4.1 draw count")
    need(f"Verified on ${kc / 1e5:.1f}\\times10^{{5}}$ random draws satisfying the corrected criterion",
         "the Proposition S7 draw count")
    need(f"the bound fails in ${100 * fnc / knc:.0f}\\%$", "the Proposition S7 failure rate")
    need(f"which with the same mapping gives $a={CO.klinman_a(1.3, 3.4):.2f}$ at $1.3$ and a range of "
         f"${CO.klinman_a(1.5, 3.4):.1f}$ to ${CO.klinman_a(1.1, 3.4):.1f}$", "the Klinman conversion at 3.4")
    ry = RA.r_at("YADH", "wild type", 25.0)
    rl = RA.r_at("LADH", "F93W", 25.0)
    need(f"$ar={6.13 * ry:.2f}$ against $\\gSC b={GSC * 0.73:.2f}$, and for horse liver F93W "
         f"${6.755 * rl:.2f}$ against ${GSC * 0.858:.2f}$", "the open-branch criterion values")
    for stale in ("by as much as $0.395$ at a total offset of $+0.050$", "$4$ gave an offset that changed sign",
                  "on about a third of draws, exactly as it does", "$ar=8.04$"):
        if stale in si:
            print(f"  FAIL equations 2026-10-22: retired value in the SI: {stale!r}")
            bad += 1
    if not bad:
        print(f"  equations 2026-10-22: {len(EC.FAILS) or 'all'} identities re-derived; series, sign, binding, "
              f"parallel ({sw}/{big}) and Klinman numbers from code")
    return bad


_NEEDS_DOCS = {
    "check_audit_2026_10_15", "check_consistency_2026_10_16",
    "check_light_example",
    "check_si_figrefs", "check_si_eqrefs", "check_si_pointers",
    "check_general_census", "check_corpus", "check_provenance_split",
    "check_cited_scripts", "check_no_ai_mentions", "check_methods_sources",
    "check_precision_limited_count", "check_ladh_zero_claim",
    "check_repo_metadata", "check_titles", "check_reference_asymmetry",
    "check_promoting_mode", "check_corpus_consequences", "check_scope_corrections",
    "check_source_fidelity", "check_readiness_scope", "check_generated_tables",
    "check_external_audit_2026_10_12",
    "check_source_review_2026_10_13", "check_text_2026_10_21", "check_tables_2026_10_22",
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
                check_ladh_zero_claim, check_ecdhfr_assignment, check_light_example,
                check_envelope, check_joint_counterexamples,
                check_repo_metadata, check_titles, check_reference_asymmetry,
                check_promoting_mode, check_corpus_consequences,
                check_scope_corrections, check_source_fidelity,
                check_readiness_scope, check_generated_tables,
                check_external_audit_2026_10_12, check_source_review_2026_10_13,
                check_figure_data, check_audit_2026_10_15,
                check_tables_2026_10_15, check_consistency_2026_10_16,
                check_bounds_2026_10_21, check_text_2026_10_21,
                check_tables_2026_10_22, check_equations_2026_10_22):
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
