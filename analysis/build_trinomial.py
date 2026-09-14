"""Extract matched competitive H/T and D/T KIE series from open-access SI PDFs.

The offset statistic F = ln K_HT - gamma_SC ln K_DT needs BOTH isotope pairs
measured against a common tritium reference at the same condition.  The v2
benchmark has 32 competitive records but zero D/T, so it cannot supply F.  This
script builds the matched set from primary sources.

Everything is parsed from the PDF text rather than retyped, so the transcription
is reproducible and auditable.  Sanity checks reject any row violating
K_HT > K_DT > 1 or with a non-positive standard error.
"""
from __future__ import annotations

import pathlib
import re
import subprocess

import pandas as pd

SI = pathlib.Path("../external_data/si")
OUT = pathlib.Path("../data")
OUT.mkdir(exist_ok=True)

# value +- error, tolerating the various unicode minus/plusminus forms
NUM = r"([0-9]+\.[0-9]+)\s*(?:±|\+-|\+/-)\s*([0-9]+\.[0-9]+)"
ROW4 = re.compile(r"^\s*([0-9]{1,2})\s+" + NUM + r"\s+" + NUM)


# Sources whose records are transcribed by hand from scanned pages and so are never
# produced by this parse.  Must agree with HAND in audit_numbers.check_provenance_split.
HAND_TRANSCRIBED = {"10.1021/bi00253a026"}


def _xlsx_cells(path):
    """Numeric cells of the first worksheet, as {(row, col): float}.

    xlsx is a zip of XML, so this needs no third-party reader and keeps the
    dependency list to numpy/scipy/pandas/sympy.
    """
    import zipfile
    with zipfile.ZipFile(path) as z:
        sheet = next(n for n in z.namelist()
                     if re.match(r"xl/worksheets/sheet\d+\.xml$", n))
        xml = z.read(sheet).decode("utf8", "replace")
    cells = {}
    for m in re.finditer(r'<c r="([A-Z]+)(\d+)"([^>]*)>(.*?)</c>', xml, re.S):
        col, row, attr, body = m.group(1), int(m.group(2)), m.group(3), m.group(4)
        v = re.search(r"<v>(.*?)</v>", body, re.S)
        if not v or 't="s"' in attr:      # skip shared strings: headers only
            continue
        try:
            cells[(row, col)] = float(v.group(1))
        except ValueError:
            pass
    return cells


def _replicates(cells, r0, r1, avg, sd, tol=5e-6):
    """The replicate set inside rows [r0, r1] reproducing `avg` and `sd`.

    Ref. Islam 2018 ships a hand-made spreadsheet: each temperature block lists
    its replicates, but the columns they sit in drift between blocks, so the
    replicates cannot be read positionally.  They are instead identified by the
    property that defines them -- their mean and population standard deviation
    must equal the values the paper publishes.  That makes the extraction
    self-validating: a wrong set does not reproduce both moments.  Returning the
    set also yields n, which the standard error needs.
    """
    import itertools, statistics as st
    lo, hi = avg - 6 * sd, avg + 6 * sd
    cand = [v for (r, c), v in sorted(cells.items()) if r0 <= r <= r1
            and lo <= v <= hi]
    # the published mean, and occasionally one stray, sit inside the window
    for k in range(0, min(3, len(cand)) + 1):
        for drop in itertools.combinations(range(len(cand)), k):
            keep = [v for i, v in enumerate(cand) if i not in drop]
            if len(keep) < 2:
                continue
            if (abs(st.mean(keep) - avg) < tol
                    and abs(st.pstdev(keep) - sd) < tol):
                return sorted(keep)
    raise ValueError(f"replicates not resolved for avg={avg} sd={sd}")


# Islam 2018, hsTSase.  s004 is hydride transfer, s005 proton abstraction; each
# holds four temperature blocks.  (row, avg column, sd column) per isotope.
_HSTS = {
    "hydride": ("PMC5929524/pone.0196506.s004.xlsx",
                [(5, 3, 31), (15, 33, 62), (25, 64, 87), (35, 89, 107)],
                ("E", "F"), ("H", "I")),
    "proton":  ("PMC5929524/pone.0196506.s005.xlsx",
                [(5, 4, 24), (15, 25, 52), (25, 53, 69), (35, 70, 94)],
                ("D", "E"), ("G", "H")),
}


def hstsase_records():
    """The eight hsTSase records, parsed from the published spreadsheets.

    The paper reports a population standard deviation; the benchmark carries a
    standard error, so each is converted with the recovered replicate count,
    se = sd * sqrt(n / (n - 1)) / sqrt(n).
    """
    import math
    out = []
    for step, (rel, blocks, ht, dt) in _HSTS.items():
        path = SI / rel
        if not path.exists():
            return []
        cells = _xlsx_cells(path)
        for T, r0, r1 in blocks:
            row = {}
            for iso, (ac, sc) in (("H", ht), ("D", dt)):
                avg, sd = cells[(r0, ac)], cells[(r0, sc)]
                n = len(_replicates(cells, r0, r1, avg, sd))
                row[iso] = (avg, sd * math.sqrt(n / (n - 1)) / math.sqrt(n))
            out.append(dict(system="Homo sapiens thymidylate synthase",
                            family="TSase", variant="WT", step=step,
                            donor_atom="C", T_C=float(T),
                            K_HT=row["H"][0], K_HT_se=row["H"][1],
                            K_DT=row["D"][0], K_DT_se=row["D"][1],
                            source_DOI="10.1371/journal.pone.0196506",
                            PMCID="PMC5929524",
                            source_table="S4 Table" if step == "hydride"
                                         else "S5 Table"))
    return out


def ectsase_records():
    """The four ecTSase Y209W proton records, from Table S2 of Abeysinghe 2015."""
    path = SI / "PMC4425018/ijms-16-07304-s001.pdf"
    if not path.exists():
        return []
    txt = text_of(path)
    i = txt.find("Table S2.")
    out = []
    for line in txt[i:i + 1200].splitlines():
        m = ROW4.match(line)
        if not m:
            continue
        T, kh, sh, kd, sd = (float(m.group(1)), float(m.group(2)),
                             float(m.group(3)), float(m.group(4)),
                             float(m.group(5)))
        out.append(dict(system="Escherichia coli thymidylate synthase",
                        family="TSase", variant="Y209W", step="proton",
                        donor_atom="C", T_C=T, K_HT=kh, K_HT_se=sh,
                        K_DT=kd, K_DT_se=sd,
                        source_DOI="10.3390/ijms16047304",
                        PMCID="PMC4425018", source_table="Table S2"))
    return out


def text_of(pdf):
    return subprocess.run(["pdftotext", "-layout", str(pdf), "-"],
                          capture_output=True, text=True).stdout


def parse_tables(txt, table_re, source):
    """Yield (label, [(T_C, kht, se_ht, kdt, se_dt), ...]) per matching table."""
    lines = txt.splitlines()
    heads = [(i, m) for i, l in enumerate(lines) for m in [table_re.match(l.strip())] if m]
    for n, (i, m) in enumerate(heads):
        stop = heads[n + 1][0] if n + 1 < len(heads) else len(lines)
        rows = []
        for l in lines[i:stop]:
            r = ROW4.match(l)
            if r:
                T, kht, sh, kdt, sd = (float(x) for x in r.groups())
                rows.append((T, kht, sh, kdt, sd))
        if rows:
            yield m.group(1).strip(), rows



# PDF source of each parsed series.  The spreadsheet source (Islam 2018) is not
# here: its values are recovered from replicate columns, not printed as value±error.
PDF_SOURCES = {
    "10.1021/ja411998h": "PMC3985941/ja411998h_si_001.pdf",
    "10.1021/acs.biochem.1c00558": "PMC8697555/bi1c00558_si_001.pdf",
    "10.1021/acscatal.9b03345": "manual/cs9b03345_si_001.pdf",
    "10.1021/ja501936d": "PMC4063187/ja501936d_si_001.pdf",
    "10.3390/ijms16047304": "PMC4425018/ijms-16-07304-s001.pdf",
}


def _num(x):
    """Shortest decimal form, which is how a trailing zero is lost in the CSV."""
    return format(float(x), ".10g")


def transcription_audit(bench=OUT / "trinomial_benchmark.csv"):
    """Every shipped PDF-derived value must appear WITH ITS ERROR in its source.

    The manuscript, the SI and SOURCES.md have long stated that an audit
    confirms every parsed number appears verbatim in its source.  Until
    2026-09-14 no deposited code performed that check; it existed only as prose.
    This is it.  It reads the SHIPPED benchmark -- the numbers the analysis
    uses -- not a fresh parse, so a hand edit to the CSV is caught even though
    the parser would not reproduce it.

    Each isotope effect is matched together with its uncertainty, e.g.
    "3.07 ± 0.02", tolerating trailing zeros the CSV drops (4.10 -> 4.1) and any
    spacing around the ±.  Matching the pair rather than the bare number means a
    value attached to the wrong error fails, which a presence test would miss.

    What it does NOT establish: that the value sits in the right ROW or COLUMN
    of the table.  A swap of the H/T and D/T columns, or of two temperatures,
    would pass.  Row and column assignment were checked against the rendered
    page images during the 2026-09 reference review, not by this function.
    """
    d = pd.read_csv(bench)
    pairs = found = 0
    misses = []
    for doi, rel in PDF_SOURCES.items():
        pdf = SI / rel
        if not pdf.exists():
            print(f"  transcription audit: {rel} not present; skipped")
            continue
        txt = text_of(pdf).replace("\u00b1", "±")
        # Collect every "number ± number" the source prints and compare
        # NUMERICALLY.  Matching formatted strings instead is brittle: the CSV
        # drops trailing zeros, so a printed "2.00 ± 0.05" does not contain "2".
        printed = {(float(a), float(b)) for a, b in
                   re.findall(r"(\d+\.?\d*)\s*±\s*(\d+\.?\d*)", txt)}
        for _, r in d[d["source_DOI"] == doi].iterrows():
            for val, err in (("K_HT", "K_HT_se"), ("K_DT", "K_DT_se")):
                v, e = float(r[val]), float(r[err])
                pairs += 1
                if any(abs(v - a) < 5e-10 and abs(e - b) < 5e-10
                       for a, b in printed):
                    found += 1
                else:
                    misses.append(f"{r['variant']} {r['T_C']:g}C "
                                  f"{val}={_num(v)}±{_num(e)}")
    print(f"\n  transcription audit: {found} of {pairs} value±error pairs "
          f"({2 * pairs} values) found verbatim in their source PDFs; "
          f"{len(misses)} missing")
    for m in misses:
        print(f"    MISSING  {m}")
    return len(misses)

def main():
    recs = []

    # ---- Singh et al. 2014, ecDHFR network variants ---------------------
    p = SI / "PMC3985941" / "ja411998h_si_001.pdf"
    if p.exists():
        rx = re.compile(r"Table S\d+:\s*Observed and intrinsic V/K KIEs for (.+)")
        for label, rows in parse_tables(text_of(p), rx, "Singh2014"):
            for T, kht, sh, kdt, sd in rows:
                recs.append(dict(system="Escherichia coli dihydrofolate reductase",
                                 family="ecDHFR", variant=label, step="hydride",
                                 donor_atom="C", T_C=T, K_HT=kht, K_HT_se=sh,
                                 K_DT=kdt, K_DT_se=sd,
                                 source_DOI="10.1021/ja411998h", PMCID="PMC3985941",
                                 source_table="Table S1-S4"))

    # ---- Reyes et al. 2021, hsDHFR loop/F32 variants ---------------------
    p = SI / "PMC8697555" / "bi1c00558_si_001.pdf"
    if p.exists():
        # the source writes "Delta" with a private-use glyph that pdftotext
        # renders as U+F004; map it to "d" to match the v2 naming convention
        rx = re.compile(r"Table S\d+:\s*(?:Intrinsic and Observed|Observed and intrinsic)"
                        r"\s*KIEs of\s+(.+?)\s+h(?:s)?DHFR\s+at pH.*")
        for label, rows in parse_tables(text_of(p), rx, "Reyes2021"):
            label = re.sub(r"[-Δ]", "d", label).strip()
            for T, kht, sh, kdt, sd in rows:
                recs.append(dict(system="Homo sapiens dihydrofolate reductase",
                                 family="hsDHFR", variant=label, step="hydride",
                                 donor_atom="C", T_C=T, K_HT=kht, K_HT_se=sh,
                                 K_DT=kdt, K_DT_se=sd,
                                 source_DOI="10.1021/acs.biochem.1c00558",
                                 PMCID="PMC8697555", source_table="Table S2-S6"))

    # ---- Pagano et al. 2019, formate dehydrogenase (third enzyme family) -
    # Table S2 is a two-panel layout: T | V123G(H/T, D/T, H/Tint) | V123A(...)
    p = SI / "manual" / "cs9b03345_si_001.pdf"
    if p.exists():
        pair = r"([0-9]+(?:\.[0-9]+)?)\s*±\s*([0-9]+(?:\.[0-9]+)?)"
        rx = re.compile(r"^\s*([0-9]{1,2})\s+" + r"\s+".join([pair] * 6))
        for line in text_of(p).splitlines():
            m = rx.match(line)
            if not m:
                continue
            g = [float(x) for x in m.groups()]
            T = g[0]
            for var, off in (("V123G", 1), ("V123A", 7)):
                recs.append(dict(system="Candida boidinii formate dehydrogenase",
                                 family="FDH", variant=var, step="hydride",
                                 donor_atom="C", T_C=T,
                                 K_HT=g[off], K_HT_se=g[off + 1],
                                 K_DT=g[off + 2], K_DT_se=g[off + 3],
                                 source_DOI="10.1021/acscatal.9b03345",
                                 PMCID="PMC8118594", source_table="Table S2"))

    # ---- Francis et al. 2014, light vs heavy (mass-modulated) ecDHFR -----
    # Two rows per temperature: the temperature-labelled row is l-DHFR, the
    # continuation row h-DHFR.  10 C has a single row; linear interpolation of
    # the neighbouring points (light 3.1->4.76 vs heavy 2.40->3.29) places the
    # value 3.69 unambiguously on the light series.
    p = SI / "PMC4063187" / "ja501936d_si_001.pdf"
    if p.exists():
        pair = r"([0-9]+(?:\.[0-9]+)?)\s*±\s*([0-9]+(?:\.[0-9]+)?)"
        head = re.compile(r"^\s*([0-9]{1,2})\s+" + pair + r"\s+" + pair)
        cont = re.compile(r"^\s{6,}" + pair + r"\s+" + pair)
        cur_T = None
        for line in text_of(p).splitlines():
            m = head.match(line)
            if m:
                cur_T = float(m.group(1))
                v = [float(x) for x in m.groups()[1:]]
                recs.append(dict(system="Escherichia coli dihydrofolate reductase",
                                 family="ecDHFR", variant="light enzyme",
                                 step="hydride", donor_atom="C", T_C=cur_T,
                                 K_HT=v[0], K_HT_se=v[1], K_DT=v[2], K_DT_se=v[3],
                                 source_DOI="10.1021/ja501936d",
                                 PMCID="PMC4063187", source_table="Table S1"))
                continue
            m = cont.match(line)
            if m and cur_T is not None:
                v = [float(x) for x in m.groups()]
                recs.append(dict(system="Escherichia coli dihydrofolate reductase",
                                 family="ecDHFR", variant="heavy enzyme",
                                 step="hydride", donor_atom="C", T_C=cur_T,
                                 K_HT=v[0], K_HT_se=v[1], K_DT=v[2], K_DT_se=v[3],
                                 source_DOI="10.1021/ja501936d",
                                 PMCID="PMC4063187", source_table="Table S1"))
                cur_T = None

    # ---- thymidylate synthase, parsed from the published supplements ------
    # Previously these twelve records were re-read from curated CSVs in this
    # repository, which is no evidence of transcription: rebuilding them
    # reproduced our own file.  They are now parsed from the publishers'
    # documents like every other source.
    recs.extend(hstsase_records())
    recs.extend(ectsase_records())

    d = pd.DataFrame(recs).drop_duplicates(
        subset=["family", "variant", "step", "T_C"]).reset_index(drop=True)

    bad = d[~((d.K_HT > d.K_DT) & (d.K_DT > 1.0) &
              (d.K_HT_se > 0) & (d.K_DT_se > 0))]
    if len(bad):
        print("REJECTED rows failing K_HT > K_DT > 1 with positive errors:")
        print(bad.to_string(index=False))
        d = d.drop(bad.index)

    d = d.sort_values(["family", "system", "variant", "step", "T_C"])

    # Guard: the repository ships the complete benchmark.  If primary sources
    # are missing (they are not redistributable, see external_data/SOURCES.md)
    # this script would otherwise silently overwrite it with a partial set.
    target = OUT / "trinomial_benchmark.csv"
    n_series = d.groupby(["family", "system", "variant", "step"]).ngroups
    have = None
    if target.exists():
        have = pd.read_csv(target)
        keys = ["family", "system", "variant", "step"]
        n_have = have.groupby(keys).ngroups
        # The shipped benchmark also holds records transcribed by hand from
        # scanned pages (SOURCES.md), which no parse can produce.  Until
        # 2026-09-14 the comparison below counted those in its baseline, so a
        # complete, exact rebuild still printed "Primary sources are missing" --
        # on every correct run.  Compare against the parseable subset instead.
        hand = have[have["source_DOI"].isin(HAND_TRANSCRIBED)]
        parseable = have[~have["source_DOI"].isin(HAND_TRANSCRIBED)]
        n_parse = parseable.groupby(keys).ngroups
        if n_series < n_have:
            # Never overwrite the complete file with a parse that lacks the
            # hand-transcribed rows, whether or not the parse itself is complete.
            target = OUT / "trinomial_benchmark_partial.csv"
        if n_series < n_parse:
            print(f"\n  WARNING: rebuilt {len(d)} records / {n_series} series, but the")
            print(f"  shipped benchmark has {len(parseable)} / {n_parse} parseable ones.")
            print("  Primary sources are missing; run external_data/fetch_sources.sh")
            print("  and see SOURCES.md.")
    d.to_csv(target, index=False)
    if have is not None and n_series < n_have:
        if n_series < n_parse:
            pass  # warned above; a partial parse is not compared
        else:
            # Compare what was WRITTEN, read back, against the shipped parseable
            # rows.  Comparing the in-memory frame instead fails on dtype
            # differences in the merge keys and reports a false difference.
            new = pd.read_csv(target)
            k = keys + ["T_C"]
            cmp = parseable.merge(new, on=k, how="outer",
                                  suffixes=("_ship", "_new"), indicator=True)
            paired = bool((cmp["_merge"] == "both").all())
            same = paired and all(
                (cmp[f"{c}_ship"] - cmp[f"{c}_new"]).abs().max() < 1e-9
                for c in ("K_HT", "K_HT_se", "K_DT", "K_DT_se"))
            print(f"\n  Rebuilt all {len(d)} parseable records / {n_series} series"
                  + (" and they reproduce the shipped values exactly." if same
                     else ", but they DIFFER from the shipped values; inspect before use."))
            print(f"  The shipped benchmark also holds {len(hand)} hand-transcribed records"
                  f" / {n_have - n_parse} series")
            print("  that no parse produces (SOURCES.md), so this is written separately.")
        print(f"  Wrote {target.name}; the complete file is left intact.\n")

    print(f"matched (H/T, D/T) competitive records: {len(d)}")
    print(f"independent series: {d.groupby(['family','system','variant','step']).ngroups}")
    print()
    print(d.groupby(["family", "system", "variant", "step"])
           .agg(n_T=("T_C", "size"), T_min=("T_C", "min"), T_max=("T_C", "max"),
                KHT_min=("K_HT", "min"), KHT_max=("K_HT", "max"))
           .to_string())
    # name the file actually written: `target` is the partial file whenever the
    # guard above fired, and reporting the canonical name there wrongly tells the
    # operator the shipped benchmark was overwritten
    print(f"\n[written] {target}")


if __name__ == "__main__":
    import sys
    main()
    sys.exit(1 if transcription_audit() else 0)
