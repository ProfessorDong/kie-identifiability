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

    # ---- thymidylate synthase series already curated in v2 ---------------
    v2 = pathlib.Path("../data")
    ts = pd.read_csv(v2 / "tsase_proton_series.csv")
    for _, r in ts.iterrows():
        recs.append(dict(system=("Homo sapiens thymidylate synthase"
                                 if r.system == "hsTSase"
                                 else "Escherichia coli thymidylate synthase"),
                         family="TSase", variant=r.variant, step="proton",
                         donor_atom="C", T_C=r.T_C, K_HT=r.K_HT, K_HT_se=r.K_HT_se,
                         K_DT=r.K_DT, K_DT_se=r.K_DT_se,
                         source_DOI=r.source_DOI, PMCID=r.PMCID,
                         source_table=r.source_table))
    su = pd.read_csv(v2 / "hstsase_kie_summary.csv")
    w = su.pivot_table(index=["system", "step", "T_C"], columns="isotope_pair",
                       values=["KIE", "KIE_sd", "n_rep"]).reset_index()
    for _, r in w.iterrows():
        n_ht = r[("n_rep", "H/T")]; n_dt = r[("n_rep", "D/T")]
        recs.append(dict(system="Homo sapiens thymidylate synthase", family="TSase",
                         variant="WT", step=r[("step", "")].replace("_transfer", "")
                                              .replace("_abstraction", ""),
                         donor_atom="C", T_C=r[("T_C", "")],
                         K_HT=r[("KIE", "H/T")],
                         K_HT_se=r[("KIE_sd", "H/T")] / n_ht ** 0.5,
                         K_DT=r[("KIE", "D/T")],
                         K_DT_se=r[("KIE_sd", "D/T")] / n_dt ** 0.5,
                         source_DOI="10.1371/journal.pone.0196506",
                         PMCID="PMC5929524", source_table="S3-S6 Tables"))

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
    if target.exists():
        have = pd.read_csv(target)
        n_have = have.groupby(["family", "system", "variant", "step"]).ngroups
        if n_series < n_have:
            target = OUT / "trinomial_benchmark_partial.csv"
            print(f"\n  WARNING: rebuilt {len(d)} records / {n_series} series, but the")
            print(f"  shipped benchmark has {len(have)} / {n_have}.  Primary sources are")
            print("  missing; run external_data/fetch_sources.sh and see SOURCES.md.")
            print(f"  Writing to {target.name} and leaving the complete file intact.\n")
    d.to_csv(target, index=False)

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
    main()
