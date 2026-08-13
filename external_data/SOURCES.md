# Primary sources

Every kinetic measurement analysed here comes from a published paper. **No
publisher PDF is redistributed in this repository.** The extracted numerical
values are facts and are redistributed with attribution in
`data/trinomial_benchmark.csv`; the documents they came from are not ours to
share.

Run `./fetch_sources.sh` to retrieve the openly accessible supplements. Two
sources are paywalled and must be downloaded manually; see below.

## Matched H/T and D/T series (the benchmark)

| Source | DOI | Contributes | Availability |
|---|---|---|---|
| Islam et al., *PLoS ONE* **13**, e0196506 (2018) | `10.1371/journal.pone.0196506` | hsTSase hydride + proton | CC-BY 4.0 |
| Wang et al., *Int. J. Mol. Sci.* **16**, 7304 (2015) | `10.3390/ijms16047304` | ecTSase Y209W proton | CC-BY, PMC4425018 |
| Singh et al., *JACS* **136**, 2575 (2014) | `10.1021/ja411998h` | 4 ecDHFR variants, 5 temperatures | PMC3985941, OA |
| Li et al., *Biochemistry* **60**, 3822 (2021) | `10.1021/acs.biochem.1c00558` | 5 hsDHFR variants | PMC8697555, OA |
| Francis et al., *JACS* **136**, 8333 (2014) | `10.1021/ja501936d` | light vs heavy ecDHFR | PMC4063187, OA |
| Pagano et al., *ACS Catal.* **9**, 11199 (2019) | `10.1021/acscatal.9b03345` | FDH V123A, V123G | **paywalled** |

## Sources examined and excluded

| Source | DOI | Why excluded |
|---|---|---|
| Stojković et al., *JACS* **134**, 1738 (2012) | `10.1021/ja209425w` | Publishes only Northrop-derived *intrinsic* effects; the observed H/T and D/T values appear nowhere. Those intrinsic pairs are pinned to `F = 0` by construction and carry no information about the discriminating invariant. |
| Austin-Kloppe et al., *Chem. Sci.* **17**, 12440 (2026) | `10.1039/D6SC01847E` | 38 solution-phase model reactions, but deuterium only. No second isotope pair, so the offset cannot be formed. |
| Kohen et al., *Nature* **399**, 496 (1999) | `10.1038/20981` | ht-ADH. Its well-known anomalous Swain–Schaad result is the **secondary** exponent, outside the model class treated here. Primary data are published as Arrhenius ratios only, which constrain the offset to about ±0.5. |

## Manual downloads

Two documents cannot be fetched programmatically. Place them in
`external_data/si/manual/` under the filenames shown, then rerun
`python analysis/build_trinomial.py`.

1. **Pagano et al. 2019** — <https://doi.org/10.1021/acscatal.9b03345>
   Supporting Information PDF → `cs9b03345_si_001.pdf`
   (PMC8118594 exists but is not in the open-access subset.)

2. *(Optional, not used in the benchmark)* Stojković et al. 2012 and
   Austin-Kloppe et al. 2026, if you wish to reproduce the exclusion analysis
   in Supplemental Material §S4.2.

Without the Pagano SI the build produces 14 series and 63 records instead of 16
and 73; every other result is unaffected.

## Transcription discipline

`analysis/build_trinomial.py` parses the supplementary tables from the PDF text
layer rather than retyping them, and runs an audit confirming that every parsed
number appears verbatim in its source. With all sources present the audit covers
244 values across four documents and passes with zero misses.

One judgment is documented rather than automated: in Francis et al. the light
and heavy enzyme rows are distinguished only by colour in the original, which is
lost in the text layer. Rows are assigned by order, the temperature-labelled row
being the light enzyme, and the unpaired 10 °C row is assigned to the light
series on interpolation (light runs 3.10 → 4.76 across that gap, heavy
2.40 → 3.29, and the value is 3.69).

## Records transcribed by hand

Three sources predate machine-readable supplements and are scanned page images
whose text layer is unreliable, so their values were transcribed by hand from the
printed tables rather than parsed. They are excluded from the automated audit
above and are held in their own files so the parsed and typed records never mix:

| Source | File | Table | Records |
|---|---|---|---|
| Grant & Klinman, *Biochemistry* **28**, 6597 (1989), [10.1021/bi00442a010](https://doi.org/10.1021/bi00442a010) | `data/bsao_grant1989.csv` | Table I averages, 25 °C | 1 |
| Jonsson, Edmondson & Klinman, *Biochemistry* **33**, 14871 (1994), [10.1021/bi00253a026](https://doi.org/10.1021/bi00253a026) | `data/trinomial_benchmark.csv` (family `MAOB`) | Tables 1 and 3 | 10 |
| Bahnson et al. 1993/1997 | `data/ladh_adh_primary.csv` | Table I / Table 1 | 9 |

Both amine oxidase sources report the primary H/T and D/T effects together with
the secondary effects on the same references; the latter fix the ratio `r` of
the two commitments used in the unequal-reference check of `analysis/curvature.py`
(`r = 1.14` at 25 °C for both).

Because these records are not produced by `build_trinomial.py`, re-running that
script rebuilds only the parsed series. It compares the rebuilt series count with
the shipped benchmark and writes to `trinomial_benchmark_partial.csv` rather than
overwriting if it finds fewer, so the hand-transcribed rows cannot be lost by a
partial rebuild.

Two source classes were examined and rejected, for reasons worth recording:

- **Deuterium-only.** Soybean lipoxygenase (Rickert & Klinman 1999; Knapp et al.
  2002) and the ADH flexibility work (Kohen & Klinman 2000) report H/D only. No
  second isotope pair exists, so the offset cannot be formed.
- **Northrop-derived intrinsic values.** Peptidylglycine α-amidating enzyme
  (Francisco et al. 1998), dopamine β-monooxygenase (Miller & Klinman 1985) and
  the ecDHFR I14 series (Stojković et al. 2012) publish intrinsic effects
  obtained by *imposing* the Swain–Schaad relation. Those values are pinned to
  `F = 0` by construction, so using them would be circular.
