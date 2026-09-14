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
| Islam et al., *PLoS ONE* **13**, e0196506 (2018) | `10.1371/journal.pone.0196506` | hsTSase hydride + proton | PMC5929524, OA (S4/S5 Tables, .xlsx) |
| Abeysinghe & Kohen, *Int. J. Mol. Sci.* **16**, 7304 (2015) | `10.3390/ijms16047304` | ecTSase Y209W proton | PMC4425018, OA (Table S2) |
| Singh et al., *JACS* **136**, 2575 (2014) | `10.1021/ja411998h` | 4 ecDHFR variants, 5 temperatures | PMC3985941, OA |
| Li et al., *Biochemistry* **60**, 3822 (2021) | `10.1021/acs.biochem.1c00558` | 5 hsDHFR variants | PMC8697555, OA |
| Wang et al., *JACS* **136**, 8333 (2014) | `10.1021/ja501936d` | light vs heavy ecDHFR | PMC4063187, OA |
| Pagano et al., *ACS Catal.* **9**, 11199 (2019) | `10.1021/acscatal.9b03345` | FDH V123A, V123G | **paywalled** |

## Sources examined and excluded

| Source | DOI | Why excluded |
|---|---|---|
| Stojković et al., *JACS* **134**, 1738 (2012) | `10.1021/ja209425w` | Publishes only Northrop-derived *intrinsic* effects; the observed H/T and D/T values appear nowhere. Those intrinsic pairs are pinned to `F = 0` by construction and carry no information about the discriminating invariant. |
| Austin-Kloppe et al., *Chem. Sci.* **17**, 12440 (2026) | `10.1039/D6SC01847E` | 38 solution-phase hydride transfers (34 from NADH/NAD+ models, analyzed quantitatively), but deuterium only. No second isotope pair, so the offset cannot be formed. |
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

Without the Pagano SI the build produces 16 series and 73 records instead of 18
and 83; every other result is unaffected.

## Transcription discipline

`analysis/build_trinomial.py` parses the supplementary tables from the PDF text
layer rather than retyping them, and runs an audit confirming that every value in
the shipped benchmark from a PDF source appears verbatim in that source together
with its stated uncertainty. It reads the shipped file, not a fresh parse, so a
hand edit to the CSV is caught. With all sources present it covers 260 values
across five documents and passes with zero misses, and the script exits non-zero
on any miss. It does not establish row or column assignment: a swap of the H/T
and D/T columns within a row would pass. Those assignments were checked against
the rendered tables.

One judgment is documented rather than automated: in Wang et al. the light
and heavy enzyme rows are distinguished only by colour in the original, which is
lost in the text layer. Rows are assigned by order, the temperature-labelled row
being the light enzyme, and the unpaired 10 °C row is assigned to the light
series on interpolation (light runs 3.10 → 4.76 across that gap, heavy
2.40 → 3.29, and the value is 3.69). The assignment has since been confirmed
directly: rendering the table as an image preserves the colour, every paired row
is blue then red, and the 10 °C row is blue. All eleven records agree.

## Records transcribed by hand

Twenty-four records come from articles without machine-readable supplements.
Their values were transcribed by hand from the printed tables, or for one record
from the text, and each was checked against the rendered page image; most of
these articles are scans, and no text layer was relied on. They are outside the
automated audit above and are held in their own files (the monoamine oxidase
records sit in the benchmark file but are flagged by DOI in `build_trinomial.py`),
so parsed and typed records never mix. `analysis/corpus.py` lists every file of
single-condition records, and every script that needs them reads that list.

| Source | File | Where in the source | Records |
|---|---|---|---|
| Cha, Murray & Klinman, *Science* **243**, 1325 (1989), [10.1126/science.2646716](https://doi.org/10.1126/science.2646716) | `data/cha1989_yadh.csv` | Table 2: three determinations and their average; the average is the record | 1 |
| Bahnson et al., *Biochemistry* **32**, 5503 (1993), [10.1021/bi00072a003](https://doi.org/10.1021/bi00072a003) | `data/ladh_adh_primary.csv` | Table I | 5 |
| Bahnson et al., *Proc. Natl. Acad. Sci. USA* **94**, 12797 (1997) | `data/ladh_adh_primary.csv` | Table 1 | 4 |
| Tsai & Klinman, *Biochemistry* **40**, 2303 (2001), [10.1021/bi002075l](https://doi.org/10.1021/bi002075l) | `data/ladh_adh_primary.csv` | Tables 2 and 4, aqueous rows at 3 °C | 2 |
| Grant & Klinman, *Biochemistry* **28**, 6597 (1989), [10.1021/bi00442a010](https://doi.org/10.1021/bi00442a010) | `data/bsao_grant1989.csv` | Table I averages, 25 °C | 1 |
| Jonsson, Edmondson & Klinman, *Biochemistry* **33**, 14871 (1994), [10.1021/bi00253a026](https://doi.org/10.1021/bi00253a026) | `data/trinomial_benchmark.csv` (family `MAOB`) | Tables 1 and 3 | 10 |
| Agrawal, Hong, Mihai & Kohen, *Biochemistry* **43**, 1998 (2004), [10.1021/bi036124g](https://doi.org/10.1021/bi036124g) | `data/agrawal2004_ectsase.csv` | Results text, 20 °C: 6.91 ± 0.05 and 1.78 ± 0.02 (the other temperatures are plotted only and are not digitized) | 1 |

The Agrawal record was first excluded on the ground that all of that paper's
isotope effects are in figures. That overlooked the 20 °C values stated in the
text, and the record was admitted on 2026-09-14 under the rule applied to every
source: an observed matched pair is admitted when both effects are stated as
numbers with their uncertainties.

Both amine oxidase sources report the primary H/T and D/T effects together with
the secondary effects on the same references; the latter fix the ratio `r` of
the two commitments used in the unequal-reference check of `analysis/curvature.py`
(`r = 1.14` at 25 °C for bovine serum amine oxidase and for monoamine oxidase B
at pH 7.5, and 1.18 at pH 6.1; the monoamine oxidase series runs from 1.02 to
1.21 over temperature). Those secondary effects are held in
`data/secondary_reference.csv`, and `analysis/reference_asymmetry.py` computes every r
from them.

Because these records are not produced by `build_trinomial.py`, re-running that
script rebuilds only the parsed series, and always writes to
`trinomial_benchmark_partial.csv` rather than overwriting the shipped file, so the
hand-transcribed rows cannot be lost. It then says which case holds: if parseable
series are missing it warns that primary sources are absent; otherwise it compares
the rebuilt records with the shipped parseable ones and reports whether they
reproduce exactly. (Before 2026-09-14 it counted the hand-transcribed series in
its baseline, so every complete run wrongly reported missing sources.)

## How each source made its two tritium references

`data/tracer_design.csv` records, for all 30 enzyme forms, which tritiated
molecule served as the reference in the H/T and in the D/T experiment. It
decides the reference asymmetry `r = c_D/c_H` of Proposition S4, and
`analysis/reference_asymmetry.py` reads it. Each entry was read from the
source's Methods (Cha 1989 from its note 27).

**Shared tracer, 17 forms: r = 1 exactly.** One tritiated isotopologue serves
both experiments and only the bulk substrate changes, so the two references are
the same molecule.

| Source | Forms | Tritium tracer (both experiments) | Bulk in H/T / in D/T |
|---|---|---|---|
| Pagano et al. 2019 | FDH V123A, V123G | [3H]formic acid | formic acid / 99.8% deuterated formic acid |
| Islam et al. 2018 | hsTSase WT hydride | (R)-6-[3H]-MTHF | (R)-6-[1H]-MTHF / (R)-6-[2H]-MTHF |
| Islam et al. 2018 | hsTSase WT proton | [5-3H]-dUMP | [2-14C]-dUMP / [5-2H,2-14C]-dUMP |
| Abeysinghe & Kohen 2015 | ecTSase Y209W proton | [5-3H] dUMP | [2-14C] dUMP / [2-14C,5-2H] dUMP |
| Agrawal et al. 2004 | ecTSase WT hydride | (R)-[6-3H]CH2H4folate | protiated / (R)-[6-2H]CH2H4folate |
| Singh et al. 2014 | ecDHFR F125M, G121V-F125M, M42W-F125M, W133F | (4R)-[4-3H]-NADPH | [Ad-14C]-NADPH / (4R)-[Ad-14C,4-2H]-NADPH |
| Wang et al. 2014 | ecDHFR light and heavy enzyme | [4R-T]-NADPH | [Ad-14C]-NADPH / [4R-D]-NADPH |
| Li et al. 2021 | five hsDHFR forms | 4R-3H-NADPH | [carbonyl-14C]-NADPH / 4R-[carbonyl-14C,4-2H]-NADPH |

Formate, C6 of the folate cofactor and C5 of dUMP each carry a single hydrogen.
NADPH carries two at C4; the transferred one is 4R, and the 4S hydrogen is
protium in both references.

**Non-transferred position labeled, 13 forms: r inferred.** The D/T reference
carries deuterium at the non-transferred position, so it differs from the H/T
reference there. This is the same design that makes the secondary effects
measurable, and these are exactly the 13 forms that report them.

| Source | Forms | H/T reference | D/T reference |
|---|---|---|---|
| Cha et al. 1989 | yeast ADH | benzyl alcohol tritiated at C1 | benzyl alcohol "containing D and T at C-1" (note 27) |
| Bahnson et al. 1993; 1997 (same protocol, "doubly labeled substrates as described") | nine horse liver ADH forms | [1-3H]benzyl alcohol, randomly tritiated | [1,1-2H2,1-3H]benzyl alcohol |
| Tsai & Klinman 2001 | horse liver WT and F93W at 3 C | [1-3H]benzyl alcohol | [1-2H,1-3H]benzyl alcohol |
| Grant & Klinman 1989 | bovine serum amine oxidase | randomly tritiated benzylamine from [1-1H]benzaldehyde | randomly tritiated benzylamine from [1-2H]benzaldehyde |
| Jonsson et al. 1994 | monoamine oxidase B, both pH | [1,1-1H2-1-3H]-p-methoxybenzylamine | [1,1-2H2-1-3H]-p-methoxybenzylamine |

Until 2026-09-14 the 17 shared-tracer forms were recorded as "assumed" and the
text said the references always differ at the non-transferred position; reading
the Methods showed the condition holds for them with equality.

## Auxiliary kinetics used to complete the identified set

These sources contribute no isotope-effect records. They supply the one extra
quantity that turns a half-line into a point, and are read by
`analysis/completion.py`:

| Source | Quantity taken | Used for |
|---|---|---|
| Klinman, *Biochemistry* **15**, 2018 (1976), [10.1021/bi00654a032](https://doi.org/10.1021/bi00654a032) | `k_-1/k_cat = 1.3–7.3` for benzyl alcohol at 25 °C, pH 8.5 (Table IV; Table III gives the underlying Michaelis-constant effects) | completes yeast ADH, whose isotope effects are Cha et al. (1989) |
| Grant & Klinman, *Biochemistry* **28**, 6597 (1989) | Table IV, pre-steady-state beside steady-state isotope effects, 15–45 °C | masking factor 0.948 ± 0.044, completes bovine serum amine oxidase |
| Sekhar & Plapp, *Biochemistry* **29**, 4289 (1990), [10.1021/bi00470a005](https://doi.org/10.1021/bi00470a005) | `K_m` = 0.03 mM (benzyl alcohol) against 0.08 mM (α,α-`d2`) | shows the alcohol Michaelis effect is inverse in horse liver ADH, which is why that system admits no completion |

These are rate ratios and Michaelis constants read from the printed tables, not
isotope-effect records, so they are held as literals in `completion.py` at the
point of use rather than in `data/`.  Running that script reproduces both
completion tables of the supplement.

Two source classes were examined and rejected, for reasons worth recording:

- **Deuterium-only.** Soybean lipoxygenase (Rickert & Klinman 1999; Knapp et al.
  2002) and the ADH flexibility work (Kohen & Klinman 2000) report H/D only. No
  second isotope pair exists, so the offset cannot be formed.
- **No matched D/T.** Peptidylglycine α-amidating enzyme (Francisco et al.
  1998) and dopamine β-monooxygenase (Miller & Klinman 1985) report H/T and H/D
  effects on `V/K`, not the matched H/T and D/T pair against a tritium reference
  that the offset requires, and derive intrinsic effects from them by Northrop's
  method. (Recorded until 2026-09-14 as publishing only intrinsic values, which
  misstated what those papers contain.)
- **Northrop-derived intrinsic values only.** The ecDHFR I14 series (Stojković et
  al. 2012) publishes intrinsic effects obtained by *imposing* the Swain–Schaad
  relation, and no observed pair. Those values are pinned to `F = 0` by
  construction, so using them would be circular.
- **Isotope effects in figures only.** Glucose oxidase (Kohen, Jonsson &
  Klinman 1997) and thermophilic ADH (Kohen et al. 1999) state no primary
  isotope effect as a number.
