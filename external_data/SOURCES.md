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
| Islam et al., *PLoS ONE* **13**, e0196506 (2018) | `10.1371/journal.pone.0196506` | hsTSase hydride + proton | PMC5929524, OA (S2/S3 Tables, files s004/s005 .xlsx) |
| Abeysinghe & Kohen, *Int. J. Mol. Sci.* **16**, 7304 (2015) | `10.3390/ijms16047304` | ecTSase Y209W proton | PMC4425018, OA (Table S2) |
| Singh et al., *JACS* **136**, 2575 (2014) | `10.1021/ja411998h` | 4 ecDHFR variants, 5 temperatures | PMC3985941, OA |
| Li et al., *Biochemistry* **60**, 3822 (2021) | `10.1021/acs.biochem.1c00558` | 5 hsDHFR variants | PMC8697555, OA |
| Wang et al., *JACS* **136**, 8333 (2014) | `10.1021/ja501936d` | heavy ecDHFR | PMC4063187, OA |
| Sikorski et al., *JACS* **126**, 4778 (2004) | `10.1021/ja031683w` | light (wild-type) ecDHFR | **paywalled** |
| Pagano et al., *ACS Catal.* **9**, 11199 (2019) | `10.1021/acscatal.9b03345` | FDH V123A, V123G | **paywalled** |

## Series admitted after the census of 2026-10-20

The census queries, their date and the stopping rule are in
`census_queries.md` beside this file.

A census of the two laboratories' publications (SI, "Source census") found
these further series. Every one is transcribed from the table named, and every
value of the PDF supplements is checked, with its error, against the PDF text
by the transcription audit (`TRANSCRIBED` in `analysis/build_trinomial.py`);
the two Word supplements and the scanned Rucker table were checked against the
document's table and the page image (`HAND_TRANSCRIBED`). Row and temperature
assignment was checked for every record: each H/T and D/T pair sits on one line
of the source beside its temperature, or was read from the rendered page.
These supplements are paywalled or distributed as Word files and are placed by
hand in `external_data/si/manual/` under the publisher's file names.

| Source | DOI | Contributes | Table | Records |
|---|---|---|---|---|
| Wang & Kohen, *JACS* **132**, 9820 (2010) | `10.1021/ja103010b` | ecTSase WT proton, 50 mM Mg2+ | SI Table S1 | 4 |
| Wang et al., *JACS* **134**, 17722 (2012) | `10.1021/ja307859m` | ecTSase Y209W hydride (no D/T at 15 °C) | SI Table S4 | 3 |
| Wang et al., *JACS* **135**, 7583 (2013) | `10.1021/ja400761x` | ecTSase WT hydride, 50 mM Mg2+ | SI Table S2 | 4 |
| Hong et al., *JACS* **128**, 5636 (2006) | `10.1021/ja060196o` | ecTSase W80M hydride (its WT 20 °C column reprints Agrawal 2004) | SI Table S1 | 4 |
| Islam et al., *JACS* **136**, 9850 (2014) | `10.1021/ja504341g` | ecTSase R166K hydride | SI Table S2 | 4 |
| Islam et al., *ACS Catal.* **5**, 6061 (2015) | `10.1021/acscatal.5b01332` | ecTSase R166K proton (its Table S2 reprints Islam 2014) | SI Table S1 | 4 |
| Gurevic et al., *ACS Catal.* **8**, 10241 (2018) | `10.1021/acscatal.8b02554` | ecTSase N177D, N177S hydride (its Table S1 reprints Wang 2013) | SI Tables S2, S3 | 8 |
| Wang et al., *Biochemistry* **45**, 1383 (2006) | `10.1021/bi0518242` | ecDHFR G121V, pH 8.0 (Word supplement) | SI Table 1S | 5 |
| Wang et al., *Phil. Trans. R. Soc. B* **361**, 1307 (2006) | `10.1098/rstb.2006.1871` | ecDHFR M42W, pH 8.0 (Word supplement) | ESM Table 1 | 5 |
| Stojković et al., *Chem. Commun.* **46**, 8974 (2010) | `10.1039/c0cc02988b` | ecDHFR I14A, pH 8.0 | ESI Table 1S | 5 |
| Singh et al., *ACS Catal.* **5**, 3067 (2015) | `10.1021/acscatal.5b00331` | ecDHFR I14A-G121V | SI Table S1 | 5 |
| Liu et al., *PNAS* **111**, 18231 (2014) | `10.1073/pnas.1415940111` | ecDHFR Y100F, D27S, D27S-Y100F (pH 9.0) and WT at pH 7.0 | SI Tables S5–S8 | 20 |
| Li et al., *Biochemistry* **58**, 3861 (2019) | `10.1021/acs.biochem.9b00460` | hsDHFR M53W, S145V, M53W-S145V (5–35 °C) | SI Table S2 | 12 |
| Bandaria et al., *JACS* **131**, 10151 (2009) | `10.1021/ja902120t` | FDH WT, commercial enzyme (a mixture of isozymes) | SI Table S2 | 5 |
| Guo et al., *Biochemistry* **55**, 2760 (2016) | `10.1021/acs.biochem.6b00181` | FDH WT, recombinant | SI Table S1 | 5 |
| Ranasinghe et al., *JACS* **139**, 17405 (2017) | `10.1021/jacs.7b08359` | FDH half-heavy and heavy enzymes (its light-enzyme Table S1 reprints Guo 2016) | SI Tables S2, S3 | 14 |
| Rucker et al., *Biochemistry* **31**, 11489 (1992) | `10.1021/bi00161a030` | yeast ADH with p-chlorobenzyl alcohol, 0 and 25 °C (scanned) | Table IV | 2 |

Single-condition records admitted at the same time are listed under "Records
transcribed by hand".

## Sources examined and excluded

| Source | DOI | Why excluded |
|---|---|---|
| Stojković et al., *JACS* **134**, 1738 (2012) | `10.1021/ja209425w` | Publishes only Northrop-derived *intrinsic* effects; the observed H/T and D/T values appear nowhere. Those intrinsic pairs are pinned to `F = 0` by construction and carry no information about the discriminating invariant. |
| Austin-Kloppe et al., *Chem. Sci.* **17**, 12440 (2026) | `10.1039/D6SC01847E` | 38 solution-phase hydride transfers (34 from NADH/NAD+ models, analyzed quantitatively), but deuterium only. No second isotope pair, so the offset cannot be formed. |
| Kohen et al., *Nature* **399**, 496 (1999) | `10.1038/20981` | ht-ADH. Its well-known anomalous Swain–Schaad result is the **secondary** exponent, outside the model class treated here. Primary data are published only as Arrhenius prefactor ratios (F_A = -0.38 ± 0.52), which is not the offset at any temperature because the activation-energy differences are not reported. |

## Manual downloads

Three documents cannot be fetched programmatically. Place them in
`external_data/si/manual/` under the filenames shown, then rerun
`python analysis/build_trinomial.py`.

1. **Pagano et al. 2019**: <https://doi.org/10.1021/acscatal.9b03345>
   Supporting Information PDF → `cs9b03345_si_001.pdf`
   (PMC8118594 exists but is not in the open-access subset.)

2. **Sikorski et al. 2004**: <https://doi.org/10.1021/ja031683w>
   Supporting Information PDF → `ja031683w_si_002.pdf`
   (its Table 1 prints the plus-minus sign as a dash; the parse reads it so)

3. *(Optional, not used in the benchmark)* Stojković et al. 2012 and
   Austin-Kloppe et al. 2026, if you wish to reproduce the exclusion analysis
   in Supplemental Material §S4.2.

Without the Pagano SI the parse produces 14 series and 63 records instead of 16
and 73; without the Sikorski SI, 15 and 67; without both, 13 and 57. The shipped
benchmark is never overwritten by a partial parse.

## Transcription discipline

`analysis/build_trinomial.py` parses the supplementary tables from the PDF text
layer rather than retyping them, and runs an audit confirming that every value in
the shipped benchmark from a PDF source appears verbatim in that source together
with its stated uncertainty. It reads the shipped file, not a fresh parse, so a
hand edit to the CSV is caught. With all sources present it covers 260 values
across six documents and passes with zero misses, and the script exits non-zero
on any miss. It does not establish row or column assignment: a swap of the H/T
and D/T columns within a row would pass. Those assignments were checked against
the rendered tables.

The ecDHFR light-enzyme series is parsed from its original report, Sikorski et
al. 2004 (SI Table 1), which lists that series alone. Wang et al. 2014 reprint it
beside their heavy enzyme, rounded (3.1 for 3.06 at 5 °C; 4.8 ± 0.1 for 4.82 ± 0.11
and 4.80 ± 0.12 at 45 and 35 °C) and distinguished only by colour. Until
2026-10-15 the benchmark took the light series from that reprint, with the 10 °C
row assigned by interpolation; only Wang's heavy rows are now parsed.

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
| Cha, Murray & Klinman, *Science* **243**, 1325 (1989), [10.1126/science.2646716](https://doi.org/10.1126/science.2646716) | `data/cha1989_yadh.csv`, `data/cha1989_yadh_dt_only.csv` | Table 2: three matched H/T and D/T determinations and their average (the average is the record), plus a fourth D/T-only determination, 1.72 ± 0.05, with no H/T partner, kept separately because it is not a matched pair | 1 |
| Bahnson et al., *Biochemistry* **32**, 5503 (1993), [10.1021/bi00072a003](https://doi.org/10.1021/bi00072a003) | `data/ladh_adh_primary.csv` | Table I | 5 |
| Bahnson et al., *Proc. Natl. Acad. Sci. USA* **94**, 12797 (1997), [10.1073/pnas.94.24.12797](https://doi.org/10.1073/pnas.94.24.12797) | `data/ladh_adh_primary.csv` | Table 1 | 4 |
| Tsai & Klinman, *Biochemistry* **40**, 2303 (2001), [10.1021/bi002075l](https://doi.org/10.1021/bi002075l) | `data/ladh_adh_primary.csv` | Tables 2 and 4, aqueous rows at 3 °C | 2 |
| Tsai & Klinman, *Biochemistry* **40**, 2303 (2001), [10.1021/bi002075l](https://doi.org/10.1021/bi002075l) | `data/trinomial_benchmark.csv` (family `LADH`, variants in 40%, 50% and 60% MeOH) | Table 2 (wild type, 50% MeOH, 3 to −50 °C) and Table 4 (F93W, 40% MeOH, 3 to −30 °C; 60% MeOH, 3 and −20 °C), all at 10 mM NAD+ | 11 |
| Grant & Klinman, *Biochemistry* **28**, 6597 (1989), [10.1021/bi00442a010](https://doi.org/10.1021/bi00442a010) | `data/bsao_grant1989.csv` | Table I averages, 25 °C | 1 |
| Jonsson, Edmondson & Klinman, *Biochemistry* **33**, 14871 (1994), [10.1021/bi00253a026](https://doi.org/10.1021/bi00253a026) | `data/trinomial_benchmark.csv` (family `MAOB`) | Tables 1 and 3 | 10 |
| Agrawal, Hong, Mihai & Kohen, *Biochemistry* **43**, 1998 (2004), [10.1021/bi036124g](https://doi.org/10.1021/bi036124g) | `data/agrawal2004_ectsase.csv` | Results text, 20 °C: 6.91 ± 0.05 and 1.78 ± 0.02 (the other temperatures are plotted only and are not digitized) | 1 |
| Rucker, Cha, Jonsson, Grant & Klinman, *Biochemistry* **31**, 11489 (1992), [10.1021/bi00161a030](https://doi.org/10.1021/bi00161a030) | `data/trinomial_benchmark.csv` (family `YADH`), `data/yadh_rucker1992.csv`, `data/bsao_rucker1992.csv` | Table IV (yeast ADH, p-chlorobenzyl alcohol, 0 and 25 °C); Table V (yeast ADH, p-methoxybenzyl alcohol, 25 °C; its H/T is an extrapolation to zero enzyme concentration, Fig. 3); Table I and text (bovine serum amine oxidase, p-methoxybenzylamine, 25 °C: H/T 20.2 ± 1.9; D/T 2.86 ± 0.05 stated in the text, Table II giving the two determinations). D/T corrected by the authors for protium | 4 |
| Chin & Klinman, *Biochemistry* **39**, 1278 (2000), [10.1021/bi9920331](https://doi.org/10.1021/bi9920331) | `data/ladh_chin2000.csv` | Table 1: F93W and the double mutants F93W:I224A, I224L, I224G, I224V at 25 °C, pH 7 (standard errors of the mean, 14 to 41 determinations) | 5 |
| Hong, Maley & Kohen, *Biochemistry* **46**, 14188 (2007), [10.1021/bi701363s](https://doi.org/10.1021/bi701363s) | `data/ectsase_hong2007.csv` | Table 1: ecTSase wild type and Y94F, proton abstraction, 25 °C, at 5 µM CH2H4folate, no Mg2+ | 2 |
| Islam et al., *PLoS ONE* **13**, e0196506 (2018), [10.1371/journal.pone.0196506](https://doi.org/10.1371/journal.pone.0196506) | `data/hstsase_islam2018_nomg.csv` | Table 2: hsTSase hydride without Mg2+ (temperature not printed; 25 °C inferred from the with-Mg2+ cell of the same table). Its two ecTSase pairs are quotations: 4.14/1.65 is Wang et al. 2013 SI Table S2 at 25 °C; 7.07/1.81 is printed in neither cited source and is not admitted | 1 |
| Sen, Yahashiri & Kohen, *Biochemistry* **50**, 6462 (2011), [10.1021/bi2003873](https://doi.org/10.1021/bi2003873) | `data/dhfr_sen.csv` | Table 1: ecDHFR with dihydrofolate, ecDHFR with dihydrobiopterin, and R67 (type II) DHFR (standard deviations; the temperature of these runs is not printed) | 3 |
| Sen, Stojković & Kohen, *Anal. Biochem.* **430**, 123 (2012), [10.1016/j.ab.2012.08.012](https://doi.org/10.1016/j.ab.2012.08.012) | `data/dhfr_sen.csv` | Table 2, [carbonyl-14C] row: ecDHFR wild type, 25 °C, pH 9.0 (its [Ad-14C] row reprints Sikorski 2004) | 1 |

The Tsai & Klinman methanol records were first excluded on the ground that a
random component in the mechanism breaks the single shared observation map.
That was wrong on both counts: the competitive runs were made at 10 mM NAD+,
for F93W far above its 0.75 mM Michaelis constant, which the authors state
ensures an ordered pathway (no Michaelis constant is reported for the wild type
in 50% methanol), and a random mechanism with no route around the chemical step
still gives the base map at fixed cosubstrate. They were admitted on
2026-10-13 under the same rule as every other source. `data/ladh_bahnson.csv`
is the 1997 re-tabulation of the Bahnson Table 1 values (L57F and F93W rounded
as reprinted there); no script reads it, and the records it holds enter the
corpus once, from `data/ladh_adh_primary.csv`.

The Agrawal record was first excluded on the ground that all of that paper's
isotope effects are in figures. That overlooked the 20 °C values stated in the
text, and the record was admitted on 2026-09-14 under the rule applied to every
source: an observed matched pair is admitted when both effects are stated as
numbers with their uncertainties.

Both amine oxidase sources report the primary H/T and D/T effects together with
the secondary effects on the same references; the ratio `R` of those observed
effects is a lower bound on the ratio `r` of the two commitments used in the
unequal-reference check of `analysis/curvature.py` (`R = 1.14` at 25 °C for
bovine serum amine oxidase and for monoamine oxidase B at pH 7.5, and 1.18 at
pH 6.1; the monoamine oxidase series runs from 1.02 to 1.21 over temperature). Those secondary effects are held in
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

`data/tracer_design.csv` records, for all 33 enzyme forms, which tritiated
molecule served as the reference in the H/T and in the D/T experiment. It
decides the reference asymmetry `r = c_D/c_H` of Proposition S4, and
`analysis/reference_asymmetry.py` reads it. Each entry was read from the
source's Methods (Cha 1989 from its note 27).

**Shared tracer, 16 forms: r = 1 exactly.** One tritiated isotopologue serves
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
| Wang et al. 2014 | ecDHFR heavy enzyme | [4R-T]-NADPH | [Ad-14C]-NADPH / [4R-D]-NADPH |
| Li et al. 2021 | five hsDHFR forms | 4R-3H-NADPH | [carbonyl-14C]-NADPH / 4R-[carbonyl-14C,4-2H]-NADPH |

Formate, C6 of the folate cofactor and C5 of dUMP each carry a single hydrogen.
NADPH carries two at C4; the transferred one is 4R, and in these protocols the
4S hydrogen is protium in both references.

**Non-transferred position labeled, 17 forms: r inferred.** The D/T reference
carries deuterium at the non-transferred position, so it differs from the H/T
reference there. This is the same design that makes the secondary effects
measurable, and these are exactly the 17 forms that report them.

| Source | Forms | H/T reference | D/T reference |
|---|---|---|---|
| Cha et al. 1989 | yeast ADH | benzyl alcohol tritiated at C1 | benzyl alcohol "containing D and T at C-1" (note 27) |
| Bahnson et al. 1993; 1997 (same protocol, "doubly labeled substrates as described") | nine horse liver ADH forms | [1-3H]benzyl alcohol, randomly tritiated | [1,1-2H2,1-3H]benzyl alcohol |
| Tsai & Klinman 2001 | horse liver WT and F93W at 3 C, and the three methanol series | [1-3H]benzyl alcohol | [1-2H,1-3H]benzyl alcohol |
| Grant & Klinman 1989 | bovine serum amine oxidase | randomly tritiated benzylamine from [1-1H]benzaldehyde | randomly tritiated benzylamine from [1-2H]benzaldehyde |
| Jonsson et al. 1994 | monoamine oxidase B, both pH | [1,1-1H2-1-3H]-p-methoxybenzylamine | [1,1-2H2-1-3H]-p-methoxybenzylamine |
| Sikorski et al. 2004 | ecDHFR light (wild-type) enzyme | (4R)-[4-3H]-NADPH vs [Ad-14C]-NADPH | (4R)-[4-3H]-NADPH with D at 4S ("4R [4,4-2H,3H]") vs [Ad-14C,4-2H2]-NADPH |

Until 2026-09-14 the 17 shared-tracer forms were recorded as "assumed" and the
text said the references always differ at the non-transferred position; reading
the Methods showed the condition holds for them with equality. Until 2026-10-15
the ecDHFR light enzyme was among them, on the strength of Wang et al.'s Methods,
which describe their own heavy-enzyme experiments; Sikorski et al.'s original
protocol puts deuterium at 4S in the D/T reference, and their secondary effects
(1.149 ± 0.005 H/T, 1.058 ± 0.003 D/T at 25 °C, in `data/secondary_reference.csv`)
give R = 1.086.

## Auxiliary kinetics used to complete the identified set

These sources contribute no isotope-effect records. They supply the one extra
quantity that turns a half-line into a point, and are read by
`analysis/completion.py`:

| Source | Quantity taken | Used for |
|---|---|---|
| Klinman, *Biochemistry* **15**, 2018 (1976), [10.1021/bi00654a032](https://doi.org/10.1021/bi00654a032) | `k_-1/k_cat = 1.3–7.3` at 25 °C, pH 8.5: Table IV's conversion of Michaelis-constant isotope effects of 1.1–1.5 (the Table IV entries, spanning the normal effects of the para-substituted alcohols whose `k_cat` Klinman attributes to hydrogen transfer), assuming an isotope effect of 4 on `k_cat` (her measured value for benzyl alcohol is 3.4) and none on binding; for benzyl alcohol itself Table III gives 1.3 on the NAD+ constant `K_A` (maps to 2.3) and an inverse 0.80 on the alcohol constant `K_B` | completes yeast ADH conditionally on that conversion; isotope effects are Cha et al. (1989) |
| Grant & Klinman, *Biochemistry* **28**, 6597 (1989) | Table IV, pre-steady-state beside steady-state isotope effects, 0–45 °C (six temperatures; the steady-state column is derived through the Swain–Schaad relation) | masking factor at 25 °C, 0.839 ± 0.097; the six temperatures are inhomogeneous (Q = 17.9, p = 0.003), so the completion is **withdrawn**: F_int = −0.038, 95% [−0.267, +0.242] |
| Sekhar & Plapp, *Biochemistry* **29**, 4289 (1990), [10.1021/bi00470a005](https://doi.org/10.1021/bi00470a005) | apparent `K_m` for the transient oxidation at pH 8.0, 0.03 mM (benzyl alcohol) against 0.08 mM (α,α-`d2`) | shows the alcohol Michaelis effect is inverse in horse liver ADH, which is why that system admits no completion |

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
