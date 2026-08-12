# One-sided identifiability in isotope tests of enzymatic hydrogen tunneling

Data, code and verification suite for the analysis of what competitive
multiple-isotope kinetics can establish about hydrogen tunneling in enzymes.

## The question and the answer

Competitive H/T and D/T measurements are the main evidence bearing on enzymatic
hydrogen tunneling. It is established that the Swain-Schaad exponent alone does
not diagnose it (Kohen & Jensen 2002; Hirschi & Singleton 2005; Shelton, Hrovat
& Borden 2007). This work asks what such an experiment does determine.

**The exponent is the wrong coordinate.** `X_H/X_D` is a ratio and is unchanged
by scaling, so it cannot see a displacement transverse to the mass-scaling line.
The quantity that can is the offset `F = ln K_HT - gamma_SC ln K_DT`, which is
zero on that line by construction.

**F is not point identified.** Inverting the observation map of the minimal
competitive scheme over all admissible commitments leaves a set. That set is a
half-line with a closed-form endpoint:

```
dF/dc = [ c (gamma b - a) - a b (gamma - 1) ] / [ c (c - a) (c - b) ],
        a = K_HT - 1,  b = K_DT - 1

L_H = a/b >= gamma :  open half-line (F_obs, inf), endpoint only as c -> inf
L_H = a/b <  gamma :  closed half-line [F_min, inf), interior minimum at
                      c* = a b (gamma - 1) / (gamma b - a)
```

Verified against direct search on 600 random admissible pairs to 9e-8.

**The direction the set is open is a property of the kinetic scheme, not of the
data.** Under a single shared additive commitment `C_f` instead of a single
shared `c`, both intrinsic effects diverge together and the set is open *below*
rather than above. The two families agree pair by pair but not jointly.

**Mechanisms enter through the offsets they can attain.** For a gated vibronic
model the attainable offsets have a *scale-dependent* envelope:

```
B(K) = ln K - (gamma/2) ln[ K^2 (t-1) / (t - d + K^2 (d-1)) ]   K < sqrt(t)
     = F0 = -0.042086                                            K >= sqrt(t)
```

with threshold `sqrt(mu_T/mu_H) = 1.2689`. There is no flat bound: `dF/dA < 0`
identically, so the unconstrained supremum is 0, approached as the isotope
effects collapse to unity. Verified against direct maximization to 7e-7.

**Neither the envelope nor the ordering is universal.** Summing over excited
vibronic channels can carry the offset above zero; a worked parameter set gives
`F = +0.00530`, stable from n=6 to n=30. So the gated family is not confined
below the semiclassical locus once excited channels contribute.

**Empirically, nothing is excluded.** Across 73 matched records in 16 series,
three enzyme families, five enzyme-organism systems and two chemical steps,
one-sided 95% confidence bounds exclude neither mechanism at any assumed
correlation between the H/T and D/T effects. The closest series on the point
estimate falls short by 0.0175, which is 0.61 standard errors.

The framework, rather than any verdict on tunneling, is the contribution:
mechanistic claims from these experiments are comparisons between an
experimentally identified set and a computed mechanism envelope.

## Layout

```
analysis/
  masses.py             isotope mass convention, imported everywhere
  qtunnel.py            exact gated-overlap model
  ridge.py              ridge limit and the offset
  partial_id.py         envelope B(K), exact identified-set endpoint,
                        two-commitment exploration
  identifiable_set.py   the half-line result and the point-estimate bounds
  bounds_uncertainty.py one-sided confidence bounds with shared-reference
                        correlation bracketed
  build_trinomial.py    parses primary SI tables -> benchmark + audit
  offset_analysis.py    profiles the offset per series (continuation method)
  offset_summary.py     goodness-of-fit stratification of those profiles
  vibronic.py           does the bound survive the full vibronic sum?
  verify_derivation.py  numerical + symbolic checks of every closed form
  audit_v3.py           adversarial audit of the reconstructed claims
  audit_referee.py      independent check of a first review's claims
  audit_r2.py           independent check of a second review's claims
  audit_tableII.py      profiling test and a non-monotonicity counterexample
  discriminate.py       model fits (superseded analysis, retained for the record)
  design.py             power calculation (superseded, retained)
  export_figs_v3.py     pgfplots tables for the figures
data/                   curated inputs, never written by any script
results/                generated reports and tables
figures/tikz/           figure sources and their data
external_data/          SOURCES.md and fetch_sources.sh (no PDFs redistributed)
```

## Reproduce

```bash
pip install -r requirements.txt

./external_data/fetch_sources.sh     # open-access supplements; see SOURCES.md
cd analysis

python masses.py             # the mass convention and its consequences
python verify_derivation.py  # must exit 0
python build_trinomial.py    # benchmark + transcription audit (244 values)
python ridge.py              # ridge limit, offset, finite-scale bound
python partial_id.py         # envelope + exact endpoint, both verified  (~5 min)
python identifiable_set.py   # the half-line result, point bounds       (~1 min)
python bounds_uncertainty.py # one-sided confidence bounds             (~25 min)
python audit_r2.py           # reconfirms the second review's claims   (~10 min)
python audit_v3.py           # adversarial audit, must exit 0          (~10 min)
python vibronic.py           # the vibronic-sum penalty                 (~5 min)
python offset_analysis.py    # per-series profiles                     (~40 min)
python offset_summary.py     # stratified interpretation
python export_figs_v3.py     # pgfplots tables
cd ../figures/tikz && ./build.sh
```

Every stochastic step is seeded, so a rerun reproduces the reported values
exactly. `data/` is inputs only: no script writes to a file it reads.

One caveat on timing: `offset_analysis.py` computes profile likelihoods by
continuation over a 301-point grid for each of 16 series and takes roughly
40 minutes on one core.

## Key numbers

Reduced masses of the X-H oscillator, CODATA/AME2020, carbon donor. This is the
convention of Kohen & Jensen, who obtain 3.34 from the reduced mass of 12C and
the hydrogen isotopes; the bare-mass value 3.2628 is also in circulation.

```
gamma_SC = 3.34887   rigid = 2.45884   F0 = -0.042086   threshold = 1.26890
```

Benchmark: 73 matched records, 16 series, 3 enzyme families, 5 enzyme-organism
systems, 3 organisms, 5-45 C. Transcription audit: 244 values, 0 misses.

Bounds: 0 of 16 series exclude either mechanism, at correlation 0, 0.5 or 0.9.
Best 95% lower bound -0.0973 (rho = 0) to -0.0935 (rho = 0.9), ecDHFR light
enzyme in both cases. Best point estimate -0.0596 (ecDHFR W133F), short of F0 by
0.0175, which is 0.61 standard errors of the simultaneous statistic and 0.43 of
the pointwise delta-method standard error at the single closest temperature.

## On the superseded analysis

`results/` retains outputs from an earlier version of this work that was
withdrawn and rebuilt, together with the audits that overturned it
(`audit_referee.txt`, `audit_tableII.txt`). Three errors were found and are
documented rather than quietly removed:

1. An optimization performed on a constraint *boundary* rather than over the
   constraint *region*, which reversed the direction of a proposed experimental
   criterion. A "cap" that decreased as the feasible set grew was an internal
   contradiction that should have been caught.
2. A supremum over a model family read as a boundary in observable space, when
   it is an asymptotic slope reached only as the isotope effects diverge.
3. A goodness-of-fit table produced by pinning the intrinsic isotope effect at a
   value itself derived from the same data through the relation under test.
   Profiling instead moved a chi-square of 845 to 6.17.

A fourth was caught internally: a monotonicity assertion that failed on 21 of 73
records, which led to the `L_H` criterion now in the paper.

A second external review then found three more, all confirmed here
(`results/audit_review2.txt`):

4. A claimed flat bound `F <= F0` over the whole parameter space. It is false:
   `dF/dA < 0` identically, so the supremum is 0, and at `A = 0.05, w -> 0` one
   finds `F = -0.0110 > F0`. The numerical check had begun at
   `K_HT = 1.5`, entirely inside the region where the flat bound does hold. It is
   replaced by the scale-dependent envelope `B(K)`.
5. An inference from `F0 + Delta > 0` that the two model classes overlap. That
   does not follow from an inequality on an upper bound. Direct evaluation is
   required, and it shows the summed model can reach `F = +0.0053`.
6. An experimental recommendation with the commitment ratio reversed. Masking
   falls as `k_off/k_chem` *increases*, not decreases.

A seventh point, that point-estimate bounds were reported without sampling
uncertainty, is addressed by `bounds_uncertainty.py`.

They are kept because the verification suite that found them is the reason to
trust what replaced them.

## Licensing

Code MIT (`LICENSE`); curated data CC BY 4.0 with attribution requirements to
the primary measurement papers (`LICENSE-DATA.md`). No publisher PDF is
redistributed; see `external_data/SOURCES.md`.

## Citation

See `CITATION.cff`. Please also cite the primary measurement papers listed in
`external_data/SOURCES.md` for any use of the benchmark values.
