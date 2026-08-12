# One-sided identifiability in isotope tests of enzymatic hydrogen tunneling

Data, code and verification suite for the analysis of what competitive
multiple-isotope kinetics can establish about hydrogen tunneling in enzymes.

## The question and the answer

Competitive H/T and D/T measurements are the evidential basis for enzymatic
hydrogen tunneling. This work asks what they determine, and finds three things.

**Gating inverts the isotope mass power.** For vibronically nonadiabatic
transfer with a thermally sampled donor–acceptor distance, the ground-channel
rate ratios can be evaluated exactly. The overlap exponent scales as `mu^(1/2)`
in a rigid barrier and as `mu^(-1/2)` under strong gating, and `mu^(-1/2)` is
precisely the zero-point-energy mass factor. The semiclassical and
gated-tunneling loci are therefore *parallel* in log-isotope-effect space, and
the Swain–Schaad exponent, being an angle from the origin, cannot see the fixed
offset between them.

**The offset is the right invariant, and it separates the mechanisms.** With
`F = ln K_HT - gamma_SC ln K_DT`, semiclassical kinetics has `F = 0` while the
gated-tunneling family satisfies `F <= F0 = -0.042086`, independent of the
isotope-effect magnitude. The two are disjoint, separated by a gap fixed by the
isotope masses alone.

**But the experiment determines `F` only from below.** Inverting the observation
map over admissible commitments gives the half-line `[F_min, infinity)`, whose
lower end is set by the classical H-reference contrast `L_H = (K_HT-1)/(K_DT-1)`
relative to `gamma_SC`. Since the semiclassical value lies *above* the tunneling
bound, any measurement refuting semiclassical kinetics also refutes gated
tunneling, and no attainable measurement refutes semiclassical kinetics alone.
The evidence is structurally one-sided.

Across 73 matched records spanning four enzyme families, five systems and three
organisms, **no series excludes either mechanism**. The closest falls short by
0.018 in `F`. The limiting factor is not precision or temperature coverage but
the kinetic commitment, which displaces the observed offset by up to eighty
times the signal.

## Layout

```
analysis/
  masses.py             isotope mass convention, imported everywhere
  qtunnel.py            exact gated-overlap model
  ridge.py              ridge limit, the offset, the finite-scale bound
  identifiable_set.py   the half-line result and the model-free bounds
  build_trinomial.py    parses primary SI tables -> benchmark + audit
  offset_analysis.py    profiles the offset per series (continuation method)
  offset_summary.py     goodness-of-fit stratification of those profiles
  vibronic.py           does the bound survive the full vibronic sum?
  verify_derivation.py  numerical + symbolic checks of every closed form
  audit_v3.py           adversarial audit of the reconstructed claims
  audit_referee.py      independent check of an external review's claims
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
python identifiable_set.py   # the half-line result, model-free bounds  (~1 min)
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

Reduced masses of the X–H oscillator, CODATA/AME2020, carbon donor:

```
gamma_SC = 3.34887    rigid = 2.45884    gated = 2.84592    F0 = -0.042086
```

Benchmark: 73 matched records, 16 series, 4 families, 5 systems, 3 organisms,
5–45 °C. Transcription audit: 244 values, 0 misses.

Model-free bounds: 0 of 16 series exclude either mechanism; best is ecDHFR
W133F at `F >= -0.0596`, short of the tunneling bound by 0.0175.

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
records and led to the `L_H` criterion now in the paper.

They are kept because the verification suite that found them is the reason to
trust what replaced them.

## Licensing

Code MIT (`LICENSE`); curated data CC BY 4.0 with attribution requirements to
the primary measurement papers (`LICENSE-DATA.md`). No publisher PDF is
redistributed; see `external_data/SOURCES.md`.

## Citation

See `CITATION.cff`. Please also cite the primary measurement papers listed in
`external_data/SOURCES.md` for any use of the benchmark values.
