# One-sided identification of quantum signatures in enzymatic hydrogen transfer

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.21913975.svg)](https://doi.org/10.5281/zenodo.21913975)


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

**The evidence is one-sided, and that is the whole result.** For any admissible
observation `K_HT > K_DT > 1`,

```
L_H = (K_HT-1)/(K_DT-1)  >  gamma_obs = ln K_HT / ln K_DT     ALWAYS
```

because `(x-1)/ln x` is strictly increasing. An exclusion is claimed exactly when
the reference exponent is below the observed one, and there `L_H > gamma`, so
Proposition 2 gives the OPEN half-line with `F_obs >= 0`: the entire identified
set lies above zero. Hence:

* an observed exponent ABOVE the semiclassical reference cannot be produced by
  kinetic complexity at any commitment -- the inference is robust with no
  assumption about complexity;
* an observed exponent BELOW it leaves a set containing the semiclassical locus
  and every mechanism above it -- uninformative at any precision, because the
  obstruction is the set's direction, not its width.

This makes exact the criterion of Saunders (1985), used by Cha et al. (1989) and
Bahnson et al. (1993). Of 83 matched primary records, 82 fall on the
uninformative side.

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

**The exponent is not even constant where the offset is.** Along the ridge
`X_L = P_L + s C_L` the Swain-Schaad exponent runs over the whole interval
[2.846, 3.349] as `s` varies, spanning essentially every value reported as
diagnostic, while `F` stays pinned at `F0` to ten decimal places because
`C_H - gamma_SC C_D` vanishes identically.

**Mechanisms enter through the offsets they can attain.** For a gated vibronic
model the attainable offsets have a *scale-dependent* envelope:

```
B(K) = ln K - (gamma/2) ln[ K^2 (t-1) / (t - d + K^2 (d-1)) ]   K < sqrt(t)
     = F0 = -0.042086                                            K >= sqrt(t)
```

with threshold `sqrt(mu_T/mu_H) = 1.2689`. There is no flat bound: `dF/dA < 0`
identically, so the unconstrained supremum is 0, approached as the isotope
effects collapse to unity. Proved in the supplement; the proof turns on
`p_L' = q_L/2`, which cancels the intercept and slope contributions to `X_D`
along the constraint and leaves `dX_D/dw < 0` fixed by `t > d > 1` alone.
Numerical agreement is 7e-7, which is the double-precision noise floor, not the
accuracy: at 50 digits the supremum approaches `F0` strictly from below.

**Neither the envelope nor the ordering is universal.** Summing over excited
vibronic channels can carry the offset above zero; a worked parameter set gives
`F = +0.00530`, stable from n=6 to n=30. So the gated family is not confined
below the semiclassical locus once excited channels contribute.

**Reversibility opens a bounded vacuity window.** Deriving the map for a
reversible chemical step (it is Northrop's equation re-referenced to tritium,
which is what makes both commitments shared) gives
`F = F_obs + gamma ln D_D - ln D_H`. The set is unbounded below exactly when the
corner `Cf* = (v-u)/d`, `Cr* = (a-b)/d`, `d = av-bu` lies in the physical
quadrant, which needs BOTH `d > 0` and `v >= u`:

```
window = [E_D*, E_D**],   E_D* = (K_HT/K_DT)^(1/(gamma-1))
E_D** = root above E_D* of  a(K_DT/E - 1) = b(K_HT/E^gamma - 1)
window is nonempty  <=>  F_obs < 0
```

It is a bounded interval, not a half-line: `d -> -(a-b) < 0` at large `E_D`, so
the set recovers a finite endpoint above the window as well as below it. And a
POSITIVE observed offset can never be made vacuous by any equilibrium isotope
effect. Across the benchmark the lower edge runs 1.099-1.694 (median 1.465),
upper edge median 1.971, median width 0.50. Verified against direct optimization
and on 3000 random pairs.

**Empirically, the temperature series exclude nothing.** Across 73 matched
records in 16 series, one-sided 95% confidence bounds exclude neither mechanism
at any assumed correlation. The closest series falls short by 0.0175 on the
point estimate and 0.0634 on its 95% bound.

**But yeast ADH does.** The primary effects of Cha, Murray & Klinman (1989)
give `F_obs = +0.129` against `F0 = -0.042`, with a 95% bound of `+0.063`
(rising to `+0.079` at rho = 0.9). Because `L_H = 8.40 > gamma_SC` the identified
set is the *open* half-line, so no commitment can explain it away -- which
vindicates their 1989 argument against the later kinetic-complexity objection.
It clears `F0`, clears the semiclassical locus, and clears `B_vib` for
lambda > 10 kcal/mol, but not at lambda = 5 with strong driving force. Because
`F_obs > 0` its vacuity window is empty, so the conclusion is unconditional in
the equilibrium isotope effect rather than conditional on a bound.

**The experiment is precision limited, not commitment limited.** 52 of 73
records already have zero commitment gap. The signal is 0.042 and the median
sampling sd of F is 0.028, so 80% power needs a 1.68x reduction in sd, about
2.8x more replication. The shared tritium reference *helps*: it enters Var(F)
with weight `(gamma-1)^2 = 5.52` instead of `1 + gamma^2 = 12.21`, a factor 2.2
for measuring H/T and D/T in one triple-label mixture. No source reports the
covariance, so this is currently discarded.

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
  reversible.py         reversible-scheme identified set and the vacuity
                        threshold E_D* = (K_HT/K_DT)^(1/(gamma-1))
  vibronic_envelope.py  summed envelope B_vib; exact gated Franck-Condon
                        factors via tilted-Gaussian moments
  design_power.py       precision, correlation and commitment requirements
  decisive_case.py      yeast ADH primary effects (Cha, Murray & Klinman 1989)
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
python bounds_uncertainty.py # one-sided confidence bounds              (~3 s)
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
Best 95% lower bound -0.0975 (rho = 0) to -0.0933 (rho = 0.9), ecDHFR light
enzyme in both cases. Best point estimate -0.0596 (ecDHFR W133F), short of F0 by
0.0175; its 95% lower bound is short by 0.0634. The shortfall is deliberately
not quoted in standard errors: the series statistic is a maximum over
temperatures and is strongly right skewed, so the sampling s.d. (0.0853) and the
Gaussian back-calculation (point - bound)/1.645 (0.0279) differ by a factor of
three. Bounds are averaged over 5 replications of 5e4 draws; largest Monte Carlo
s.d. over all 48 cells is 7.4e-4.

## Relation to recent work

Williams (*J. Phys. Chem. B* **129**, 3604, 2025) observes that many parameter
combinations reproduce the same apparent KIE, and concludes derived intrinsic
values are of doubtful validity. That is the non-identification problem stated
qualitatively. This work characterizes the identified set exactly instead, and
finds it one-sided, so a robust inference survives in one direction.

Smedarchina & Siebrand (*Chem. Phys. Lett.* **410**, 370, 2005) reached the
direction of the envelope result in 2005 by instanton methods on a 2D model of
vibrationally assisted tunneling: assistance drives the tritium-referenced
exponent *toward* the semiclassical value from below, opposite to the then-
conventional assumption. Their Eq. (16) gives `e2` in [2.3032, 3.25725], and the
upper limit is exactly the integer-mass semiclassical exponent 3.25725 that
`masses.py` returns independently -- a five-decimal cross-check on both.

Their result is also the sharpest illustration of why the offset is the right
coordinate. The exponent converges to the semiclassical value along the ray,
which in exponent coordinates looks like the gated family merging with the
semiclassical locus. It does not merge: the offset is pinned at `F0 = -0.042086`
for every point on that ray. A ratio of two diverging logarithms cannot resolve
the bounded difference between them.

Siebrand & Smedarchina (*J. Phys. Chem. B* **108**, 4185, 2004) is the
full-length treatment of the same promoting-mode model class.

Francis & Kohen (*Perspect. Sci.* **1**, 110, 2014) is the community reporting
standard. It specifies error propagation but does not call for the covariance
between H/T and D/T effects sharing a tritium reference, which `design_power.py`
shows is worth a factor of 2.2 in that variance component.

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
