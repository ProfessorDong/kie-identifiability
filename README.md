# Network topology fixes the geometry of partial identification in enzymatic hydrogen transfer

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.21913975.svg)](https://doi.org/10.5281/zenodo.21913975)


Data, code and verification suite for the analysis of what competitive
multiple-isotope kinetics can establish about hydrogen tunneling in enzymes.

## The question and the answer

Competitive H/T and D/T measurements are the main evidence bearing on enzymatic
hydrogen tunneling. It is established that the Swain-Schaad exponent alone does
not diagnose it (Kohen & Jensen 2002; Hirschi & Singleton 2005; Shelton, Hrovat
& Borden 2007). This work asks what such an experiment does determine.

**The exponent is the wrong coordinate.** `X_H/X_D` is a ratio and is unchanged
when both logarithms are rescaled together, so it reports a displacement from the
mass-scaling line only relative to the size of the effects: it equals
`gamma + F/X_D` and therefore confounds the displacement with `X_D`. The quantity
that separates them is the offset `F = ln K_HT - gamma_SC ln K_DT`, which is zero
on that line by construction.

**The evidence is one-sided, and that is the whole result.** For any admissible
observation `K_HT > K_DT > 1`,

```
L_H = (K_HT-1)/(K_DT-1)  >  gamma_obs = ln K_HT / ln K_DT     ALWAYS
```

because `(x-1)/ln x` is strictly increasing. An exclusion is claimed exactly when
the reference exponent is below the observed one, and there `L_H > gamma`, so
Proposition S2 gives the OPEN half-line with `F_obs >= 0`: the entire identified
set lies above zero. Hence:

* an observed exponent ABOVE the semiclassical reference cannot be produced by
  series masking at any commitment -- the inference is robust with no assumption
  about how much masking there is;
* an observed exponent BELOW it leaves a set containing the semiclassical locus,
  which no precision can then exclude, because the obstruction is the set's
  direction, not its width. Such an observation is NOT empty: the set is still
  bounded below, so a mechanism whose attainable offsets lie entirely beneath
  that bound is still excluded.

This makes exact the criterion of Saunders (1985), used by Cha et al. (1989) and
Bahnson et al. (1993). Across 94 matched primary records comprising 29
independent systems, 28 systems fall on the uninformative side.

**Why it is one-sided: curvature against homogeneity.** The theorem above is a
computation; this is the reason behind it, and it is more general than the
scheme. Work in log-rate coordinates and write `h` for the map carrying an
intrinsic log-effect to the observed one, so `ln K_LT = h(X_L)` with the SAME `h`
for both isotope pairs. Then

* *mass scaling is homogeneity* -- `X_H = gamma X_D` is a ray through the origin
  and `F` is the displacement from it;
* *masking is curvature* -- for the competitive scheme
  `h(t) = ln[e^t (1+c)/(e^t + c)]`, with `h(0)=0`, `h' = c/(e^t+c)` in `(0,1)`,
  and `h'' < 0`.

A concave map fixing the origin is subhomogeneous, so it removes a larger
*fraction* of the larger log-effect -- which is the light isotope's -- and drags
every mass-scaled pair below the ray. Writing `psi = t - h(t)`, the whole proof is

```
F_int - F_obs = psi(X_H) - gamma psi(X_D)  >=  0
  psi convex, psi(0)=0  =>  psi(gamma t) >= gamma psi(t)     (superhomogeneity)
  F_int >= 0            =>  X_H >= gamma X_D, and psi is nondecreasing
```

needing only four transparent axioms on `h`: (A1) `h(0)=0`, (A2) nondecreasing,
(A3) `h(t) <= t` (masking cannot amplify), (A4) concave.

**Three consequences.** *Composition*: concavity through the origin survives
composition, so a chain of series bottlenecks of ANY length leaves the conclusion
intact -- this is the precise content of "kinetic complexity cannot manufacture
an above-reference observation". *Curvature fixes the direction*: `V/K` is a
conductance and a commitment puts the isotope-sensitive step in SERIES with an
isotope-blind one (`1/(V/K) = 1/k_on + k_off/(k_on k_i)`), which saturates and so
is concave; an isotope-blind route in PARALLEL adds conductances, makes `h`
convex, and pushes observations ABOVE the ray instead. *Sharing is a hypothesis*:
if `h_H != h_D` nothing follows.

**Sharing is approximate even in the primary protocol, and harmlessly so.** The
H/T and D/T effects use tritium references that differ at the non-transferred
position, so `c_D = r c_H` with `r` the secondary H/D effect on the reference.
What matters is the sign of `r-1`. Secondary effects here are normal (`r ~ 1.14`
in every source used), so the H comparison is the MORE masked one and the
imbalance pushes `F_obs` down: over 3e5 draws with `r` in [1,2] there is no case
with `F_obs > F_int`. An INVERSE secondary effect would reverse this and break
the argument (3e5 draws with `r` in [0.2,1] give violations in a third of cases,
by up to +0.43). The condition to check is therefore not that the commitments
match, which they never do, but that the H/T reference is the faster molecule.

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

**The direction the set is open follows the sign of `h''`.** Under a single
shared ADDITIVE commitment `C_f` instead of a single shared `c`,
`h(t) = ln[(e^t + C_f)/(1 + C_f)]` is convex rather than concave, both intrinsic
effects diverge together, and the set is open *below* rather than above. The two
families agree pair by pair but not jointly; the curvature is what distinguishes
them.

**The exponent is not even constant where the offset is.** Along the ridge
`X_L = P_L + s C_L` the Swain-Schaad exponent sweeps the OPEN interval
(2.846, 3.349) as `s` varies, approaching but not attaining its endpoints and
covering the range in which primary exponents are ordinarily read as diagnostic,
while `F` stays pinned at `F0` to ten decimal places because
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
below the semiclassical locus once excited channels contribute. The relaxation is
not uniform either: what Jensen's inequality establishes is that the BOUND moves
up, which is weaker than a statement about the attained maximum, and the sampled
grid contains parameter sets where the summed maximum lies BELOW the
single-channel value `F0` (-0.2645 at lambda=10, dG=-15, 3600 cm-1). Tabulated
`B_vib` values are maxima over a finite grid, not certified suprema.

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

**Empirically, the temperature series exclude nothing.** Across 83 matched
records in 18 series, one-sided 95% confidence bounds exclude neither mechanism
at any assumed correlation, now bracketed over `rho` in {-1, 0, 0.5, 0.9} with
`rho = -1` the maximally adverse case. The closest series falls short by 0.0175
on the point estimate and 0.0634 on its 95% bound.

**But yeast ADH does.** The primary effects of Cha, Murray & Klinman (1989)
give `F_obs = +0.129` against `F0 = -0.042`, with a 95% bound of `+0.063`
(rising to `+0.079` at rho = 0.9, falling to `+0.049` at the maximally adverse
rho = -1). Singled out from 29 systems it also deserves a multiplicity
adjustment: one-sided Bonferroni across 29 gives `+0.012` at rho = 0 and
`-0.013` at rho = -1, so the exclusion of the SEMICLASSICAL LOCUS survives
multiplicity and survives adverse correlation but not both at once, while the
exclusion of the GATED ENVELOPE survives every combination since all four bounds
exceed `F0`. Because `L_H = 8.40 > gamma_SC` the identified
set is the *open* half-line, so no commitment can explain it away -- which
vindicates their 1989 argument against the later kinetic-complexity objection.
It clears `F0`, clears the semiclassical locus, and clears `B_vib` for
lambda > 10 kcal/mol, but not at lambda = 5 with strong driving force. Because
`F_obs > 0` its vacuity window is empty, so the conclusion is unconditional in
the equilibrium isotope effect rather than conditional on a bound.

**Two obstacles, and only one yields to more measurement.** Precision is the
tractable one: the signal is 0.042 and the median sampling sd of F is 0.028, so
80% power needs a 1.68x reduction in sd, about 2.8x more replication. Masking is
not tractable. Replication locates `F_min` more precisely rather than higher, so
an endpoint that masking has pushed down stays down: a system exactly
semiclassical intrinsically but carrying `c = 1` has `F_min = -0.447` and cannot
be made decisive at any sample size. The frequent case `L_H >= gamma` (62 of 83
records) means the endpoint EQUALS `F_obs`, so computing it costs no assumption
about `c`; it does NOT mean masking is absent. The attainable margin is
`F_min - B`, fixed by the enzyme rather than the experiment. The shared tritium reference *helps*: it enters Var(F)
with weight `(gamma-1)^2 = 5.52` instead of `1 + gamma^2 = 12.21`, a factor 2.2
for measuring H/T and D/T in one triple-label mixture. No source reports the
covariance, so this is currently discarded.

**The argument does NOT reach the mixed-labeling secondary record.** Theorem 1
needs one commitment shared by both isotope pairs. That holds for the primary
experiment. It fails for the mixed-labeling design behind most secondary
measurements, where (Kohen & Jensen 2002) the secondary H/T effect is measured
for C-H cleavage while the secondary D/T effect accompanies C-D cleavage, so the
two ratios are referred to different molecules whose rates differ by the PRIMARY
isotope effect. With `c_H != c_D` the inequality fails outright: an intrinsic
pair sitting exactly on their 4.8 locus (`x_D = 1.10`, `x_H = x_D^4.8`) returns
an observed exponent of **7.85** at `c_H = 5`, `c_D = 1`. We therefore draw no
identification inference from that record in either direction; the observed
exponents are tabulated without inference. A single-label secondary design, with
both ratios against a common primary background, does satisfy the condition.

**Two amine oxidases test the direction, and the classical reading fails both.**
Bovine serum amine oxidase (Grant & Klinman 1989) is close to commitment free --
the authors establish by stopped-flow comparison that `V/K` is rate limited by
C-H cleavage -- and tunneling there is established by a criterion independent of
Swain-Schaad, namely `A_H/A_T = 0.12 +/- 0.04` against a semiclassical floor of
0.6. Its exponent is nonetheless BELOW `gamma_SC` at all six measured
temperatures (by 10.4 standard errors at 5 C), giving `F_min = -0.195` at 25 C.
The authors recorded the anomaly themselves, noting their failure to see the
inequality Saunders predicted. Monoamine oxidase B (Jonsson, Edmondson & Klinman
1994) adds 10 records over 2-43 C at two pH values, every exponent below
`gamma_SC`, with the lowest at the lowest temperature where the authors
independently identify a change in rate-limiting step. Masking and tunneling push
the observation the same way, which is why neither can be read off the exponent.

**What does NOT work.** Two recombinations of the existing data were tried and
both are rejected. Pooling across temperatures within a series would gain ~sqrt(n)
if the offset were constant; it is not (median chi2/dof = 7.8, up to 88), and the
weighted means are LOWER than the maxima. Pooling the nine horse liver forms
under a shared intrinsic pair with nine free commitments gives chi2 = 55.8 on 7
dof, p = 1e-9: the mutations change the intrinsic chemistry, not merely the
commitment, so the forms are not replicates. The obstacle is the systems, not the
estimator.

**Network topology fixes the geometry.** For any steady-state mechanism with one
isotope-sensitive step, each King-Altman term is a product of distinct rate
constants, so that constant enters `V/K` to at most first power above and below
the line. The observation map is therefore Moebius, `K(x) = (Ax+B)/(Cx+D)`
fixing `x=1`, with

```
h'' proportional to -(AD - BC)(AC x^2 - BD)
```

Two coefficients decide everything, and both are read off the topology rather
than fitted. `B != 0` only if some route reaches the monitored product without
passing through the isotope-sensitive step; `C != 0` only if that step
partitions against another, which is what a commitment is. So:

| competing branch | isotope-blind bypass | identified set |
|---|---|---|
| yes | no  | upper half-line |
| no  | yes | lower half-line |
| yes | yes | upper above `x* = sqrt(BD/(AC))`, lower below |

Both schemes in this work have `B = 0`; the reversible result is now a corollary
rather than a separate calculation.

What a bypass does to the identified set follows from an exact reduction, not
from curvature. With `phi = k_b/k_T` the bypass fraction and `q = k_2/k_T` the
competing rate, the general map factors as the series map acting on a contracted
effect,

```
K(x; phi, q) = K_series(Y; c),   Y = (x + phi)/(1 + phi),   c = q/(1 + phi)
```

so a bypass does not change the FORM of the masking, it pulls the intrinsic
effect toward unity first. Because `Y - 1 = (x-1)/(1+phi)`, the `L_H` case split
is unchanged. The endpoint then falls smoothly with `phi`, with no threshold,
and the two tritium references do not scale together: `phi_D = r phi_H`, so
bypass and reference asymmetry must be profiled jointly. With the measured
`r = 1.31`, the one exclusion in the record survives an isotope-blind route
below **18%** of the isotope-sensitive flux (13% against the semiclassical
locus). An earlier version of this file claimed a bound of "a fifth to a half"
derived from the curvature switch evaluated at the SMALLEST admissible competing
rate; that is the infimum of the switch, not a bound above it, and it established
nothing. (`analysis/network_geometry.py`)

**The half-line is one measurement from being a number, and one system has that
measurement.** Inverting the observation map at known commitment gives
the intrinsic effect in closed form,

```
x = K c / (1 + c - K)          finite above unity exactly when c > K - 1
```

and the offset follows as a point rather than a bound, once the reference
asymmetry `r = c_D/c_H` is known. One system in the record completes from
kinetics already published, with no new experiment; a second is attempted and
fails, instructively:

* *Yeast ADH completes.* Klinman (1976) measured `k_-1/k_cat = 1.3-7.3` for
  benzyl alcohol on the same enzyme, in the same laboratory, at the conditions
  Cha, Murray & Klinman (1989) used thirteen years later. Referenced to tritium
  that ratio is `c_H = (K_HT - 1) + K_HT a`, giving `c_H = 15.4-58.2`, well above
  the singular value `K_HT - 1 = 6.13`. With `c_D = 1.31 c_H` the identified set
  collapses from `F > +0.129` to `F_int = +0.208 to +0.513`, an intrinsic
  exponent of 3.72-4.23 against a semiclassical 3.349. Propagating every measured
  input (both effects, `r`, and `a`) by Monte Carlo leaves `P(F_int > 0) = 1.0000`
  at every `a`: the conclusion is unchanged in direction and strengthened in size.
* *Bovine serum amine oxidase does not.* Grant & Klinman's Table IV (1989) sets
  the pre-steady-state isotope effect beside the steady-state one at **six**
  temperatures, and their ratio is a masking factor directly. But the six are not
  consistent with one value (`Q = 17.9` on 5 dof, `p = 0.003`), and masking
  cannot exceed unity, which the 0 C ratio does by 3.5 sigma. The competitive
  effects exist only at 25 C, where the matched factor is `0.839 +/- 0.097`, the
  most extreme of the six. Full propagation gives `F_int = -0.038`, 95%
  `[-0.267, +0.242]`, straddling both thresholds: `P(F_int < F0) = 0.49`. An
  earlier version of this file reported `-0.188 to -0.107` from a factor pooled
  over 15-45 C with only that factor's uncertainty propagated; the pooling is
  what the data reject.

A commitment closes the half-line only when it is precise enough AND measured at
the condition of the isotope pair. Yeast ADH meets both; the amine oxidase meets
neither. Two further systems resist completion for reasons that are diagnostic
rather than accidental.
Monoamine oxidase B does not complete because the stopped-flow comparison Jonsson
et al. (1994) invoke is attributed to an earlier source and its value not
reproduced, while their own `D(V/K)` is a second measurement of the same masked
quantity. Horse liver ADH does not complete because its alcohol Michaelis-constant
isotope effect is *inverse* -- `K_m` = 0.03 mM against 0.08 mM for the
`alpha,alpha-d2` analogue in Sekhar & Plapp (1990), a ratio of 0.375 -- and
solving the two-step relation for `k_off` then returns a negative number. Both
remain bounds.

The same inversion is what keeps the yeast range from being narrowed further:
Klinman's Table III reports 1.3 on the coenzyme constant, which maps to `a = 2.3`,
but 0.80 on the alcohol constant, which the two-step partition model cannot
produce at all, since `K_m = (k_off + k)/k_on` with `k_H > k_D` forces that ratio
above unity. The inverse alcohol effect is systematic in these enzymes rather
than an outlier, and the mapping inherits the model risk. What survives it is the
direction, which survives completely: for *any* admissible commitment the
completed yeast offset exceeds the half-line endpoint.

The precision required of `c` is mild. Writing `c = m(K_HT - 1)`, recovering the
yeast offset to within half the mechanistic signal needs `c_H` to only **5.3%**
at `m = 3`, relaxing to 11.2% at `m = 5` (with `c_D = r c_H`). Commitments are routinely measured to that
precision. The bound is not merely loose; it is loose in the direction that hides
the signal.

**The obvious alternative fails, instructively.** Since `k_off` is diffusional
and the chemical step is not, Stokes-Einstein makes `c` scale as `1/eta`, so one
system at two viscosities over-determines the intrinsic pair and returns it in
closed form, with the two isotope pairs giving the same `c` as a consistency
check. Exact recovery on 2e4 random triples. But the inversion amplifies
relative error by 3x at `c = 10` and 30x at `c = 1`, demanding the isotope
effects to 0.7%-0.07% against roughly 1% in practice, and it is most demanding
exactly where masking is heaviest.

**What a closed offset would bear on: the promoting mode.** In the gated model
`F` depends on the gating width `w = 2 kappa sqrt(mu_H) sigma^2(T)` AND on the
scale parameter `A`, and for a harmonic promoting mode the exact thermal
variance is

```
w(T) = w0 coth(hbar omega_g / 2 k_B T)
```

a quantum expression in which zero-point motion holds `w` fixed below
`theta = hbar omega_g / k_B`. The offset inherits that structure. Across
250-350 K the predicted span of `F` at `A = 5`, `w0 = 0.5` is 0.12|F0| at
50 cm-1, **1.26|F0| at 200**, **1.47|F0| at 400**, and 0.51|F0| at 800, while
the Swain-Schaad exponent moves 0.01-0.02 over the same interval. Sensitivity peaks where the mode sits at the
quantum-classical crossover `hbar omega_g ~ 2 k_B T`, which at 300 K is
417 cm-1; sweeping `A` over [1,50] and `w0` over [0.1,2] moves that peak only
between 140 and 480 cm-1, so it tracks the crossover rather than the
parameterization. Promoting modes invoked in this literature lie in that window.
(`analysis/completion.py`)

Sensitivity is not identifiability, and only the first is claimed. `F` depends on
`A`, `w0` and `theta` jointly and the three are badly confounded: for an
11-point series over 250-350 K the Jacobian in `(ln A, ln w0, ln theta)` has
singular values 0.60, 2.7e-2, 1.4e-4, so `cond(J'J) ~ 1.7e7`. Below `theta` the
data carry almost nothing about `theta`; above it only `w0/theta` is determined.
No system supplies a completed `F(T)` series: yeast ADH gives one temperature,
the amine oxidase one isotope pair. Recovering `omega_g` is therefore a design
target, not a result here. An experiment long read as a yes-or-no tunneling test,
completed and read in the offset coordinate, returns a number whose temperature
dependence carries the gating variance; that is where a quantitative test of
promoting-mode models would begin.

The framework, rather than any verdict on tunneling, is the contribution:
mechanistic claims from these experiments are comparisons between an
experimentally identified set and a computed mechanism envelope.

## Layout

```
analysis/
  masses.py             isotope mass convention, imported everywhere
  audit_numbers.py      sweeps every derived number against both manuscript
                        files; catches fabricated and stale values that a
                        whitelist check cannot
  network_geometry.py   Moebius structure of V/K; curvature class and the
                        identified-set direction from network topology
  completion.py         closing the half-line with one commitment measurement,
                        and what the closed offset measures about the gating mode
  curvature.py          WHY the evidence is one-sided: masking is concave
                        through the origin in log-rate coordinates, mass scaling
                        is homogeneity; axioms, composition, direction, and the
                        unequal-reference check of the primary protocol
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
  export_fig_asym.py    identified-set figure data (29 systems)
  export_si_tables.py   supplementary tables; writes to the manuscript tree if
                        present, otherwise to results/
  export_benchmark_table.py  complete 94-record benchmark table
  export_profiles_table.py   per-series profile table
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
python bounds_uncertainty.py # one-sided confidence bounds, rho in {-1,0,.5,.9}
python curvature.py          # the one-sidedness theorem + controls    (~2 min)
python network_geometry.py   # topology -> curvature -> set direction
python completion.py         # closing the set; the gating-mode readout
python audit_r2.py           # reconfirms the second review's claims   (~10 min)
python audit_v3.py           # adversarial audit, must exit 0          (~10 min)
python vibronic.py           # the vibronic-sum penalty                 (~5 min)
python offset_analysis.py    # per-series profiles                     (~50 min)
python offset_summary.py     # stratified interpretation
python export_figs_v3.py     # pgfplots tables
cd ../figures/tikz && ./build.sh
```

Every stochastic step is seeded, so a rerun reproduces the reported values
exactly. `data/` is inputs only: no script writes to a file it reads.

One caveat on timing: `offset_analysis.py` computes profile likelihoods by
continuation over a 301-point grid for each of 18 series and takes roughly
50 minutes on one core.

## Key numbers

Reduced masses of the X-H oscillator, CODATA/AME2020, carbon donor. This is the
convention of Kohen & Jensen, who obtain 3.34 from the reduced mass of 12C and
the hydrogen isotopes; the bare-mass value 3.2628 is also in circulation.

```
gamma_SC = 3.34887   rigid = 2.45884   F0 = -0.042086   threshold = 1.26890
```

Temperature-series benchmark: 83 matched records, 18 series, 5 enzyme families,
6 enzyme-organism systems, 2-45 C. Transcription audit: 244 values across the
four PDF-derived sources, 0 misses. The wider matched-primary record used for the
identified-set figure is 94 records over 29 systems, of which 28 fall on the
uninformative side.

Bounds: 0 of 18 series exclude either mechanism, at correlation -1, 0, 0.5 or
0.9. Best 95% lower bound -0.0987 (rho = -1) to -0.0933 (rho = 0.9), ecDHFR light
enzyme in every case. Best point estimate -0.0596 (ecDHFR W133F), short of F0 by
0.0175; its own 95% lower bound is short by 0.0634 (rho = 0) and 0.0682 at the
maximally adverse rho = -1. The shortfall is deliberately not quoted in standard
errors: the series statistic is a maximum over temperatures and is strongly right
skewed, so the sampling s.d. (0.0853) and the Gaussian back-calculation
(point - bound)/1.645 (0.0279) differ by a factor of three. Bounds are averaged
over 5 replications of 5e4 draws; largest Monte Carlo s.d. over all 72 cells is
7.4e-4.

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

A fourth was caught internally: a monotonicity assertion that failed on 21 of the
73 records then in the benchmark, which led to the `L_H` criterion now in the paper.

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
