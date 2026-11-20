# Limits on the Intrinsic Swain-Schaad Exponent from Competitive Isotope Effects When the Commitment to Catalysis Is Unknown

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.21913975.svg)](https://doi.org/10.5281/zenodo.21913975)


Data, code and verification suite for the analysis of competitive H/T and D/T effects on enzyme specificity constants.

For the data, start with `results/matched_measurements.csv`, the ledger of 247 distinct published pairs. The dataset analyzed is its 244 pairs in 76 analysis units (48 temperature series and 28 single-condition records; `analysis/corpus.py`). The other three pairs are the additional concentrations of a cofactor titration that its source used as a control; `analysis/cofactor_benchmark.py` tests the kinetic scheme on that titration and repeats the counts with the three conditions treated as units (`results/nad_condition_bounds.csv`, `results/nad_scheme_fit.csv`, `results/all_condition_units.csv`, `results/condition_screen_sensitivity.csv`). The Reproduce section below gives instructions for the other calculations.

## The question and the answer

Competitive H/T and D/T measurements are the main evidence bearing on enzymatic
hydrogen tunneling. It is established that the Swain-Schaad exponent alone does
not diagnose it (Kohen & Jensen 2002; Hirschi & Singleton 2005; Shelton, Hrovat
& Borden 2007). Here the observed isotope effects constrain the corresponding intrinsic effects on the chemical step.

**The bound on the intrinsic exponent.** For normally ordered effects in the
minimal competitive scheme with a shared tritium tracer, isotope-insensitive
binding, one required irreversible hydrogen-transfer step and positive commitment, the intrinsic Swain-Schaad exponent can take any value above the observed
exponent and none at or below it: the allowed range is `(gamma_obs, infinity)`.
The measured pair bounds the intrinsic relation from below and not from above.

**The offset is the coordinate for model tests.** `X_H/X_D` is a ratio, equal to
`gamma + F/X_D`, so it reports a displacement from the mass-scaling line only
relative to the size of the effects. The offset
`F = ln K_HT - gamma_SC ln K_DT` is zero on that line by construction, and the
attainable offsets of a model can be compared directly with the allowed range.
It adds no information to the pair.

**The inequality behind the lower bound.** For any admissible
observation `K_HT > K_DT > 1`,

```
L_H = (K_HT-1)/(K_DT-1)  >  gamma_obs = ln K_HT / ln K_DT     ALWAYS
```

because `(x-1)/ln x` is strictly increasing. An exclusion is claimed exactly when
the reference exponent is below the observed one, and there `L_H > gamma`, so
Proposition S2 gives the OPEN half-line with `F_obs >= 0`: the entire identified
set lies above zero. Hence:

* an observed exponent ABOVE the chosen mass reference requires an intrinsic
  exponent above that reference for every admissible commitment in this scheme;
* an observed exponent BELOW it leaves a set containing the semiclassical reference,
  which no precision can then exclude, because the obstruction is the set's
  direction, not its width. Such an observation is NOT empty: the set is still
  bounded below, so a mechanism whose attainable offsets lie entirely beneath
  that bound is still excluded.

This makes exact the kinetic part of the criterion of Saunders (1985) as used by
Cha et al. (1989) and Bahnson et al. (1993). The dataset contains 244 matched primary records in 76 analysis units
from a targeted literature search. At their reported central values,
68 units lie below the semiclassical reference at every admissible temperature, so no gain in precision at those values can show an
intrinsic departure from it; eight are observed above it, four of them at the
nominal 95% level and one after adjustment across all units
(`reference_counts.py`). Against the ground-channel model with unrestricted
parameters, 61 are compatible with it at their reported values and nine more
are limited by precision alone; six exclude it
nominally (yeast ADH with benzyl and with p-chlorobenzyl alcohol, horse liver
ADH F93W in 40% and in 60% methanol, and its F93W-I224G and I224V double
mutants), and only the two yeast units survive adjustment across all 76 at
rho = -1. Only the p-chlorobenzyl alcohol unit also survives removal of its
isotopic-purity correction at the same time (`isotopic_purity.exclusions`).

**A cofactor titration tests the scheme.** [Tsai and Klinman](https://doi.org/10.1021/bi002075l)
(Table 5) measured the F93W horse liver enzyme at 3 °C in 40% methanol at four
NAD⁺ concentrations, as a control on the order of substrate binding; the 10 mM
condition is the 3 °C record of their temperature series. The four conditions
share one intrinsic pair, so under the minimal scheme they can differ only in
commitment. That requires a D/T effect above 1.62 and 1.63 at 0.02 and 1 mM,
where 1.03 ± 0.04 and 1.2 ± 0.1 are measured (12.8 and 4.2 standard errors
below); no single intrinsic pair fits the four conditions (chi-square 276 at
r = 1, worse at larger r). The minimal scheme is therefore contradicted at the
two subsaturating concentrations. The two saturating conditions attain the
largest D/T effect and carry no test; their pairs agree with each other, as
the scheme requires when the cofactor binds first. A random order of binding
does not explain the pattern: with one ternary complex it changes only the
commitment and would raise both effects at low cofactor, whereas the H/T and
D/T effects fall by 5.1 and 15.4 standard errors between 10 and 0.02 mM
(`random_order_check`, `saturating_fall`). The three additional
concentrations are not counted as analysis units. The yeast pairs were
measured at 10 mM NAD⁺, twice the Michaelis constant their source quotes, and
have not been checked in this way; the titration fails at 1.3 times the
Michaelis constant of the titrated enzyme (`cofactor_saturation`). If they are
counted, 68 of 79 units are below the semiclassical reference, six of the
eleven above it clear zero nominally and one after adjustment, and the
ground-channel counts are 61 compatible, ten limited by precision and eight
nominal exclusions, with the same two yeast units after adjustment under
lognormal errors. With Student t tails at ten degrees of freedom both lose
adjusted significance, with 76 or with 79 units
(`analysis/cofactor_benchmark.py`).

Two conditions go with these counts. The comparison is with the model's
parameters unrestricted: F0 is approached only as the action and the gating
width grow without bound, so an exclusion based on it holds under any
restriction, while the 61 units that cannot exclude the unrestricted model fall
to 56, 42 and 35 as the parameters are bounded more tightly (`model_domain.py`).
All eight nominal ground-channel exclusions in the expanded screen use a mixed-label design, whose D/T effect
is that of the doubly labeled substrate; comparing it with a single-coordinate
relation or model assumes the rule of the geometric mean
(`reference_asymmetry.geminal_ratio`).

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
isotope-insensitive one (`1/(V/K) = 1/k_on + k_off/(k_on k_i)`), which saturates and so
is concave; an isotope-insensitive route in PARALLEL adds conductances, makes `h`
convex, and pushes observations ABOVE the ray instead. *Sharing, or ordering, is a
hypothesis*: if `h_H` and `h_D` are unrelated nothing follows, but ordered
commitments (`c_D >= c_H`) suffice.

**Sharing is approximate even in the primary protocol, and harmlessly so.** The
H/T and D/T effects of 27 of the 74 forms use tritium references that differ at
the non-transferred position, so `c_D = r c_H` with `r` the secondary H/D effect
on the reference (25 of them report the secondary effects that check its sign);
the other 47 use one tritiated molecule for both, so `r = 1`
(`data/tracer_design.csv`).
What matters is the sign of `r-1`. Observed secondary ratios here are normal
(`R = 1.02-1.38`, a lower bound on `r`), so the H comparison is the MORE masked one and the
imbalance pushes `F_obs` down: over 2e5 draws with `r` in [1,2] there is no case
with `F_obs > F_int`. An INVERSE secondary effect would reverse this and break
the argument (2e5 draws with `r` in [0.2,1] give violations in a third of cases,
by up to +0.43). The condition to check is therefore not that the two tracer
commitments match (they do exactly for the 47 shared-tracer forms, and need not
otherwise), but that the H/T tracer is the faster molecule.

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

**Mechanisms enter through the offsets they can attain.** For the ground-channel
tunneling model the attainable offsets have a *scale-dependent* envelope:

```
B(K) = ln K - (gamma/2) ln[ K^2 (t-1) / (t - d + K^2 (d-1)) ]   K < sqrt(t)
     = F0 = -0.042086                                            K >= sqrt(t)
```

with threshold `sqrt(t) = (mu_T/mu_H)^(1/4) = 1.2689`. There is no flat bound: `dF/dA < 0`
identically, so the unconstrained supremum is 0, approached as the isotope
effects collapse to unity. Proved in the supplement; the proof turns on
`p_L' = q_L/2`, which cancels the intercept and slope contributions to `X_D`
along the constraint and leaves `dX_D/dw < 0` fixed by `t > d > 1` alone.
Numerical agreement is 7e-7, which is the double-precision noise floor, not the
numerical precision: at 50 digits the supremum approaches `F0` strictly from below.

**Excited vibrational channels.** Summing over channels can produce a positive
offset. A moderate-parameter example gives `F = +0.005118` (H/T = 1.52050,
D/T = 1.13156) at channel cutoffs 6, 30 and 60. Directly multiplying rates at
large action caused floating-point underflow in the earlier extended-grid
calculation. Logarithmic evaluation gives sampled maxima +0.52503, +0.19068,
+0.34910 and +0.59232 for the four parameter sets in SI Note 14, each at the
chosen upper action limit. These grid values leave the global supremum open.
The experimental comparisons use only the separate ground-channel envelope.

**Reversibility opens a bounded vacuity window.** Deriving the map for a
reversible chemical step (it is Northrop's equation re-referenced to tritium,
which is what makes both commitments shared) gives
`F = F_obs + gamma ln D_D - ln D_H`. The set is unbounded below exactly when the
corner `Cf* = (v-u)/d`, `Cr* = (a-b)/d`, `d = av-bu` lies in the physical
quadrant, which needs BOTH `d > 0` and `v >= u`:

```
window = [E_D*, E_D**),   E_D* = (K_HT/K_DT)^(1/(gamma-1))
E_D** = root above E_D* of  a(K_DT/E - 1) = b(K_HT/E^gamma - 1)
window is nonempty  <=>  F_obs < 0
```

It is a bounded interval, not a half-line: `d -> -(a-b) < 0` at large `E_D`, so
the set recovers a finite endpoint above the window as well as below it. And a
positive observed offset implies an intrinsic offset at least as large, for the
shared-tracer reversible model with mass-scaled equilibrium isotope effects.
The algebraic crossing near 1.83 for the yeast benzyl pair is outside the
physical domain. The result needs the shared tracer: with two tracers whose
commitments differ, a reference-obeying pair can be observed above the reference
(`reversible.ordered_tracer_counterexample`), so applying it to the mixed-label
yeast assays needs an additional condition. Across the 209 temperature-resolved
records that have a window, its lower edge runs 1.097-2.663 (median 1.451), with
upper edge median 1.961 and median width 0.408. Verified against direct
optimization and on 3000 random pairs.

**Empirically, few temperature series exclude anything.** Across the 216 temperature-resolved
records in 48 series, one-sided 95% confidence bounds bracketed over `rho` in
{-1, 0, 0.5, 0.9} (`rho = -1` the maximally adverse case) exclude both mechanisms
at each of the four tested correlations in three series: horse liver ADH F93W in 40% and 60%
methanol (Tsai & Klinman 2001), bounds `+0.032` and `+0.070` at rho = -1, and
yeast ADH with p-chlorobenzyl alcohol (Rucker et al. 1992), `+0.170`. Of the
rest, only the thymidylate synthase Y209W hydride transfer clears `F0` at any
correlation (`-0.039` at rho = 0.9); its point endpoint, `-0.009`, and those of
W80M, `-0.001`, and R67 DHFR, `-0.017`, lie above `F0` (see Bounds below).

**Yeast alcohol dehydrogenase.** The primary effects of Cha, Murray & Klinman (1989)
give `F_obs = +0.129` against `F0 = -0.042`, with a 95% bound of `+0.063`
(rising to `+0.079` at rho = 0.9, falling to `+0.049` at the maximally adverse
rho = -1). Singled out from 76 analysis units it also deserves a
multiplicity adjustment: one-sided Bonferroni across 76 gives `+0.0004` at rho = 0
and `-0.027` at rho = -1, so the exclusion of the SEMICLASSICAL REFERENCE survives
multiplicity and survives adverse correlation but not both at once, while the
exclusion of the GROUND-CHANNEL MODEL survives every combination since all four bounds
exceed `F0` -- though with little room: under both at once the published
errors need only be understated by 9.6% to overturn it. The same enzyme with
p-chlorobenzyl alcohol (Rucker et al. 1992, D/T corrected for protium) excludes
more firmly: `+0.170` at rho = -1, and `+0.050` after the adjustment across 76
units (`yadh_robustness.pchloro_series`). These bounds use a
normal approximation that treats the published errors as standard errors. Cha
et al.'s Table 2 lists four D/T determinations, three paired with H/T (the
fourth is in `data/cha1989_yadh_dt_only.csv`). Because `L_H = 8.40 > gamma_SC`
the identified set is the *open* half-line, so no forward commitment can explain
it away, which is the exact form of the argument Cha et al. made from Northrop's
expressions. No exclusion of the summed vibronic model is claimed at any
reorganization energy. For the shared-tracer reversible map with mass-scaled equilibrium effects, the
window is empty at positive observed offset. Applying this shared-tracer
result to the actual mixed-label experiment requires a separate tracer argument.

**Two obstacles, and only one yields to more measurement.** Precision is the
tractable one: the signal is 0.042 and the median sampling sd of F is 0.038 (rho = 0), so
80% power needs a 2.24-fold reduction in sd, about 5.0 times more replicates
(`results/design_power.txt`). Masking is
not tractable. Replication locates `F_min` more precisely rather than higher, so
an endpoint that masking has pushed down stays down: a system exactly
semiclassical intrinsically but carrying `c = 1` has `F_min = -0.447` and cannot
be made decisive at any sample size. The case `L_H >= gamma` means the endpoint EQUALS `F_obs`, so computing it costs no assumption
about `c`; it does NOT mean masking is absent. The attainable margin is
`F_min - B`, fixed by the enzyme rather than the experiment. If both measurements contained the same reference-rate error, its contribution to Var(F) would enter
with weight `(gamma-1)^2 = 5.52` instead of `1 + gamma^2 = 12.21`, a factor 2.2
for an idealized triple-label mixture. The sources do not report the actual
covariance, and the mixed-label experiments can use chemically different tracers.

**Forward masking in the mixed-labeling secondary record is ordered, and it
deflates.** In that design (Kohen & Jensen 2002) the secondary H/T effect is
measured for C-H cleavage and the secondary D/T effect for C-D cleavage, so the
two reference molecules (HT and DT, transferred isotope first) differ at the
transferred position and `c_D/c_H = k_HT/k_DT`, the primary H/D effect: above
one for a normal primary effect. That is the ordering of Proposition S4, whose
proof needs only a reference exponent above one, so forward masking can only
lower an observed secondary exponent. An intrinsic pair on the 4.8 locus
(`x_D = 1.10`) returns 2.57 at `c_H = 1`, `c_D = 5`; only the anti-ordered
`c_H = 5`, `c_D = 1` gives 7.85, and that needs an inverse primary effect of
five. No ordered map in 4e5 random draws raises the offset (`mixed_label.py`).
Kohen and Jensen's inflation comes from a reverse commitment acting with a
secondary equilibrium isotope effect, outside forward masking. No envelope
inference is drawn from that record, because the reverse commitment and
equilibrium effect are unbounded and the ground-channel envelope has no
secondary coordinate. For the 18 published secondary exponent comparisons,
recalculation from the printed ratios differs by up to 0.198 in the exponent.
All differences lie within the quoted uncertainties. (Until 2026-09-15 this file presented the anti-ordered
example as what masking does in that design.)

**Additional evidence in two amine oxidases.**
Grant & Klinman (1989) interpreted bovine serum amine oxidase as having a
small commitment on the basis of their stopped-flow comparison. The steady-state
effect used for that comparison was itself inferred through a mass relation.
They also interpreted the Arrhenius prefactor ratio `A_H/A_T = 0.12 +/- 0.04`
against their semiclassical lower limit of 0.6 as evidence for tunneling. Its exponent is nonetheless BELOW `gamma_SC` at all six measured
temperatures (by 10.4 standard errors at 5 C), giving `F_min = -0.195` at 25 C.
The authors recorded the anomaly themselves, noting their failure to see the
inequality Saunders predicted. Monoamine oxidase B (Jonsson, Edmondson & Klinman
1994) adds 10 records over 2-43 C at two pH values, every exponent below
`gamma_SC`, with the lowest at the lowest temperature where the authors
independently identify a change in rate-limiting step. Masking and tunneling push
the observation the same way, which is why neither can be read off the exponent.

**What does NOT work.** Two recombinations of the existing data were tried and
both are rejected. Pooling across temperatures within a series would gain ~sqrt(n)
if the observed offset were constant; it is not (median chi2/dof = 3.5, up to
71), and the weighted means are LOWER than the maxima. With one shared tracer per pair, pooling enzyme variants
that share one intrinsic pair but differ in unknown commitment cannot identify
the scale even with exact data: with `z = 1 - 1/K`, `u = 1 - 1/x` and `s = c/(1+c)`, the
map is exactly `z_Lj = u_L s_j`, so `(u, s) -> (lam u, s/lam)` leaves every
observation unchanged (Proposition S8). With a known common unequal tracer
ratio, two distinct masking levels instead determine both intrinsic effects.
The intercept and slope of 1/z_D as a function of 1/z_H identify the intrinsic
pair. A ratio merely assumed to exceed unity does not fix that intercept.
Of 11 nominal parameters for the nine
horse liver forms only 10 are identified; the identifiable refit gives
chi2 = 55.8 on 8 dof, p = 3e-9, so the common-pair model is rejected, without
saying which of its assumptions fails (`analysis/pooling_check.py`). The
obstacle is the systems, not the estimator.

**Network topology classifies the geometry.** For any steady-state mechanism with one
isotope-sensitive step, each King-Altman term is a product of distinct rate
constants, so that constant enters `V/K` to at most first power above and below
the line. The observation map is therefore Moebius, `K(x) = (Ax+B)/(Cx+D)`
fixing `x=1`, with

```
h'' proportional to -(AD - BC)(AC x^2 - BD)
```

Two coefficients fix which behaviors are possible, and both are read off the
topology rather than fitted. `B != 0` only if some route reaches the monitored product without
passing through the isotope-sensitive step; `C != 0` only if that step
partitions against another, which is what a commitment is. So:

| commitment | bypass | class | identified set of the offset |
|---|---|---|---|
| yes | no  | concave | upper half-line (a lower limit, no upper limit) |
| negligible (a limit) | yes | convex | lower half-line (an upper limit, no lower limit) |
| yes | yes | conditional | all of R; upper half-line once the bypass ratio is bounded |
| label after the first irreversible step | any | blind | `K = 1`: nothing is learned |

The convex row is a limit, not a scheme: it is approached when dissociation is
much faster than both routes to product, and without any competing exit the
competitive effect is unity. A symbolic enumeration of 194 chain mechanisms
with one to four bound states gives 62 concave, 62 conditional and 70 blind
cases and none that is convex; in all 62 conditional mechanisms the behavior
of the offset depends on the rate constants (`network_atlas.py 4`).

The identified sets are for unknown rate constants. With `K(1) = 1`, `B = 0`
forces the one-parameter map `x(1+c)/(x+c)` and `C = 0` the one-parameter map
`(x+phi)/(1+phi)`, whatever the network. With both present the map has two free
parameters, commitment and bypass, and nothing is identified until one is
bounded. The curvature switches at `x* = sqrt(BD/(AC))`; at given rate constants
the offset along the mass-scaling ray takes both signs exactly when `x* > 1`,
but its sign change is not at `x*`. Both schemes in this work have
`B = 0` for each map, but the reversible scheme inherits the half-line only when
its equilibrium isotope effects are unity. When they depend on isotope, the H/T
and D/T comparisons carry different maps, and the set is the vacuity-window
result above: inside the window it is unbounded below (`reversible.py`).

What a bypass does to the identified set follows from an exact reduction, not
from curvature. With `phi = k_b/k_T` the bypass fraction and `q = k_2/k_T` the
competing rate, the general map factors as the series map acting on a contracted
effect,

```
K(x; phi, q) = K_series(Y; c),   Y = (x + phi)/(1 + phi),   c = q/(1 + phi)
```

so a bypass does not change the FORM of the masking, it pulls the intrinsic
effect toward unity first. Because `Y - 1 = (x-1)/(1+phi)`, the `L_H` case split
is unchanged. For the yeast pairs the endpoint then falls smoothly with `phi`, with no
threshold (this is not general: for observed effects near unity it can first
rise), and the two tritiated tracers do not scale together: `phi_D = r phi_H`, so
bypass and reference asymmetry must be profiled jointly. The tolerance falls as
`r` grows, and the observed secondary ratio `1.31` is only a lower bound on `r`
(masking lowers it). Solving `r` jointly with the published commitment gives
1.36-1.57, and the yeast exclusion then survives an isotope-insensitive
route below **13-17%** of the isotope-sensitive flux (10-12% against the
semiclassical reference); propagated, median 16% with 95% interval 8-25%. Because
that `r` was solved without a bypass, these are sensitivities at stipulated `r`,
not a joint fit under the enlarged model. The closed form
`E_r(phi) = ln[K_HT + (K_HT-1)phi] - gamma ln[K_DT + (K_DT-1) r phi]` is the
endpoint only when the commitment infimum sits at `c -> infinity`; otherwise the
minimum is interior and lower, and `network_geometry.endpoint_exact` gives it
from a stationarity quadratic (it matches direct profiling to 1e-9 on 3000
random triples; the yeast system is in the `c -> infinity` branch throughout).
Through v1.8.1, `network_geometry.endpoint` profiled the D side at the H-side
commitment rather than at `r(1+phi)/(1+r phi)` times it, wrong for `r != 1` in
the interior branch only; no reported number was in that branch. With binding
effects in [0.90, 1.10] as well, the joint condition
`F_bind <= 0.164 - 1.465 phi` holds at every stipulated `r` from 1.31 to 1.57
(`joint_nuisance.check_envelope_r`); the line `0.1656 - 1.0866 phi` holds at
`r = 1.31` only. No bound on the bypass follows from the curvature switch
evaluated at the SMALLEST admissible competing rate: that is the infimum of the
switch, not a bound above it. (`analysis/network_geometry.py`)

**How an independent commitment measurement would recover the intrinsic pair.** Inverting the observation map at a known reciprocal commitment `c` gives
the intrinsic effect in closed form,

```
x = K c / (1 + c - K)          finite above unity exactly when c > K - 1
```

and the offset follows as a point rather than a bound, once the reference
asymmetry `r = c_D/c_H` is known. One system in the record completes, conditionally, from
a published kinetic estimate, with no new experiment; a second is attempted and
fails, instructively:

* *Yeast ADH completes, conditionally.* Klinman (1976) inferred
  `k_-1/k_cat = 1.3-7.3` across substituted alcohols using a historical
  kinetic reconstruction, under different conditions from the competitive
  measurements of Cha, Murray & Klinman (1989). Specifically, her Table IV converts Michaelis-constant isotope effects into
  `k_-1/k_cat` assuming an isotope effect of 4 on `k_cat` and none on binding,
  and 1.3-7.3 is the span over the effects (1.1-1.5) seen across her
  para-substituted alcohols, not an interval for benzyl alcohol. Referenced to tritium
  that ratio is `c_H = (K_HT - 1) + K_HT a`, giving `c_H = 15.4-58.2`, well above
  the singular value `K_HT - 1 = 6.13`. With `c_D = r c_H` and `r` solved jointly
  from the primary and secondary effects at the same commitment (1.36-1.57), the
  identified set narrows from `F > +0.129` to the interval
  `F_int = +0.209 to +0.534`, an intrinsic exponent of 3.72-4.27 against a
  semiclassical 3.349. Propagating every measured input (both primary and both
  secondary effects, and `a` uniform on its span as a sensitivity choice) by
  Monte Carlo gives no draw of 4e5 with `F_int <= 0` at any `a`.
* *Bovine serum amine oxidase does not.* Grant & Klinman's Table IV (1989) sets
  the pre-steady-state isotope effect beside the steady-state one at **six**
  temperatures. Their ratio would be a masking factor if the steady-state column
  were measured independently; it was derived from the competitive effects
  through the Swain-Schaad relation, so it is not an independent constraint.
  Read naively anyway, the six are not
  consistent with one value (`Q = 17.9` on 5 dof, `p = 0.003`), and masking
  cannot exceed unity, which the 0 C ratio does by 3.5 sigma. The competitive
  effects exist only at 25 C, where the matched factor is `0.839 +/- 0.097`, the
  most extreme of the six. Full propagation gives `F_int = -0.038`, 95%
  `[-0.267, +0.242]`, straddling both thresholds: `P(F_int < F0) = 0.49`. A factor
  pooled over 15-45 C, with only that factor's uncertainty propagated, would
  give `-0.188 to -0.107`; the pooling is what the data reject.

A commitment closes the half-line only when it is precise enough AND measured at
the condition of the isotope pair. The historical yeast estimate illustrates
the calculation conditional on transferring its assumptions across conditions.
The amine oxidase comparison cannot provide a matched commitment. Two further systems resist completion for reasons that are diagnostic
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
than an outlier, and the mapping inherits the model risk. What does not depend on the estimate is
the direction, within the base map: for *any* admissible commitment the
completed yeast offset exceeds the half-line endpoint.

The commitment's own precision budget is mild. Writing the reciprocal commitment as `c = m(K_HT - 1)`, its
contribution to the yeast offset stays within half the mechanistic signal if `c_H`
is known to **5.3%** at `m = 3` (11.2% at `m = 5`). That is not the whole budget:
at known `c` the published isotope errors alone give `F` a standard error of
0.044 (rho = 0) to 0.057 (rho = -1), above the 0.021 target, so the isotope pair
must be remeasured too (`completion.isotope_error_budget`). The bound is not merely loose; it is loose in the direction that hides
the signal.

**The obvious alternative fails, instructively.** Since `k_off` is diffusional
and the chemical step is not, Stokes-Einstein makes `c` scale as `1/eta`, so one
system at two viscosities over-determines the intrinsic pair and returns it in
closed form, with the two isotope pairs giving the same `c` as a consistency
check. Exact recovery on 2e4 random triples. But the inversion is poorly
conditioned: with equal relative errors on the effects measured at the two
viscosities it amplifies relative error 3.8-fold at `c = 10` and 76-fold at
`c = 1`, demanding the isotope effects to 0.55%-0.03% against roughly 1% in
practice, and it is most demanding exactly where masking is heaviest
(`completion.viscosity_conditioning`).

**What a closed offset would bear on: the promoting mode.** In the ground-channel model
`F` depends on the gating width `w = 2 kappa sqrt(mu_H) sigma^2(T)` AND on the
scale parameter `A`, and for a harmonic promoting mode the exact thermal
variance is

```
w(T) = w0 coth(hbar omega_g / 2 k_B T)
```

a quantum expression in which zero-point motion holds `w` fixed below
`theta = hbar omega_g / k_B`. The offset inherits that structure. Across
250-350 K the predicted span of `F` at `A = 5`, `w0 = 0.5` is 0.07|F0| at
50 cm-1, **0.85|F0| at 200**, **1.00|F0| at 400**, and 0.35|F0| at 800, while
the intrinsic Swain-Schaad exponent moves 0.009-0.035 over the same interval.
These use the reduced masses in the tunneling model as well as in `gamma_SC`;
until v1.8.6 `qtunnel.py` used integer masses 1:2:3, which made the spans and
the Fisher example below too large. Sensitivity peaks where the mode sits at the
quantum-classical crossover `hbar omega_g ~ 2 k_B T`, which at 300 K is
417 cm-1; sweeping `A` over [1,50] and `w0` over [0.1,2] moves that peak only
between 140 and 470 cm-1, so it tracks the crossover rather than the
parameterization. The promoting-mode frequencies invoked in this literature lie
mostly below the 417 cm-1 crossover, in the lower part of that window or beneath
it (50 to 300 cm-1).
(`analysis/completion.py`)

Sensitivity is not identifiability, and only the first is claimed. `F` depends on
`A`, `w0` and `theta` jointly and the three are badly confounded: for an
11-point series over 250-350 K the Jacobian in `(ln A, ln w0, ln theta)` has
singular values 0.40, 1.8e-2, 1.2e-4, so `cond(J'J) ~ 1.0e7`. Below `theta` the
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
    corpus.py             the dataset: temperature series and single-condition
                        records (244 pairs in 76 analysis units)
  cofactor_benchmark.py the NAD+ titration of Tsai and Klinman (Table 5): test
                        of the minimal scheme, random-order check, 247-pair
                        ledger, counts with the three extra conditions as units
    reference_asymmetry.py  the tracer ratio r (the "reference asymmetry" of the
                        code and of some passages below): exactly 1 where one
                        tritiated isotopologue served both experiments
                        (data/tracer_design.csv), otherwise inferred from the
                        sources' secondary effects
  manuscript_root.py    locates a manuscript for the text cross-checks
                        (KIE_MANUSCRIPT, or an untracked manuscript_root.local);
                        without one those checks are skipped
  holdout.py            leave-one-temperature-out stability; Theorem 1
                        consistency check (algebraic, cannot fail)
  yadh_robustness.py    the yeast result under mass convention, error
                        inflation, pooling, reversibility and multiplicity
  audit_numbers.py      sweeps every derived number against both manuscript
                        files; catches fabricated and stale values that a
                        whitelist check cannot
  network_geometry.py   Moebius structure of V/K; curvature class and the
                        identified-set direction from network topology
  completion.py         closing the half-line with one commitment measurement,
                        and what the closed offset measures about the gating mode
    curvature.py          why the allowed range is one-sided: masking is concave
                        through the origin in log-rate coordinates, mass scaling
                        is homogeneity; axioms, composition, direction, and the
                        ordered-tracer check of the primary protocol
  mixed_label.py        forward masking in the mixed-label secondary design is
                        ordered by the primary effect, so it only deflates
  network_atlas.py      symbolic enumeration and classification of the chain
                        mechanisms (194 with one to four bound states); the
                        full check is `python network_atlas.py 4`
  atlas_general.py      the same classification for branched schemes with one
                        to three bound states, each step labeled in turn
  atlas_general_census.py  those schemes in exact rational arithmetic at
                        sampled rate constants (764 schemes, 5538 labeled steps)
  joint_nuisance.py     binding isotope effects and a bypass acting together:
                        where the yeast exclusion survives, and a finite model
                        that reproduces the pair where it does not
  secondary.py          published secondary H/T and D/T effects against the
                        exponents printed in the same papers
  qtunnel.py            exact ground-channel overlap model (reduced masses,
                        from masses.py)
  ridge.py              ridge limit and the offset
  partial_id.py         envelope B(K), exact identified-set endpoint,
                        two-commitment exploration
  identifiable_set.py   the half-line result and the point-estimate bounds
  bounds_uncertainty.py one-sided confidence bounds with shared-reference
                        correlation bracketed; exact_bound gives the exact
                        rho = -1 quantile of the endpoint (single records,
                        purity scenarios, multiplicity)
  isotopic_purity.py    protium in the deuterated comparators: fractional-
                        conversion model by labeling design (exact for one
                        shared tracer, a lower limit for dideuterated
                        comparators), corrections applied draw by draw, each
                        published correction undone with its own constants,
                        and the six exclusions per unit and across all units
                        (data/deuterium_purity.csv; results/sm_table_exclusions.tex)
  repeat_determinations.py  independent determinations of one condition:
                        do they detect an offset discrepancy?
  model_domain.py       the model comparison when the action and gating width
                        are bounded: largest offset, and the units that then
                        cannot exclude the model (results/model_domain.csv)
  reference_counts.py   the analysis units against the semiclassical reference
                        itself: observed below it at every admissible
                        temperature, above it, above at 95%, and above after
                        adjustment across units (results/reference_counts.csv)
  export_fig_example.py the worked example: the intrinsic exponent of two
                        published pairs at each assumed commitment (figure
                        data and table rows)
  synthetic_validation.py  the procedure run forward and back on known
                        truths: false exclusions, power, coverage, and what
                        a commitment measurement recovers
  equation_checks.py    every displayed identity of the article and the
                        supplement re-derived with sympy (48 checks)
  reversible.py         reversible-scheme identified set and the vacuity
                        threshold E_D* = (K_HT/K_DT)^(1/(gamma-1))
  vibronic_envelope.py  summed envelope B_vib; exact gated Franck-Condon
                        factors via tilted-Gaussian moments
  design_power.py       precision, correlation and commitment requirements
  decisive_case.py      yeast ADH primary effects (Cha, Murray & Klinman 1989)
  adh_primary.py        alcohol dehydrogenase primary records; horse liver
                        correlation bracket (results/sm_table_ladhprim.tex)
    build_trinomial.py    parses primary SI tables -> series data file + audit
  offset_analysis.py    profiles the offset per series (continuation method)
    offset_summary.py     goodness-of-fit stratification of those profiles
  profile_edges.py      continues the profiles whose grid optimum sat at
                        the edge F = 1.5: none has a finite minimum
  pooling_check.py      within-series pooling, and the variant-pooling scale
                        symmetry with the identifiable horse liver refit
  vibronic.py           does the bound survive the full vibronic sum? and the
                        penalty table (results/sm_table_vibronic.tex)
  published_intrinsic.py  published intrinsic pairs solved with an exponent
                        (Stojkovic 2012): offsets fixed by construction
                        (results/sm_table_stoj.tex)
  verify_derivation.py  numerical + symbolic checks of every closed form
    audit_v3.py           adversarial checks of the main claims, each probed
                        off its easy case
  audit_model_limits.py limiting behavior of the model: which constraints cap
                        the exponent, the ridge limit, dimensions, mass
                        convention (results/audit_model_limits.txt)
  audit_envelope.py     the envelope of attainable offsets, the direction of
                        masking, the sampling error of one record
                        (results/audit_envelope.txt)
  audit_profiling.py    a fit at fixed intrinsic magnitude against a profiled
                        one, and a non-monotonic K_HT(T) example
                        (results/audit_profiling.txt)
    discriminate.py       two-model fits of a temperature series at fixed
                        intrinsic magnitude (read by audit_profiling.py; no
                        reported number depends on them)
  design.py             a power calculation for separating the two model
                        classes (read by export_figures.py; no reported
                        number depends on it)
  export_figs_v3.py     pgfplots tables for the figures
  export_fig_asym.py    identified-set figure data (76 analysis units, two panels)
  export_fig_overview.py   overview figure: the compatible set as the range of
                        the commitment curve, for two observations
  export_fig_envelope.py   envelope figure: the region the ground-channel model
                        fills (edges checked on a dense grid), and the
                        promoting-mode readout
  export_fig_design.py     design figure: joint bypass-binding region, the yeast
                        completion, the pooling symmetry, the error budget
    export_fig_atlas.py   scheme-class figure: one mechanism of each class at
                        stated rate constants
  export_fig_geom.py    geometry figure: reference ray, model ridge, the two
                        worked examples
    export_figures.py     pgfplots tables for fig1_separation and fig2_bounds
                        (figures not used in the article)
  export_script_index.py   index of the scripts the supplement names, from
                        script_index.csv (results/sm_table_scripts.tex)
  export_si_tables.py   supplementary tables (results/)
  export_bounds_table.py     per-series bounds table (results/)
  export_benchmark_table.py  complete 216-record temperature-series table
  export_profiles_table.py   per-series profile table
data/                   curated inputs; no script overwrites a shipped file
results/                generated reports and tables
figures/tikz/           figure sources and their data
external_data/          SOURCES.md, census_queries.md and fetch_sources.sh
                        (no PDFs redistributed)
```

## Reproduce

```bash
pip install -r requirements.txt

./external_data/fetch_sources.sh     # open-access supplements; see SOURCES.md
cd analysis

python masses.py             # the mass convention and its consequences
python verify_derivation.py  # must exit 0
python build_trinomial.py    # series data + transcription audit (700 values)
python ridge.py              # ridge limit, offset, finite-scale bound
python partial_id.py         # envelope + exact endpoint, both verified
python identifiable_set.py   # the half-line result, point bounds
python bounds_uncertainty.py # one-sided confidence bounds, rho in {-1,0,.5,.9}
python curvature.py          # the one-sidedness theorem + controls
python network_geometry.py   # topology -> curvature -> set direction
python completion.py         # closing the set; the gating-mode readout
python audit_envelope.py > ../results/audit_envelope.txt          # envelope, masking direction
python audit_model_limits.py > ../results/audit_model_limits.txt  # limits of the model
python audit_profiling.py > ../results/audit_profiling.txt        # fixed against profiled fit (~3 min)
python audit_v3.py           # adversarial checks, must exit 0
python vibronic.py           # the vibronic-sum penalty and its table
python adh_primary.py        # alcohol dehydrogenase records, correlation bracket
python published_intrinsic.py  # intrinsic pairs that cannot test the reference
python offset_analysis.py    # per-series profiles   (about 3 min per series, ~2 h)
python holdout.py            # deletion sensitivity: point and confidence (~1 min)
python design_power.py       # precision requirements; dataset statistics
python repeat_determinations.py  # independent determinations of one condition
python reference_counts.py   # the dataset against the semiclassical reference
python cofactor_benchmark.py # the NAD+ titration: scheme test, table, sensitivity
python reversible.py         # reversible step: windows, endpoints, counterexample
python network_atlas.py 4    # enumeration and classification of mechanisms (~3 min)
python model_domain.py       # the model comparison under bounded parameters (~40 s)
python isotopic_purity.py    # purity scenarios and the exclusions table (~1 min)
python synthetic_validation.py   # the procedure on synthetic data with a known truth
python equation_checks.py    # every displayed identity re-derived symbolically
python export_fig_asym.py    # data for the asymmetry figure (exact single-record bounds)
python offset_summary.py     # stratified interpretation
python export_figs_v3.py     # pgfplots tables
python export_fig_overview.py; python export_fig_envelope.py; python export_fig_design.py
python export_fig_example.py # the worked example
python export_script_index.py   # index of the scripts the supplement names
cd ../figures/tikz && ./build.sh
```

Every stochastic step is seeded, so a rerun reproduces the reported values
exactly. `data/` is inputs only: no script writes to a file it reads.

Verified from a clean checkout and a fresh virtualenv on 2026-08-18, against
newer libraries than the pins require (numpy 2.5.2, scipy 1.18.0, pandas 3.0.5):
all twenty scripts exit 0 and every regenerated file in `results/` matches,
except that `offset_profiles*.csv` move by 1e-5 in `chi2_min` and 1e-11 in
`p_fit`. Those are optimizer-convergence differences between scipy releases.
Every reported column -- `F_hat`, `F_lo`, `F_hi`, `width`, `open_interval`,
`excludes_semiclassical`, `excludes_ridge`, `adequate` -- is unchanged, as is
every verdict. Only the timings above are machine dependent; most steps finish
in seconds, and `offset_analysis.py` then took about fifty minutes (about two
hours with the series added since).

Re-verified for v1.8.0 on 2026-09-14 from a clean clone (Python 3.12.13, numpy
2.2.6, scipy 1.17.1, pandas 3.0.2, sympy 1.14.0), with the fetched supplements in
place: `verify_derivation.py`, `corpus.py`, `reference_asymmetry.py`,
`holdout.py`, `build_trinomial.py` (transcription audit 130 of 130 pairs),
`network_geometry.py`, `identifiable_set.py`, `bounds_uncertainty.py`,
`completion.py`, `decisive_case.py`, `yadh_robustness.py`, every `export_*`
script and `audit_numbers.py` exit 0, the six figures build, and no tracked file
changes. `offset_analysis.py` and `network_atlas.py`, whose computations this
release does not touch, were not rerun.

Re-verified for v1.8.1 on 2026-09-15 from a clean clone, same environment: the
same scripts plus `export_fig_atlas.py`, `export_figures.py` and
`export_profiles_table.py` exit 0 (`bounds_uncertainty.py` about two minutes),
the six figures build, `audit_numbers.py` passes, and no tracked file changes.
The changes in `network_atlas.py` are to comments and docstrings only, and
`offset_analysis.py` is untouched; neither was rerun.

Re-verified for v1.8.13 on 2026-11-20 from a clean clone, same environment: the
scripts this release changes, adds or depends on (`cofactor_benchmark.py`,
`bounds_uncertainty.py`, `completion.py`, `reversible.py`, `curvature.py`,
`network_geometry.py`, `vibronic.py`, `vibronic_envelope.py`, `secondary.py`,
`decisive_case.py`, `adh_primary.py`, `isotopic_purity.py`,
`yadh_robustness.py`, `pooling_check.py`, `reference_counts.py`,
`repeat_determinations.py`, `holdout.py`, `model_domain.py`,
`joint_nuisance.py`, `synthetic_validation.py`, `equation_checks.py`,
`verify_derivation.py`, `build_trinomial.py`, `identifiable_set.py`,
`partial_id.py`, `export_fig_asym.py`, `export_fig_overview.py`,
`export_fig_envelope.py`, `export_fig_design.py`, `export_fig_example.py`,
`export_fig_atlas.py`, `export_fig_geom.py`, `export_si_tables.py`,
`export_bounds_table.py`, `export_benchmark_table.py`,
`export_profiles_table.py`, `export_script_index.py` and `audit_v3.py`) exit 0;
`audit_model_limits.py`, `audit_envelope.py` and `audit_profiling.py`
reproduce their stored outputs exactly; `network_atlas.py 4` passes; the ten
figures build; `audit_numbers.py` passes with the manuscript and exits 0
without it; and no tracked file changes. The data and the inputs of
`offset_analysis.py` are unchanged since v1.8.11, so it was not rerun.

Re-verified for v1.8.12 on 2026-11-14 from a clean clone, same environment: the
scripts this release changes, adds or depends on (`partial_id.py`,
`reference_asymmetry.py`, `reference_counts.py`, `model_domain.py`,
`adh_primary.py`, `completion.py`, `export_fig_example.py`,
`export_fig_atlas.py`, `export_fig_overview.py`, `export_fig_envelope.py`,
`export_fig_design.py`, `export_fig_asym.py`, `equation_checks.py`,
`verify_derivation.py`, `curvature.py`) exit 0, `network_atlas.py 4` passes, the
figures build, `audit_numbers.py` passes with the manuscript and exits 0
without it, and no tracked file changes. The data and the inputs of
`offset_analysis.py` are unchanged since v1.8.11, so it was not rerun.

Re-verified for v1.8.10 on 2026-10-22 from a clean clone, same environment,
with the fetched supplements in place: the scripts this release changes or
depends on (`partial_id.py`, `isotopic_purity.py`, `equation_checks.py`,
`export_si_tables.py`, `export_fig_asym.py`, `synthetic_validation.py`,
`design_power.py`, `repeat_determinations.py`) exit 0, the figures build, the
audit passes, and no tracked file changes. The inputs of `offset_analysis.py`
and `network_atlas.py` are unchanged since v1.8.9, so they were not rerun.

Re-verified for v1.8.9 on 2026-10-22 from clean clones, same environment, with
the fetched supplements in place: 46 scripts exit 0 (the v1.8.8 list plus
`design_power.py`, `repeat_determinations.py`, `isotopic_purity.py` and
`synthetic_validation.py`), `verify_derivation.py` reports no failures, the
transcription audit finds all 648 values verbatim, the figures build, the audit
passes, and no tracked file changes. The census changes the inputs of
`offset_analysis.py`, so it was rerun from a clean clone (about 131 minutes): it
exits 0 and reproduces its tracked outputs exactly. The scripts changed after
that run (`equation_checks.py`, `curvature.py`, `completion.py`,
`atlas_general.py`, `export_si_tables.py`) were rerun from a second clean clone
of the release commit, with the same result. `network_atlas.py`, whose
computation this release does not touch, was not rerun.

Re-verified for v1.8.8 on 2026-10-16 from a clean clone, same environment, with
the fetched supplements in place: 42 scripts exit 0 (the v1.8.7 list plus
`adh_primary.py`, `published_intrinsic.py`, `export_fig_overview.py`,
`export_fig_envelope.py`, `export_fig_design.py`, `joint_nuisance.py` and
`reversible.py`), `verify_derivation.py` reports no failures in 16 sections, the
figures build, the audit passes, and no tracked file changes. The light ecDHFR
series changes the inputs of `offset_analysis.py`, so it was rerun from the clean
clone (about 55 minutes): it exits 0 and reproduces its tracked outputs exactly.
`network_atlas.py`, whose computation this release does not touch, was not
rerun.

Re-verified for v1.8.7 on 2026-10-13 from a clean clone, same environment, with
the fetched supplements in place: the same 35 scripts exit 0, the six figures
build, `verify_derivation.py` reports no failures in 16 sections, the audit
passes, and no tracked file changes. Because the admitted methanol records
change its inputs, `offset_analysis.py` was rerun from the clean clone as well
(about 56 minutes): it exits 0 and reproduces its tracked outputs exactly.
`network_atlas.py`, whose computation this release does not touch, was not
rerun.

Re-verified for v1.8.6 on 2026-10-13 from a clean clone, same environment, with
the fetched supplements in place: 35 scripts exit 0, namely every script in the
reproduce list except `offset_analysis.py`, every `export_*` script, and
`corpus.py`, `reference_asymmetry.py`, `holdout.py`, `decisive_case.py`,
`yadh_robustness.py`, `mixed_label.py`, `pooling_check.py`, `profile_edges.py`,
`vibronic_envelope.py`, `discriminate.py`, `design.py` and `audit_numbers.py`.
`verify_derivation.py` reports no failures in 16 sections, the six figures
build, the audit passes, and no tracked file changes. `offset_analysis.py` and
`network_atlas.py`, whose computations this release does not touch, were not
rerun.

Re-verified for v1.8.5 on 2026-09-25 from a clean clone, same environment, with
the fetched supplements in place: 31 scripts exit 0, namely every script in the
reproduce list except `offset_analysis.py`, every `export_*` script, and
`corpus.py`, `reference_asymmetry.py`, `holdout.py`, `decisive_case.py`,
`yadh_robustness.py`, `mixed_label.py` and `audit_numbers.py`.
`verify_derivation.py` reports no failures in 16 sections, the six figures
build, the audit passes, and no tracked file changes. This release changes the
title and two text checks in `audit_numbers.py` only; `offset_analysis.py` and
`network_atlas.py`, whose computations it does not touch, were not rerun.

Re-verified for v1.8.4 on 2026-09-16 from a clean clone, same environment: every
script in the reproduce list exits 0, `verify_derivation.py` reports no failures
in 16 sections, the six figures build, `audit_numbers.py` passes, no tracked file
changes, and all four generated table inputs match the copies used in the
manuscript.

Re-verified for v1.8.3 on 2026-09-15 from a clean clone, same environment: the
same scripts exit 0, `verify_derivation.py` reports no failures in 16 sections,
the six figures build, `audit_numbers.py` passes, and no tracked file changes.
The asymmetry figure the clone builds is identical, page for page and character
for character, to the one used in the manuscript.

Re-verified for v1.8.2 on 2026-09-15 from a clean clone, same environment: the
v1.8.1 list plus `joint_nuisance.py` and `mixed_label.py` exit 0,
`verify_derivation.py` reports no failures in 16 sections, the transcription
audit finds 130 of 130 pairs, the six figures build, `audit_numbers.py` passes,
and no tracked file changes. `offset_analysis.py` and `network_atlas.py` are
untouched and were not rerun.

One caveat on timing: `offset_analysis.py` computes profile likelihoods by
continuation over a 301-point grid for each of the 46 series with at least
three temperatures, about three minutes per series on one core (roughly 55
minutes for the 20 series before the census of 2026-10-20).

## Key numbers

Reduced masses of the X-H oscillator, CODATA/AME2020, carbon donor. This is the
convention of Kohen & Jensen, who obtain 3.34 from the reduced mass of 12C and
the hydrogen isotopes; the bare-mass value 3.2628 is also in circulation.

```
gamma_SC = 3.34887   rigid = 2.45884   F0 = -0.042086   threshold = 1.26890
```

Temperature-series data: 216 temperature-resolved records, 48 series, -50 to 45 C.
Transcription audit: 700 values across the twenty-two PDF supplements (six parsed, sixteen
transcribed), each matched with its uncertainty, 0 misses. The wider matched-primary record used for the
identified-set figure comprises 244 matched primary records in 76 analysis units.
The 247 unique observations in `results/matched_measurements.csv` include the
three additional cofactor concentrations. `results/nad_condition_bounds.csv`
gives all four titration rows, the row that also belongs to a temperature
series, and the test of the minimal scheme across them.
`results/all_condition_units.csv` gives the bounds for each of the 79 units
when the three are counted. `results/condition_screen_sensitivity.csv`
summarizes the counts with and without them, and with and without the
heterogeneous commercial *Candida boidinii* FDH preparation.

Bounds: 3 of 48 series exclude both mechanisms, at every correlation -1, 0, 0.5
and 0.9 (horse liver ADH F93W in 40% and 60% methanol, yeast ADH with
p-chlorobenzyl alcohol); of the other 45 only the Y209W thymidylate synthase
hydride transfer clears F0 at any correlation, at rho = 0.9. A series bound is on the largest endpoint over its temperatures: the
per-temperature one-sided bounds at level 0.05/n_T, then their maximum, which
covers that maximum with probability at least 95% by the union bound whenever
each per-temperature bound has its nominal coverage. Each is the 0.05/n_T
quantile of the endpoint over lognormal draws with median at the measurement and
log-scale standard deviations `s/K` taken at the observed effects. Because the
endpoint is monotone in each effect, that would be exact at rho = -1 if the
log-scale standard deviations were known; estimated at the observation it is a
plug-in bound with approximate coverage: in repeated sampling the nominal 5%
bound fails 4.8% of the time for the yeast errors and 7.6% for a hypothetical
`2.0 +- 0.2`, `1.1 +- 0.001` pair (`bounds_uncertainty.plugin_coverage`). Through
v1.8.0 it was the 5% quantile of the resampled maximum, which is not a valid
bound for a maximum: taking ecDHFR W133F as truth, it exceeded the true maximum
in 20% of 2e4 simulated trials, against 2.4% for the corrected bound
(`bounds_uncertainty.coverage_check`); every
corrected bound is lower, and no verdict changes. Among the series that do not
exclude, the best 95% lower bound is -0.059 (rho = -1) to -0.039 (rho = 0.9),
thymidylate synthase Y209W hydride transfer in every case, and the best point
estimate is -0.001 (W80M), above F0, with 95% lower bounds of -0.086 (rho = 0)
and -0.101 (rho = -1). Bounds at rho = -1 are closed-form quantiles. Those at
the other correlations are averaged over 5 replications of 1e6 draws; the
largest Monte Carlo s.d. over the 144 simulated cells is 1.1e-3 (the wild type
in 50% methanol, whose inadmissible record leaves few finite tail draws), and
at rho = -1 the simulated values agree with the closed form within 0.0004.

## Relation to recent work

Williams (*J. Phys. Chem. B* **129**, 3604, 2025) shows by numerical example
that the same apparent KIE can arise from many combinations of step isotope
effects and weighting factors, and proposes comparing computed effects with
observed ones rather than with derived intrinsic values "of possibly dubious
validity". That is the non-identification problem, illustrated. This work characterizes the identified set exactly instead, and
finds it one-sided, so a robust inference survives in one direction.

Smedarchina & Siebrand (*Chem. Phys. Lett.* **410**, 370, 2005) reached the
direction of the envelope result in 2005 by instanton methods on a 2D model of
vibrationally assisted tunneling: assistance drives the tritium-referenced
exponent *toward* the semiclassical value from below, opposite to the then-
conventional assumption. Their Eq. (16) gives `e2` in [2.3032, 3.25725], and the
upper limit is the integer-mass semiclassical exponent 3.25725 by construction
(their Eq. 17 lower limit coincides with their Eq. 1); `masses.py` returns the
same value.

Their result is also the sharpest illustration of why the offset is the
coordinate used for model tests. The exponent converges to the semiclassical
value along the ray, which in exponent coordinates looks like the model family
merging with the semiclassical reference. It does not merge: the offset is pinned at `F0 = -0.042086`
for every point on that ray. A ratio of two diverging logarithms cannot resolve
the bounded difference between them.

Siebrand & Smedarchina (*J. Phys. Chem. B* **108**, 4185, 2004) give the
gating-averaged overlap of Proposition S1 (their Eq. 7) and the quantum thermal
width (their Eq. 4).

Francis & Kohen (*Perspect. Sci.* **1**, 110, 2014) is the community reporting
standard. It specifies error propagation but does not call for the covariance
between H/T and D/T effects sharing a tritium reference, which `design_power.py`
shows is worth a factor of 2.2 in that variance component.

## Pitfalls the suite checks for

Several mistakes are easy to make in this problem. Each has a script that
tests it directly, and the outputs are kept in `results/`:

1. *Optimizing on a constraint boundary instead of over the constraint
   region.* A cap on the exponent that decreases as the feasible set grows is a
   contradiction: the supremum over the region is the semiclassical value for
   every upper bound on the gating ratio (`audit_model_limits.py`).
2. *Reading a supremum over a model family as a boundary in observable space.*
   It is an asymptotic slope, reached only as the isotope effects diverge
   (`audit_model_limits.py`).
3. *Fixing the intrinsic isotope effect at an assumed value and fitting only
   the commitment.* At an exponent of 2.90, profiling the intrinsic magnitude
   instead moves the chi-square from 842 to 6.25 (`audit_profiling.py`).
4. *Assuming that an observed effect is monotonic in temperature.* A
   semiclassical pair with an Arrhenius commitment can give a turning point
   (`audit_profiling.py`); the `L_H` criterion needs no such assumption.
5. *A flat bound `F <= F0` over the whole parameter space.* It is false:
   `dF/dA < 0` identically, so the supremum is 0, and at `A = 0.05, w -> 0` one
   finds `F = -0.0110 > F0`. A numerical check that begins at `K_HT = 1.5` sits
   entirely inside the region where the flat bound holds. The correct statement
   is the scale-dependent envelope `B(K)` (`audit_envelope.py`).
6. *Inferring from `F0 + Delta > 0` that two model classes overlap.* An
   inequality on an upper bound does not show it. Direct evaluation is
   required, and it shows that the summed model can reach `F = +0.005118` at
   the worked parameter values (`vibronic_envelope.py`).
7. *Reversing the commitment ratio.* Masking falls as `k_off/k_chem`
   increases, not as it decreases (`audit_envelope.py`).

Point bounds are never reported without their sampling uncertainty
(`bounds_uncertainty.py`).

## Licensing

Code MIT (`LICENSE`); curated data CC BY 4.0 with attribution requirements to
the primary measurement papers (`LICENSE-DATA.md`). No publisher PDF is
redistributed; see `external_data/SOURCES.md`.

## Citation

See `CITATION.cff`. Please also cite the primary measurement papers listed in
`external_data/SOURCES.md` for any use of the curated values.
