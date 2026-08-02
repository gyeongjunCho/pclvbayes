# pcLVbayes Roadmap

## Vision

pcLVbayes focuses on reliable Bayesian inference of directed pair-to-rest log-ratio dynamics from longitudinal compositional data. These coefficients support interaction screening but are not generally absolute direct-gLV coefficients.

## Scientific principle

The objective of pcLVbayes is not to force every pair to yield an interaction estimate.

The objective is to identify directed pair-to-rest dynamics only when they are supported by the posterior, without presenting them as absolute direct-gLV interactions.

Directions that remain non-identifiable after convergence diagnostics are reported as **indeterminate** rather than being interpreted as zero interaction or forced toward an arbitrary estimate.

Scientific conclusions should become progressively more complete as the amount,
quality, and diversity of longitudinal information increase.

Small datasets may support only a subset of interactions, whereas larger and
more informative datasets should identify additional interactions while
reducing uncertainty, rather than changing interactions that are already well
identified.

## Core principles

- One canonical scientific model
- Deterministic preprocessing
- Irregular-time OU residuals
- Fail-fast validation
- Structured runtime failures
- No silent scientific fallback
- Conservative Bayesian summarization
- Interaction identifiability is evaluated separately from residual-parameter identifiability
- Diagnostic failure is not interpreted as evidence of zero interaction
- Reproducible benchmarking

---

## v0.2 — Core reduction and release preparation

### Completed

- [x] Characterization tests
- [x] Canonical Core settings
- [x] Canonical preprocessing
- [x] Training-only global scaling
- [x] Predictor-variation validation
- [x] Failure-aware K-fold evaluation
- [x] Posterior Student-t ν estimation
- [x] Structured failure schema
- [x] Silent fallback removal
- [x] Centralized input validation
- [x] Unified sequential/parallel execution
- [x] MTIST adapter
- [x] Three-species smoke benchmark
- [x] Reference benchmark
- [x] OU geometry investigation
- [x] Multi-chain residual-identifiability investigation
- [x] Decision to retain the canonical OU parameterization in v0.2

### OU diagnostic conclusion

- Increasing warmup, sampling iterations, and chain count did not resolve the
  residual-geometry problem.
- Reparameterizing bounded persistence as log-decay did not provide a robust
  publication-scale improvement.
- A total-residual-scale plus OU-variance-fraction parameterization degraded
  sampling and destabilized a previously successful direction.
- Multi-chain fitting revealed stable interaction signs in some directions
  despite diagnostically invalid residual posteriors.
- Other directions showed distinct observation-noise and OU-allocation regimes,
  multimodality, or nonstationarity.
- Stable pooled interaction signs alone are insufficient for declaring a
  successful fit.
- The canonical Stan model remains unchanged in v0.2.

### Identifiability interpretation policy

pcLVbayes distinguishes between interaction identifiability and residual-parameter identifiability.

- If `a_ij` is stable across chains while residual parameters remain poorly
  identified, classify the direction as **interaction-stable but
  residual-unstable**.
- Such a direction is not automatically accepted as a completed fit; its
  residual diagnostic failure must remain visible.
- If residual regimes change the magnitude or sign of `a_ij`, classify the
  direction as **interaction-unstable / indeterminate**.
- Indeterminate directions are excluded from significant interaction lists.
- Diagnostic failure is not interpreted as zero interaction or evidence of
  no interaction.
- PSP/LFSR is interpreted only when interaction-level identification and the
  required diagnostic gates are satisfied.
- v0.3 should aim to reduce residual non-identifiability while preserving
  interactions whose `a_ij` and PSP are already stable across chains.

### Remaining

- [x] Strengthen multi-chain diagnostic reporting — deterministic interaction
  and residual identifiability classes now gate interpretation separately.
- [x] Record chain-specific sign agreement and residual-allocation disagreement
  — indeterminate directions remain explicit and are never treated as zero.
- [x] Audit temporary-file and CmdStan output lifecycle — temporary roots, workers, CmdStan processes, and executable reuse were verified.
Maintainer commit of the completed MTIST benchmark infrastructure is an external release action, not a v0.2.0 functional task.
- [x] Profile R preprocessing, K-fold, memory, serialization, and sampling —
  MTIST 361 baseline separates directly measurable R components from combined
  CmdStan/process boundaries and ranks measured optimization candidates.
- [x] Optimize measured posterior-summary bottleneck — sampler-only retry
  diagnostics now precede one retained scientific summary bundle; measured
  parent-R time fell without changing Bayesian evidence or downstream ELPD.
- [x] Optimize measured outer-worker orchestration — explicit immutable
  exports and pair-level load balancing reduced the representative two-worker
  MTIST profile without changing scientific signatures or one-worker runtime.
Further optimization is an optional post-baseline backlog item and is undertaken only when profiling demonstrates a material release-relevant benefit.
- [x] Run representative 10-species MTIST benchmarks — dataset 37 completed
  all 45 pair tasks and 90 directed fits in smoke and multi-chain reference
  configurations with coverage-aware scoring and explicit indeterminacy.
- [x] Documentation
- [x] Reduce and freeze the public `fit_pclv_bayes()` API
Maintainer-managed packaging, tagging, and source release operations are outside the v0.2.0 functional completion list.

### Scientific estimand hardening and conservative sign interpretation

The targeted four-chain confirmation remained converged for 6/6 prespecified directions, and posterior signs matched the observed canonical deterministic transformed-data projection for 6/6. Posterior signs matched absolute MTIST `A[target, source]` for 3/6. This is an estimand distinction, not a software defect or recovery of the absolute direct-gLV coefficient: 44/90 instantaneous mechanistic pair-to-rest contrasts changed sign across states, and 16/90 exact-zero absolute coefficients induced nonzero transformed effects. No focused discrepancy was caused by MCMC instability, finite-chain error, observation noise, finite-interval approximation, or a posterior-versus-observed-oracle discrepancy. `species_7 -> species_4` is the verified special case where absolute `A` was zero, the finite-interval transformed projection was positive, and canonical smoothing changed the projection and posterior sign to negative.

The v0.2 objective is to retain the canonical posterior model while preventing stable pair-to-rest coefficients from being misrepresented as context-independent absolute biological interactions, separating statistical convergence from scientific sign robustness, and completing this bounded remediation before release.

#### Estimand contract

- [x] Define `a_ij` consistently as the posterior coefficient of the lagged source pair-to-rest log-ratio predictor in the target pair-to-rest log-ratio rate model, conditional on the self predictor and canonical preprocessing.
- [x] State prominently that `a_ij` is not generally identical to absolute joint-gLV coefficient `A[i,j]`.
- [x] Distinguish absolute direct-gLV `A[i,j]`, state-dependent mechanistic contrast `C_ij(x)`, deterministic canonical transformed-data projection, and Bayesian posterior coefficient `a_ij`.
- [x] Explain that the sign of `a_ij` describes a fitted pair-to-rest dynamic direction, not alone biological facilitation, inhibition, or absolute interaction strength.
- [x] Preserve indeterminate-not-zero semantics.
- [x] Explain that a nonzero transformed coefficient may occur when `A[i,j]` is zero.

The mechanistic contrast is

    C_ij(x) = A[i,j] - sum_{k in rest} w_k(x) A[k,j]

and is state dependent; it is not identical to the fitted constant posterior coefficient.

#### Terminology, metrics, and interpretation

- [x] Audit user-facing documentation, examples, reports, plots, tables, and summaries for unqualified direct-interaction, facilitation, inhibition, or absolute-strength claims.
- [x] Prefer qualified pair-to-rest terminology; preserve historical names only with definitions.
- [ ] Do not change the public API solely for terminology cleanup.
- [x] Make transformed-estimand validation primary: posterior versus canonical transformed oracle = 6/6.
- [x] Report absolute-A comparison separately: posterior versus absolute-A cross-estimand agreement = 3/6; do not call it primary pcLV recovery.
- [x] Prefer explicit metrics `transformed_oracle_sign_agreement`, `absolute_A_sign_agreement`, `absolute_A_zero_to_nonzero`, and `state_dependent_contrast`.
- [x] Preserve explicit denominators and the prespecified six-direction rule.

#### Interpretation and preprocessing robustness

- [ ] Define candidate interpretation states separate from `diagnostic_class`: `pair_to_rest_direction_supported`, `preprocessing_sensitive`, `empirically_context_sensitive`, `interpretation_indeterminate`, and `insufficient_information`; implement them only after v0.2.1 validation.
- [ ] Allow statistical convergence with scientific interpretation sensitivity; convergence alone must not create an unqualified ecological sign.
- [ ] Preserve missing and indeterminate values rather than converting them to zero.
- [ ] Compare canonical smoothed projections with unsmoothed finite-interval projections using QR/SVD rank and conditioning diagnostics.
- [ ] Mark identifiable disagreements as preprocessing-sensitive without selecting the variant that best matches absolute truth.
- [ ] Preserve the canonical posterior and document `species_7 -> species_4` as the first verified smoothing-sensitive case.
- [ ] Measure smoothing-sign sensitivity on additional datasets before any calibrated hard exclusion threshold.

#### Empirical sign-heterogeneity diagnostics

- [x] Document benchmark-side deterministic robustness checks and their limits; production interpretation classification and calibrated gates are deferred to v0.2.1.


#### Empirical pair-to-rest sign-reversal susceptibility

Define a numerical measure of how readily the canonical posterior sign changes
under predefined, scientifically defensible deterministic robustness analyses.
`empirical_sign_reversal_susceptibility` is an empirical sign-sensitivity
measure calculable without knowing absolute `A[i,j]`; it is not the probability
that the absolute joint-gLV coefficient has the opposite sign and is separate
from posterior sign probability, PSP, and LFSR.

For each predefined robustness axis `g`, define

    r_g = (n_opposite,g + 0.5 * n_ambiguous,g) /
          (n_same,g + n_opposite,g + n_ambiguous,g)

where `same` and `opposite` are identifiable comparisons with the same or
opposite sign to the canonical posterior coefficient. `ambiguous` includes
near-zero, rank-deficient, or otherwise unsupported comparisons. Unavailable
comparisons remain missing and are never converted to zero risk. The initial
descriptive summary is

    empirical_sign_reversal_susceptibility = mean_g(r_g)

and reporting also requires

    worst_axis_sign_reversal_susceptibility = max_g(r_g)
    worst_axis
    number_of_valid_axes
    number_of_valid_comparisons

The arithmetic mean is a transparent initial descriptive summary, not a
calibrated probability or optimized weighted score; v0.2 introduces no
data-dependent weights.

Candidate axes are canonical smoothing versus unsmoothed finite-interval
projection, per-subject projection signs, leave-one-subject-out projection
signs, predefined early/late or time-window projections, and predefined
denominator or balance variants when implemented. v0.2 defines the contract and
uses benchmark-side comparisons where available; expanded subject, LOSO,
time-window, and denominator components remain v0.2.1 validation work unless
already implemented and validated.

- [x] Define `empirical_sign_reversal_susceptibility` on a 0–1 scale from
  predefined deterministic robustness comparisons.
- [x] Report component-level `r_g` values, `worst_axis_sign_reversal_susceptibility`,
  `worst_axis`, `number_of_valid_axes`, and `number_of_valid_comparisons`.
- [x] Keep posterior sign probability, PSP/LFSR, sampler diagnostics, and
  empirical sign-reversal susceptibility as separate outputs and concepts.
- [x] Preserve missing, rank-deficient, and indeterminate comparisons rather
  than treating them as agreement or zero risk.
- [x] Use only prospectively defined robustness analyses; do not search for
  transformations that best match absolute truth.
- [x] Treat the v0.2 score as an interpretation warning and descriptive measure,
  not as a calibrated hard exclusion threshold.
- [x] Do not call the score `absolute_sign_reversal_risk` or interpret it as a
  universal ecological probability.

#### Conservative output and release gate

- [ ] Distinguish posterior sign support from scientific sign robustness in documentation and benchmark reports; do not add a new public result field or calibrated hard gate in v0.2.
- [x] Document the 44/90 state-dependent contrast, 16/90 absolute-zero induced-effect, and smoothing-sensitive focused findings.
- [ ] Document the distinction among posterior sign support, empirical sign-reversal susceptibility, and absolute-gLV sign interpretation; state that low empirical susceptibility does not guarantee agreement with absolute `A[i,j]`.
- [x] Keep the existing model, Stan code, priors, posterior thresholds, public API, and result schema unchanged unless a separately reviewed defect is found.

Explicit v0.2 exclusions:

- recovering absolute `A[i,j]` from relative abundance alone;
- forcing posterior signs to agree with absolute MTIST truth;
- replacing the canonical spline model;
- latent biomass, nonlinear/state-dependent coefficients, or a full joint absolute interaction matrix.

### Final public API reduction

Perform only after benchmarking, profiling, and measured optimization are
complete.

- Classify every public argument as:
  - required scientific input;
  - justified user-facing configuration;
  - computational resource control;
  - internal implementation detail;
  - obsolete experimental option.
- Remove internal numerical constants, fallback controls, retry internals,
  Pathfinder tuning parameters, and mutually incompatible preprocessing
  options from the public signature.
- Do not hide obsolete options inside a large `control` object.
- Freeze one canonical zero-handling policy.
- Freeze one canonical smoothing policy.
- Preserve only reproducibility, sampling effort, K-fold, progress, and
  resource controls that users genuinely need.
- Add migration notes for every removed argument.
- Update documentation and examples.
- Require package tests and representative benchmark equivalence before
  release.

---

## v0.2.0 — Frozen functional baseline

v0.2.0 is the first frozen functional baseline. Packaging and source-tarball
operations are maintained separately and are not part of this functional
roadmap.

### Implemented functionality

- phyloseq-based longitudinal input with selected-taxa pairwise analysis;
- deterministic unordered-pair and directed source-to-target enumeration;
- pair-to-rest compositional modeling and directional posterior fitting;
- posterior interaction and residual diagnostics with separate Bayesian
  eligibility;
- significance evaluation, conditional K-fold execution, ELPD calculation,
  and stacking when common evidence is available;
- complete directed matrices with explicit masks for unavailable or
  non-reportable values;
- pair-specific posterior self-effect aggregation for diagonal output;
- bounded outer parallel execution, deterministic direction seeds, executable
  reuse, and explicit retention of failed or incomplete states;
- one canonical preprocessing path and the frozen public
  `fit_pclv_bayes()` interface.

The obsolete `fit_pclv_bayes2()` prototype, duplicate public controls, public
canonical preprocessing constants, retry controls, and Pathfinder tuning
controls were removed from the public surface. Those policies remain private;
they are not public configurability.

### Stability contract

Public stability includes the exact exported functions, the frozen
`fit_pclv_bayes()` formals, public result structure, source-to-target and matrix
orientation conventions, diagnostic/status vocabulary, unavailable-value
semantics, and zero-placeholder semantics. A zero in a complete external
matrix is a placeholder for an unavailable or non-reportable direction, not an
inferred scientific zero.

Computational stability includes deterministic task ordering and direction
indices, unique deterministic seeds, bounded workers, one canonical
preprocessing path, executable reuse, no worker-side model compilation,
explicit retry outcomes, and explicit completed, incomplete, failed, skipped,
and unavailable states. Unstable directions are never silently promoted and
indeterminate estimates are never converted to scientific zero.

Inferential stability includes separate interaction and residual assessment.
Residual instability does not automatically invalidate an internal interaction
coefficient, but it prevents promotion to the fully eligible class.
Interaction-indeterminate, sampler-failed, insignificant, residual-unstable,
and unavailable states remain distinct. Bayesian evidence is finalized before
K-fold and stacking; stacking does not rewrite finalized Bayesian evidence;
and opposite directions are assessed independently.

Only resource behavior actually exercised by the benchmark and lifecycle tests
is treated as established: worker and CmdStan cleanup, package-owned
successful temporary-root return, and compiled-executable reuse. No broader
resource guarantee is implied.

### Scientific philosophy

The model estimates compositional pair-to-rest effects and directional effects
under the current pairwise conditional formulation, together with posterior
uncertainty and diagnostic eligibility. It does not directly estimate an
absolute-abundance joint-gLV interaction matrix, causal ecological interactions
from observational data alone, a truth-equivalent MTIST absolute-gLV
coefficient, or an inferred zero whenever a value is unavailable, unstable, or
indeterminate.

The v0.2.0 principles are:

1. Keep uncertainty visible.
2. Keep failure, instability, indeterminacy, insignificance, and zero as
   different scientific states.
3. Do not report a coefficient merely because a finite posterior exists.
4. Assess directions independently.
5. Prefer conservative omission to confident reporting of an unreliable sign.
6. Never flip an estimated sign using benchmark truth.
7. Use external truth to evaluate diagnostics only; never use it in inference.
8. Report accuracy with coverage and explicit denominators.
9. Interpret MTIST results in light of the absolute-gLV versus compositional
   pair-to-rest estimand mismatch.
10. Document limitations rather than hiding them with placeholders or
    optimistic summaries.

### Validated 10-species MTIST baseline

The validated dataset contains 10 taxa, 10 independent series, 15 observations
per series, 45 unordered pairs, and 90 directed interactions. Orientation is
`source -> target = A[target, source]`.

Stage B classified 6 directions as `converged`, 61 as
`interaction_stable_residual_unstable`, 23 as `interaction_indeterminate`, and
0 as `sampler_diagnostics_failed`. Posterior records were retained for 90/90
directions; 6/90 were Bayesian-eligible; K-fold completed for 6/6 eligible
directions; finite ELPD was available for 6/6; stacking was available for 5/6;
and 6 directions were final-significant.

Conditional sign accuracy was 3/6. Significant-direction coverage was 6/90;
69/74 nonzero truth directions were omitted from significant output; 21/74
nonzero truth directions were interaction-indeterminate; and 1/16 truth-zero
directions was called nonzero.

This validates operational scalability to the representative benchmark, not
broad recovery of the MTIST absolute interaction matrix. Six reportable
directions are insufficient to calibrate a reliable sign-reversal diagnostic;
the larger benchmark in v0.2.1 is required to study when a statistically
significant sign should still be withheld.

### v0.2.0 functional completion gate

The v0.2.0 functional baseline is complete when the public API and canonical
behavior are frozen, supporting tests pass, the 10-species operational and
scientific baseline is documented, the estimand and limitations are explicit,
and no unresolved functional or stability defect is known to invalidate that
baseline. Maintainer-managed packaging operations are separate.

## v0.2.1 — Sign-reversal risk diagnosis and calibration

v0.2.1 develops and validates a conservative sign-reversal diagnostic for the
frozen v0.2.0 Core. The primary real-data-capable measure is
`empirical_sign_reversal_susceptibility`, an uncalibrated within-estimand
robustness score. The optional benchmark-calibrated layer reports
`absolute_same_sign_probability`, `absolute_sign_reversal_risk`,
`absolute_zero_probability`, and `direct_effect_misinterpretation_risk` only
when held-out validation and calibration-domain checks pass. These quantities
are not interchangeable, and the empirical score is never an absolute-A
probability.

Zero-heavy, extinction, irregular-time, synthetic-recovery, residual-model,
and other scientific-development studies belong to v0.3.0 unless explicitly
used as prespecified strata of the v0.2.1 calibration corpus. They are not
v0.2.1 feature work.

### Deferred interpretation validation

- [ ] Validate per-subject sign projections, leave-one-subject-out projections,
  predefined time-window sign sensitivity, and multi-dataset smoothing-sign
  calibration.
- [ ] Consider production `interpretation_class`, calibrated hard interpretation
  gates, and public significant-edge schema extensions only after that validation.
- [ ] Preserve indeterminate-not-zero semantics and require a separately
  reviewed promotion decision for any production schema change.


### Benchmark-calibrated absolute-sign interpretation risk

Use multiple truth-known MTIST datasets to estimate separate outcome
probabilities relative to the canonical pcLV posterior sign:

- same nonzero sign: `absolute_same_sign_probability`;
- opposite nonzero sign: `absolute_sign_reversal_risk` = P(absolute `A` has
  the opposite nonzero sign);
- absolute coefficient equals zero: `absolute_zero_probability`.

Also define

    direct_effect_misinterpretation_risk =
      1 - P(absolute A has the same nonzero sign)

`absolute_zero_probability` is not an opposite-sign event, while
`direct_effect_misinterpretation_risk` includes both opposite-sign and
absolute-zero possibilities. These are benchmark-calibrated probabilities,
distinct from the uncalibrated v0.2 empirical susceptibility score. Multi-dataset MTIST calibration belongs to v0.2.1, not v0.3.

Candidate calibration inputs include empirical susceptibility components and
worst-axis values, posterior sign probability and LFSR, diagnostic and
residual-identifiability status, preprocessing disagreement, subject and LOSO
sign agreement, time-window sign agreement, denominator sensitivity, predictor
correlation and design conditioning, numbers of subjects and time points,
sampling irregularity, zero prevalence, community dominance, and information
content. State-dependent contrasts and absolute-zero truth are calibration
outcomes or truth-known annotations, not production inputs.

- [ ] Construct a multi-dataset truth-known calibration corpus spanning
  interaction matrices, initial states, noise levels, sample sizes, time-point
  counts, interval structures, and community-dominance regimes.
- [ ] Require interaction-matrix-level train/test separation and held-out
  simulation regimes.
- [ ] Evaluate multiclass log loss, classwise Brier scores, calibration error,
  opposite-sign discrimination, and comparisons with class-frequency and
  simple-information baselines.
- [ ] Report uncertainty intervals and the calibration domain for every
  estimated probability.
- [ ] Estimate same-sign, opposite-sign, and absolute-zero probabilities as
  separate outcomes, with uncertainty intervals.
- [ ] Report the simulation and data domain over which calibration is supported.
- [ ] Evaluate calibration, discrimination, class imbalance, and held-out
  benchmark performance.
- [ ] Use prospectively fixed training and held-out MTIST datasets.
- [ ] Evaluate transport across interaction matrices and simulation regimes,
  not only random rows from one dataset.
- [ ] Do not present the result as a universal ecological probability outside
  the validated calibration domain.
- [ ] Do not add a production hard gate until thresholds are prospectively fixed
  and validated on held-out benchmark datasets.
- [ ] Preserve uncalibrated susceptibility when a real dataset lies outside the
  supported calibration domain.
- [ ] Return an explicit out-of-calibration-domain or insufficient-calibration
  state rather than fabricating a probability.

A direction must not be accepted solely because its pooled interaction sign is
stable. However, interaction stability and residual instability should be
reported separately rather than collapsed into one undifferentiated failure.
Diagnostic classification must also consider:

- R-hat
- Bulk and tail ESS
- Divergence count
- Maximum-treedepth saturation
- E-BFMI
- Chain-specific sign agreement
- Chain-specific residual-allocation regimes



## Sign interpretation reporting hierarchy

    statistical evidence
      - posterior sign
      - posterior sign probability
      - PSP/LFSR
      - diagnostic class

    empirical interpretation robustness
      - empirical_sign_reversal_susceptibility
      - worst_axis_sign_reversal_susceptibility
      - component-level robustness values
      - valid comparison counts

    benchmark-calibrated absolute-sign interpretation
      - absolute_same_sign_probability
      - absolute_sign_reversal_risk
      - absolute_zero_probability
      - direct_effect_misinterpretation_risk
      - uncertainty interval
      - calibration-domain status

The benchmark-calibrated third level is optional and may remain unavailable when held-out MTIST validation does not demonstrate sufficient predictive information.

High posterior sign probability means the fitted pair-to-rest coefficient has a
stable posterior sign. High empirical susceptibility means that sign is
sensitive to predefined analysis perturbations. High absolute-sign reversal
risk means that, within the validated MTIST calibration domain, the absolute
direct-gLV coefficient is estimated to have the opposite nonzero sign. High
direct-effect misinterpretation risk means the pcLV sign is unlikely to
represent the same nonzero absolute direct-gLV sign.

Version boundary: v0.2 defines and reports the empirical 0–1 susceptibility
without claiming an absolute-A probability; v0.2.1 calibrates the three
truth-known outcome probabilities with multiple MTIST datasets and held-out
validation; v0.3 investigates model changes that may reduce smoothing,
interval, denominator, omitted-community, and state-dependence distortions
while retaining this calibrated framework as a validation tool. Absolute-scale
modeling is not moved into v0.2 or v0.2.1.

Safeguards: stable posterior signs do not guarantee absolute direct-effect
signs; low empirical susceptibility does not prove equality with absolute `A`;
`A[i,j] = 0` is not an opposite-sign event; indeterminate and missing
comparisons are not zero risk; the v0.2 index is descriptive rather than a
validated probabilistic model; and no threshold is selected retrospectively to
maximize MTIST truth agreement.

---

## v0.3.0 — Scientific and statistical development

Core v0.3 work is the release-bound scope below. Exploratory pair-to-rest research tracks are not all mandatory for v0.3 release and require separate promotion decisions. External absolute-scale confirmation is not mandatory for v0.3 release.

### Core v0.3 work

### Shared time infrastructure

- Common time utilities for `fit_pclv_bayes()` and `cor_meta_resid()`
- Numeric time with explicit units
- `Date` and `POSIXct`
- Subject-local ordering
- Subject-local predecessor construction
- Continuous-time OU utilities
- OU half-life reporting

### Continuous-time OU

- Estimate one OU decay posterior from the complete retained data for each
  directed microbial interaction model.
- Apply the inferred decay to each subject-specific irregular interval using
  `rho_k = exp(-lambda * dt_k)`.
- Use a common decay rate within each directed fit rather than estimating an
  independent decay parameter for every interval.
- Preserve the current scientific OU process while improving time-unit handling,
  interval construction, and reporting.
- Keep time-axis modernization separate from sampler-coordinate changes.
- Do not assume that persistence-to-decay reparameterization alone improves
  posterior geometry.
- Do not silently switch to a different residual model.

### OU and observation-scale identifiability

Treat residual identifiability as a separate scientific problem.

Investigate:

- Observation-noise versus OU-residual allocation
- OU scale–decay ridges
- Student-t ν and observation-scale confounding
- Self-effect and microbial-partner collinearity
- Multiple residual-allocation regimes
- Multimodality and nonstationarity across chains
- Scientifically justified priors based on response scale and ecological time
  scale

Validation requires:

- Multiple chains
- Previously successful and failed directions
- Observation-dominant, balanced, and OU-dominant synthetic cases
- Stable interaction inference
- Acceptable R-hat, ESS, divergence, treedepth, and E-BFMI
- No prior-driven forced convergence

No new residual parameterization should replace the canonical model until it
performs consistently across real and synthetic cases.

A successful v0.3 change should reduce smoothing-, interval-, denominator-,
omitted-community-, and state-dependence distortions without destabilizing
directions whose pair-to-rest coefficient and PSP/LFSR are already identified.
It must not promise recovery of absolute `A[i,j]` from relative abundance alone.

### Latent and integrated interval dynamics

- [ ] Replace fixed external spline smoothing with a latent continuous-time trajectory and observation model.
- [ ] Propagate smoothing uncertainty into posterior inference.
- [ ] Reduce boundary-derivative artifacts and directly support irregular intervals while preserving subject resets, structured failures, and conservative identifiability semantics.
- [ ] Benchmark against exact noiseless transformed trajectories and the v0.2 smoothing-sensitive cases.
- [ ] Do not promise that this alone recovers absolute `A`.

### Integrated interval dynamics

- [ ] Distinguish instantaneous and finite-interval estimands.
- [ ] Evaluate interval-averaged or integrated source predictors.
- [ ] Model irregular `dt` directly; retain one decay posterior per directed fit and preserve `rho(dt) = exp(-lambda * dt)`.
- [ ] Quantify sign sensitivity to observation interval length.

### Exploratory v0.3 research tracks

### Context-dependent pair-to-rest coefficients

- [ ] Investigate low-dimensional state-dependent or regime-dependent `a_ij` coefficients, including early/late coefficients, predefined regimes, varying-coefficient splines, Gaussian-process varying coefficients, and low-dimensional neural varying coefficients.
- [ ] Report state-dependent sign changes, require sufficient information, preserve indeterminate results when context dependence is not identifiable, and protect pairwise scalability.

### Rest-composition conditioning and denominator robustness

- [ ] Evaluate low-dimensional rest-composition covariates (ILR/principal balances, latent community factors or observed total-load covariates) with shrinkage and collinearity monitoring; do not call the result an unconditional direct effect.
- [ ] Predefine defensible denominator/balance variants, measure denominator-sensitive directions, and never select a denominator by truth agreement. Preserve the canonical denominator until a new contract is formally adopted.

### External scale-aware confirmation workflow

- [ ] Document a handoff workflow from pcLVbayes pair-to-rest screening to an
  external absolute-scale joint dynamical model when sequencing counts and
  total-load qPCR, spike-in, flow-cytometry, biomass, or equivalent scale
  information are available.
- [ ] Use MDSINE2 as an example of an external scale-aware Bayesian microbial
  dynamics framework, without claiming guaranteed one-to-one API integration
  or edge-level confirmation.
- [ ] Keep pcLVbayes and external-model estimands, coefficient signs,
  uncertainty summaries, and identifiability claims explicitly separate.
- [ ] Evaluate whether posterior-supported pcLVbayes directions can be used as
  a candidate-screening or prioritization layer for external confirmation.
- [ ] Do not treat directions omitted by pcLVbayes screening as confirmed zero
  interactions.
- [ ] Do not require an external model to preserve the exact pcLV pairwise
  parameterization or coefficient meaning.
- [ ] Document the data requirements and interpretation boundary for the
  external confirmation workflow.
- [ ] Do not implement a duplicate absolute-scale joint-gLV engine inside
  pcLVbayes.

### Executable v0.3.0 scientific-development plan

v0.3.0 is the first version allowed to change model assumptions, priors,
residual model, diagnostics, estimand, preprocessing, eligibility thresholds,
significance policy, or interaction-recovery behavior. Research themes are
residual instability and indeterminacy, sign-reversal mechanisms,
estimand-aligned truth, compositional versus absolute relationships, prior and
regularization alternatives, residual-process alternatives, timing/interval,
series-count, observation-count, read-depth, noise, zeros/extinction,
coverage/precision trade-offs, and diagonal estimands.

Every task must define a scientific hypothesis, proposed change, estimand
impact, benchmark datasets, primary metrics, failure criteria, comparison with
frozen v0.2.x baselines, and any API/schema migration. No v0.3 implementation
starts before v0.2.2 is complete.

1. **V030-01 — Hypothesis and estimand register.** Approve hypotheses and
   baseline comparisons before coding. No implementation. Acceptance: each
   candidate has explicit failure criteria.
2. **V030-02 — Benchmark matrix.** Define compositional and absolute-scale
   truth relationships across timing, series count, observations, read depth,
   noise, zeros, and extinction. Deterministic simulation only after approval.
3. **V030-03 — Residual and indeterminacy research.** Test proposed residual,
   prior, and regularization alternatives against frozen baselines. Sampling
   permitted only under an approved benchmark plan.
4. **V030-04 — Estimand/preprocessing research.** Evaluate interval, smoothing,
   denominator, and diagonal alternatives with explicit estimand labels and
   no truth-driven sign flipping.
5. **V030-05 — Decision and migration.** Promote only changes meeting declared
   metrics; document rejected hypotheses, schema/API impact, and migration.



### Core validation and calibrated screening

### Multi-dataset oracle validation

This work may use the v0.2.1 calibrated risk framework to evaluate model improvements; it does not own the initial MTIST calibration task.

- [ ] Extend oracle audits across matrices, initial states, noise levels, sampling intervals, series lengths, dominance, extinction/entry scenarios, and denominators.
- [ ] Report transformed-oracle agreement, absolute-A agreement, state-dependent contrasts, smoothing/interval/noise sign changes, denominator sensitivity, posterior-oracle disagreement, and indeterminate rates with explicit denominators.
- [ ] Acceptance criteria must follow the declared estimand rather than force agreement with absolute `A`.

### Backward compatibility

- [ ] Retain the v0.2 pair-to-rest model as an explicitly supported model.
- [ ] Do not silently change the meaning of `a_ij`; version schemas when a distinct estimand is introduced and label relative- and absolute-scale coefficients separately.
- [ ] Preserve conservative failure and indeterminate semantics.

### Calibrated Laplace screening

Pipeline:

```text
Candidate generation
→ optional cor_meta_resid screening
→ Laplace approximation
→ reject or escalate
→ full NUTS + K-fold
→ conservative summarization
```

Requirements:

- Approximate the same canonical posterior as full NUTS
- Explicit screening decisions
- Failed Laplace approximations escalate to NUTS
- Never interpret approximation failure as evidence of no interaction
- Calibrate thresholds using synthetic data and MTIST
- Measure false-negative rate
- Preserve strong microbial interactions
- Reduce runtime without changing scientific conclusions

### Adaptive escalation (optional)

- Tier 0: reliable Laplace rejection
- Tier 1: short multi-chain diagnostic NUTS
- Tier 2: full NUTS + K-fold

---

## Post-v0.2 executable development plans

### A. v0.2.1 sign-reversal risk diagnostic and 100-species study

v0.2.1 aims to reduce confidently reported incorrect signs while preserving as
much useful coverage as possible. The diagnostic uses observed-data, posterior,
sampler, residual, K-fold, ELPD, and stacking information only. It never uses
MTIST truth for inference or feature generation, never flips coefficients, and
never overwrites original posterior, significance, or diagnostic outputs.

1. **V021-01 — Truth-isolation contract.** Define pair-level splits, inference
   versus truth-labeling boundaries, and the three outcome labels. Likely files:
   benchmark configuration and audit documentation. Add leakage tests. No
   sampling. Acceptance: truth is unavailable to inference and both directions
   of an unordered pair remain in one split. Commit boundary: contract only.
2. **V021-02 — Resource and scheduling policy.** Encode 12 logical threads with
   2 reserved, at most 10 active CmdStan chains, one thread per chain, fixed
   numerical-library/OpenMP threading. Safe outer concurrency is derived by
   preflight from observed chain behavior; ten outer workers are not assumed
   safe when each fit can launch four chains. Add preflight capacity tests. Sampling permitted only in
   preflight. Acceptance: observed peak never exceeds 10 chains. Commit
   boundary: scheduler/resource policy.
3. **V021-03 — Checkpoint and manifest architecture.** Add deterministic pair,
   direction, task, and chain seeds; atomic writes; restart; completed/failed/
   skipped/incomplete states; elapsed times; retry/Pathfinder records; and
   cleanup checks. Add interruption/restart fixture tests. No full benchmark.
   Acceptance: restart does not duplicate or alter completed results.
4. **V021-04 — Diagnostic-feature schema.** Define reproducible feature records
   for PSP/LFSR, posterior distance and spread, chain agreement, R-hat/ESS,
   divergences, treedepth/E-BFMI, interaction/residual classes, predictive
   completion, abundance/sparsity, and temporal sufficiency. Add schema and
   missingness tests. No truth-derived features.
5. **V021-05 — Four-chain preflight.** Run a bounded final-configuration
   preflight at 4 chains, 2000 warmup, 2000 sampling, validating scheduler,
   checkpointing, manifests, cleanup, and executable reuse. Acceptance: all
   states and resource measurements are retained; no K-fold truth leakage.
6. **V021-06 — Full 100-species pairwise coverage.** Run 4,950 unordered
   pairs and 9,900 directions (not a joint 100-species NUTS model),
   with 8,000 retained draws per direction under the global ceiling. Record all
   explicit states and denominators. Benchmark/sampling permitted. Acceptance:
   complete task manifest and reproducible restart.
7. **V021-07 — Post-inference truth labeling.** Join truth only after inference
   completion, keeping both directions in pair-level splits. No inference
   result may be changed. Acceptance: leakage audit passes and all labels have
   explicit denominators.
8. **V021-08 — Calibration and locked evaluation.** Freeze development,
   threshold-calibration, and held-out evaluation sets. Evaluate error rate,
   retained coverage, incorrect signs withheld, correct signs withheld,
   truth-zero behavior, calibration, and diagnostic-class strata. Acceptance:
   no threshold selected on locked data and no near-zero-coverage solution is
   approved.
9. **V021-09 — Conservative output integration.** Add independent concepts
   `conservative_sign_withholding_risk`, `sign_risk_class`, `sign_reportable`,
   `sign_withheld`, and `sign_withhold_reason` while retaining original output.
   Add schema/regression tests. Acceptance: original posterior and significance
   are byte/schema-equivalent and withholding is reversible and explicit.
10. **V021-10 — Documentation and regression release gate.** Document exact
    denominators, coverage/error trade-offs, resource ceilings, calibration
    domain, and limitations. Acceptance: all v0.2.1 completion gates below
    pass. Recommended commit boundary: documentation and tests only.

v0.2.1 completion gates: every planned task has an explicit state; incomplete
and failed directions remain denominators; truth leakage is excluded and
verified; development and locked evaluation are separated; features and
thresholds are reproducible and frozen before evaluation; original results
remain unchanged; no sign is flipped; withholding reasons are explicit;
error reduction, coverage loss, over-withholding, uncertainty, and denominators
are reported; checkpoint recovery and the 10-chain ceiling are validated.

### v0.2.2 — Behavior-preserving internal simplification

v0.2.2 begins only after v0.2.1 policy is frozen. It may remove dead private
code, obsolete compatibility paths, duplicate helpers, nested control flow,
private naming inconsistencies, and coupling among inference, diagnostics,
assembly, and benchmark code. It may centralize private constants and simplify
worker/checkpoint/cleanup logic. It may not change public formals/exports,
result schema or masks, diagnostic vocabulary, sign-risk policy, estimand,
preprocessing, priors, thresholds, Stan, eligibility, significance, K-fold,
ELPD, stacking, ordering, seeds, or unavailable/placeholder semantics.

1. **V022-01 — Characterization inventory.** Capture public names/classes/
   dimensions/masks, fixed preprocessing outputs, Stan data, task ordering,
   seed maps, statuses, and fixed-fixture outputs. Add tests; no sampling.
2. **V022-02 — Private dead-code cleanup.** Remove only code proven unused by
   call-graph and tests. Acceptance: characterization suite unchanged.
3. **V022-03 — Helper and control-flow consolidation.** Simplify private helpers
   and normalize errors/statuses in isolated commits. Acceptance: fixed-fixture
   outputs and failure semantics unchanged.
4. **V022-04 — Worker/checkpoint/resource simplification.** Refactor only after
   lifecycle and restart tests exist. Acceptance: ordering, seeds, manifests,
   ceilings, and cleanup unchanged.
5. **V022-05 — Architecture documentation and final regression.** Document
   concrete duplication/complexity reductions and run the complete suite.

v0.2.2 completion requires unchanged public API, exports, schema, masks,
seeds, generated Stan data, diagnostics, withholding classifications, and all
tests passing. Unfinished cleanup does not leak into v0.3.0.

### Transition gates

- **v0.2.0 functional completion:** frozen API and canonical behavior; passing
  characterization tests; documented 10-species operational/scientific
  baseline; explicit estimand and limitations; no known functional or
  stability defect invalidating the baseline.
- **v0.2.1 start:** v0.2.0 baseline and diagnostic semantics frozen, truth
  isolation approved, and the 100-species execution plan documented.
- **v0.2.2 start:** v0.2.1 diagnostic policy complete and frozen, with
  characterization targets identified and correctness separated from cleanup.
- **v0.3.0 start:** v0.2.2 cleanup complete, frozen v0.2.x baselines
  reproducible, and scientific hypotheses/success criteria approved before
  implementation.

## Post-v0.3 — CI and reproducible benchmark automation

Add automation only after the time-input contract, continuous-time OU
implementation, OU-identifiability policy, and Laplace-screening workflow are
stable.

Planned work:

- GitHub Actions package and adapter smoke tests
- Calibrated Laplace-screening regression tests
- Representative full-NUTS reference benchmarks
- Manual publication-scale benchmark workflows
- Durable benchmark history
- Release validation automation

---

## v0.4 — Optional cross-kingdom extension

Separate from the microbial Core.

Possible scope:

- Bacteria–fungi
- Bacteria–pathogen
- Cross-kingdom validation
- Colletotrichum suppression-candidate inference

---

## Future research — not currently planned

These are intentionally **outside the current roadmap**.

- Sequencing-count observation models
- Read-depth-aware latent composition
- Environmental predictors
- Joint microbial/environmental models
- Hierarchical OU
- Experiment-specific OU
- AI surrogate approximation
- Full joint 100+ taxa multivariate NUTS inference

---

## v1.0 — Stable pair-to-rest microbial dynamics release

- Canonical pair-to-rest microbial dynamics Core for qualified compositional interaction screening
- Unit-aware irregular-time continuous-time OU
- Multi-chain convergence and residual-identifiability diagnostics
- Structured directional and fold failures
- Conservative Bayesian summaries
- Optional calibrated Laplace screening
- Safe temporary-file and resource handling
- Reproducible benchmarks
- Durable benchmark history
- GitHub Actions
- Comprehensive documentation
