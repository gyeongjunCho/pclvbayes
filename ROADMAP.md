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
- [ ] Commit MTIST benchmark infrastructure — representative 10-species benchmark, targeted four-chain confirmation, and oracle estimand audit are complete but remain uncommitted.
- [x] Profile R preprocessing, K-fold, memory, serialization, and sampling —
  MTIST 361 baseline separates directly measurable R components from combined
  CmdStan/process boundaries and ranks measured optimization candidates.
- [x] Optimize measured posterior-summary bottleneck — sampler-only retry
  diagnostics now precede one retained scientific summary bundle; measured
  parent-R time fell without changing Bayesian evidence or downstream ELPD.
- [x] Optimize measured outer-worker orchestration — explicit immutable
  exports and pair-level load balancing reduced the representative two-worker
  MTIST profile without changing scientific signatures or one-worker runtime.
- [ ] Optimize other measured bottlenecks only when profiling demonstrates a material release-relevant benefit.
- [x] Run representative 10-species MTIST benchmarks — dataset 37 completed
  all 45 pair tasks and 90 directed fits in smoke and multi-chain reference
  configurations with coverage-aware scoring and explicit indeterminacy.
- [ ] Documentation
- [ ] Reduce and freeze the public `fit_pclv_bayes()` API
- [ ] Release preparation

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

- [ ] Define `empirical_sign_reversal_susceptibility` on a 0–1 scale from
  predefined deterministic robustness comparisons.
- [ ] Report component-level `r_g` values, `worst_axis_sign_reversal_susceptibility`,
  `worst_axis`, `number_of_valid_axes`, and `number_of_valid_comparisons`.
- [ ] Keep posterior sign probability, PSP/LFSR, sampler diagnostics, and
  empirical sign-reversal susceptibility as separate outputs and concepts.
- [ ] Preserve missing, rank-deficient, and indeterminate comparisons rather
  than treating them as agreement or zero risk.
- [ ] Use only prospectively defined robustness analyses; do not search for
  transformations that best match absolute truth.
- [ ] Treat the v0.2 score as an interpretation warning and descriptive measure,
  not as a calibrated hard exclusion threshold.
- [ ] Do not call the score `absolute_sign_reversal_risk` or interpret it as a
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

## v0.2.1 — Frozen Core release validation

Validate the frozen Core without reopening removed model options.

- Posterior ν calibration
- Multi-chain diagnostic gating
- Chain-specific sign agreement
- OU and observation-scale identifiability reporting
- CV spline robustness
- Zero-heavy datasets
- Extinction and external-entry scenarios
- Irregular sampling
- Retry policy validation
- Synthetic parameter recovery
- Synthetic sign recovery
- MTIST validation
- Publication-grade sampling recommendations

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
distinct from the uncalibrated v0.2 empirical susceptibility score.

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

## v0.3 — Time infrastructure, OU identifiability, and acceleration

Core v0.3 work is the release-bound scope below. Exploratory research tracks are not all mandatory for v0.3 release and require separate promotion decisions.

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

A successful v0.3 change should reduce the number of residual-indeterminate
directions without destabilizing directions whose interaction coefficient and
PSP/LFSR are already identified.

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

### Optional absolute-scale information and joint confirmation

- [ ] Support uncertain optional total-load or absolute-abundance information from qPCR, spike-ins, flow cytometry, biomass, or calibrated sequencing.
- [ ] Separate measured from prior-imputed scale and retain pair-to-rest inference when scale is unavailable; absolute direct-gLV claims require an explicitly scale-aware model with sufficient information.
- [ ] Preserve pcLV as scalable screening and add optional sparse joint confirmation for selected candidates, keeping estimands separate and never converting omitted candidates to confirmed zeros.

### Core validation and calibrated screening

### Multi-dataset oracle validation

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
- Full unrestricted 100+ taxa NUTS inference

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
