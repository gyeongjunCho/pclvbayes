# pcLVbayes Roadmap

## Vision

pcLVbayes focuses on reliable Bayesian inference of **microbial–microbial interactions** from longitudinal compositional data.

## Scientific principle

The objective of pcLVbayes is not to force every pair to yield an interaction estimate.

The objective is to identify interactions only when they are supported by the posterior.

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
- [ ] Audit temporary-file and CmdStan output lifecycle
- [ ] Commit MTIST benchmark infrastructure
- [x] Profile R preprocessing, K-fold, memory, serialization, and sampling —
  MTIST 361 baseline separates directly measurable R components from combined
  CmdStan/process boundaries and ranks measured optimization candidates.
- [x] Optimize measured posterior-summary bottleneck — sampler-only retry
  diagnostics now precede one retained scientific summary bundle; measured
  parent-R time fell without changing Bayesian evidence or downstream ELPD.
- [x] Optimize measured outer-worker orchestration — explicit immutable
  exports and pair-level load balancing reduced the representative two-worker
  MTIST profile without changing scientific signatures or one-worker runtime.
- [ ] Optimize other measured bottlenecks only
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

- [ ] Define `a_ij` consistently as the posterior coefficient of the lagged source pair-to-rest log-ratio predictor in the target pair-to-rest log-ratio rate model, conditional on the self predictor and canonical preprocessing.
- [ ] State prominently that `a_ij` is not generally identical to absolute joint-gLV coefficient `A[i,j]`.
- [ ] Distinguish absolute direct-gLV `A[i,j]`, state-dependent mechanistic contrast `C_ij(x)`, deterministic canonical transformed-data projection, and Bayesian posterior coefficient `a_ij`.
- [ ] Explain that the sign of `a_ij` describes a fitted pair-to-rest dynamic direction, not alone biological facilitation, inhibition, or absolute interaction strength.
- [ ] Preserve indeterminate-not-zero semantics.
- [ ] Explain that a nonzero transformed coefficient may occur when `A[i,j]` is zero.

The mechanistic contrast is

    C_ij(x) = A[i,j] - sum_{k in rest} w_k(x) A[k,j]

and is state dependent; it is not identical to the fitted constant posterior coefficient.

#### Terminology, metrics, and interpretation

- [ ] Audit user-facing documentation, examples, reports, plots, tables, and summaries for unqualified direct-interaction, facilitation, inhibition, or absolute-strength claims.
- [ ] Prefer qualified pair-to-rest terminology; preserve historical names only with definitions.
- [ ] Do not change the public API solely for terminology cleanup.
- [ ] Make transformed-estimand validation primary: posterior versus canonical transformed oracle = 6/6.
- [ ] Report absolute-A comparison separately: posterior versus absolute-A cross-estimand agreement = 3/6; do not call it primary pcLV recovery.
- [ ] Prefer explicit metrics `transformed_oracle_sign_agreement`, `absolute_A_sign_agreement`, `absolute_A_zero_to_nonzero`, and `state_dependent_contrast`.
- [ ] Preserve explicit denominators and the prespecified six-direction rule.

#### Interpretation and preprocessing robustness

- [ ] Add interpretation states separate from `diagnostic_class`: `pair_to_rest_direction_supported`, `preprocessing_sensitive`, `empirically_context_sensitive`, `interpretation_indeterminate`, and `insufficient_information`.
- [ ] Allow statistical convergence with scientific interpretation sensitivity; convergence alone must not create an unqualified ecological sign.
- [ ] Preserve missing and indeterminate values rather than converting them to zero.
- [ ] Compare canonical smoothed projections with unsmoothed finite-interval projections using QR/SVD rank and conditioning diagnostics.
- [ ] Mark identifiable disagreements as preprocessing-sensitive without selecting the variant that best matches absolute truth.
- [ ] Preserve the canonical posterior and document `species_7 -> species_4` as the first verified smoothing-sensitive case.
- [ ] Measure smoothing-sign sensitivity on additional datasets before any calibrated hard exclusion threshold.

#### Empirical sign-heterogeneity diagnostics

- [ ] Add lightweight diagnostics for per-subject projections, leave-one-subject-out projections, predefined time-window signs, canonical-versus-unsmoothed signs, rank-deficient subset reasons, and the fraction agreeing with the canonical posterior sign.
- [ ] Document that these diagnostics cannot recover a mechanistic state-dependent contrast without absolute-model information.
- [ ] Do not define an arbitrary hard gate before multi-dataset calibration.

#### Conservative output and release gate

- [ ] Distinguish posterior sign support from scientific sign robustness and expose posterior coefficient, PSP, LFSR, diagnostic and interpretation classes, preprocessing sensitivity, and heterogeneity where available.
- [ ] Do not give preprocessing-sensitive or interpretation-indeterminate directions unqualified facilitation/inhibition labels.
- [ ] Document the 44/90 state-dependent contrast, 16/90 absolute-zero induced-effect, and smoothing-sensitive focused findings.
- [ ] Keep the existing model, Stan code, priors, posterior thresholds, and public API unchanged unless a separately reviewed defect is found.

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

---

## v0.3 — Time infrastructure, OU identifiability, and acceleration

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

### Latent continuous-time trajectory model

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

### Context-dependent pair-to-rest coefficients

- [ ] Investigate low-dimensional state-dependent or regime-dependent `a_ij` coefficients, including early/late coefficients, predefined regimes, varying-coefficient splines, Gaussian-process varying coefficients, and low-dimensional neural varying coefficients.
- [ ] Report state-dependent sign changes, require sufficient information, preserve indeterminate results when context dependence is not identifiable, and protect pairwise scalability.

### Rest-composition conditioning and denominator robustness

- [ ] Evaluate low-dimensional rest-composition covariates (ILR/principal balances, latent community factors, environmental or total-load covariates) with shrinkage and collinearity monitoring; do not call the result an unconditional direct effect.
- [ ] Predefine defensible denominator/balance variants, measure denominator-sensitive directions, and never select a denominator by truth agreement. Preserve the canonical denominator until a new contract is formally adopted.

### Optional absolute-scale information and joint confirmation

- [ ] Support uncertain optional total-load or absolute-abundance information from qPCR, spike-ins, flow cytometry, biomass, or calibrated sequencing.
- [ ] Separate measured from prior-imputed scale and retain pair-to-rest inference when scale is unavailable; absolute direct-gLV claims require an explicitly scale-aware model with sufficient information.
- [ ] Preserve pcLV as scalable screening and add optional sparse joint confirmation for selected candidates, keeping estimands separate and never converting omitted candidates to confirmed zeros.

### Multi-dataset oracle validation

- [ ] Extend oracle audits across matrices, initial states, noise levels, sampling intervals, series lengths, dominance, extinction/entry scenarios, and denominators.
- [ ] Report transformed-oracle agreement, absolute-A agreement, state-dependent contrasts, smoothing/interval/noise sign changes, denominator sensitivity, posterior-oracle disagreement, and indeterminate rates with explicit denominators.
- [ ] Acceptance criteria must follow the declared estimand rather than force agreement with absolute `A`.

### Backward compatibility

- [ ] Retain the v0.2 pair-to-rest model as an explicitly supported model.
- [ ] Do not silently change the meaning of `a_ij`; version schemas when a distinct estimand is introduced and label relative- and absolute-scale coefficients separately.
- [ ] Preserve conservative failure and indeterminate semantics.

### Laplace screening

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

## v1.0 — Stable microbial interaction release

- Canonical microbial–microbial Core
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
