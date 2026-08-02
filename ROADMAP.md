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
- [ ] Run representative 10-species MTIST benchmarks
- [ ] Documentation
- [ ] Reduce and freeze the public `fit_pclv_bayes()` API
- [ ] Release preparation

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
