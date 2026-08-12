# pcLVbayes Roadmap

## Vision

**pcLVbayes** develops a reproducible Bayesian framework for estimating directed **pair-to-rest compositional Lotka–Volterra dynamics** from longitudinal microbiome data.

The package is designed for conservative interaction screening and inference under sparse, irregularly sampled compositional time-series data. It preserves uncertainty explicitly and distinguishes unsupported, indeterminate, failed, and non-significant results rather than converting them to scientific zero.

The fitted interaction coefficient is a **triplet-ALR / pair-to-rest estimand**. It should not automatically be interpreted as an absolute direct-gLV coefficient or as proof of biological causality.

---

## Scientific principles

- One canonical scientific model
- Deterministic compositional preprocessing
- Directed target–source–rest triplet-ALR representation
- Student-t observation model
- Irregular-time OU residual structure
- Conservative Bayesian inference and reporting
- Explicit interaction-level and residual-level diagnostics
- No silent scientific fallback
- `indeterminate != zero`
- Screening must never redefine the compositional rest universe
- Benchmark truth is used only for post-inference validation, never to guide inference
- Reproducible seeds, execution manifests, restart behavior, and benchmark provenance

Conceptual inference target:

```text
longitudinal compositional data
        -> pair-to-rest transformation
        -> Bayesian pcLV inference
        -> posterior interaction distribution p(a_ij | D)
```

Orientation convention:

```text
source j -> target i
```

---

## v0.2 — Core reduction and baseline hardening

### Completed

- [x] Canonical pcLV preprocessing and model path
- [x] Directed pair-to-rest enumeration
- [x] Student-t likelihood and irregular-time OU residuals
- [x] Multi-chain posterior diagnostics
- [x] Separate interaction and residual identifiability assessment
- [x] Failure-aware K-fold evaluation
- [x] Structured failure and unavailable-state handling
- [x] Deterministic direction-specific seeds
- [x] Removal of silent scientific fallback behavior
- [x] Centralized input validation
- [x] Unified sequential/parallel execution behavior
- [x] MTIST benchmark adapter
- [x] 3-species and 10-species benchmark validation
- [x] Public `fit_pclv_bayes()` API reduction
- [x] Explicit documentation that pcLV estimates a pair-to-rest compositional interaction rather than a direct absolute-gLV coefficient

### Current interpretation policy

pcLVbayes separates statistical convergence from scientific interpretability.

- Stable interaction posterior + unstable residual allocation may be reported as interaction-stable / residual-unstable.
- Residual instability must remain visible in diagnostics.
- If interaction magnitude or sign is unstable across chains or posterior regimes, classify the direction as indeterminate.
- Indeterminate directions are excluded from significant-interaction lists.
- Posterior sign probability, PSP/LFSR, and predictive summaries are interpreted only after the required diagnostic gates pass.

---

## v0.2.1 — Publication-scale benchmark and execution freeze

Primary goal: freeze a reproducible, resource-safe pcLVbayes baseline at publication-scale network size.

### Execution architecture

- [x] Use one unordered pair as the unit of outer parallel work
- [x] Let one worker own both directions and all pair-specific work
- [x] Run four CmdStan chains serially within each pair worker (`parallel_chains = 1`)
- [x] Remove nested process-level K-fold worker pools
- [x] Precompute reusable deterministic pair state before worker dispatch
- [x] Bound BLAS/OpenMP/numerical-library threading inside workers
- [x] Validate fresh-session execution to avoid stale installed namespaces

Current benchmark policy on a 16-logical-CPU host:

```text
15 outer pair workers
x 1 active CmdStan chain per worker
= bounded process-level utilization
```

### Reproducibility and restart

- [x] Immutable execution manifest
- [x] Atomic status and completion artifacts
- [x] Deterministic resume
- [x] Explicit interrupted-work reconciliation
- [x] No reuse of invalidated benchmark roots

### 100-species MTIST benchmark

Target:

```text
100 species
4,950 unordered pairs
9,900 directed fits
4 chains
2,000 warmup + 1,000 retained iterations per chain
```

- [x] Launch benchmark under the pair-owned scheduler
- [ ] Complete all 4,950 unordered pairs / 9,900 directed fits
- [ ] Record completed, failed, indeterminate, and unavailable directions
- [ ] Preserve posterior summaries, PSP/LFSR, diagnostics, K-fold outputs, seeds, and denominators
- [ ] Freeze the completed benchmark corpus and provenance
- [ ] Summarize accuracy, coverage, calibration, failure modes, and runtime

### Release gate

- [ ] Full package regression suite
- [ ] `R CMD check`
- [ ] Core-lock verification
- [ ] Benchmark coverage and failure documentation
- [ ] Resource-policy documentation
- [ ] Reproducible benchmark provenance
- [ ] Updated user-facing limitations and interpretation guidance

---

## v0.2.2 — Internal simplification

Primary goal: simplify internals **without changing the frozen scientific behavior**.

Hotfix (2026-08-12): corrected held-out-subject `sd_r0` integration in predictive scoring and prevented main-fit sampler/retry state from leaking into K-fold fits. Canonical triplet preprocessing and Stan outputs were aligned without changing the pcLV estimand.

Completed early:

- [x] Pair-owned outer execution
- [x] Removal of nested K-fold process scheduling
- [x] Vectorized reusable pair-state preprocessing
- [x] Shared compositional helper foundation
- [x] `cor_meta_resid()` migration onto shared compositional helpers
- [x] Regression tests for shared triplet-ALR semantics
- [x] Effective-sample-size and time-aggregation helper hardening

Planned:

- [x] Characterize all public API, result schemas, masks, seeds, and fixed fixtures
- [ ] Remove only proven-dead private code
- [x] Migrate remaining duplicated triplet construction onto the shared helper foundation
- [x] Prove numerical equivalence before deleting legacy implementations
- [ ] Consolidate checkpoint and lifecycle helpers only after dedicated tests exist
- [ ] Preserve all public scientific behavior unless a later version explicitly declares a model/API change

---

## v0.3.0 — Laplace staged inference

Primary goal: make large candidate sets computationally practical while preserving the canonical pcLV estimand.

Target workflow:

```text
all candidate pairs
    -> cor_meta_resid() screening
    -> Laplace approximation
    -> selective full NUTS
```

Required work:

- [ ] Freeze `cor_meta_resid()` as the low-cost candidate-screening stage
- [ ] Reuse the canonical observed-support and triplet-ALR preprocessing contract
- [ ] Implement Laplace approximation for the canonical pcLV model
- [ ] Validate Laplace coefficient estimates against full NUTS
- [ ] Compare sign probability, uncertainty intervals, and failure diagnostics
- [ ] Define conservative promotion criteria from screening -> Laplace
- [ ] Define conservative promotion criteria from Laplace -> full NUTS
- [ ] Preserve `not-run`, `failed`, `indeterminate`, and `unavailable` states explicitly
- [ ] Keep the full analysis universe fixed for rest construction regardless of screening
- [ ] Benchmark runtime and memory scaling across increasing taxa counts

Scientific contract:

```text
correlation screening selects candidate pairs;
it does not redefine the compositional denominator.
```

---

## v0.3.1 — Automatic staged pcLV pipeline

Primary goal: expose the staged inference system through a stable one-call workflow.

Target user flow:

```text
longitudinal microbiome data
    -> observed-support filtering
    -> cor_meta_resid()
    -> candidate selection
    -> Laplace pcLV
    -> selective full NUTS
    -> unified pcLVbayes result object
```

Required work:

- [ ] Add a high-level automatic pipeline wrapper
- [ ] Keep expert lower-level functions independently callable
- [ ] Preserve the original full-community rest definition throughout all stages
- [ ] Preserve stage-specific diagnostics and withholding reasons
- [ ] Propagate missingness and unavailable states without numerical placeholders
- [ ] Provide reproducible resource controls and deterministic seeds
- [ ] Add unified summaries, matrices, and network-ready outputs
- [ ] Freeze a stable machine-readable result schema

---

## v0.4.0 — Ratio-robust cross-kingdom pipeline

Primary target: bacteria–fungi interaction inference when the absolute bacterial:fungal biomass ratio is not known precisely.

Initial sensitivity design:

```text
within-kingdom bacterial RA + fungal RA
        |
        +--> assumed joint mixture 5:5  -> primary pcLV analysis
        |
        +--> assumed joint mixture 7:3  -> sensitivity analysis
        |
        `--> assumed joint mixture 3:7  -> sensitivity analysis
```

The objective is not to claim mathematical independence from the bacterial:fungal ratio. Instead, pcLVbayes will identify interactions that remain sufficiently stable across a predefined plausible range of kingdom-scale mixture assumptions.

Required work:

- [ ] Define principled scaling of separately closed bacterial and fungal relative-abundance tables into a joint compositional universe
- [ ] Use a prespecified primary kingdom mixture
- [ ] Re-evaluate reportable interactions across prespecified alternative mixtures
- [ ] Define ratio-robustness using posterior sign, effect magnitude/uncertainty, diagnostics, and reportability
- [ ] Flag ratio-sensitive interactions explicitly
- [ ] Preserve the full joint taxonomic universe for rest construction within each assumed mixture
- [ ] Build the cross-kingdom workflow on top of the v0.3.1 staged pipeline
- [ ] Validate on bacteria–fungi longitudinal data
- [ ] Prioritize biologically testable candidate interactions for downstream experimental validation

---

## v1.0 — Stable pcLVbayes release

Target: a stable, documented, reproducible Bayesian pair-to-rest microbial-dynamics package.

Expected v1.0 capabilities:

- conservative Bayesian pair-to-rest inference;
- sparse and irregular longitudinal sampling support;
- Student-t observation model and irregular-time OU residuals;
- explicit diagnostic and indeterminate semantics;
- automatic correlation -> Laplace -> selective-NUTS inference;
- resource-safe and restartable large-scale execution;
- complete provenance and deterministic seeds;
- stable public API and result schema;
- network-ready posterior summaries;
- optional ratio-robust cross-kingdom workflow;
- publication-scale synthetic and biological validation;
- clear interpretation limits for compositional interaction coefficients.

---

## Known limitations and future scientific questions

pcLVbayes estimates compositional pair-to-rest interaction dynamics from relative-abundance time series. Therefore:

- it does not directly recover absolute gLV coefficients from relative abundance alone;
- a statistically supported pcLV sign is not by itself proof of absolute biological facilitation or inhibition;
- observational interaction inference is not equivalent to causal experimental confirmation;
- high-priority candidates should ideally be followed by microbial isolation, controlled co-culture, perturbation, or host-associated validation where feasible;
- improved procedures for interpreting the relationship between compositional interaction signs and underlying absolute biological interactions remain a separate future methodological question.

These limitations do not change the pcLVbayes inferential target: the package is intended to provide a reproducible posterior estimate of the supported triplet-ALR / pair-to-rest interaction from the observed longitudinal compositional data.

---

## Explicit exclusions

The pcLVbayes roadmap does **not** currently include:

- recovery of a full absolute interaction matrix from relative abundance alone;
- automatic reinterpretation or flipping of posterior interaction signs;
- simulation-trained machine-learning correction inside the Bayesian Core;
- full joint 100-species NUTS fitting;
- treating unavailable or indeterminate values as zero;
- choosing preprocessing rules according to agreement with simulation truth;
- claiming observational posterior interactions as causal proof.
