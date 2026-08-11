# pcLVbayes

> **Status: Beta testing (experimental)**
>
> pcLVbayes is under active scientific and software validation.
> APIs, defaults, preprocessing rules, diagnostic thresholds, scoring methods,
> and output schemas may change before the stable release.
>
> Results should not yet be treated as production-ready or used for biological
> conclusions without independent validation.

Bayesian **pairwise compositional Lotka–Volterra** regressions for amplicon
time-series data.

pcLVbayes includes:

- correlation-based pair screening;
- pair-to-rest additive log-ratio transformation;
- irregular-time Ornstein–Uhlenbeck residuals;
- Student-t observation models;
- posterior sign probabilities and LFSR;
- repeated subject-level K-fold evaluation;
- Student-t scale-mixture Kalman OU predictive scoring;
- diagnostics-aware sampling and retry policies;
- tidy pairwise and directed-edge outputs.

The package is built on **R**, **cmdstanr**, and **CmdStan**.

[![R](https://img.shields.io/badge/R-%3E%3D%204.2-276DC3.svg)](https://cran.r-project.org/)
[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
[![Status: beta testing](https://img.shields.io/badge/status-beta%20testing-yellow.svg)](#beta-testing-status)

---

## Beta-testing status

pcLVbayes is currently in **beta testing**.

The current version is intended for:

- scientific method development;
- simulation studies;
- benchmark evaluation;
- reproducibility testing;
- diagnostic validation;
- limited technical review.

The current version is not yet intended for:

- production pipelines;
- clinical or regulatory decisions;
- unattended large-scale deployment;
- definitive ecological interaction claims;
- direct interpretation as an absolute gLV interaction network.

During beta testing, previously generated results may become incompatible with
new versions when the scientific preprocessing or scoring contract changes.

In particular, changes to any of the following require regeneration of affected
results:

- compositional closure;
- zero replacement;
- spline smoothing;
- lag construction;
- time-unit handling;
- predictor standardization;
- posterior likelihood;
- predictive scoring;
- diagnostic classification;
- checkpoint or result schemas.

Users should record the exact package commit, CmdStan version, configuration,
random seeds, preprocessing schema, and output schema for every analysis.

---

## Licensing during beta testing

This beta version is provided for limited technical review.

Unless otherwise stated in `LICENSE`, the beta repository is not licensed for
redistribution, modification, public deployment, or production use.

A stable open-source release is planned under **GPL-3** after scientific and
software validation is complete.

© 2025–2026 Rural Development Administration (RDA)
© 2025–2026 National Institute of Agricultural Sciences (NAS)

---

## What pcLVbayes does

### Correlation screening — `cor_meta_resid()`

`cor_meta_resid()` computes within-subject pairwise associations using Pearson
or Spearman correlations and combines them across subjects using
meta-analytic estimators such as DL, with optional KNHA adjustment.

The screening stage can include effective-sample-size correction for temporal
autocorrelation and uneven sampling.

Correlation screening is intended to reduce the number of candidate pairs
before Bayesian fitting. It is not itself interpreted as a causal interaction
analysis.

### Bayesian pairwise pcLV Core — `fit_pclv_bayes()`

For each unordered taxon pair `{i, j}`, the Core fits two directed regressions:

- `j → i`;
- `i → j`.

Each direction uses:

- lag-1 target pair-to-rest ALR;
- lag-1 source pair-to-rest ALR;
- the target response `ΔALR_i / Δt`;
- irregular-time OU residual dependence;
- Student-t observation noise;
- one posterior-estimated degrees-of-freedom parameter `nu > 2`;
- posterior sign probabilities and LFSR;
- repeated subject-level K-fold predictive evaluation.

---

## Scientific estimand and interpretation of `a_ij`

pcLVbayes estimates a **directed pair-to-rest dynamic coefficient**.

For the fitted direction `j → i`, `a_ij` is the posterior coefficient of the
lagged source pair-to-rest log-ratio predictor in the target pair-to-rest
log-ratio rate model, conditional on:

- the target self predictor;
- the selected pair;
- the remainder of the observed community;
- the compositional preprocessing;
- the observed state distribution;
- the residual model;
- the prior and likelihood.

The fitted model can be written schematically as:

```text
ΔALR_i / Δt
    =
r0
+ a_ii × lag(ALR_i)
+ a_ij × lag(ALR_j)
+ residual process
```

where:

```text
ALR_i = log(x_i / x_rest)
ALR_j = log(x_j / x_rest)
```

and:

```text
x_rest = sum of all taxa other than i and j
```

after the declared smoothing, reclosure, and zero-handling procedure.

### What `a_ij` is

`a_ij` is:

- a fitted model-level compositional coefficient;
- a directed pair-to-rest dynamic association;
- conditional on the target self term;
- specific to the declared preprocessing and residual model.

### What `a_ij` is not

`a_ij` is not generally:

- the absolute joint-gLV coefficient `A[i,j]`;
- a context-independent direct ecological interaction;
- an absolute-abundance effect size;
- proof of facilitation or inhibition;
- proof of causality.

A positive posterior sign means that the fitted source pair-to-rest coordinate
has a positive conditional association with the target pair-to-rest rate under
the model.

A negative posterior sign means that the fitted association is negative under
the model.

Neither sign alone proves a biological mechanism.

---

## Relationship to absolute gLV coefficients

The following quantities are related but are not interchangeable:

1. the absolute gLV coefficient `A[i,j]`;
2. the state-dependent compositional contrast

   ```text
   C_ij(x) = A[i,j] - sum_{k in rest} w_k(x) A[k,j]
   ```

3. the deterministic coefficient obtained after transformation;
4. the Bayesian posterior coefficient `a_ij`.

The pair-to-rest coefficient may differ in sign from the absolute coefficient
because it represents the source effect relative to the aggregated remainder
of the community.

A stable posterior sign validates the fitted transformed coefficient. It does
not by itself validate the absolute direct-gLV sign.

A nonzero transformed coefficient may occur even when:

```text
A[i,j] = 0
```

An indeterminate posterior direction should not be interpreted as a confirmed
zero interaction.

Absolute direct-interaction claims require additional assumptions,
absolute-abundance measurements, perturbation data, or external mechanistic
evidence.

---

## Core compositional transformation

The Core uses zero-aware pair-to-rest ALR coordinates for the triplet:

```text
(i, j, rest)
```

where `rest` contains all taxa other than `i` and `j`.

### Smoothing and reclosure

Relative abundance trajectories are smoothed taxon-by-taxon on the log scale
within each subject.

Independent taxon-wise smoothing does not automatically preserve the
compositional unit-sum constraint. Therefore, the complete smoothed community
is reclosed sample-by-sample before any pair is extracted.

For each sample:

```text
x_k,reclosed
    =
x_k,smoothed / sum_l(x_l,smoothed)
```

After reclosure:

```text
sum_k(x_k,reclosed) = 1
```

Consequently:

```text
rest = 1 - x_i - x_j
```

is numerically equivalent to:

```text
rest = sum of all reclosed taxa other than i and j
```

### Zero-aware triplet construction

The `(i, j, rest)` triplet is processed using the configured zero-replacement
and rest-floor policy.

The triplet is then closed again before ALR calculation.

The resulting coordinates are:

```text
ALR_i = log(x_i / x_rest)
ALR_j = log(x_j / x_rest)
```

ALR magnitudes are protected by the canonical cap policy.

### Predictors and response

Predictors:

```text
xi = lag-1 ALR_i
xj = lag-1 ALR_j
```

Response:

```text
y = ΔALR_i / Δt
```

The response is not standardized.

Predictors are centered and scaled globally for the main posterior fit.

Within K-fold evaluation, predictor centering and scaling are estimated using
training subjects only.

---

## Time-axis policy

pcLVbayes supports irregularly spaced observations.

Use one time unit consistently across the complete analysis, such as:

- days;
- weeks;
- hours.

Decimal values are allowed, for example:

```text
3.5 days
1.25 weeks
```

`Date` or `POSIXct` values may be converted using `difftime()` before fitting.

Do not mix time units across subjects or samples.

Changing the time unit changes the numerical scale and interpretation of:

- `ΔALR_i / Δt`;
- `a_ij`;
- `a_ii`;
- the intercept;
- OU decay;
- OU persistence;
- residual scale parameters.

The selected unit should be reported explicitly in publications and analysis
records.

Example:

```text
All Δt values were expressed in weeks.
```

Within each subject, observation times must be finite and strictly increasing.

---

## Observation and residual model

The main posterior uses an irregular-time Ornstein–Uhlenbeck residual process.

A schematic state model is:

```text
e_t = rho_t × e_(t-1) + eta_t
```

where persistence depends on the observed time interval.

For a decay parameter `lambda`:

```text
rho_t = exp(-lambda × Δt)
```

The observation model uses Student-t noise:

```text
y_t ~ Student-t(nu, mu_t + e_t, sigma)
```

with:

```text
nu = 2 + exp(log_nu_minus_two)
```

Therefore:

```text
nu > 2
```

for every posterior draw.

The Student-t likelihood allows heavier tails than a Gaussian observation
model and reduces the influence of unusually large residual observations.

---

## Posterior sign support

For each directed coefficient, posterior draws are used to compute:

```text
P(a_ij > 0)
P(a_ij < 0)
```

The posterior sign probability can be summarized as:

```text
PSP = max(P(a_ij > 0), P(a_ij < 0))
```

The local false sign rate is:

```text
LFSR = min(P(a_ij > 0), P(a_ij < 0))
```

Equivalent two-sided posterior tail summaries may also be returned.

Posterior sign support should only be interpreted together with:

- R-hat;
- bulk ESS;
- tail ESS;
- divergences;
- maximum-treedepth hits;
- E-BFMI;
- chain-specific sign behavior;
- chain-specific residual behavior.

A high pooled sign probability is not sufficient when chains disagree or the
sampler diagnostics are poor.

---

## Predictive validation and ELPD

Repeated subject-level K-fold predictive evaluation is used to assess
out-of-sample predictive support.

Subjects, rather than individual observations, are assigned to folds.

This prevents observations from the same longitudinal subject from appearing
in both training and test data.

### Fold-level leakage guards

Each fold uses:

- training-subject-only predictor centering;
- training-subject-only predictor scaling;
- independent sampler adaptation;
- no reuse of the main-fit step size;
- no reuse of the main-fit inverse metric;
- no reuse of fold-incompatible initialization artifacts.

Failed fold evaluations remain explicit failures.

They do not contribute artificial zero-valued ELPD observations.

Model comparisons should use only common successfully evaluated evidence.

---

## Student-t scale-mixture Kalman OU scoring

The latent OU state transition is Gaussian and remains Kalman-based.

The observation model is Student-t and is not replaced by the former
variance-matched Gaussian approximation.

The Student-t distribution is represented as a Gaussian scale mixture:

```text
epsilon_t | omega_t
    ~ Normal(0, sigma^2 / omega_t)

omega_t
    ~ Gamma(nu / 2, nu / 2)
```

Conditional on the local scale `omega_t`, the observation update is Gaussian
and can be evaluated using Kalman recursion.

The local scale is integrated using deterministic quadrature.

The current beta scorer is identified as:

```text
student-t-scale-mixture-kalman-ou-q16
```

This is a deterministic numerical approximation to the Student-t predictive
integral. It should not be described as an exact closed-form Student-t Kalman
filter.

The scorer uses posterior draws of:

- regression coefficients;
- Student-t `nu`;
- observation scale;
- OU scale;
- OU persistence or decay;
- optional subject-level random-intercept scale.

Predictive log densities are aggregated over posterior draws using stable
log-mean-exp calculations.

---

## ELPD aggregation

Repeated K-fold evaluation produces subject-level out-of-fold predictive
scores.

Aggregation is **subject-uniform** by default so that subjects with more
observations do not automatically dominate the final summary.

The output may include:

- total subject-level ELPD;
- per-observation ELPD;
- successful evaluation counts;
- failed evaluation counts;
- fold diagnostics;
- fold-specific posterior `nu` summaries;
- repeated split manifests.

ELPD is supporting evidence for model comparison. It should not replace
posterior diagnostics or scientific interpretation.

---

## Sampling robustness

The Core uses a diagnostics-aware sampling policy.

The retry policy monitors:

- divergent transitions;
- maximum-treedepth saturation;
- E-BFMI;
- sampler execution failures.

Retries may adjust computational sampling settings according to the private
canonical policy, such as:

- `adapt_delta`;
- warmup length;
- metric type;
- step size;
- initialization.

Retries do not change:

- the scientific response;
- the pair-to-rest transformation;
- the Student-t likelihood;
- the prior on `nu`;
- the directed estimand.

All retry attempts and final sampling settings should be retained in the
output provenance.

A completed fit with poor diagnostics may still be returned for explicit
classification, but it should not automatically be treated as a converged
scientific result.

---

## Diagnostic interpretation

Directed fits may be classified into categories such as:

- `converged`;
- `interaction_indeterminate`;
- `interaction_stable_residual_unstable`;
- `sampler_diagnostics_failed`.

The exact classification contract may change during beta testing.

A coefficient should normally be treated as interpretable only when:

- chains agree on the interaction sign;
- interaction magnitudes are compatible across chains;
- R-hat is acceptable;
- ESS is sufficient;
- divergences are absent or below the declared policy;
- treedepth saturation is absent or below the declared policy;
- E-BFMI is acceptable.

Residual non-identifiability should be reported even when the interaction sign
appears stable.

---

## Correlation screening and ALR

For compositional correlation screening, pair-to-rest ALR is recommended for
the triplet:

```text
(i, j, rest)
```

The default ALR coordinates are:

```text
log(x_i / rest)
log(x_j / rest)
```

ILR was evaluated during method development but is not the default
correlation-screening transformation because it produced an undesirable
increase in negative associations in the evaluated setting.

This observation is empirical and dataset-dependent. It should not be treated
as a general claim that ILR is inappropriate for microbiome data.

---

## Parallelism

The main computational levels are:

1. outer pair-level workers;
2. CmdStan chains;
3. optional fold-level workers.

Approximate active chain concurrency is:

```text
n_workers_outer × chains
```

Nested parallelism should be avoided.

When outer pair-level parallelism is active, K-fold parallelism is normally
restricted to prevent oversubscription.

Recommended policy:

```text
n_workers_outer × chains
    <= available logical CPU threads
```

Additional concurrency may be possible on lightly loaded machines, but memory,
thermals, disk I/O, and process counts should be monitored.

BLAS, OpenMP, and other implicit thread pools should normally be restricted to
one thread per R or CmdStan worker.

---

## Progress reporting

When `progressr` is installed and:

```r
progress = "bar"
```

progress bars may be shown for pair-level and fold-level execution.

Otherwise, progress reporting falls back to quiet or text-based execution
according to the selected options.

Progress reporting does not alter the scientific model.

---

## Recommended interpretation workflow

For each directed edge:

1. verify the preprocessing and time-unit contract;
2. inspect MCMC diagnostics;
3. inspect chain-specific coefficient behavior;
4. inspect posterior sign probabilities and LFSR;
5. inspect credible intervals;
6. inspect residual identifiability;
7. use subject-level ELPD as supporting predictive evidence;
8. compare the transformed estimand with the intended biological question;
9. avoid direct absolute-gLV claims without additional information.

The primary reported interpretation should normally be:

```text
posterior-supported pair-to-rest dynamic direction
```

rather than:

```text
confirmed direct ecological interaction
```

---

## Reproducibility requirements

For every beta analysis, record:

- pcLVbayes Git commit;
- R version;
- cmdstanr version;
- CmdStan version;
- operating system;
- compiler version;
- random seeds;
- chain count;
- warmup iterations;
- sampling iterations;
- retry policy version;
- preprocessing schema;
- predictive scorer identifier;
- selected time unit;
- zero-replacement settings;
- ALR cap;
- spline configuration;
- K-fold split manifest;
- dataset checksum or immutable identifier.

Results generated under different preprocessing or scoring schemas should not
be pooled without explicit compatibility validation.

---

## Installation

### Beta-testing period

Public installation may be disabled during the private beta-testing period.

Access to source code does not necessarily grant permission to install,
redistribute, modify, or deploy the package. Refer to `LICENSE` for the
applicable beta terms.

Beta testers should install only from an authorized repository or provided
source archive and should preserve the exact commit identifier.

### Planned stable release

After the stable open-source release, installation is expected to use the
following dependencies.

```r
install.packages(c(
  "dplyr",
  "tidyr",
  "tibble",
  "purrr",
  "metafor",
  "posterior"
))

install.packages(
  "cmdstanr",
  repos = c(
    "https://mc-stan.org/r-packages/",
    getOption("repos")
  )
)

cmdstanr::install_cmdstan()
```

Bioconductor dependency:

```r
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

BiocManager::install("phyloseq")
```

Planned post-release development installation:

```r
if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}

remotes::install_github("gyeongjunCho/pcLVbayes")
```

The installation interface and dependency list may change before the stable
release.

---

## Reporting issues during beta testing

Beta testers should report:

- the exact Git commit;
- a minimal reproducible example;
- sanitized input dimensions and metadata structure;
- the complete error message;
- failure payloads;
- sampling diagnostics;
- relevant CmdStan output;
- operating-system and toolchain information.

Do not upload confidential sample metadata or restricted biological data to a
public issue tracker.

---

## Citation

A formal citation will be provided with the stable release.

Until then, analyses should cite the exact software commit and identify
pcLVbayes as beta software.

Example:

```text
pcLVbayes beta, Git commit <commit>, accessed <date>.
```

---

## Disclaimer

pcLVbayes is experimental research software.

The beta version is provided without guarantees of:

- API stability;
- numerical equivalence across versions;
- backward-compatible checkpoints;
- production reliability;
- biological validity for every dataset;
- fitness for clinical, regulatory, or commercial use.

All scientific conclusions remain the responsibility of the analyst.
