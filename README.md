# pcLVbayes

> **Status: Beta (experimental)** — APIs and defaults may change; validate on your data and report issues.

Bayesian **pairwise compositional Lotka-Volterra** regressions for amplicon time-series.
Includes **correlation-based screening**, **Δt-aware OU residuals**, K×R **repeated K-fold**
with **Kalman ELPD**, tidy outputs, and diagnostics. Built on **cmdstanr/CmdStan**.

[![R](https://img.shields.io/badge/R-%3E%3D%204.2-276DC3.svg)](https://cran.r-project.org/)
[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)

---

## Licensing (Pre-release)

This pre-release is **view-only via the official GitFront link** and **not licensed for redistribution or use**.
A stable open-source release is planned under **GPL-3**.  
See `LICENSE` for preview-only terms.

This repository is temporarily viewable via GitFront for a conference demo at the
**2025 KSPP Fall International Conference (Oct 21–24, 2025; The-K Hotel Gyeongju)**.

© 2025 Rural Development Administration (RDA) & National Institute of Agricultural Sciences (NAS).

---

## ✨ What it does

- **Correlation screening — `cor_meta_resid()`**  
  Computes within-subject correlations (Pearson/Spearman) and meta-analyzes them across subjects (e.g., **DL**, optional **KNHA**).  
  Includes **ACF/AR(1) effective-n** correction (`effn_*`) so uneven sampling and autocorrelation are penalized.  
  Use this to shortlist partners before MCMC.

- **Bayesian pairwise pcLV Core — `fit_pclv_bayes()`**
  For each unordered pair {i, j}, fits two ordered regressions (j→i and i→j). Predictors are lag-1 pair-to-rest ALR values and the response is ΔALR\_i / Δt. The model uses irregular-time OU residuals, Student-t observations with one posterior-estimated `nu > 2` per directed fit, and repeated subject-level K-fold Kalman ELPD.

---

## Scientific estimand and interpretation of `a_ij`

pcLVbayes estimates a **directed pair-to-rest dynamic coefficient**: the posterior coefficient of the lagged source pair-to-rest log-ratio predictor in the target pair-to-rest log-ratio rate model, conditional on the target self predictor and the canonical preprocessing and residual model.

This is a model-level compositional coefficient, not generally the absolute joint-gLV coefficient `A[i,j]`, a context-independent direct ecological interaction, or an absolute-abundance effect size. A positive `a_ij` is a positive fitted pair-to-rest dynamic direction under the declared model and observed state distribution; a negative value is a negative fitted direction. Neither sign alone proves biological facilitation or inhibition.

The absolute coefficient `A[i,j]`, the state-dependent mechanistic contrast `C_ij(x) = A[i,j] - sum_{k in rest} w_k(x) A[k,j]`, the canonical deterministic transformed-data projection, and the Bayesian posterior coefficient `a_ij` are related but not interchangeable. A stable posterior sign validates the fitted transformed coefficient, not an absolute direct-gLV sign. A nonzero transformed coefficient may occur when `A[i,j] = 0`; an indeterminate direction is not a zero interaction. Absolute direct-interaction claims require additional assumptions or absolute-scale information.

## Core transformation

The Core uses zero-aware pair-to-rest ALR for the `(i, j, rest)` triplet. Predictors are lag-1 ALR values; the unscaled response is ΔALR\_i / Δt.

---

## Validation & ELPD

Repeated subject-level K-fold ELPD is always computed. Each fold uses training-only global predictor scaling and Kalman OU scoring. Failed evaluations remain unavailable (`NA`) rather than contributing zero ELPD; model weights use only common successfully evaluated evidence.

---

## Noise model & scoring

The Core estimates one shared Student-t degrees-of-freedom parameter per directed fit, constrained by `nu = 2 + exp(log_nu_minus_two)`. Irregular-time OU prediction is scored through the Kalman filter, including the Student-t variance expansion.

---

## Sampling robustness

A diagnostics-aware retry policy monitors divergences, tree depth, and E-BFMI without changing the Student-t model or its `nu` prior. Retry and Pathfinder settings follow the private canonical v0.2 policy; diagnostics are returned for downstream filtering.

---

## Residual model

The Core uses an irregular-time Ornstein–Uhlenbeck residual process with Δt-aware persistence.

---

## Time axis (Δt) policy

- Use **one unit** project-wide (e.g., **days** or **weeks**). **Decimals are allowed** (e.g., 3.5 days, 1.25 weeks).  
- **POSIXct/Date** is recommended but optional; if provided, values are converted via `difftime(..., units = "<unit>")`.  
- **Do not mix units** across subjects/samples—this changes the scale/interpretation of **ΔALR/Δt**, coefficients (`a_ij`, `a_ii`), and OU decay (`λ`/`φ`).  
- Document the chosen unit (e.g., *All Δt are in weeks*).

---

## Compositional handling (ALR)

- Recommended: **ALR (pair-to-rest)** on the **(i, j, rest)** triplet with **zero-aware replacement** (optionally capped).  
- Predictors: **lag-1 pair-to-rest ALR**. Response: **ΔALR\_i / Δt**.
- **ILR** was evaluated but **not used** for correlation screening here, as it tended to **inflate negative correlations** in this setting; **ALR** is the default.

---

## K-fold, leakage guards & aggregation

- **Repeated K-fold** at the **subject** level (K×R) yields **per-subject OOF ELPD**.  
- **No leakage**:
  - **Train-only global predictor scaling/centering** inside each fold.
  - Fold runs do **not** reuse main-run tuning artifacts.  
- ELPD uses the `"kalman-ou"` scorer, and **aggregation is subject-uniform** by default.

---

## Parallelism & progress

- Effective concurrency ≈ **`n_workers_outer × chains`**.  
  - Outer pair loop: `n_workers_outer` (PSOCK via **future/furrr**).  
  - Within each pair: Stan **chains** (`cmdstanr::sample()` with `parallel_chains`).  
  - If `n_workers_outer > 1`, **K-fold parallelism is disabled** to avoid nested parallelism.
- **Oversubscription guidance**  
  - Aim for `n_workers_outer × chains` ≤ **logical CPU threads**.  
  - In practice chains under-utilize cores; on idle machines, up to **~2×** can be acceptable—monitor thermals.
- If **progressr** is installed and `progress = "bar"`, progress bars are shown for both outer pairs and K-fold;
  otherwise it falls back silently.

---

## Interpretation

For edges, prioritize **posterior-supported pair-to-rest directions** (sign probabilities / LFSR from draws) and **MCMC diagnostics**;
use **per-subject ELPD** (from K×R) as **supporting evidence**.

---

## Installation

### Preview period
Installation is **disabled during the preview** (the repository is private and view-only via GitFront).

### After the stable open-source release
```r
# CRAN deps (examples)
install.packages(c(
  "dplyr","tidyr","tibble","purrr",
  "metafor","posterior","cmdstanr"
))

# CmdStan toolchain
install.packages("cmdstanr", repos = c("https://mc-stan.org/r-packages/", getOption("repos")))
cmdstanr::install_cmdstan()  # requires a C++ toolchain (Rtools/CLT/build-essential)

# Bioconductor
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install("phyloseq")

# Dev install (post-release; example)
if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")
remotes::install_github("gyeongjunCho/pcLVbayes")
```
