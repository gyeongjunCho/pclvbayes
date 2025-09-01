# pGLVbayes

> **Status: Beta (experimental)** — APIs and defaults may change; validate on your data and report issues.

Bayesian **pairwise gLV-inspired** regressions for amplicon time-series.
Includes **correlation-based screening**, **Δt-aware OU residuals**, K×R **repeated K-fold**
with **Kalman ELPD**, tidy outputs, and diagnostics. Built on **cmdstanr/CmdStan**.

[![R](https://img.shields.io/badge/R-%3E%3D%204.2-276DC3.svg)](https://cran.r-project.org/)
[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)

---

## Licensing (Pre-release)

This pre-release is **view-only via the official GitFront link** and **not licensed for redistribution or use**.  
A stable open-source release is planned under **MIT** or **GPL-3**.  
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

- **Bayesian pairwise gLV (core) — `fit_glv_pairwise()`**  
  For each unordered pair {i, j}, fits **two ordered** regressions (**j→i**, **i→j**).  
  Response: **ΔALR\_i / Δt** (standardized upstream as needed). Predictors: **lagged** states of i and j (ALR or raw RA).  
  Residuals: **OU (≈ AR(1))** with **Δt-aware** persistence **or** **white-noise**; optional **Student-t** observation noise.  
  Returns directional effects (`a_ij`, `a_ii`/`a_jj`), **posterior sign probabilities** and an **approx. local false sign rate**  
  *(LFSR ≈ 1 − max{Pr(a>0), Pr(a<0)} from posterior draws)*, **MCMC diagnostics**, and **K×R ELPD per subject**.

---

## Transforms

- `transform = "alr"` (**recommended**): zero-aware replacement on the triplet **(i, j, rest)**;  
  predictors are **lagged ALR**, response is **ΔALR\_i / Δt** (then standardized as configured).
- `transform = "raw"`: predictors use **lagged relative abundances**; the response remains **ΔALR\_i / Δt**.

---

## Validation & ELPD (automatic modes)

When `compute_elpd = TRUE`, the package runs **Repeated K-fold (K×R)** at the **subject level** and computes
**per-subject projection log-likelihood** under the **fold-train posterior**, returning the **mean ELPD per subject**.
**Train-only scaling** is used to avoid out-of-fold leakage.

- If `elpd_mode` is **omitted/NULL**, the scorer is chosen automatically:  
  - `resid_mode = "ou"` ⟶ **`kalman`** (exact OU state-space marginal log-lik via Kalman filter).  
  - `resid_mode = "wn"` ⟶ **`indep`** (i.i.d. Gaussian likelihood).
- You can override with `elpd_mode = "kalman"` or `"indep"` explicitly.
- **`kalman` (OU only)**: irregular-interval OU innovations are integrated as a **linear Gaussian state-space model**;  
  measurement variance uses `sigma^2` (or a t→Gaussian variance expansion when `use_student_t = TRUE`).  
  If required OU parameters are unavailable, it **falls back** to `indep`.
- **`indep`**: independent Gaussian scoring.  
  Predictive SD is  
  \[
    \text{sd}_{\text{pred}} = \begin{cases}
      \sqrt{\sigma^2 + \mathrm{sd\_ou}^2} & \text{(OU-like)} \\
      \sqrt{\sigma^2} & \text{(white-noise)}
    \end{cases}
  \]
  (or uses `sigma_pred` if provided).

---

## Noise model & scoring

- Optional **Student-t** observation noise via `use_student_t = TRUE` (also honored during K-fold scoring).  
- For `elpd_mode="kalman"` (OU), OU transition \(a_t=\exp(-\lambda \Delta t)\) and process noise \(q_t=\mathrm{sd\_ou}^2(1-a_t^2)\)  
  are handled **inside** the Kalman filter; measurement variance uses `sigma^2` (or its t-variance expansion).

---

## Sampling robustness

A diagnostics-aware **retry** policy (up to `max_retries`) monitors **divergences**, **tree depth**, and **E-BFMI**.  
If needed, it **escalates `nu_fixed` within [4, 7]** while keeping other sampler controls stable for fold runs.  
Diagnostics are returned for downstream filtering.

---

## Residual model — `resid_mode`

- `"ou"` (**recommended**): Ornstein–Uhlenbeck (≈ AR(1)) residuals with **Δt-aware** persistence; best when residual autocorrelation exists or sampling is uneven.  
  For ELPD, prefer the **Kalman** marginal-likelihood path.  
- `"wn"`: white-noise (independent) residuals — useful as a fast baseline or for very short per-subject series.

---

## Time axis (Δt) policy

- Use **one unit** project-wide (e.g., **days** or **weeks**). **Decimals are allowed** (e.g., 3.5 days, 1.25 weeks).  
- **POSIXct/Date** is recommended but optional; if provided, values are converted via `difftime(..., units = "<unit>")`.  
- **Do not mix units** across subjects/samples—this changes the scale/interpretation of **ΔALR/Δt**, coefficients (`a_ij`, `a_ii`), and OU decay (`λ`/`φ`).  
- Document the chosen unit (e.g., *All Δt are in weeks*).

---

## Compositional handling (ALR)

- Recommended: **ALR (pair-to-rest)** on the **(i, j, rest)** triplet with **zero-aware replacement** (optionally capped).  
- Predictors: **lagged** ALR (or raw RA with `transform = "raw"`). Response: **ΔALR\_i / Δt**.  
- **ILR** was evaluated but **not used** for correlation screening here, as it tended to **inflate negative correlations** in this setting; **ALR** is the default.

---

## K-fold, leakage guards & aggregation

- **Repeated K-fold** at the **subject** level (K×R) yields **per-subject OOF ELPD**.  
- **No leakage**:
  - **Train-only scaling/centering** (optionally **per subject**) inside each fold.  
  - Fold runs do **not** reuse main-run tuning artifacts.  
- ELPD is reported with the scorer used (e.g., `"kalman-ou"` or `"indep"`), and **aggregation is subject-uniform** by default.

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

For edges, prioritize **signs** (posterior sign probabilities / LFSR from draws) and **MCMC diagnostics**;  
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
remotes::install_github("gyeongjunCho/pGLVbayes")
```
