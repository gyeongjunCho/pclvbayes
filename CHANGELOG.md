# CHANGELOG

## v0.2 (In Progress)

### Scientific estimand clarification
- Clarify that `a_ij` is a directed pair-to-rest dynamic coefficient, distinct from absolute direct-gLV `A[i,j]`.
- Record the oracle audit findings: transformed-oracle agreement 6/6, absolute-A cross-estimand agreement 3/6, 44/90 state-dependent contrasts, 16/90 absolute-zero induced transformed effects, and the `species_7 -> species_4` smoothing-sensitive case.
- Separate transformed-oracle validation from absolute-A comparison in benchmark reporting.
- No posterior-model, computation, public-API, or result-schema change.


### v0.2 public API freeze
- Retain only scientific inputs, sampling/reproducibility controls, progress,
  resource controls, and K-fold controls in `fit_pclv_bayes()`.
- Remove fixed canonical preprocessing controls and exposed numerical constants:
  `zero_mode_alr`, `minpos_alpha`, `minpos_base`, `smooth_scale`,
  `alr_spline_df`, `alr_spline_spar`, `alr_spline_cv`, `eps`, `eps_fixed`,
  `lib_eps_c`, and `rest_floor_frac`.
- Remove `metric`, `quiet`, `progress_every`, `silent_sampler`, `max_retries`,
  and all Pathfinder tuning arguments; canonical retry and Pathfinder behavior
  is internal, and `progress` is the sole public progress control.
- Remove the obsolete `fit_pclv_bayes2()` cross-kingdom prototype and its
  generated help; no unique canonical Core behavior was migrated from it.
- Removed arguments now fail as unused arguments; migration uses canonical v0.2
  preprocessing, internal retry/Pathfinder settings, and retained sampling/K-fold
  controls. No posterior result schema or scientific estimand changed.

### Planned
- Simplify public API
- Remove obsolete alpha parameters
- Eliminate silent fallback logic
- Reduce conditional branching
- Introduce validated intermediate objects
- Preserve Core OU model
