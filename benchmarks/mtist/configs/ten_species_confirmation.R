list(
  stage = "targeted_four_chain_confirmation",
  dataset_id = 37L,
  stage_b_result_label = "ten_species_37_reference_250_500",
  result_label = "ten_species_37_confirmation_4x2000_2000",
  seed = 20260802L,
  chains = 4L,
  iter_warmup = 2000L,
  iter_sampling = 2000L,
  n_workers_outer = 3L,
  max_retries = 0L,
  use_pathfinder_init = FALSE,
  run_kfold = FALSE,
  calculate_elpd = FALSE,
  calculate_stacking = FALSE,
  maximum_projected_seconds = 7200,
  maximum_projected_temp_bytes = 10 * 1024^3
)
