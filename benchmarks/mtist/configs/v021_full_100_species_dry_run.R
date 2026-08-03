dry_run_root <- Sys.getenv("PCLV_V021_DRY_RUN_ROOT", "")
if (!nzchar(dry_run_root))
  stop("PCLV_V021_DRY_RUN_ROOT is required for the no-sampling dry run.")
build_v021_full_config(dry_run_root)
