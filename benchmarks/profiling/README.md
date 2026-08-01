# pcLVbayes runtime profiling

This benchmark-only harness measures deterministic R components and two small
end-to-end MTIST 361 fits. Run it from the repository root:

```sh
Rscript benchmarks/profiling/run_profile.R
```

The default configuration uses one chain, 50 warmup and 50 retained draws,
two subject-level folds, no retries, no Pathfinder, and compares one versus two
outer workers. Model compilation is measured separately and excluded from fit
timings. Generated machine-specific output is ignored under `results/`.

`end_to_end_fit` is necessarily a combined wall-time phase: the current
CmdStanR boundary does not reliably separate JSON writing, process startup,
warmup, sampling, CSV loading, and cleanup after the fit object is deliberately
released. The direct JSON probe and R call-stack profile provide supporting
overhead measurements but are not subtracted from end-to-end time.

`resources.csv` reports the parent R process high-water RSS and temporary
directory snapshots before and after each run. It does not claim a process-tree
memory peak or an in-run disk peak; those fields are explicitly `NA` in metadata
until a portable external monitor is added.
