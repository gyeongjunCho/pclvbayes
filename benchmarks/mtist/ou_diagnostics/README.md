# MTIST OU diagnostics

This experiment isolates the `species_2 -> species_0` direction from MTIST
dataset 361. It compares the unchanged canonical model with an exactly
equivalent continuous-decay coordinate and a minimally regularized version.
All variants use the same canonical preprocessed data, seed, sampler settings,
and Student-t observation model. Generated output is ignored under `results/`.

Run from the package root:

```sh
Rscript benchmarks/mtist/ou_diagnostics/run_diagnostics.R
```
