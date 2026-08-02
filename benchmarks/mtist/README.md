# MTIST smoke benchmark

This directory contains benchmark-only infrastructure for converting one MTIST
absolute-density study to a `phyloseq` object and evaluating pcLVbayes output.
It is not sourced by package tests or package builds.

The smoke configuration selects dataset 361: three taxa, 10 independent time
series, 15 time points per series, even sampling, noise 0.01, and truth
`3_sp_gt_1`. Absolute densities are strictly validated and converted row-wise
to relative abundance; no read counts are simulated.

Run from the package repository with:

```sh
Rscript benchmarks/mtist/run_smoke.R
```

Set `MTIST_ROOT` to override the default `~/mtist` checkout.

Matrices follow MTIST orientation: row = target and column = source. Failed
directions are numeric zero only in required matrices and remain explicit in
`failure_mask.tsv`. The diagonal is the **median across pair-specific posterior
self-effect means**. This is not a posterior median. The same diagonal is used
in both matrices because no new self-effect uncertainty aggregation is
invented. Conservative filtering is benchmark-specific and off-diagonal only:
retain a successful posterior mean when LFSR (`p_sign2 / 2`) is at most 0.05.

The primary scientific result is the off-diagonal, diagonal-excluded score.
MTIST truth contains joint absolute-abundance gLV coefficients, whereas
pcLVbayes estimates compositional pair-to-rest effects and repeated
pair-specific self effects. Scores including the diagonal are retained only
for benchmark compatibility.

## Ten-species benchmark

The deterministic 10-species benchmark selects dataset 37 and runs fixed smoke
and multi-chain reference stages. Commands, selection criteria, orientation,
configurations, and reviewed results are in `TEN_SPECIES_BASELINE.md`.
Generated outputs remain ignored under `results/`.

Run the targeted four-chain main-posterior confirmation with:

    Rscript benchmarks/mtist/run_ten_species_confirmation.R

This command deterministically selects the six Stage B converged directions and
does not run K-fold, Kalman scoring, ELPD, or stacking. The public canonical fit
path remains unchanged and continues through predictive evaluation.
