# MTIST 10-species benchmark baseline

## Selection and truth contract

Measured 2026-08-02 with pcLVbayes `c14e9f1` and MTIST
`b0da1e5c41e9562d6113c14aa4de9d553dfa79a6`. The local catalog contains 324
ten-species studies. The deterministic rule retains compatible studies with
noise 0.01, even sampling, at least 10 series, and at least 15 time points,
then chooses the lowest dataset ID. Twenty-seven studies meet the filter;
dataset 37 is first.

Dataset 37 uses truth `10_sp_gt_1`, taxa `species_0` through `species_9`, 10
independent series, 15 observations per series, and times from 0 to 30. MTIST
labels the scheme `even`; its floating grid gives 13 intervals of 2.121212 and
one of 2.424242 per subject. All 1,500 abundance cells are finite and positive
(minimum 0.004279422), with no observed zeros or extinctions.

MTIST source defines interaction rows as focal/target species and columns as
interacting/source species; its simulator uses `dot(A, y)`. The benchmark maps
`source -> target` to `A[target, source]`. The truth is non-symmetric in 74
cells and has 16 zero off-diagonal entries. Tests use unequal transposed cells,
so orientation is not inferred from prediction signs.

## Configurations and runtime

Stage A command:

```sh
Rscript benchmarks/mtist/run_ten_species.R
```

It uses one chain, 50 warmup, 50 retained draws, K=2/R=1, six outer workers,
one K-fold worker, seed 20260802, no retries, and no Pathfinder.

Stage B command:

```sh
PCLV_MTIST_10_CONFIG=benchmarks/mtist/configs/ten_species_reference.R \
  Rscript benchmarks/mtist/run_ten_species.R
```

It uses two chains, 250 warmup, 500 retained draws per chain, with the same
K-fold, workers, seed, retry, Pathfinder, and canonical sampler controls. This
is a development reference, not publication-grade inference. Stage A took
128.046 seconds. The pre-run projection was about 16 minutes for main fits and
48 minutes if every direction also ran two folds, below the fixed two-hour
ceiling. Stage B was therefore run and took 396.332 seconds.

| Stage | Pair tasks | Directed records | Timed wall seconds |
|---|---:|---:|---:|
| A: 1 x 50/50 | 45 | 90 | 128.046 |
| B: 2 x 250/500 | 45 | 90 | 396.332 |

Compilation and isolated package installation were excluded. One outer pool
used explicit inputs, one future per pair, canonical order restoration, and
one effective K-fold worker. Up to 11 simultaneous CmdStan chains were
observed (theoretical maximum 12). Per-task elapsed time is not exposed, so
named slowest tasks and min/median/max task times remain unavailable. The
observed taper from 11 to 6 chains is the task-imbalance proxy. No retry ran.
All 90 seeds were unique and direction indices were exactly 1 through 90.

## Diagnostics and downstream stages

All Stage A directions were `sampler_diagnostics_failed`, as expected for the
execution-only short chains. None reached evidence, K-fold, ELPD, stacking, or
significance. Its ES values were all 0.5 because external matrices were zero;
this is MTIST's neutral-zero baseline, not successful inference.

Stage B produced:

| Diagnostic class | Count | Fraction |
|---|---:|---:|
| converged | 6 | 6.67% |
| interaction_stable_residual_unstable | 61 | 67.78% |
| interaction_indeterminate | 23 | 25.56% |
| sampler_diagnostics_failed | 0 | 0% |

| Stage | Count |
|---|---:|
| main attempted / posterior retained | 90 / 90 |
| converged / Bayesian eligible | 6 / 6 |
| K-fold attempted / completed | 6 / 6 |
| finite ELPD / stacking available | 6 / 5 |
| final significant | 6 |

Residual-unstable directions remain explicit and are not promoted. The 23
indeterminate directions retain finite coefficients internally and are not
interpreted as zero. Only complete external matrices use zero placeholders,
with separate execution, diagnostic, indeterminate, residual, significance,
and ELPD masks.

Of 45 unordered pairs, 39 had neither direction converge, six had one, and
none had both. Twenty-two pairs had different directional classes, confirming
opposite-direction independence.

## Scores and coverage

Off-diagonal results are primary because MTIST truth is joint absolute-gLV,
whereas pcLVbayes estimates compositional pair-to-rest effects.

| Matrix | ES with diagonal | ES without diagonal |
|---|---:|---:|
| posterior mean | 0.523810 | 0.506757 |
| conservative | 0.523810 | 0.506757 |

All six converged directions passed LFSR <= 0.05. Conditional sign accuracy
was 3/6 (50.0%) at 6/90 coverage. Eligible signs covered 1/33 positive and
4/41 negative truth edges. Of 74 nonzero truth directions, 69 (93.24%) were
omitted from significant output; 21/74 (28.38%) were specifically interaction
indeterminate and none were sampler failures. One of 16 truth-zero cells was
called nonzero. Accuracy must always be reported with these denominators.

All six eligible directions completed ELPD. Five obtained stacking weights;
one target lacked sufficient common model evidence. Bayesian evidence was
finalized before K-fold and stacking, and stacking did not mutate it.

## Diagonal

The diagonal is the **median across pair-specific posterior self-effect
means**, not a posterior median. Contributors existed for species 4 (2,
-0.0270475), species 7 (3, -0.487793), and species 9 (1, -0.0373434). Seven
diagonals were unavailable, zero only in required matrices, and explicitly
marked unavailable. Both diagonal-included and excluded scores are retained.

## Memory, disk, and cleanup

Stage B parent RSS rose from 871,068 to 879,720 kB; parent HWM was 879,720 kB.
The fit object was 1,354,936 bytes. Full process-tree peak RSS was unavailable
after removing an interfering forked monitor and is not inferred.

The historical CmdStan cache baseline/final snapshot was 191 files, 67
directories, and 387,936 kB. A live peak was 227 files and 404,296 kB: 36
active files and about 16 MB above baseline. It returned to baseline. The
dedicated parent root was 0 files/0 bytes before and after. No worker or
CmdStan process remained. Session/future files are not called leaks, and
SIGKILL cleanup is not claimed. The installed and source executables were
unchanged; no worker compiled a model.

## Release implication

Operationally, the frozen Core completed all 45 tasks and 90 directions in a
practical runtime with deterministic ordering, bounded outputs, executable
reuse, and complete failures. This supports proceeding to packaging and
documentation work.

Scientifically, 6.67% full identifiability and 50% conditional sign accuracy
do not support broad MTIST recovery claims. Release preparation must preserve
explicit indeterminate/residual-unstable reporting and document the
compositional-versus-absolute estimand mismatch. No model, scoring,
preprocessing, K-fold, or production optimization change was made.

## Targeted four-chain confirmation

The six directions provisionally classified as converged in Stage B were
selected deterministically from the retained Stage B direction_results.tsv
using only diagnostic_class == "converged". Selection did not use MTIST truth,
sign correctness, coefficient magnitude, or taxon names. The six records belong
to six distinct unordered-pair tasks.

The confirmation used the private canonical main-posterior phase introduced by
the execution-boundary refactor. The public fit_pclv_bayes() path still calls
that phase followed by the unchanged eligible-direction K-fold evaluation; ELPD
and stacking remain downstream summarization inputs. The confirmation stopped
after sampling/retry, retained-attempt diagnostics, the posterior summary
bundle, sign probability/LFSR, and diagnostic classification. It invoked no
K-fold, Kalman scoring, ELPD, or stacking.

The fixed configuration was four chains, 2,000 warmup iterations, 2,000 retained
draws per chain, seed 20260802 with unchanged direction-seed derivation,
Pathfinder disabled, zero retries, canonical adapt_delta = 0.98, canonical
max_treedepth = 14, canonical diag_e metric, and three outer workers. Thus at
most 12 CmdStan chains ran simultaneously.

Stage B consumed 153,000 chain-iterations when its 90 main fits and six
two-fold predictive evaluations are counted. Scaling its 396.332 seconds to
the confirmation's 96,000 main-fit chain-iterations projected 248.679 seconds;
the conservative threefold projection was 746.036 seconds. Projected temporary
disk was 85 MiB. Actual main-fit wall time was 27.021 seconds.

Stage B medians are NA because the retained Stage B schema explicitly recorded
that pooled posterior medians were not exposed. Long-run medians are retained
by the private main-result contract.

| Source -> target | Task / direction | Seed | Stage B mean | Long mean / median | Long P(+)/P(-), LFSR | Stage B -> long class | Outcome | Truth | Correct |
|---|---:|---:|---:|---:|---:|---|---|---:|---|
| species_1 -> species_7 | 15 / 30 | 20468804 | 0.168680 | 0.170775 / 0.169259 | 1 / 0, 0 | converged -> converged | confirmed_converged_same_sign | -0.017456 (-) | no |
| species_2 -> species_9 | 24 / 48 | 20570804 | 0.026579 | 0.026718 / 0.026488 | 0.992375 / 0.007625, 0.007625 | converged -> converged | confirmed_converged_same_sign | 0.073251 (+) | yes |
| species_3 -> species_7 | 28 / 56 | 20668804 | -0.097970 | -0.098411 / -0.097429 | 0.001125 / 0.998875, 0.001125 | converged -> converged | confirmed_converged_same_sign | -0.227619 (-) | yes |
| species_7 -> species_4 | 33 / 65 | 20768803 | -0.028977 | -0.029051 / -0.028895 | 0 / 1, 0 | converged -> converged | confirmed_converged_same_sign | 0 (zero) | no |
| species_8 -> species_4 | 34 / 67 | 20769803 | -0.018702 | -0.018587 / -0.018556 | 0 / 1, 0 | converged -> converged | confirmed_converged_same_sign | -0.268544 (-) | yes |
| species_5 -> species_7 | 37 / 74 | 20868804 | 0.329838 | 0.332788 / 0.333871 | 1 / 0, 0 | converged -> converged | confirmed_converged_same_sign | -0.348405 (-) | no |

All 24 chain-specific interaction medians retained their Stage B dominant sign,
and every direction had chain-sign agreement. R-hat ranged from 1.00058 to
1.00158, minimum bulk ESS from 2,948.87 to 4,708.52, and minimum tail ESS from
3,444.68 to 3,901.57. Every direction had zero divergences and zero
maximum-treedepth saturations. All 24 chain-specific E-BFMI values exceeded
0.899. Chain-specific sigma, sd_ou, phi, lambda, and nu intervals overlapped for
every direction; residual-regime disagreement was false in all six, and no
structured diagnostic reason was emitted. Detailed chain summaries remain in
the ignored confirmation_comparison.tsv.

Outcome counts were 6 confirmed_converged_same_sign, 0
confirmed_converged_sign_changed, 0 downgraded_residual_unstable, 0
downgraded_indeterminate, and 0 long_run_sampler_failed. Six of six retained
the Stage B sign and zero of six changed sign. Six of six remained converged and
zero of six were downgraded. Long-run sign accuracy was 3/6 (50%) over the full
prespecified denominator and 3/6 (50%) among directions still converged.
Therefore the provisional Stage B 3/6 result was confirmed, not improved or
made non-interpretable.

Parent RSS increased from 783,576 to 848,816 kB and parent HWM was 848,816 kB.
The returned six-result object was 138,912 bytes. The package-owned temporary
root was 0 files, 0 directories, and 0 bytes before and after. The historical
CmdStan cache was unchanged at 191 files, 67 directories, and 396,443,655
bytes. The installed executable modification time was unchanged, so no
recompilation occurred. No benchmark worker or CmdStan process remained.

Same-sign convergence supports posterior sign stability for these six
directions. It does not remove the mismatch between MTIST's joint
absolute-abundance gLV coefficients and pcLVbayes's compositional pair-to-rest
effects. A future indeterminate or failed direction must remain explicit and
must not be interpreted as a zero interaction.
