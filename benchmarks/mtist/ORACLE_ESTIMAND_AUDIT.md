# MTIST oracle estimand audit

## Scientific question

This benchmark-side audit diagnoses the three Stage B/four-chain sign
disagreements without treating the MTIST absolute-gLV coefficient and the
pcLVbayes posterior coefficient as the same estimand. It uses dataset 37
(10_sp_gt_1) and retained four-chain results. It launches no Stan sampling,
K-fold, ELPD, or stacking.

## Verified estimands and source locations

MTIST defines the absolute joint gLV system as dx_i/dt = x_i g_i(x), with
g_i(x) = r_i + sum_j A[i,j] x_j. Rows are focal/target species and columns
are source species. This is implemented by mtist_utils.py:create_lv_dicts()
and lvutils.py:run_lv(), whose RHS uses dot(A, y) + mu and multiplies by the
state. The adapter preserves this orientation in mtist_adapter.R.

For target i and source j, the canonical rest set is taxa minus i and j. The
canonical denominator is exactly the sum of all other relative abundances,
implemented by R/pclv_helpers.R:.pair_to_rest_abundance() and
.make_pair_inputs_glv(). The transformed coordinates are log(x_i/R) and
log(x_j/R), followed by within-subject lagging, Delta ALR / Delta t, global
predictor standardization, and the Stan mean equation in inst/stan/pclv.stan:
mu[n] = r0 + r0_sub[sid[n]] + a_ii * xi[n] + a_ij * xj[n].

Differentiating the denominator gives
d log(x_i/R) / dt = g_i(x) - sum_k w_k(x) g_k(x), where w_k = x_k/R.
The source-j contrast is C_ij(x) = A[i,j] - sum_k w_k A[k,j], with
instantaneous source contribution x_j C_ij(x). The audit helper computes both
the complete derivative and this source-specific contrast.

## Dataset and trajectory source

Dataset 37 has 10 taxa, 10 subjects, 15 retained observations per subject,
times 0 to 30, and noise metadata 0.01. The dataset stores absolute abundances;
the adapter converts them row-wise to relative abundance.

The retained official noiseless trajectories in
results/oracle_estimand_audit/official_noiseless_100.csv and
official_noiseless_dense.csv were regenerated through MTIST's own Python
simulator using stored A, growth rates, subject seeds, initial states, tend=30,
dt=0.1, and sample_freq=100. All ten initial states reproduce dataset 37
exactly. A forward difference against the dense solver output has median
absolute discrepancy 0.00743 and maximum 0.0950 at the solver's finite
internal step; this is a numerical finite-step check, not the algebraic
derivative used for Oracle 3.

## Oracle construction

oracle_estimand_audit.R implements seven levels: absolute A truth; statewise
instantaneous C_ij summaries; exact instantaneous QR projection; exact sampled
finite-interval projection; canonical noiseless preprocessing projection;
observed-data preprocessing projection; and retained four-chain posterior
values. Projections use intercept, subject effects where identifiable, lagged
target ALR, and lagged partner ALR. They report rank, condition number, and
structured rank failures. No regularization or zero substitution is used.

All 90 directions were full-rank at every deterministic projection level.
The exact derivative and finite-interval designs use the same positive global
predictor scaling as production preprocessing. Subject boundaries reset exactly
as in .delta_alr_over_dt().

## Focused six-direction comparison

The six directions were selected only from Stage B records with diagnostic_class
equal to converged.

| Source -> target | A | Contrast | Instant | Finite | Noiseless | Observed | Posterior | Primary category |
|---|---:|---|---:|---:|---:|---:|---:|---|
| species_1 -> species_7 | - | changes sign | + | + | + | + | + | state-dependent_no_single_sign |
| species_2 -> species_9 | + | changes sign | + | + | + | + | + | state-dependent_no_single_sign |
| species_3 -> species_7 | - | changes sign | - | - | - | - | - | state-dependent_no_single_sign |
| species_7 -> species_4 | 0 | changes sign | + | + | - | - | - | truth_zero_induced_compositional_effect |
| species_8 -> species_4 | - | constant negative | - | - | - | - | - | unresolved |
| species_5 -> species_7 | - | changes sign | + | + | + | + | + | state-dependent_no_single_sign |

The three A-versus-posterior disagreements are species_1 -> species_7,
species_7 -> species_4 (A exactly zero), and species_5 -> species_7. Every
deterministic projection has the posterior sign for these directions. The
first and last have state-dependent contrasts; the zero-truth direction has a
nonzero induced compositional contrast and a preprocessing sign change from
finite projection + to canonical/observed projection -.

The three signs agreeing with A do not establish A recovery. Species_2 ->
species_9 and species_3 -> species_7 also have state-dependent contrasts;
species_8 -> species_4 has a constant-negative contrast. There were no focused
finite-interval sign flips, no noiseless-to-observed sign flips, and no
posterior-versus-observed deterministic sign flips.

## All-direction summaries

All denominators are 90 directed interactions.

| Comparison with A | Agreement | Positive coverage | Negative coverage |
|---|---:|---:|---:|
| instantaneous contrast | 49/90 | 22/33 | 27/41 |
| finite-interval projection | 46/90 | 21/33 | 25/41 |
| canonical noiseless projection | 42/90 | 18/33 | 24/41 |
| observed deterministic projection | 42/90 | 18/33 | 24/41 |

There are 44/90 state-dependent sign-changing contrasts. Sixteen exact-zero A
entries produce nonzero transformed signs at every oracle level. All 90
projections are full rank. Complete 3x3 sign cross-tabs are in the ignored
results/oracle_estimand_audit/sign_cross_tabs.tsv. Official MTIST ES is not
used as the primary audit metric.

## Sensitivity checks

For the six focused directions, exact derivative and finite-interval signs
agree for all six. No-spline and canonical noiseless signs agree except for
species_7 -> species_4, where smoothing changes + to -. Observed pooled signs
remain unchanged under subject-wise and leave-one-subject-out checks, although
species_2 -> species_9 has some subject-specific sign instability. Contrast
weight ranges span multiple states for every state-dependent focused direction.

## Interpretation and release recommendation

Recommendation A is supported for the confirmed directions: posterior signs
agree with canonical noiseless and observed deterministic transformed-data
projections. They should be described as pair-to-rest compositional effects or
transformed-data projection coefficients, not absolute direct MTIST gLV
coefficients.

The audit does not justify broad recovery claims. State-dependent contrasts
mean even positive/negative pair-to-rest labels summarize projections over
heterogeneous states. The A==0 direction demonstrates induced compositional/
projection effect, not an absolute false-positive interaction. The estimand
mismatch must remain explicit in release documentation.

Limitations include solver finite-step numerical validation, pooled linear
projections of nonlinear/state-dependent dynamics, and no full posterior oracle
for the 84 directions without valid retained long-run classifications. No
production R, Stan, public API, diagnostic threshold, scoring policy, K-fold,
ELPD, stacking, or scientific model changed.
