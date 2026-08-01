# Profiling baseline: MTIST 361

Measured 2026-08-02 at pcLVbayes `a352313` on Ubuntu 24.04, R 4.6.1,
CmdStan 2.36.0, an AMD Ryzen 7 8700G (16 logical cores), and 31,935,200 kB RAM.
The fixture had three taxa and used one chain, 50 warmup, 50 sampling draws,
K=2, R=1, no retries, no Pathfinder, and one versus two outer workers. The
unchanged package was installed into an isolated temporary library; its 83.879
second installation cost and model-load cost were excluded from fit timings.

## Wall time

| Measurement | Calls | Total seconds | Per call |
|---|---:|---:|---:|
| Sequential end-to-end fit | 1 | 83.465 | 83.465 |
| Two-worker end-to-end fit | 1 | 56.959 | 56.959 |
| One-direction, two-fold K-fold probe | 1 | 8.121 | 8.121 |
| Global CV spline smoothing | 5 | 0.160 | 0.0320 |
| Canonical pair preprocessing | 5 | 0.025 | 0.0050 |
| K-fold split construction | 100 | 0.009 | 0.00009 |
| CmdStan JSON serialization probe | 100 | 0.073 | 0.00073 |

Sequential Rprof sampled 5.68 seconds of active parent-R execution. The other
77.785 seconds (93.2% of wall time) were spent while R was waiting at external
or non-profiled boundaries. This is a **combined CmdStan JSON/process startup,
warmup, sampling, and output boundary**, not a measurement of pure NUTS time.
The current cleaned-fit lifecycle does not preserve sufficiently reliable
CmdStan timing metadata to split those components after fitting, so no pure
sampling percentage is claimed.

## R-side hotspots

Within the 5.68 seconds sampled while sequential parent R was active:

| Inclusive stack | Seconds | Parent-R active time | End-to-end wall time |
|---|---:|---:|---:|
| `.sample_with_retry()` | 5.39 | 94.9% | 6.5% |
| `.summarise_diag()` | 4.57 | 80.5% | 5.5% |
| `posterior::summarise_draws()` | 4.43 | 78.0% | 5.3% |
| Canonical deterministic pair preprocessing (direct) | 0.005/call | — | <0.01%/call |

Call-stack percentages are inclusive and must not be summed. Diagnostic draw
summarization, rather than ALR/index loops, is the measured R-side hotspot.

## K-fold attribution

The two-fold one-direction probe took 8.121 seconds (~4.06 seconds/fold).
Split construction (0.00009 seconds/call) and direct JSON serialization
(0.00073 seconds/call) were negligible at this size, so repeated CmdStan model
execution dominates the measured fold cost. Both intentionally short fold fits
failed sampler diagnostic gates; therefore successful Kalman scoring time was
not separately measured. Fold preparation, process startup, sampling, failed
diagnostic extraction, and aggregation remain a combined phase.

## Sequential versus parallel

Two outer workers reduced wall time from 83.465 to 56.959 seconds: 26.506
seconds (31.8%) and a 1.47x speedup. Scientific result signatures were
identical. Relative to the unattainable ideal of half the sequential time,
15.227 seconds is an upper bound on combined worker startup, serialization,
scheduling, task imbalance, and worker-local R overhead. These contributions
were not individually observable, and worker CPU is absent from parent Rprof.

## Memory and temporary disk

The parent process rose from 823,984 kB RSS/HWM before sequential fitting to
1,028,996 kB afterward and 1,033,788 kB after parallel fitting: a measured
parent-process high-water increase of 209,804 kB. Worker RSS was not included.
Major retained objects were small: study 89,840 bytes, smoothed matrix 17,512,
pair input 10,408, Stan list 6,384, and each final fit result 112,504 bytes.

The process temporary directory changed from 123 files / 6,237,427 bytes to
129 files / 6,302,711 bytes after sequential fitting, then 130 files /
6,302,717 bytes after parallel fitting. This snapshot includes R session and
future infrastructure, so the 65,290-byte net increase is not attributed to
package leakage. Existing lifecycle tests verify package-owned CmdStan outputs
return to baseline. In-run peak disk use and full process-tree peak RSS were not
captured and are explicitly recorded as unavailable (`NA`).

## Ranked candidates

1. **Number and duration of CmdStan fits.** Combined external boundary was
   77.785 seconds (93.2% of sequential wall time); reducing scientifically
   unnecessary fits or sampler cost has the largest plausible impact, but the
   correctness risk is high and no change is made here.
2. **Posterior diagnostic summarization.** `posterior::summarise_draws()` used
   4.43 seconds (5.3% end-to-end; 78.0% of active parent R). Restricting or
   consolidating repeated summaries is a moderate-impact, moderate-risk R
   candidate. Rcpp is not the first mechanism because the hotspot is library
   summary machinery, not a demonstrated package loop.
3. **Outer scheduling/startup.** Parallel execution saved 31.8%, but the
   15.227-second gap from ideal includes task imbalance and startup. Reusing
   workers or reducing exports is plausible, subject to measurement with
   worker-level instrumentation.
4. **CV spline and pair preprocessing.** At 0.032 and 0.005 seconds/call they
   are low-impact for this three-species fixture. Rcpp is not supported by the
   current end-to-end evidence.
5. **Split creation and JSON writing.** Sub-millisecond direct costs are not
   optimization priorities here. Process startup cannot be inferred from JSON
   writing alone.

## Limitations

This is an attribution smoke profile, not a scientific fit. Short chains were
diagnostically poor, prevented canonical main-fit K-fold entry, and made the
separate K-fold probe fail before successful scoring. Rprof records active
parent-R time, not child CmdStan or worker CPU. Parent HWM excludes workers;
temporary disk is before/after rather than peak. The measurements support a
Stan/process-dominated conclusion but do not distinguish pure NUTS from process
startup or output handling.

## Posterior-summary optimization (2026-08-02)

The retained scientific order remains sampling/retry, sampler-only attempt
diagnostics, final-attempt selection, one retained posterior summary bundle,
Bayesian sign evidence, K-fold ELPD, stacking, and final assembly. Attempt-level
retry diagnostics no longer load posterior variables or call
`posterior::summarise_draws()`; they read sampler diagnostics once. The retained
attempt materializes selected draws once, adds R-hat/ESS once, and reuses an
immutable bundle for coefficients, quantiles, ν, sign probabilities, and
multi-chain residual summaries.

The repeatable 4-chain × 2,000-draw component benchmark (15 repetitions) was:

| Path | Median seconds | `summarise_draws()`/invocation | posterior draw loads/invocation |
|---|---:|---:|---:|
| Legacy attempt diagnostics | 0.091 | 1 | 3 |
| Optimized attempt diagnostics | 0.001 | 0 | 0 |
| Optimized retained final summary | 0.099 | 1 | 1 |

Thus rejected attempts avoid essentially all scientific-summary work. A
no-retry retained fit still performs one exact posterior convergence summary,
as scientifically required, but eliminates duplicate CmdStan draw loading.
The private bundle was 1,168 bytes in the synthetic benchmark and does not
retain the CmdStan fit or a duplicate draw matrix.

Using the identical MTIST 361 configuration, sequential wall time decreased
from 83.465 to 78.282 seconds (-5.183 seconds, -6.2%). Active parent-R time fell
from 5.68 to 1.20 seconds (-78.9%); `summarise_draws()` disappeared from the top
20 Rprof stacks (previously 4.43 seconds). The two-worker run decreased from
56.959 to 55.121 seconds (-3.2%), and the two-fold failed-diagnostic probe from
8.121 to 7.252 seconds (-10.7%). Sequential and parallel scientific signatures
remained identical.

For this profile, `summarise_draws()` calls fell from eight to six: six retained
main attempts still required exact convergence summaries, while two rejected
fold attempts no longer did. Posterior draw loads/conversions fell from 30
(three attempt loads plus one retained load per main direction, and three per
rejected fold) to six (one retained load per main direction), an 80% reduction.

Parent-process HWM after the parallel run decreased from 1,033,788 to 1,028,504
kB (-5,284 kB); final fit objects decreased from 112,504 to 109,464 bytes. These
single smoke measurements include normal run-to-run process noise, but show no
material memory regression. Numerical tests preserve means, medians, exact R
quantiles, sign probabilities, LFSR, R-hat, ESS, all diagnostic classes and
reasons, one-chain `NA` agreement, and downstream weight immutability.

The maximum end-to-end impact remains modest because CmdStan/process execution
dominates. No Rcpp or new dependency was introduced.
