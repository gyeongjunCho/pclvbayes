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
