# pcLVbayes Roadmap

**Canonical roadmap revision:** 2026-08-31  
**Current development line:** `v0.3.0-dev`  
**Frozen baseline:** `v0.2.2`

---

# 0. Current project state

```text
Frozen baseline
    v0.2.2
        canonical Bayesian pcLV Core
        + hardened diagnostic semantics
        + pair-specific predictive contract
        + frozen 100-species HMC reference

Active development
    v0.3.0-dev
        Laplace-only audit
        -> frozen-HMC calibration
        -> conservative Laplace-to-HMC promotion
        -> staged inference

Next
    v0.3.1
        automatic staged pipeline
        + stable machine-readable H_pcLV handoff

Later v0.3.x
        canonical local-risk q_local calibration

v0.4.0
        ratio-robust cross-kingdom workflow
```

이 문서는 v0.2.2 frozen baseline부터 현재/향후 개발을 관리한다. v0.2.0–v0.2.1의 세부 hotfix, scheduler 변천, benchmark recovery history는 Git history와 benchmark provenance에 보존하며 여기서 반복하지 않는다.

---

# 1. Vision

**pcLVbayes**는 sparse, irregular longitudinal microbiome relative-abundance data에서 directed **pair-to-rest compositional Lotka–Volterra dynamics**를 reproducibly 추론하는 Bayesian package다.

orientation:

```text
source j -> target i
```

canonical local target:

```text
longitudinal compositional data
        -> deterministic full-community preprocessing
        -> target/source/rest triplet construction
        -> canonical pcLV posterior p(a_ij | D_Omega)
```

pcLV coefficient는 **triplet-ALR / pair-to-rest compositional dynamic estimand**이며 physical absolute \(A_{ij}\)와 동일하지 않다.

장기 staged architecture:

```text
full analysis universe
        |
        v
observed-support / preprocessing
        |
        v
cor_meta_resid candidate evidence
        |
        v
rough canonical Laplace approximation
        |
        v
selective full NUTS-HMC
        |
        v
diagnostics + posterior sign evidence
        |
        v
pair-specific predictive evidence where scheduled
        |
        v
stage-aware machine-readable H_pcLV
```

Laplace와 correlation screening은 **계산량을 줄이는 inference staging**이며 scientific estimand를 변경하지 않는다.

---

# 2. Scientific invariants

## 2.1 One canonical scientific model

- one canonical deterministic preprocessing contract;
- one canonical triplet-ALR / pair-to-rest estimand;
- Student-t observation likelihood;
- irregular-time OU residual process;
- subject-level random intercept;
- explicit posterior uncertainty;
- deterministic seeds;
- no silent scientific fallback.

Approximate inference stages는 이 estimand의 rough evidence만 제공하며 별도 scientific model을 정의하지 않는다.

## 2.2 Full-community rest universe

screening은 compositional denominator를 재정의하지 않는다.

```text
full analysis universe
    -> deterministic preprocessing / zero convention
    -> full-community smoothed state
    -> target/source/rest construction
```

corr/Laplace/HMC가 어떤 pair를 실행하는지는 계산 schedule을 바꿀 뿐 rest universe를 바꾸지 않는다.

## 2.3 Observed-support before numerical zero repair

원칙:

```text
observed support / detectability
        -> support eligibility
        -> zero handling as finite-log safeguard
        -> reclosure
        -> target/source/rest ALR
```

zero replacement는 존재하지 않던 biological evidence를 만드는 절차가 아니다.

## 2.4 Explicit missingness and failure

```text
indeterminate != zero
screen-out != zero
not-run != zero
engine failure != zero
unavailable != zero
```

모든 state는 machine-readable하게 분리한다.

## 2.5 Truth is post-inference validation only

simulation truth는 inference rule, preprocessing, screen threshold를 truth agreement에 맞춰 선택하는 데 사용하지 않는다.

---

# 3. Frozen baseline — v0.2.2

## Status: FROZEN

v0.2.2는 canonical HMC Core, internal simplification, hardened diagnostic semantics, publication-scale frozen HMC reference를 통합한 baseline이다.

## 3.1 Frozen execution architecture

- [x] one unordered pair = outer work unit
- [x] one worker owns both directed fits
- [x] MCMC chains serial within pair worker
- [x] K-fold work serial within pair worker
- [x] nested process-level scheduling removed
- [x] deterministic reusable pair state precomputed before dispatch
- [x] vectorized pair ALR/lag/difference preparation
- [x] bounded BLAS/OpenMP threading
- [x] deterministic direction-specific seeds
- [x] structured retry/failure handling
- [x] restart/provenance hardening

## 3.2 Frozen preprocessing / scientific Core

- [x] shared compositional helper foundation
- [x] deterministic support and zero-handling semantics
- [x] full-community rest universe preserved
- [x] `cor_meta_resid()` aligned with shared compositional helpers
- [x] canonical target/source/rest triplet construction
- [x] Student-t + irregular-time OU Stan model
- [x] subject-level random intercept
- [x] predictor/response standardization inside fit path
- [x] public output retains explicit scientific-scale semantics

## 3.3 Diagnostic semantics

pcLVbayes separates four questions.

### Sampler adequacy

HMC sampling quality:

```text
R-hat
ESS bulk
ESS tail
divergence
treedepth saturation
E-BFMI
```

### Between-chain direction stability

```text
all chain posterior-median signs agree
    -> chain_direction_stable = TRUE

chain posterior-median signs disagree
    -> interaction_indeterminate
```

historical per-chain `>=0.95` sign-mass rule은 hard gate가 아니다.

historical custom q05/q95 magnitude overlap도 hard interaction gate가 아니다.

### Pooled posterior sign evidence

edge-wise local false sign rate:

```text
LFSR <= 0.05          stringent
0.05 < LFSR <= 0.20  balanced
0.20 < LFSR <= 0.40  permissive
LFSR > 0.40           unsupported
```

LFSR은 convergence criterion이 아니다.

### Residual identifiability

residual allocation stability는 interaction reportability와 별도 diagnostic dimension이다.

current contract:

```text
interaction_reportable
    = sampler_diag_ok
      AND chain_direction_stable
```

stable interaction + unstable residual은 explicit warning과 함께 interaction-reportable일 수 있다.

## 3.4 Diagnostic engineering invariants

```text
missing diagnostics fail closed
diagnostics_unavailable explicit
non-finite coefficient sign = NA
divergence / treedepth counts must be valid finite integers
directed cross identity join is one-to-one
self identity preserves taxon + partner context
```

## 3.5 Predictive evidence contract

pcLVbayes는 pair-specific held-out predictive adequacy를 repeated K-fold fingerprint로 보존한다.

\[
\boxed{
E^{pred}_{ij}
=
\text{pair-specific repeated-K-fold predictive fingerprint}
}
\]

서로 다른 pair-to-rest model은 서로 다른 transformed outcome을 가질 수 있으므로:

```text
cross-pair stacking
cross-pair pseudo-BMA
cross-pair pseudo-BMA+
```

는 canonical pcLVbayes handoff/output에 포함하지 않는다.

가능한 predictive semantics:

```text
subject-level pointwise ELPD
ELPD total
ELPD per test observation
subject-level ELPD mean / SD
repeat/fold variability
n_test
n predictive subjects
n successful folds/repeats
predictive availability/status
```

가능하면 aggregation 전 long-form predictive artifact를 보존한다.

## 3.6 Frozen 100-species HMC reference

```text
MTIST dataset             1
taxa                      100
samples                   150
unordered pairs           4,950
directed fits             9,900

chains                    4
warmup / chain            2,000
sampling / chain          1,000
adapt_delta               0.98
max_treedepth             14
seed                      20260802
```

total runtime:

```text
236.396 h
~ 9 d 20 h 24 min
```

current diagnostic classification:

```text
converged                                  5,486
interaction_stable_residual_unstable           5
interaction_indeterminate                    228
sampler_diagnostics_failed                 4,181
------------------------------------------------
total                                      9,900
```

interaction-reportable:

```text
5,491 / 9,900 = 55.46%
```

reportable HMC LFSR:

```text
<= 0.05 stringent              3,620
0.05-0.20 balanced               850
0.20-0.40 permissive             745
> 0.40 unsupported               276
-------------------------------------
total                           5,491
```

cumulative:

```text
<= 0.05   3,620 / 5,491 = 65.9%
<= 0.20   4,470 / 5,491 = 81.4%
<= 0.40   5,215 / 5,491 = 95.0%
```

이 frozen HMC posterior가 v0.3.0 Laplace approximation의 primary reference다.

## 3.7 Frozen-reference caveats

### Predictive K-fold caveat

100-species benchmark predictive evaluation은 당시 active diagnostic gate의 영향을 받았다. later diagnostic hardening으로 rescued된 direction은 retrospective K-fold output이 없을 수 있다.

이 predictive missingness는 explicit availability state로 보존하며 Laplace calibration 때문에 HMC를 재실행하지 않는다.

### MTIST absolute A caveat

absolute MTIST `A[target, source]`는 pcLV posterior recovery의 primary truth가 아니다.

```text
pcLV posterior
    -> pair-to-rest compositional estimand

absolute A
    -> secondary cross-estimand comparison
```

---

# 4. v0.3.0-dev — Laplace staged inference

## Status: ACTIVE

primary goal:

> expensive full NUTS-HMC를 모든 candidate direction에 실행하지 않고도, HMC-important directions의 높은 recall을 유지하는 conservative pre-HMC Laplace gate를 만든다.

target workflow:

```text
candidate pair
    |
    v
canonical pair preprocessing
    |
    v
MAP + Laplace
    |
    +-------------------------+
    |                         |
failure                    succeeds
    |                         |
    v                         v
fail-open to HMC      calibrated promotion rule
                         |               |
                         v               v
                       HMC       HMC missing-by-design
```

Laplace는 HMC posterior replacement가 아니다.

---

# 5. v0.3.0 Phase A — Laplace-only implementation

## 5.1 Current state

- [x] standalone `fit_pclv_laplace.R` implementation drafted
- [ ] package documentation/export check
- [ ] CmdStanR `optimize()` / `laplace()` compatibility check
- [ ] 2–3 taxon smoke test
- [ ] deterministic seed audit
- [ ] temporary CmdStan output cleanup audit
- [ ] output schema audit
- [ ] explicit confirmation that HMC/Pathfinder/retries/K-fold are not invoked

## 5.2 Canonical preprocessing reuse

Laplace audit must reuse:

```text
same input validation
same observed-support rules
same zero-handling convention
same full-community smoothing
same pair state
same triplet-ALR definition
same lag/difference construction
same prev / dt construction
same canonical pclv.stan
same coefficient scale semantics
```

only inference engine differs.

## 5.3 Laplace inference engine

```text
canonical Stan data
    -> MAP optimize(jacobian = TRUE)
    -> Laplace approximation around MAP
    -> approximate a_ij marginal draws
    -> approximate coefficient summaries
    -> sign probability / LFSR
```

## 5.4 Laplace-native diagnostics only

Laplace must not fabricate:

```text
R-hat
ESS
divergence
treedepth
E-BFMI
chain-direction stability
```

Laplace-native audit quantities include:

```text
MAP success / return code
Laplace success / return code
finite draw fraction
mode vs mean displacement
approximate variance / interval
timing
failure stage / message
```

---

# 6. v0.3.0 Phase B — 100-species Laplace audit

- [ ] run all 4,950 unordered pairs
- [ ] obtain all 9,900 directed Laplace rows
- [ ] preserve failed/unavailable rows
- [ ] record MAP success/failure
- [ ] record Laplace/Hessian success/failure
- [ ] record approximate coefficient summaries
- [ ] record positive/negative sign probabilities
- [ ] record LFSR
- [ ] record optimize/laplace/total timing
- [ ] record worker/runtime configuration
- [ ] freeze audit corpus and provenance

minimum directed result schema:

```text
from
to
direction
n_pairs
N
S

map_ok
laplace_ok

failure_stage
failure_reason
failure_message

map_return_code
laplace_return_code

map_aij
a_mean
a_median
a_sd
a_q2.5
a_q97.5

positive_sign_probability
negative_sign_probability
dominant_sign_probability
lfsr
p_sign2

laplace_draws_requested
laplace_draws_finite
finite_draw_fraction
mode_mean_shift_sd

optimize_sec
laplace_sec
runtime_sec
```

---

# 7. v0.3.0 Phase C — Frozen-HMC calibration

join key:

```text
from + to
```

reference:

```text
frozen v0.2.2 current-policy HMC summary
```

required comparisons:

- [ ] Laplace vs HMC coefficient mean
- [ ] Laplace vs HMC coefficient median
- [ ] sign agreement
- [ ] sign-probability agreement
- [ ] LFSR association/calibration
- [ ] uncertainty-width comparison
- [ ] MAP/Laplace failure vs HMC diagnostic class
- [ ] runtime distribution

## 7.1 HMC LFSR tier preservation

explicitly evaluate:

```text
HMC <= 0.05
HMC <= 0.20
HMC <= 0.40
```

outputs:

- [ ] exact tier confusion
- [ ] cumulative cutoff retention
- [ ] false promotion
- [ ] false skip
- [ ] sign disagreement
- [ ] calibration near cutoff boundaries

## 7.2 HMC-important reference sets

initial reference sets:

```text
interaction_reportable & HMC LFSR <= 0.05
interaction_reportable & HMC LFSR <= 0.20
interaction_reportable & HMC LFSR <= 0.40
```

absolute MTIST A is not used as the primary Laplace-calibration truth.

---

# 8. v0.3.0 Phase D — Laplace-to-HMC promotion rule

Laplace screen은 recall-oriented다.

```text
false promotion
    -> extra HMC computation

false skip
    -> potentially lost scientifically useful HMC direction
```

false skip이 더 중요한 error mode다.

## 8.1 Threshold sweep

- [ ] sweep Laplace LFSR cutoff over broad range
- [ ] HMC-important recall
- [ ] HMC workload retained
- [ ] HMC workload reduction
- [ ] false-skip rate
- [ ] sign disagreement
- [ ] robustness near threshold
- [ ] Laplace failure-open behavior
- [ ] freeze rule only after empirical calibration

production decision variable:

```text
send_to_hmc
    = laplace_failure
      OR calibrated_promotion_condition
```

## 8.2 No prespecified final threshold

0.40, 0.45, 0.50 등은 calibration candidate이지 scientific constant가 아니다.

primary optimization target:

```text
maximize recall of HMC-important directions
while minimizing retained HMC workload
```

---

# 9. v0.3.0 Phase E — Production staged inference

calibration 후에만 production path에 통합한다.

- [ ] explicit Laplace gate
- [ ] fail-open on Laplace/MAP failure
- [ ] preserve missing-by-design
- [ ] preserve full rest universe
- [ ] HMC remains canonical final posterior
- [ ] deterministic staged seeds
- [ ] restart/checkpoint support
- [ ] stage-specific provenance
- [ ] runtime/memory benchmark
- [ ] compare against frozen 236.396 h HMC baseline

---

# 10. v0.3.1 — Automatic staged pipeline + stable Handoff

## Status: PLANNED

primary goal:

> `precheck -> corr -> Laplace -> selective HMC -> diagnostics/predictive -> H_pcLV`를 stable one-call workflow와 machine-readable schema로 노출한다.

target user flow:

```text
longitudinal microbiome data
    -> validation / observed support
    -> cor_meta_resid
    -> candidate selection
    -> Laplace pcLV
    -> selective full NUTS-HMC
    -> diagnostic-aware posterior summary
    -> pair-specific predictive evidence where scheduled
    -> unified directed-pair result universe
    -> H_pcLV
```

## 10.1 Canonical stage grammar

```text
PRECHECK_UNAVAILABLE

CORR_SCREEN_OUT

LAPLACE_SCREEN_OUT

MCMC_UNREPORTABLE

MCMC_REPORTABLE

MCMC_ELIGIBLE_BUT_BUDGET_SUBSAMPLED
```

stage state는 interaction absence를 의미하지 않는다.

## 10.2 Engineering / validity state

```text
corr_engine_failure
laplace_engine_failure
sampler_failure
posterior_invalid
diagnostics_unavailable
predictive_failure
```

## 10.3 Missing-by-design semantics

```text
CORR_SCREEN_OUT
    corr evidence exists
    Laplace/MCMC = missing-by-design

LAPLACE_SCREEN_OUT
    corr + Laplace evidence exists
    MCMC = missing-by-design

MCMC_ELIGIBLE_BUT_BUDGET_SUBSAMPLED
    not a failure
    not a biological zero
    explicit acquisition provenance
```

## 10.4 Stable machine-readable handoff

canonical conceptual handoff:

\[
\boxed{
H_{\mathrm{pcLV},ij}
=
[
E^{corr}_{ij},
E^{Lap}_{ij},
E^{MCMC}_{ij},
E^{pred}_{ij},
PSP/LFSR_{ij},
D^{diag}_{ij},
M^{stage}_{ij},
q^{local}_{ij}
]
}
\]

v0.3.1 may initially ship before q_local is calibrated; in that case q_local must be explicit `unavailable/not-yet-calibrated`, not fabricated.

## 10.5 Coefficient-scale contract

- [ ] audit predictor/response scaling
- [ ] audit back-transform
- [ ] freeze exposed coordinate-scale semantics
- [ ] add `coefficient_scale_version`
- [ ] guarantee Laplace/HMC comparable exposed scale
- [ ] preserve zero-handling / ALR convention version

## 10.6 Predictive evidence contract

- [ ] preserve pair-specific repeated-K-fold ELPD
- [ ] preserve support metadata
- [ ] preserve predictive availability
- [ ] preserve pointwise long-form artifact where feasible
- [ ] prohibit cross-pair stacking/pseudo-BMA in canonical output

## 10.7 Stable result identity

- [ ] cross directed identity one-to-one
- [ ] self identity includes taxon + partner
- [ ] full directed-pair result universe preserved
- [ ] stage/failure/provenance machine-readable
- [ ] stable output-schema version

---

# 11. Later v0.3.x — Canonical local-risk calibration

## Status: PLANNED

LFSR와 q_local을 구분한다.

```text
LFSR
    pooled posterior sign uncertainty

q_local
    calibrated risk that the reportable local pcLV claim is wrong
    relative to canonical a*_Omega truth
```

local truth hierarchy:

\[
A_{ij}
\rightarrow
a^*_{ij;\Omega}
\rightarrow
\widehat a_{ij}.
\]

local claim-failure target:

\[
Y^{local}_{ij;\Omega}
=
I[
|\widetilde a^*_{ij;\Omega}|\le\delta_*
\ \lor\
\operatorname{sign}(a^*_{ij;\Omega})
\neq
\operatorname{sign}(\widehat a_{ij})
].
\]

calibrator:

\[
\boxed{
q^{local}_{ij}
=
P(
Y^{local}_{ij}=1
\mid
E^{local}_{ij},
R^{local}_{ij}=1
)
}
\]

required work:

- [ ] freeze canonical \(a^*_\Omega\) oracle
- [ ] freeze oracle scale/back-transform contract
- [ ] define local claim-failure labels
- [ ] ecosystem-level leakage-safe split
- [ ] OOF local-risk features for downstream training
- [ ] calibration split / frozen mapping
- [ ] version q_local artifact
- [ ] validate risk-coverage
- [ ] add q_local to stable H_pcLV schema

TASR \(A\to a^*_\Omega\) is not pcLVbayes local-risk target.

---

# 12. Relationship to pcLVsp / pcLV-RA

pcLVbayes responsibility ends at **qualified local evidence**.

```text
pcLVbayes
    Local Bayesian Guard
        |
        v
    H_pcLV
        |
        v
pcLV-RA / pcLVsp
    global proportional reconstruction
```

pcLVbayes does not:

```text
recover physical absolute A
reconstruct ecosystem-wide G ≈ X A
predict structural C/R/I/L/N state
compute q_prop
perform RA OOD inference
auto-flip pcLV signs
```

downstream pcLV-RA target:

\[
\boxed{
\widehat G^{RA}\approx X_eA
}
\]

이며 physical absolute \(A\) recovery가 아니다.

pcLVbayes는 downstream RA가 사용할 **derived local structural/reliability view**를 제공하지만 independent biological information을 추가한다고 주장하지 않는다.

---

# 13. v0.4.0 — Ratio-robust cross-kingdom workflow

## Status: PLANNED

primary target: absolute bacterial:fungal biomass ratio가 정확히 알려지지 않은 bacteria–fungi longitudinal interaction inference.

initial sensitivity design:

```text
within-kingdom bacterial RA + fungal RA
        |
        +--> joint mixture 5:5  primary
        +--> joint mixture 7:3  sensitivity
        `--> joint mixture 3:7  sensitivity
```

required work:

- [ ] principled scaling into joint composition
- [ ] prespecified primary mixture
- [ ] prespecified alternative mixtures
- [ ] interaction stability across mixture assumptions
- [ ] ratio-sensitive flag
- [ ] full joint rest universe preserved within each scenario
- [ ] build on validated staged pipeline
- [ ] bacteria–fungi longitudinal validation

ratio robustness는 mathematical independence claim이 아니다.

---

# 14. v1.0 — Stable pcLVbayes release

expected capabilities:

- conservative Bayesian pair-to-rest inference
- sparse/irregular longitudinal support
- Student-t observation model
- irregular-time OU residual
- explicit sampler/direction/residual diagnostic semantics
- correlation → Laplace → selective NUTS staging
- deterministic seeds and provenance
- restartable resource-safe execution
- stable result schema
- pair-specific predictive evidence
- calibrated q_local
- stable H_pcLV handoff
- network-ready local posterior summaries
- optional ratio-robust cross-kingdom workflow
- clear claim boundary

---

# 15. Known limitations

pcLVbayes estimates pair-to-rest compositional interaction dynamics.

따라서:

- relative abundance alone에서 physical absolute \(A_{ij}\)를 직접 recover하지 않는다;
- supported pcLV sign은 absolute biological facilitation/inhibition의 직접 증명이 아니다;
- observational posterior는 causal experimental proof가 아니다;
- observation domain에 따라 pseudo-true local effect가 달라질 수 있다;
- screening은 computation을 줄일 뿐 missing biological truth를 자동으로 복원하지 않는다;
- high-priority local interactions은 가능한 경우 isolate/co-culture/perturbation/host validation이 필요하다.

---

# 16. Explicit exclusions from the Bayesian Core

- full absolute interaction matrix recovery
- automatic posterior sign flipping
- simulation-trained reconstruction inside Bayesian Core
- full joint 100-species NUTS model
- converting unavailable/indeterminate/not-run to zero
- preprocessing chosen by agreement with simulation truth
- cross-pair stacking/pseudo-BMA as canonical predictive evidence
- claiming observational pcLV interaction as causal proof

---

# 17. Versioned contracts

future compatibility를 위해 최소 다음을 version한다.

```text
pcLVbayes package version
scientific estimand version
handoff schema version
coefficient-scale / back-transform version
zero-handling / ALR convention version
predictive-evidence version
diagnostic-policy version
q_local calibration artifact version
```

sampler backend가 바뀌어도 estimand/preprocessing/handoff semantics가 유지되면 downstream pcLVsp compatibility를 최대한 유지한다.

---

# 18. Release / provenance contract

release/frozen benchmark manifest에는 최소:

```text
package git commit / tag
input hash / dataset hash
Stan model hash/version
zero-handling / preprocessing version
coefficient-scale version
diagnostic-policy version
predictive-evidence version
seed policy
CmdStan / CmdStanR version
worker/resource policy
benchmark artifact checksum
result-schema version
```

을 남긴다.

---

# 19. Roadmap maintenance rule

다음이 바뀔 때 ROADMAP을 갱신한다.

1. canonical triplet-ALR estimand
2. observed-support / zero / rest-universe preprocessing contract
3. canonical Stan scientific model
4. diagnostic/reportability semantics
5. LFSR/PSP semantic contract
6. Laplace→HMC staged inference policy
7. stage/failure/missing-by-design grammar
8. pair-specific predictive ELPD contract
9. coefficient-scale/back-transform contract
10. q_local truth/calibration contract
11. machine-readable H_pcLV schema
12. cross-kingdom scientific interpretation contract

일반 scheduler tuning, worker 수, minor private refactor는 maintenance trigger가 아니다.

---

# 20. Current execution checklist

## v0.2.2

```text
[x] frozen canonical HMC baseline
[x] frozen 100-species HMC reference
[x] hardened current diagnostic semantics
[x] predictive K-fold caveat recorded
[x] provenance frozen
```

## v0.3.0-dev

```text
[x] fit_pclv_laplace.R initial draft

[ ] package/API smoke
[ ] 2–3 taxon Laplace smoke test
[ ] 100-species / 9,900-direction Laplace audit
[ ] frozen HMC join
[ ] coefficient/sign/LFSR comparison
[ ] threshold sweep
[ ] HMC-important recall vs workload
[ ] failure-open policy
[ ] promotion rule freeze
[ ] production staged integration
```

## v0.3.1

```text
[ ] stage-aware one-call pipeline
[ ] full directed-pair result universe
[ ] predictive fingerprint
[ ] coefficient-scale version
[ ] stable H_pcLV schema
```

## later v0.3.x

```text
[ ] canonical oracle
[ ] q_local
[ ] OOF/calibration contract
```

---

# Canonical summary — 2026-08-31

```text
pcLVbayes scientific target
    pair-to-rest compositional local dynamic estimand

v0.2.2
    frozen canonical HMC Core + 100-species reference

v0.3.0-dev
    approximate local evidence for HMC staging
    NOT HMC replacement

v0.3.1
    stable stage-aware H_pcLV handoff

later v0.3.x
    local claim-risk q_local

pcLV-RA
    downstream proportional G ≈ X A reconstruction
    outside pcLVbayes Bayesian Core
```

\[
\boxed{
\text{Decompose / preprocess}
\rightarrow
\text{corr evidence}
\rightarrow
\text{Laplace evidence}
\rightarrow
\text{selective HMC}
\rightarrow
\text{diagnostics + predictive evidence}
\rightarrow
H_{\mathrm{pcLV}}
}
\]

