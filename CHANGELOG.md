# CHANGELOG

## v0.2 (In Progress)

### Scientific estimand clarification
- Clarify that `a_ij` is a directed pair-to-rest dynamic coefficient, distinct from absolute direct-gLV `A[i,j]`.
- Record the oracle audit findings: transformed-oracle agreement 6/6, absolute-A cross-estimand agreement 3/6, 44/90 state-dependent contrasts, 16/90 absolute-zero induced transformed effects, and the `species_7 -> species_4` smoothing-sensitive case.
- Separate transformed-oracle validation from absolute-A comparison in benchmark reporting.
- No posterior-model, computation, public-API, or result-schema change.

### Planned
- Simplify public API
- Remove obsolete alpha parameters
- Eliminate silent fallback logic
- Reduce conditional branching
- Introduce validated intermediate objects
- Preserve Core OU model
