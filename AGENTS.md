# AGENTS.md

## Project purpose

pcLVbayes is an experimental R package for Bayesian pairwise compositional
Lotka–Volterra modeling of longitudinal microbiome data.

The current task is limited to the v0.2 Core Reduction and Logical
Refactoring described in `TODO.md`.

## Canonical files

- `R/fit_pclv_bayes.R`
  - canonical fitting workflow
  - preserve as the main fitting entry point

- `R/pclv_helpers.R`
  - internal computational helpers for the canonical fitting workflow

- `R/summarize_bayes_pclv.R`
  - conservative post-inference summarization
  - keep separate from fitting

- summary helper files
  - internal support for Bayesian sign-error filtering and model weighting

## Excluded files

- `R/fit_pclv_bayes2.R`
  - unfinished cross-kingdom prototype
  - do not refactor, import from, delete, or use as a design reference
  - ignore unless explicitly instructed otherwise

## Scientific invariants

Do not change the scientific meaning of the Core model without explicit
approval.

The Core model includes:

- pair-to-rest ALR transformation
- delta ALR divided by sampling interval as the response
- lagged self predictor
- lagged partner predictor
- irregular-time OU residual structure
- Bayesian posterior inference
- posterior sign probability
- conservative post-inference filtering and model weighting

## Refactoring principles

- Delete obsolete alpha-stage options rather than hiding all of them in a
  control object.
- Remove silent fallbacks that substitute a different scientific model.
- Prefer fail-fast validation over runtime repair.
- Make invalid internal states unrepresentable where practical.
- Reduce nested conditional logic.
- Do not remove essential boundary validation.
- Do not silently skip malformed data or failed fits.
- Return structured failure reasons.
- Preserve fitting and conservative summarization as separate layers.
- Do not add new statistical features during v0.2.

## Change policy

For each task:

1. inspect the relevant code;
2. explain the proposed change;
3. identify behavior that will be deleted or preserved;
4. add or update tests;
5. implement the smallest coherent change;
6. run checks;
7. summarize changed files and remaining risks.

Do not perform broad repository-wide rewrites in a single task.

## Git policy

- Never commit directly to `main`.
- Work on a dedicated branch.
- Keep commits small and logically scoped.
- Do not rewrite published Git history.
- Do not merge pull requests.
- Open a PR for owner review.
- Stop when the assigned TODO item is complete.

## Verification

At minimum, run the applicable package checks and tests after changes.

Preferred commands:

```bash
Rscript -e 'devtools::document()'
Rscript -e 'devtools::test()'
R CMD build .
