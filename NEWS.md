# SAEforestLearner 2.0.0

Renamed from SAEforestArrow (2026-10-20): the package is a fork of SAEforest 1.0.1 with the scalable
engine and the learner interface below.

Scalable fork of SAEforest 1.0.1 for very large (raster) populations.

## New

* `SAEpopulation()`: population from a raster (all layers or a subset as covariates, a mask
  layer, domains from polygons rasterized by cell centre). Valid cells are extracted once,
  in parallel and block-wise, into an on-disk `fst` cache that is streamed in chunks.
  Also wraps in-memory data.frames and re-opens existing caches.
* `SAEforest_model()` accepts `SAEpop` objects or rasters (`pop_options`) as `pop_data`,
  and gains `total_factor` (domain totals and their MSE), `cores`, `chunk_rows`,
  `forest_batch`, `seed`, `verbose` and `engine`.
* Scalable engine for domain means (`meanOnly = TRUE`, unit-level data), used by default:
  memory no longer grows with population size x `B`. Point estimates are identical to the
  original implementation; the nonparametric MSE bootstrap follows the same distribution.
  - only domain sums of forest predictions are kept, the population is streamed;
  - random effects via `ranef()` look-up (identical to `predict.merMod` on new data);
  - bootstrap unit errors are drawn in C++ one replicate at a time and summed per domain;
  - bootstrap MERFs are fitted in parallel on the (small) bootstrap samples only.
* `flatten_forest()`, `predict_flat_forest()`: multi-threaded C++ prediction of ranger
  regression forests, bit-identical to `ranger::predict()` and ~10-15x faster with many
  covariates (ranger sorts all covariates at each prediction call).
* `summarize_indicators()` reports `Total` besides `Mean` when available.
* `SAEaggregate()`: means and totals of any aggregation of domains (e.g. municipalities to
  provinces, regions, country) with exact bootstrap MSEs, computed from the bootstrap
  replicates stored in the model (covariances between domains included). It is identical to
  re-running the model with `aggregate_to`.
* Box-Cox transformed MERF (scalable engine): arguments `transform = "boxcox"`, `bc_lambda`,
  `bc_shift`, `bt_nodes`, `bt_grid`.
  - The MERF is fitted to the transformed target.
  - Domain means are empirical-best-type predictions on the original scale: expected
    back-transformed values given the predicted random effect, integrating the unit error
    and the remaining random-effect uncertainty by Gauss-Hermite quadrature (Molina & Rao,
    2010).
  - The nonparametric bootstrap runs on the transformed scale; true and estimated
    bootstrap domain means are back-transformed unit by unit.
  - With `bc_lambda = 1`, `bc_shift = 1` (identity) the results equal `transform = "none"`
    with the same bootstrap error pool, as long as no transformed value falls below
    `-bc_shift`, where the back-transformation is truncated.
* `boot_error_pool` (`"auto"`, `"domain"`, `"all"`): pool of the unit errors of the
  bootstrap populations.
  - `"domain"` (the original scheme) draws the errors of sampled domains from their own,
    within-domain centred residuals; `"all"` draws them from all residuals.
  - `"auto"` uses `"domain"` without a transformation (results unchanged) and `"all"` with
    Box-Cox. After a back-transformation the true domain means depend on the whole error
    distribution, which the few own residuals of a sampled domain do not represent.
* `boot_hetero` (`"none"`, `"unit"`, `"unit_re"`): optional heteroscedastic variants of the
  bootstrap (unit errors with SD proportional to a power of the level; optionally random
  effects proportional to the domain level). `"none"` (default) is the original bootstrap.

* Learner interface (scalable engine). The fixed part of the mixed model can be any supervised
  learner: `SAEforest_model(learner = )` takes an object created by `sae_learner()` (functions
  `fit`, `predict`, optional `oof`, `pack`, `unpack`, `info`); `learner_ranger()` (default) is
  the ranger forest of the original package. Learners without native out-of-sample predictions
  (out-of-bag) get out-of-fold predictions by K-fold cross-fitting (`oof_folds`) inside the MERF
  iterations and in the bias correction of the error variance. `MERFranger()` is now a wrapper of
  the generic algorithm and returns the out-of-sample predictions (`OOBpred`) and the learner.
  Point estimates with `learner_ranger()` are identical to those of version 2.0.0; the bootstrap
  random streams differ by two draws (`predict.ranger`, which drew a seed from the R random
  number generator in two places, is no longer called on the sample), so MSE estimates with a
  given seed are not reproduced bit by bit (same distribution). Predefined learners:
  `learner_xgboost()`, `learner_lightgbm()`, `learner_glmnet()`, `learner_gam()` (suggested
  packages).

* Convergence of the MERF iterations: `em_criterion = "robust"` (default with learners other
  than ranger) compares the means of consecutive windows of the log-likelihood, stops after
  `patience` iterations without improvement and keeps the best iteration; `"original"` (default
  with `learner_ranger()`) is the criterion of the original package. `learner_lm()`: linear
  learner with exact leave-one-out predictions, which makes the MERF the nested-error linear
  mixed model. `tune_learner()` and `spatial_folds()`: tuning of a learner's parameters by
  spatial-block cross-validation of the mixed model, on the original scale with Box-Cox.
  New vignette "SAEforestLearner: differences from SAEforest and the learner interface" (eusilc
  example with ranger, xgboost and a user-defined support vector learner).

## Differences from SAEforest 1.0.1

* `MERFranger()`: a unit that is in-bag in every tree (possible with few trees) has no
  out-of-bag prediction. `ranger` returns `NaN`, `lmer` then drops the unit, and the fit stopped
  with "Missing data in dependent variable". Such units now get the prediction of the whole
  forest. Results are unchanged whenever all out-of-bag predictions exist (always in practice
  with a few hundred trees).

* Block bootstrap of errors: in the original code, a domain with a single sample unit
  calls `sample(x, ...)` with a length-one numeric `x`, so that its errors are drawn from
  `1:x` instead of `x`. The scalable engine draws the errors in C++ and is not affected.
  `engine = "legacy"` keeps the original code unchanged, including this behaviour, for
  comparison.

# SAEforest 1.0.1

# - remove dependencies to rgeos, maptools and sp
