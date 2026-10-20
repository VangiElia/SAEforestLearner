# SAEforestLearner

Scalable Mixed Effects Random Forests (MERFs) for Small Area Estimation.

`SAEforestLearner` is a fork of [SAEforest](https://github.com/krennpa/SAEforest) 1.0.1
(Krennmair & Frink) that makes MERF-based small area estimation of domain **means and
totals with nonparametric bootstrap MSE** feasible for very large unit-level populations,
e.g. tens of millions of raster cells with hundreds of covariates, on a single workstation.

Main additions with respect to SAEforest:

* `SAEpopulation()` builds the population directly from a raster: a mask layer
  selects the population cells, domains are given as polygons (rasterized by cell centre)
  and all or selected layers are the covariates. Valid cells are extracted once, in
  parallel and block-wise, into an on-disk `fst` cache that is streamed in chunks.
* A scalable estimation engine, used by default by `SAEforest_model()` for domain means.
  Its memory use does not grow with population size × bootstrap replications. Point
  estimates are identical to the original implementation, and the bootstrap follows the
  same scheme (Krennmair & Schmid, 2022).
* Domain totals (`total_factor`) and their MSE; `SAEaggregate()` for exact bootstrap MSEs
  of any aggregation of domains (e.g. municipalities to provinces, regions, country).
* `flatten_forest()` / `predict_flat_forest()`: multi-threaded C++ prediction of
  `ranger` regression forests, bit-identical to `ranger::predict()` but much faster
  with many covariates.
* Optional Box-Cox transformed MERF (`transform = "boxcox"`, `bc_lambda`, `bc_shift`) for
  skewed, heteroscedastic targets. Domain means are empirical-best-type back-transformations
  (expected back-transformed value given the predicted random effect, Gauss-Hermite
  quadrature), and the bootstrap runs on the transformed scale with unit-level
  back-transformed true values. With a transformation, bootstrap unit errors are drawn from
  all residuals (`boot_error_pool`).
* Optional heteroscedastic bootstrap variants (`boot_hetero`). The default reproduces the
  bootstrap of the original package.
* Any supervised learner in place of the random forest (`learner`). The MERF algorithm
  alternates the learner and the random-effects model exactly as with the forest; learners
  without out-of-bag predictions get out-of-fold predictions by cross-fitting (`oof_folds`).
  Predefined: `learner_ranger()` (default, identical results to the forest-only package),
  `learner_xgboost()`, `learner_lightgbm()`, `learner_glmnet()`, `learner_gam()`; user-defined
  learners through `sae_learner()`.

## Installation

The package contains C++ code: on Windows [Rtools](https://cran.r-project.org/bin/windows/Rtools/)
matching your R version is required, on macOS the Xcode command line tools.

```r
# from GitHub
install.packages("remotes")
remotes::install_github("VangiElia/SAEforestLearner")

# or from a local copy
remotes::install_local("path/to/SAEforestLearner")
# or, for development
devtools::load_all("path/to/SAEforestLearner")
```

## Example: forest growing stock per municipality from a raster

```r
library(SAEforestLearner)
library(sf)

# population: forest cells (layer "forest_mask" == 1) of a multi-layer covariate raster,
# municipalities by cell centre; built once, then re-opened with SAEpopulation(cache_dir)
pop <- SAEpopulation("predictors_70m.tif",
                     domains = st_read("municipalities.shp"), domain_field = "municipality",
                     mask = "forest_mask", mask_value = 1, exclude = c("id_layer"),
                     extra_fields = c("region", "province"),
                     cache_dir = "cache/population", cores = 12)

# sample: plot-level target (m3/ha) and the covariates of the cell containing each plot
model <- SAEforest_model(Y = smp$volume, X = smp[, pop$predictors], dName = "municipality",
                         smp_data = smp, pop_data = pop,
                         MSE = "nonparametric", B = 100, B_adj = 100,
                         total_factor = 0.49,          # 70 m cells = 0.49 ha -> totals in m3
                         num.trees = 500, seed = 1, cores = 40)

summarize_indicators(model, MSE = TRUE, CV = TRUE)     # Mean and Total per municipality
SAEaggregate(model, unique(smp[, c("municipality", "region")]), name = "region")

# the same with a Box-Cox transformed MERF (lambda chosen beforehand, e.g. by profile likelihood)
model_bc <- SAEforest_model(Y = smp$volume, X = smp[, pop$predictors], dName = "municipality",
                            smp_data = smp, pop_data = pop,
                            MSE = "nonparametric", B = 100, B_adj = 100, total_factor = 0.49,
                            transform = "boxcox", bc_lambda = 0.35, bc_shift = 1,
                            num.trees = 500, seed = 1, cores = 40)

# the same with another learner of the fixed part (boosting; the ... are its parameters)
model_xgb <- SAEforest_model(Y = smp$volume, X = smp[, pop$predictors], dName = "municipality",
                             smp_data = smp, pop_data = pop,
                             MSE = "nonparametric", B = 100, B_adj = 100, total_factor = 0.49,
                             transform = "boxcox", bc_lambda = 0.35, bc_shift = 1,
                             learner = learner_xgboost(nrounds = 300, eta = 0.05),
                             seed = 1, cores = 40)

# a user-defined learner: self-contained fit / predict functions (parameters in args)
lrn_lm <- sae_learner("lm",
  fit = function(x, y, args, threads) { d <- as.data.frame(x); d$.y <- y; stats::lm(.y ~ ., d) },
  predict = function(model, x, threads) as.numeric(stats::predict(model, as.data.frame(x))))
```

The original SAEforest functionality (nonlinear indicators, aggregated covariates,
`engine = "legacy"`) is unchanged for data.frame populations. The vignette shipped with the
package is the original SAEforest vignette: it describes this original functionality, not
the scalable extensions documented above and in the help pages.

## References

Krennmair, P., & Schmid, T. (2022). Flexible domain prediction using mixed effects random
forests. *Journal of the Royal Statistical Society: Series C*, 71(5), 1865-1894.

Krennmair, P., Würz, N., & Schmid, T. (2022). Analysing opportunity cost of care work using
mixed effects random forests under aggregated census data. arXiv:2204.10736.

Box, G. E. P., & Cox, D. R. (1964). An analysis of transformations. *Journal of the Royal
Statistical Society: Series B*, 26(2), 211-252.

Molina, I., & Rao, J. N. K. (2010). Small area estimation of poverty indicators. *Canadian
Journal of Statistics*, 38(3), 369-385.

Mendez, G., & Lohr, S. (2011). Estimating residual variance in random forest regression.
*Computational Statistics & Data Analysis*, 55(11), 2937-2950.

## License

GPL (>= 2), as SAEforest. Original authors of SAEforest: Patrick Krennmair and Nicolas
Frink. Scalable extensions: Elia Vangi.
