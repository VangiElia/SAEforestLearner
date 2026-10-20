#' 'SAEforestLearner' - Scalable Mixed Effects Random Forests for Small Area Estimation
#'
#' \pkg{SAEforestLearner} is a scalable fork of the package \pkg{SAEforest} 1.0.1 (Krennmair &
#' Frink), whose functionality it keeps. The package \pkg{SAEforest} promotes the use of Mixed Effects Random Forests (MERFs) for
#' applications of Small Area Estimation (SAE). The package effectively combines functions for
#' the estimation of regionally disaggregated linear and nonlinear economic and inequality
#' indicators using survey sample data. Estimated models increase the precision of direct
#' estimates from survey data, combining unit-level and aggregated population level covariate
#' information from census or register data. Apart from point estimates, MSE estimates for requested
#' indicators can be easily obtained. The package provides procedures to facilitate the analysis
#' of model performance of MERFs and visualizes predictive relations from covariates and variable
#' importance. Additionally, users can summarize and map indicators and corresponding measures of
#' uncertainty. Methodological details for the functions in this package are found in Krennmair & Schmid (2022),
#' Krennmair et al. (2022a) and Krennmair et al. (2022b).
#'
#' @details
#' This package includes a main function \code{\link{MERFranger}} that is wrapped in
#' \code{\link{SAEforest_model}} for an improved SAE workflow.
#' Each function produces an object inheriting requested results of regionally disaggregated point
#' and uncertainty estimates. Additionally, statistical information on model fit and variable
#' importance is accessible through generic functions such as a summary (\code{\link{summary.SAEforest}})
#' or a class-specific plot function (\code{\link{plot.SAEforest}}). For a full documentation of
#' objects of class \code{SAEforest} see \code{\link{SAEforestObject}}. An overview of all currently
#' provided functions within this package can be seen with \code{help(package="SAEforestLearner")}.
#'
#' @references
#' Krennmair, P., & Schmid, T. (2022). Flexible Domain Prediction Using Mixed Effects
#' Random Forests. Journal of Royal Statistical Society: Series C (Applied Statistics) (forthcoming).
#'
#' Krennmair, P., & Würz, N. & Schmid, T. (2022a). Analysing Opportunity Cost of Care Work using
#' Mixed Effects Random Forests under Aggregated Census Data.
#'
#' Krennmair, P., & Schmid, T & Tzavidis, Nikos. (2022b). The Estimation of Poverty Indicators Using
#' Mixed Effects Random Forests. Working Paper.
#'
#' Box, G. E. P., & Cox, D. R. (1964). An Analysis of Transformations. Journal of the Royal
#' Statistical Society: Series B (Methodological), 26(2), 211–252.
#'
#' Molina, I., & Rao, J. N. K. (2010). Small Area Estimation of Poverty Indicators. Canadian
#' Journal of Statistics, 38(3), 369–385.
#'
#' @section Scalable extension (SAEforestLearner):
#' Domain means and totals with nonparametric bootstrap MSE for very large, e.g. raster-based,
#' populations: see \code{\link{SAEpopulation}}, \code{\link{SAEforest_model}} (argument
#' \code{total_factor}), \code{\link{SAEaggregate}} and \code{\link{predict_flat_forest}}.
#' Optional extensions of \code{\link{SAEforest_model}} for domain means: a Box-Cox transformed
#' MERF with empirical-best-type back-transformation (arguments \code{transform},
#' \code{bc_lambda}, \code{bc_shift}; Box & Cox, 1964; Molina & Rao, 2010), the pool of the
#' bootstrap unit errors (\code{boot_error_pool}) and heteroscedastic bootstrap variants
#' (\code{boot_hetero}).
#'
#' @aliases SAEforest
#' @useDynLib SAEforestLearner, .registration = TRUE
#' @importFrom Rcpp sourceCpp
"_PACKAGE"
