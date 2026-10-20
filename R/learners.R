# Learner interface: the fixed part of the mixed model can be any supervised learner -------
#
# The MERF algorithm (Hajjem et al., 2014; Krennmair & Schmid, 2022) alternates a learner fitted
# to Y minus the random effects and a linear mixed model of the random effects with the learner's
# out-of-sample predictions as offset. Nothing in the scalable engine depends on the learner
# being a random forest: it only needs fit(), predict(), out-of-sample predictions of the training
# units (for the random-effect model and the bootstrap residuals), and a compact serialisable
# form of the fitted model for the parallel bootstrap.
#
# Learner functions are shipped to PSOCK workers with detach_env() (fit, oof, pack and, for
# cross-fitting, predict), so they must be
# self-contained: no free variables (parameters go in `args`), package functions called with
# `pkg::fun`.

#' Define a learner for the fixed part of the mixed-effects model
#'
#' A learner replaces the random forest of the MERF by any supervised regression method.
#' The functions must be self-contained (they are sent to parallel workers without their
#' enclosing environment): refer to package functions as \code{pkg::fun} and keep tuning
#' parameters in \code{args}, which is passed to \code{fit} and \code{oof}.
#'
#' @param name Character, the name of the learner.
#' @param fit \code{function(x, y, args, threads)}: fits the learner to the covariates \code{x}
#' (data.frame, one row per sample unit) and target \code{y}; returns the fitted model.
#' \code{args} is the list of learner parameters (\code{args} below, completed by the
#' \code{...} of \code{\link{SAEforest_model}} and by \code{importance}); \code{threads} is the
#' number of threads the learner may use.
#' @param predict \code{function(model, x, threads)}: predictions for the rows of \code{x}
#' (a numeric matrix with the covariates in the order of the sample covariates, or a
#' data.frame); returns a numeric vector.
#' @param oof Optional \code{function(model, x, y, args, threads)}: out-of-sample predictions of
#' the training units available from the fit itself, e.g. out-of-bag predictions of a random
#' forest. If \code{NULL}, out-of-fold predictions are obtained by K-fold cross-fitting
#' (\code{\link{SAEforest_model}}, argument \code{oof_folds}).
#' @param pack Optional \code{function(model)}: compact, serialisable form of the model used
#' to store bootstrap models on disk (default: the model itself).
#' @param unpack Optional \code{function(packed)}: restores an object usable by
#' \code{predict} (default: identity). Runs in the main process only.
#' @param info Optional \code{function(model)}: named list describing the fit (shown by
#' \code{summary}), e.g. the number of trees and an in-fit R-squared.
#' @param packages Character vector of packages loaded in the parallel workers.
#' @param args Named list of learner parameters.
#' @param bounded Logical: are the predictions bounded by the range of the training target
#' (tree ensembles)? For unbounded learners (linear models, splines) the tables of the
#' Box-Cox back-transformation are built on a wider range of the linear predictor.
#'
#' @return An object of class \code{sae_learner}.
#' @seealso \code{\link{learner_ranger}}
#' @export
sae_learner <- function(name, fit, predict, oof = NULL, pack = NULL, unpack = NULL, info = NULL,
                        packages = character(), args = list(), bounded = FALSE) {
  if (!is.character(name) || length(name) != 1) stop("name must be a single character string.")
  if (!is.function(fit) || !is.function(predict)) stop("fit and predict must be functions.")
  for (f in list(oof, pack, unpack, info)) if (!is.null(f) && !is.function(f)) stop("oof, pack, unpack and info must be functions or NULL.")
  if (!is.list(args) || (length(args) > 0 && is.null(names(args)))) stop("args must be a named list.")
  structure(list(name = name, fit = fit, predict = predict, oof = oof,
                 pack = pack, unpack = unpack, info = info,
                 packages = packages, args = args, bounded = isTRUE(bounded)),
            class = "sae_learner")
}

#' @export
print.sae_learner <- function(x, ...) {
  cat("SAE learner '", x$name, "'", if (is.null(x$oof)) " (out-of-fold predictions by cross-fitting)" else "", "\n", sep = "")
  if (length(x$args)) cat("  parameters:", paste(names(x$args), vapply(x$args, function(a) paste(format(a), collapse = ","), ""), sep = " = ", collapse = ", "), "\n")
  invisible(x)
}

# learner functions detached from the package namespace, for the workers
learner_detach <- function(l) {
  for (n in c("fit", "predict", "oof", "pack", "unpack", "info")) {
    if (is.function(l[[n]])) l[[n]] <- detach_env(l[[n]])
  }
  l
}

# fills the default methods and merges the learner parameters with those of the call
learner_prepare <- function(l, params = list(), importance = NULL) {
  if (!inherits(l, "sae_learner")) stop("learner must be an object created by sae_learner().")
  if (is.null(l$pack)) l$pack <- function(model) model
  if (is.null(l$unpack)) l$unpack <- function(packed) packed
  args <- l$args
  for (n in names(params)) args[[n]] <- params[[n]]
  if (!is.null(importance)) args$importance <- importance
  l$args <- args
  # the functions keep their environment here (package namespace for the predefined learners);
  # only the copy shipped to the bootstrap workers is detached (learner_detach)
  l
}


#' Random forest learner (ranger)
#'
#' The default learner of \code{\link{SAEforest_model}}: a \pkg{ranger} regression forest with
#' out-of-bag predictions as out-of-sample predictions, flattened for fast C++ prediction
#' (\code{\link{flatten_forest}}). With this learner the results are identical to those of the
#' package without the learner interface.
#'
#' @param ... Parameters passed to \code{\link[ranger]{ranger}} (e.g. \code{num.trees},
#' \code{mtry}, \code{min.node.size}); those given to \code{SAEforest_model} take precedence.
#' @return An object of class \code{sae_learner}.
#' @export
learner_ranger <- function(...) {
  sae_learner(
    name = "ranger",
    packages = "ranger",
    bounded = TRUE,
    args = list(...),
    fit = function(x, y, args, threads) {
      imp <- if (is.null(args$importance)) "none" else args$importance
      args$importance <- NULL
      rf <- do.call(ranger::ranger, c(list(x = x, y = y, importance = imp, num.threads = threads), args))
      # units that are in-bag in every tree (possible with few trees) have no OOB prediction:
      # use the prediction of the whole forest for them (otherwise lmer drops them)
      nan_oob <- is.na(rf$predictions)
      if (any(nan_oob)) {
        rf$predictions[nan_oob] <- stats::predict(rf, x[nan_oob, , drop = FALSE], num.threads = threads)$predictions
      }
      rf
    },
    oof = function(model, x, y, args, threads) model$predictions,
    # flat C++ prediction in the main process; plain ranger prediction where the package
    # functions are not available (parallel workers)
    predict = function(model, x, threads) {
      if (is.null(threads)) threads <- parallel::detectCores()
      if (exists("predict_flat_forest", mode = "function")) {
        return(predict_flat_forest(model, x, num_threads = threads))
      }
      stats::predict(model, data = x, num.threads = threads)$predictions
    },
    pack = flatten_forest,
    unpack = function(packed) {
      if (is.null(packed$ranger) && is.null(packed$depth)) {
        packed$depth <- cpp_tree_depth(packed$left, packed$right, packed$tree_start)
      }
      packed
    },
    info = function(model) {
      list(type = model$treetype, num.trees = model$num.trees, num.independent.variables = model$num.independent.variables,
           mtry = model$mtry, min.node.size = model$min.node.size, importance.mode = model$importance.mode,
           splitrule = model$splitrule, r.squared = model$r.squared)
    }
  )
}
