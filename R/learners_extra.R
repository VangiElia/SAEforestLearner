# Predefined learners other than the random forest -----------------------------------------
#
# All functions are self-contained (see sae_learner): parameters travel in `args`, package
# functions are called with pkg::fun. None of these learners has native out-of-sample
# predictions, so the engine uses K-fold cross-fitting (oof_folds) for the random-effects
# offset and the bootstrap residuals. The packages are suggested, not imported.

# numeric matrix with the columns in the order used for fitting
lrn_matrix <- function(x) {
  m <- if (is.matrix(x)) x else data.matrix(x)
  if (!is.double(m)) storage.mode(m) <- "double"
  m
}

#' Gradient boosting learner (xgboost)
#'
#' Boosted regression trees fitted with \pkg{xgboost}. Out-of-fold predictions of the sample
#' units are obtained by cross-fitting (\code{oof_folds} of \code{\link{SAEforest_model}}).
#' Subsampling uses the R random number generator, so results are reproducible with a seed.
#'
#' @param nrounds Number of boosting rounds.
#' @param eta,max_depth,subsample,colsample_bytree,min_child_weight Booster parameters
#' (see \code{\link[xgboost]{xgb.train}}).
#' @param ... Further booster parameters.
#' @return An object of class \code{sae_learner}.
#' @export
learner_xgboost <- function(nrounds = 300, eta = 0.05, max_depth = 6, subsample = 0.8,
                            colsample_bytree = 0.8, min_child_weight = 5, ...) {
  sae_learner(
    name = "xgboost",
    packages = "xgboost",
    bounded = TRUE,
    args = list(nrounds = nrounds, eta = eta, max_depth = max_depth, subsample = subsample,
                colsample_bytree = colsample_bytree, min_child_weight = min_child_weight, ...),
    fit = function(x, y, args, threads) {
      args$importance <- NULL
      nrounds <- args$nrounds
      args$nrounds <- NULL
      m <- if (is.matrix(x)) x else data.matrix(x)
      if (!is.double(m)) storage.mode(m) <- "double"
      params <- c(list(objective = "reg:squarederror", nthread = if (is.null(threads)) 1L else threads), args)
      dm <- xgboost::xgb.DMatrix(m, label = y)
      xgboost::xgb.train(params = params, data = dm, nrounds = nrounds, verbose = 0)
    },
    predict = function(model, x, threads) {
      m <- if (is.matrix(x)) x else data.matrix(x)
      if (!is.double(m)) storage.mode(m) <- "double"
      stats::predict(model, m)
    },
    pack = function(model) xgboost::xgb.save.raw(model),
    unpack = function(packed) xgboost::xgb.load.raw(packed, as_booster = TRUE),
    info = function(model) list(type = "xgboost regression booster", nrounds = model$niter,
                                num.independent.variables = model$nfeatures)
  )
}


#' Gradient boosting learner (lightgbm)
#'
#' Boosted regression trees fitted with \pkg{lightgbm}; out-of-fold predictions by
#' cross-fitting. The booster seed is drawn from the R random number generator.
#'
#' @param nrounds Number of boosting rounds.
#' @param learning_rate,num_leaves,min_data_in_leaf,feature_fraction,bagging_fraction
#' Booster parameters (see \code{\link[lightgbm]{lgb.train}}).
#' @param ... Further booster parameters.
#' @return An object of class \code{sae_learner}.
#' @export
learner_lightgbm <- function(nrounds = 300, learning_rate = 0.05, num_leaves = 15, min_data_in_leaf = 10,
                             feature_fraction = 0.8, bagging_fraction = 0.8, ...) {
  sae_learner(
    name = "lightgbm",
    packages = "lightgbm",
    bounded = TRUE,
    args = list(nrounds = nrounds, learning_rate = learning_rate, num_leaves = num_leaves,
                min_data_in_leaf = min_data_in_leaf, feature_fraction = feature_fraction,
                bagging_fraction = bagging_fraction, ...),
    fit = function(x, y, args, threads) {
      args$importance <- NULL
      nrounds <- args$nrounds
      args$nrounds <- NULL
      m <- if (is.matrix(x)) x else data.matrix(x)
      if (!is.double(m)) storage.mode(m) <- "double"
      params <- c(list(objective = "regression", num_threads = if (is.null(threads)) 1L else threads,
                       bagging_freq = 1L, verbose = -1L, seed = sample.int(.Machine$integer.max, 1)), args)
      ds <- lightgbm::lgb.Dataset(m, label = y, params = list(verbose = -1L))
      lightgbm::lgb.train(params = params, data = ds, nrounds = nrounds, verbose = -1L)
    },
    predict = function(model, x, threads) {
      m <- if (is.matrix(x)) x else data.matrix(x)
      if (!is.double(m)) storage.mode(m) <- "double"
      stats::predict(model, m, params = list(num_threads = if (is.null(threads)) 1L else threads))
    },
    pack = function(model) model$save_model_to_string(),
    unpack = function(packed) lightgbm::lgb.load(model_str = packed),
    info = function(model) list(type = "lightgbm regression booster", nrounds = model$current_iter())
  )
}


#' Penalised linear learner (glmnet)
#'
#' Elastic-net linear regression (\pkg{glmnet}) with the penalty chosen by internal
#' cross-validation (\code{lambda.min}); out-of-fold predictions by cross-fitting.
#'
#' @param alpha Elastic-net mixing parameter (1 = lasso, 0 = ridge).
#' @param nfolds Folds of the internal cross-validation of the penalty.
#' @param ... Further arguments of \code{\link[glmnet]{cv.glmnet}}.
#' @return An object of class \code{sae_learner}.
#' @export
learner_glmnet <- function(alpha = 1, nfolds = 5, ...) {
  sae_learner(
    name = "glmnet",
    packages = "glmnet",
    args = list(alpha = alpha, nfolds = nfolds, ...),
    fit = function(x, y, args, threads) {
      args$importance <- NULL
      m <- if (is.matrix(x)) x else data.matrix(x)
      if (!is.double(m)) storage.mode(m) <- "double"
      do.call(glmnet::cv.glmnet, c(list(x = m, y = y), args))
    },
    predict = function(model, x, threads) {
      m <- if (is.matrix(x)) x else data.matrix(x)
      if (!is.double(m)) storage.mode(m) <- "double"
      as.numeric(stats::predict(model, newx = m, s = "lambda.min"))
    },
    info = function(model) list(type = "elastic-net linear model (glmnet)", lambda.min = model$lambda.min,
                                nonzero = model$nzero[model$index["min", ]])
  )
}


#' Generalised additive learner (mgcv)
#'
#' Additive model with a thin-plate spline for every covariate (\code{\link[mgcv]{gam}},
#' REML); covariates with few distinct values enter linearly. Out-of-fold predictions by
#' cross-fitting.
#'
#' @param k Basis dimension of the smooth terms.
#' @param min_unique Minimum number of distinct values for a smooth term.
#' @param ... Further arguments of \code{\link[mgcv]{gam}}.
#' @return An object of class \code{sae_learner}.
#' @export
learner_gam <- function(k = 5, min_unique = 10, ...) {
  sae_learner(
    name = "gam",
    packages = "mgcv",
    args = list(k = k, min_unique = min_unique, ...),
    fit = function(x, y, args, threads) {
      args$importance <- NULL
      k <- args$k; min_unique <- args$min_unique
      args$k <- NULL; args$min_unique <- NULL
      d <- as.data.frame(if (is.matrix(x)) x else data.matrix(x))
      vn <- paste0("v", seq_along(d))
      names(d) <- vn
      d$.y <- y
      nu <- vapply(d[vn], function(v) length(unique(v)), integer(1))
      terms <- ifelse(nu >= max(min_unique, k + 1), paste0("s(", vn, ", k = ", k, ")"), vn)
      f <- stats::as.formula(paste(".y ~", paste(terms, collapse = " + ")))
      m <- do.call(mgcv::gam, c(list(formula = f, data = d, method = "REML"), args))
      m$.var_names <- colnames(x)
      m
    },
    predict = function(model, x, threads) {
      d <- as.data.frame(if (is.matrix(x)) x else data.matrix(x))
      names(d) <- paste0("v", seq_along(d))
      as.numeric(mgcv::predict.gam(model, newdata = d))
    },
    pack = function(model) {
      # drop what prediction does not need
      model$model <- NULL; model$y <- NULL; model$residuals <- NULL; model$fitted.values <- NULL
      model$linear.predictors <- NULL; model$weights <- NULL; model$prior.weights <- NULL
      model
    },
    info = function(model) list(type = "generalised additive model (mgcv)", edf = sum(model$edf),
                                r.squared = summary(model)$r.sq)
  )
}


#' Linear learner (lm) with exact leave-one-out predictions
#'
#' Ordinary least squares on all covariates. The out-of-sample predictions of the sample
#' units are the exact leave-one-out predictions, \code{y - e / (1 - h)} with \code{h} the
#' leverages, so no cross-fitting is needed. With this learner the MERF iterations converge to
#' the nested-error linear mixed model (Battese, Harter & Fuller, 1988), i.e. the EBLUP of the
#' domain means; with \code{transform = "boxcox"} the back-transformed means are the empirical
#' best predictions of Molina & Rao (2010) for the mean.
#'
#' @return An object of class \code{sae_learner}.
#' @export
learner_lm <- function() {
  sae_learner(
    name = "lm",
    fit = function(x, y, args, threads) {
      m <- if (is.matrix(x)) x else data.matrix(x)
      if (!is.double(m)) storage.mode(m) <- "double"
      fit <- stats::lm.fit(cbind(1, m), y)
      q <- fit$qr
      h <- rowSums(qr.Q(q)[, seq_len(q$rank), drop = FALSE]^2)
      coef <- fit$coefficients
      coef[is.na(coef)] <- 0
      list(coef = coef, loo = y - fit$residuals / (1 - h))
    },
    oof = function(model, x, y, args, threads) model$loo,
    predict = function(model, x, threads) {
      m <- if (is.matrix(x)) x else data.matrix(x)
      if (!is.double(m)) storage.mode(m) <- "double"
      as.numeric(model$coef[1] + m %*% model$coef[-1])
    },
    pack = function(model) list(coef = model$coef),
    info = function(model) list(type = "linear model (OLS, exact leave-one-out)", coefficients = length(model$coef))
  )
}
