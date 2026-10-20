#' Main function for unit-level MERF
#'
#' This function enables the use of Mixed Effects Random Forests (MERFs) by effectively
#' combining a random forest from \pkg{ranger} with a model capturing random effects from
#' \pkg{lme4}. The MERF algorithm is an algorithmic procedure reminiscent of an EM-algorithm
#' (see Details). The function is the base-function for the wrapping function (\code{\link{SAEforest_model}}
#' and should not be directly used by the ordinary user. Recommended exceptions are applications exceeding
#' the scope of existing wrapper functions or further research. The function \code{MERFranger}
#' allows to model complex patterns of structural relations (see Examples). The function returns
#' an object of class \code{MERFranger}, which can be used to produce unit-level predictions. In contrast to
#' the wrapping functions, this function does not directly provide SAE estimates on domain-specific indicators.
#'
#' The random forest can be replaced by any supervised learner through the argument
#' \code{learner} (see \code{\link{sae_learner}}): the algorithm then alternates that learner and
#' the random-effects model in the same way.
#'
#' @param Y Continuous input value of target variable.
#' @param X Matrix of predictive covariates.
#' @param random Specification of random effects terms following the syntax of \link[lme4]{lmer}.
#' Random effect terms are specified by vertical bars \code{(|)} separating expressions for design matrices
#' from grouping factors. For further details see \link[lme4]{lmer} and the example below.
#' @param data data.frame of sample data including the specified elements of \code{Y} and
#' \code{X}.
#' @param initialRandomEffects Numeric value or vector of initial estimate of random effects.
#' Defaults to 0.
#' @param ErrorTolerance Numeric value to monitor the MERF algorithm's convergence. Defaults to 1e-04.
#' @param MaxIterations Numeric value specifying the maximal amount of iterations for the
#' MERF algorithm. Defaults to 25.
#' @param importance Variable importance mode processed by the
#' random forest from the \pkg{ranger}. Must be 'none', 'impurity', 'impurity_corrected',
#' 'permutation'. For further details see \link[ranger]{ranger}.
#' @param na.rm Logical. Whether missing values should be removed. Defaults to \code{TRUE}.
#' @param learner Learner of the fixed part (see \code{\link{sae_learner}}); \code{NULL} (default) uses
#' \code{\link{learner_ranger}}, i.e. a \pkg{ranger} forest as in the original package.
#' @param oof_folds Number of folds of the cross-fitting used for the out-of-sample predictions
#' of learners without native ones (see \code{\link{SAEforest_model}}). Defaults to 5.
#' @param em_criterion Convergence criterion of the MERF iterations: \code{"original"} (relative
#' change of the log-likelihood below \code{ErrorTolerance}, as in the original package),
#' \code{"robust"} (for learners whose log-likelihood oscillates between refits: the means of two
#' consecutive windows of iterations are compared, the iterations stop when no improvement has
#' been seen for \code{patience} iterations, and the best iteration is kept) or \code{"auto"}
#' (default: \code{"original"} with \code{learner_ranger}, \code{"robust"} otherwise).
#' @param em_control List of settings of the robust criterion: \code{min_iter} (3),
#' \code{window} (3), \code{patience} (5), \code{keep_best} (\code{TRUE}).
#' @param ... Additional parameters are directly passed to the random forest \link[ranger]{ranger}
#' (or to the \code{fit} function of \code{learner}).
#' Most important parameters are for instance \code{mtry} (number of variables to possibly split at
#' in each node), or \code{num.trees} (number of trees). For further details on possible parameters
#' see \link[ranger]{ranger} and the example below.
#'
#' @return An object of class MERFranger includes the following elements:
#'
#' \item{\code{Forest}}{A random forest of class \link[ranger]{ranger} modelling fixed effects
#' of the model (or the fitted model of \code{learner}).}
#' \item{\code{EffectModel}}{A model of random effects of class \code{\link[lme4]{merMod}} capturing
#' structural components of MERFs and modeling random components.}
#' \item{\code{RandomEffects}}{List element containing the values of random intercepts from \code{EffectModel}.}
#' \item{\code{RanEffSD}}{Numeric value of the standard deviation of random intercepts.}
#' \item{\code{ErrorSD}}{Numeric value of standard deviation of unit-level errors.}
#' \item{\code{VarianceCovariance}}{VarCorr matrix from \code{EffectModel}.}
#' \item{\code{LogLik}}{Vector with numerical entries showing the loglikelihood of the MERF algorithm.}
#' \item{\code{IterationsUsed}}{Numeric number of iterations used until convergence of the MERF algorithm.}
#' \item{\code{BestIteration}}{The iteration whose fit is returned (equal to \code{IterationsUsed}
#' with the original criterion).}
#' \item{\code{OOBresiduals}}{Vector of OOB-residuals.}
#' \item{\code{OOBpred}}{Vector of out-of-sample (OOB or out-of-fold) predictions of the learner
#' for the sample units, at the last iteration.}
#' \item{\code{Random}}{Character specifying the random intercept in the random effects model.}
#' \item{\code{ErrorTolerance}}{Numerical value to monitor the MERF algorithm's convergence.}
#' \item{\code{initialRandomEffects}}{Numeric value or vector of initial specification of random effects.}
#' \item{\code{MaxIterations}}{Numeric value specifying the maximal amount of iterations for the
#' MERF algorithm.}
#' \item{\code{learner}}{The learner used (class \code{sae_learner}).}
#'
#' @details
#' There exists a generic function for \code{predict} for objects obtained by \code{MERFranger}.
#'
#' The MERF algorithm iteratively optimizes two separate steps: a) the random forest
#' function, assuming the random effects term to be correct and b) estimates the random
#' effects part, assuming the OOB-predictions from the forest to be correct. Overall convergence
#' of the algorithm is monitored by the log-likelihood of a joint model of both components. For
#' further details see Krennmair & Schmid (2022) or Hajjem et al. (2014).
#'
#' Note that the \code{MERFranger} object is a composition of elements from a random forest of class
#' \code{ranger} and a random effects model of class \code{\link[lme4]{merMod}}. Thus, all generic functions are
#' applicable to corresponding objects. For further details on generic functions see \code{\link[ranger]{ranger}}
#' and \code{\link[lme4]{lmer}} as well as the examples below.
#'
#' @references
#' Hajjem, A., Bellavance, F., & Larocque, D. (2014). Mixed-Effects Random Forest for Clustered
#' Data. Journal of Statistical Computation and Simulation, 84 (6), 1313–1328.
#'
#' Krennmair, P., & Schmid, T. (2022). Flexible Domain Prediction Using Mixed Effects
#' Random Forests. Journal of Royal Statistical Society: Series C (Applied Statistics) (forthcoming).
#'
#' @seealso \code{\link{SAEforest}}, \code{\link[ranger]{ranger}}, \code{\link[lme4]{lmer}},
#' \code{\link{SAEforest_model}}, \code{\link{sae_learner}}
#'
#' @export
#' @import ranger
#' @import lme4
#'
#' @examples
#' \donttest{
#' # Load Data
#' data("eusilcA_pop")
#' data("eusilcA_smp")
#'
#' income <- eusilcA_smp$eqIncome
#' X_covar <- eusilcA_smp[, -c(1, 16, 17, 18)]
#'
#' # Example 1:
#' # Calculating general model used in wrapper functions
#'
#' model1 <- MERFranger(Y = income, X = X_covar, random = "(1|district)",
#'                      data = eusilcA_smp, num.trees=50)
#'
#' # get individual predictions:
#'
#' ind_pred <- predict(model1, eusilcA_pop)
#'}
MERFranger <- function(Y,
                       X,
                       random,
                       data,
                       importance = "none",
                       initialRandomEffects = 0,
                       ErrorTolerance = 0.0001,
                       MaxIterations = 25,
                       na.rm = TRUE,
                       learner = NULL,
                       oof_folds = 5,
                       em_criterion = c("auto", "original", "robust"),
                       em_control = list(),
                       ...) {
  dots <- list(...)
  threads <- dots$num.threads
  dots$num.threads <- NULL
  if (is.null(learner)) learner <- learner_ranger()
  learner <- learner_prepare(learner, params = dots, importance = importance)
  ctrl <- em_control_prepare(match.arg(em_criterion), em_control, learner)
  MERF_learner(Y = Y, X = X, random = random, data = data, learner = learner, threads = threads,
               initialRandomEffects = initialRandomEffects, ErrorTolerance = ErrorTolerance,
               MaxIterations = MaxIterations, na.rm = na.rm, oof_folds = oof_folds, em_control = ctrl)
}

# Settings of the EM convergence criterion. "auto": "original" with learner_ranger (results
# identical to the original package), "robust" with any other learner.
em_control_prepare <- function(criterion, control, learner) {
  if (criterion == "auto") criterion <- if (identical(learner$name, "ranger")) "original" else "robust"
  ctrl <- list(criterion = criterion, min_iter = 3L, window = 3L, patience = 5L, keep_best = TRUE)
  for (n in names(control)) ctrl[[n]] <- control[[n]]
  if (!ctrl$criterion %in% c("original", "robust")) stop("em_criterion must be 'auto', 'original' or 'robust'.")
  ctrl
}


# MERF algorithm with a generic learner (internal; also run in the bootstrap workers) --------
# The learner is a prepared sae_learner (learner_prepare): its functions are self-contained.
# Out-of-sample predictions of the training units: the learner's own (oof, e.g. OOB of a forest),
# otherwise K-fold cross-fitting (fit on K-1 folds, predict the held-out fold; the model fitted on
# all units is kept for prediction). The random effects model uses them as offset and the
# residuals returned (OOBresiduals) are based on them.
MERF_learner <- function(Y, X, random, data, learner, threads = NULL, initialRandomEffects = 0,
                         ErrorTolerance = 0.0001, MaxIterations = 25, na.rm = TRUE, oof_folds = 5,
                         em_control = list(criterion = "original")) {
  if (na.rm == TRUE) {
    comp_smp <- stats::complete.cases(data)
    data <- data[comp_smp, ]
    Y <- Y[comp_smp]
    X <- X[comp_smp, ]
  }
  Target <- Y
  ContinueCondition <- TRUE
  iterations <- 0
  AdjustedTarget <- Target - initialRandomEffects
  oldLogLik <- 0
  args <- learner$args
  native_oof <- is.function(learner$oof)
  folds <- NULL
  if (!native_oof) {
    # fixed folds for all iterations (one RNG draw, taken before any fit)
    folds <- sample(rep(seq_len(oof_folds), length.out = length(Target)))
  }
  # convergence: "original" (relative change of the log-likelihood below ErrorTolerance, as in
  # SAEforest) or "robust" for stochastic learners whose log-likelihood oscillates from one
  # refit to the next: compare the means of two consecutive windows of iterations, stop when no
  # improvement has been seen for `patience` iterations, and keep the best iteration
  ctrl <- em_control
  robust <- identical(ctrl$criterion, "robust")
  ll_hist <- numeric(0)
  best <- list(ll = -Inf, iter = 0)
  while (ContinueCondition) {
    iterations <- iterations + 1
    model <- learner$fit(X, AdjustedTarget, args, threads)
    if (native_oof) {
      forest_preds <- learner$oof(model, X, AdjustedTarget, args, threads)
    } else {
      forest_preds <- numeric(length(Target))
      for (k in seq_len(oof_folds)) {
        te <- folds == k
        m_k <- learner$fit(X[!te, , drop = FALSE], AdjustedTarget[!te], args, threads)
        forest_preds[te] <- learner$predict(m_k, X[te, , drop = FALSE], threads)
      }
    }
    if (anyNA(forest_preds)) stop("The learner returned missing out-of-sample predictions.")
    f0 <- stats::as.formula(paste0("Target ~ -1+", random))
    lmefit <- lme4::lmer(f0, data = data, REML = FALSE, offset = forest_preds)

    newLogLik <- as.numeric(stats::logLik(lmefit))

    if (robust) {
      ll_hist <- c(ll_hist, newLogLik)
      if (newLogLik > best$ll) {
        best <- list(ll = newLogLik, iter = iterations, model = model, lmefit = lmefit, forest_preds = forest_preds)
      }
      w <- ctrl$window
      converged <- FALSE
      if (iterations >= max(ctrl$min_iter, 2 * w)) {
        m_new <- mean(ll_hist[(iterations - w + 1):iterations])
        m_old <- mean(ll_hist[(iterations - 2 * w + 1):(iterations - w)])
        converged <- abs(m_new - m_old) / abs(m_old) < ErrorTolerance
      }
      stalled <- iterations - best$iter >= ctrl$patience
      ContinueCondition <- !converged && !stalled && iterations < MaxIterations
    } else {
      ContinueCondition <- (abs(abs(newLogLik - oldLogLik[iterations]) / oldLogLik[iterations]) > ErrorTolerance &
        iterations < MaxIterations)
    }
    oldLogLik <- c(oldLogLik, newLogLik)
    AllEffects <- stats::predict(lmefit)
    AdjustedTarget <- Target - (AllEffects - forest_preds)
  }
  best_iter <- iterations
  if (robust && isTRUE(ctrl$keep_best) && best$iter > 0) {
    model <- best$model
    lmefit <- best$lmefit
    forest_preds <- best$forest_preds
    best_iter <- best$iter
  }

  data$forest_preds <- NULL
  residuals <- Target - stats::predict(lmefit)

  result <- list(
    Forest = model,
    EffectModel = lmefit,
    RandomEffects = lme4::ranef(lmefit),
    RanEffSD = as.data.frame(lme4::VarCorr(lmefit))$sdcor[1],
    ErrorSD = stats::sigma(lmefit),
    VarianceCovariance = lme4::VarCorr(lmefit),
    LogLik = oldLogLik,
    IterationsUsed = iterations,
    BestIteration = best_iter,
    OOBresiduals = residuals,
    OOBpred = forest_preds,
    Random = random,
    ErrorTolerance = ErrorTolerance,
    initialRandomEffects = initialRandomEffects,
    MaxIterations = MaxIterations,
    em_control = ctrl,
    learner = learner,
    threads = threads
  )

  class(result) <- "MERFranger"

  return(result)
}

# Generic predict function for MERFranger ------------------------------------------------
#' @export
predict.MERFranger <- function(object, newdata, ...) {
  fixed <- object$learner$predict(object$Forest, newdata, object$threads)
  fixed + predict(object$EffectModel, newdata = newdata, allow.new.levels = TRUE)
}
