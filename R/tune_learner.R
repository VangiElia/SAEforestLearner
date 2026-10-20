# Tuning of a learner's parameters by spatial cross-validation of the MERF ------------------
#
# The learner is tuned inside the mixed model (MERF iterations with random intercepts), on
# held-out plots of spatial blocks, with predictions on the original scale when a Box-Cox
# transformation is used (same back-transformation as the estimator). Functions are
# self-contained for the PSOCK workers (see learners.R).

#' Spatial cross-validation folds from square blocks
#'
#' Assigns the units to \code{k} folds by square spatial blocks of side \code{block_size}
#' (in the units of the coordinates), balancing the number of units per fold: blocks are
#' taken in random order, largest first, each to the currently smallest fold.
#'
#' @param x,y Coordinates of the units (projected, e.g. metres).
#' @param block_size Side of the square blocks, in the units of the coordinates.
#' @param k Number of folds.
#' @param seed Optional seed.
#' @return An integer vector of fold labels (1 to \code{k}), with the block identifier as
#' attribute \code{"block"}.
#' @export
spatial_folds <- function(x, y, block_size, k = 5, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  block <- paste(floor(x / block_size), floor(y / block_size), sep = "_")
  tb <- table(block)
  tb <- tb[sample(length(tb))]
  tb <- tb[order(-as.numeric(tb))]
  load <- numeric(k)
  fold_of <- integer(length(tb))
  for (j in seq_along(tb)) {
    f <- which(load == min(load))
    if (length(f) > 1) f <- sample(f, 1)
    fold_of[j] <- f
    load[f] <- load[f] + tb[j]
  }
  out <- as.integer(fold_of[match(block, names(tb))])
  attr(out, "block") <- block
  out
}


#' Tune a learner by spatial cross-validation of the mixed-effects model
#'
#' Fits the MERF (\code{\link{MERFranger}}) with the learner returned by
#' \code{learner_fun(...)} for every row of \code{grid} and every fold, predicts the held-out
#' units (fixed part plus the random effect of their domain when it is in the training part)
#' and summarises the accuracy on the original scale. With \code{transform = "boxcox"} the
#' model is fitted to the transformed target and the predictions are the empirical-best-type
#' back-transformations used by \code{\link{SAEforest_model}}.
#'
#' @param Y,X,dName,smp_data As in \code{\link{SAEforest_model}}.
#' @param learner_fun A learner constructor, e.g. \code{\link{learner_xgboost}}; its arguments
#' are the columns of \code{grid}.
#' @param grid A data.frame, one row per parameter combination (columns = arguments of
#' \code{learner_fun}; \code{NULL} for the learner defaults only); list columns are allowed.
#' @param folds Integer vector of fold labels, e.g. from \code{\link{spatial_folds}}.
#' @param transform,bc_lambda,bc_shift Transformation of the target, as in
#' \code{\link{SAEforest_model}}.
#' @param cores Number of parallel processes (one fit per process).
#' @param threads Threads per fit (passed to the learner).
#' @param seed Seed (one per fit is derived from it).
#' @param oof_folds,em_criterion,em_control,MaxIterations As in \code{\link{MERFranger}}.
#' @param ... Further fixed arguments of \code{learner_fun}.
#'
#' @return A list with \code{results} (the grid with mean RMSE, bias, R-squared, slope of
#' the observed on the predicted values and the fold-to-fold SE of the RMSE), \code{best}
#' (the row with the lowest RMSE), \code{by_fold} (all fold results) and the folds used.
#' @export
tune_learner <- function(Y, X, dName, smp_data, learner_fun, grid, folds,
                         transform = c("none", "boxcox"), bc_lambda = NULL, bc_shift = 1,
                         cores = 1, threads = 1, seed = NULL, oof_folds = 5,
                         em_criterion = "auto", em_control = list(), MaxIterations = 25, ...) {
  transform <- match.arg(transform)
  if (transform == "boxcox" && is.null(bc_lambda)) stop("bc_lambda is required with transform = 'boxcox'.")
  if (length(folds) != length(Y)) stop("folds must have one label per sample unit.")
  grid <- if (is.null(grid)) data.frame(row.names = 1) else as.data.frame(grid)   # NULL: the learner defaults
  fixed <- list(...)
  if (!is.null(seed)) set.seed(seed)
  tasks <- expand.grid(g = seq_len(nrow(grid)), fold = sort(unique(folds)))
  tasks$seed <- sample.int(.Machine$integer.max, nrow(tasks))

  X <- as.data.frame(X)
  smp_data <- as.data.frame(smp_data)
  dom <- dom_key(smp_data[[dName]])
  yt <- if (transform == "boxcox") bc_fwd(Y, bc_lambda, bc_shift) else Y
  learners <- lapply(seq_len(nrow(grid)), function(g) {
    a <- c(lapply(grid[g, , drop = FALSE], function(v) if (is.list(v)) v[[1]] else v), fixed)
    l <- do.call(learner_fun, a)
    learner_detach(learner_prepare(l))
  })
  ctrl <- em_control_prepare(em_criterion, em_control, learners[[1]])
  job <- list(X = X, yt = yt, dom = dom, dName = dName, folds = folds, learners = learners, tasks = tasks,
              transform = transform, lambda = bc_lambda, shift = bc_shift, threads = threads,
              oof_folds = oof_folds, ctrl = ctrl, MaxIterations = MaxIterations, obs = Y,
              gh = gh_normal(20))

  run <- function(i) {
    t <- job$tasks[i, ]
    te <- job$folds == t$fold
    lrn <- job$learners[[t$g]]
    d <- job$X[!te, , drop = FALSE]
    d[[job$dName]] <- factor(job$dom[!te], levels = unique(job$dom[!te]))
    set.seed(t$seed)
    m <- job$merf(Y = job$yt[!te], X = job$X[!te, , drop = FALSE], random = paste0("(1|", job$dName, ")"), data = d,
                  learner = lrn, threads = job$threads, na.rm = FALSE, oof_folds = job$oof_folds,
                  em_control = job$ctrl, MaxIterations = job$MaxIterations)
    eta <- lrn$predict(m$Forest, job$X[te, , drop = FALSE], job$threads)
    re <- lme4::ranef(m$EffectModel)[[job$dName]]
    u <- re[match(job$dom[te], rownames(re)), 1]
    known <- !is.na(u); u[!known] <- 0
    if (job$transform == "boxcox") {
      ni <- table(job$dom[!te])[job$dom[te]]; ni[is.na(ni)] <- 0
      su2 <- m$RanEffSD^2; se2 <- m$ErrorSD^2
      gam <- ifelse(known, su2 / (su2 + se2 / pmax(as.numeric(ni), 1)), 0)
      tau <- sqrt(se2 + su2 * (1 - gam))
      gh <- job$gh
      lam <- job$lambda; sh <- job$shift
      bci <- function(t) if (lam == 0) exp(t) - sh else pmax(lam * t + 1, 0)^(1 / lam) - sh
      pred <- 0
      for (q in seq_along(gh$nodes)) pred <- pred + gh$weights[q] * bci(eta + u + gh$nodes[q] * tau)
    } else {
      pred <- eta + u
    }
    obs <- job$obs[te]
    data.frame(g = t$g, fold = t$fold, n = sum(te), RMSE = sqrt(mean((pred - obs)^2)), bias = mean(pred - obs),
               R2 = 1 - sum((obs - pred)^2) / sum((obs - mean(obs))^2),
               slope = unname(stats::coef(stats::lm(obs ~ pred))[2]), iterations = m$IterationsUsed,
               sigma_u = m$RanEffSD, sigma_e = m$ErrorSD)
  }
  job$merf <- detach_env(MERF_learner)

  if (cores > 1) {
    cl <- parallel::makePSOCKcluster(min(cores, nrow(tasks)))
    on.exit(parallel::stopCluster(cl), add = TRUE)
    pk <- unique(unlist(lapply(learners, `[[`, "packages")))
    parallel::clusterCall(cl, detach_env(function(job, pk) {
      suppressPackageStartupMessages({ requireNamespace("lme4"); for (p in pk) requireNamespace(p) })
      assign("job", job, envir = globalenv()); invisible(TRUE)
    }), job, pk)
    res <- pbapply::pblapply(seq_len(nrow(tasks)), detach_env(run), cl = cl)
  } else {
    res <- lapply(seq_len(nrow(tasks)), run)
  }
  by_fold <- do.call(rbind, res)
  agg <- do.call(rbind, lapply(split(by_fold, by_fold$g), function(z) {
    data.frame(g = z$g[1], RMSE = mean(z$RMSE), RMSE_se = stats::sd(z$RMSE) / sqrt(nrow(z)), bias = mean(z$bias),
               R2 = mean(z$R2), slope = mean(z$slope), iterations = mean(z$iterations),
               sigma_u = mean(z$sigma_u), sigma_e = mean(z$sigma_e))
  }))
  results <- cbind(grid[agg$g, , drop = FALSE], agg[, -1])
  rownames(results) <- NULL
  list(results = results, best = results[which.min(results$RMSE), , drop = FALSE], by_fold = by_fold, folds = folds)
}
