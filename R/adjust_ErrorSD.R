# Auxiliary function for MSEs to perform bias correction ---------------------------------

adjust_ErrorSD <- function(Y, X, smp_data, mod, B = 100, ...) {
  pred_OOB <- matrix(mod$Forest$predictions, ncol = B,
                     nrow = length(mod$Forest$predictions), byrow = FALSE)

  e_ij <- Y - predict(mod$Forest, smp_data)$predictions
  e_ij <- e_ij - mean(e_ij)

  y_star_OOB <- pred_OOB + sample(e_ij, size = length(pred_OOB), replace = TRUE)

  my_estim_f2 <- function(x) {
    ranger::ranger(y = x, x = X, data = smp_data, ...)
  }
  my_f_n2 <- pbapply::pbapply(y_star_OOB, 2, my_estim_f2)

  my_pred_f <- function(x) {
    x$predictions
  }
  pred_OOB_star <- sapply(my_f_n2, my_pred_f)

  Adjustment <- (pred_OOB_star - pred_OOB)^2

  outvar <- sqrt(mod$ErrorSD^2 - mean(rowMeans(Adjustment)))

  return(outvar)
}


# Bias correction of the error SD (Mendez & Lohr, 2011) with a generic learner --------------
# Scalable engine. Same computation as adjust_ErrorSD(): the out-of-sample predictions of the
# sample units are perturbed with re-sampled in-sample residuals, the learner is refitted B
# times and the mean squared difference of its out-of-sample predictions from the original ones
# is subtracted from the residual variance. Learners without native out-of-sample predictions
# use the same K-fold cross-fitting as MERF_learner().
adjust_ErrorSD_learner <- function(Y, X, mod, B, learner, threads, pred_in, oof_folds = 5) {
  oob <- mod$OOBpred
  pred_OOB <- matrix(oob, ncol = B, nrow = length(oob), byrow = FALSE)

  e_ij <- Y - pred_in
  e_ij <- e_ij - mean(e_ij)

  y_star_OOB <- pred_OOB + sample(e_ij, size = length(pred_OOB), replace = TRUE)

  args <- learner$args
  args$importance <- "none"
  native_oof <- is.function(learner$oof)
  folds <- if (native_oof) NULL else sample(rep(seq_len(oof_folds), length.out = length(Y)))
  refit_oof <- function(ys) {
    if (native_oof) {
      m <- learner$fit(X, ys, args, threads)
      return(learner$oof(m, X, ys, args, threads))
    }
    p <- numeric(length(ys))
    for (k in seq_len(oof_folds)) {
      te <- folds == k
      m_k <- learner$fit(X[!te, , drop = FALSE], ys[!te], args, threads)
      p[te] <- learner$predict(m_k, X[te, , drop = FALSE], threads)
    }
    p
  }
  pred_OOB_star <- pbapply::pbapply(y_star_OOB, 2, refit_oof)

  Adjustment <- (pred_OOB_star - pred_OOB)^2

  sqrt(mod$ErrorSD^2 - mean(rowMeans(Adjustment)))
}
