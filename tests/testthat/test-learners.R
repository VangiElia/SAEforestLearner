# Tests of the learner interface (scalable engine)

fit_lrn <- function(learner = NULL, B = 3, ...) {
  data("eusilcA_pop", package = "SAEforestLearner", envir = environment())
  data("eusilcA_smp", package = "SAEforestLearner", envir = environment())
  smp <- eusilcA_smp
  X <- smp[, -c(1, 2, 16, 17, 18)]   # numeric covariates only
  # the ... of SAEforest_model are learner parameters: num.trees only for the forest
  args <- list(Y = smp$eqIncome, X = X, dName = "district", smp_data = smp, pop_data = eusilcA_pop,
               MSE = if (B > 0) "nonparametric" else "none", B = B, B_adj = B,
               seed = 11, cores = 2, verbose = FALSE, learner = learner, ...)
  if (is.null(learner) || identical(learner$name, "ranger")) args$num.trees <- 30
  do.call(SAEforest_model, args)
}

test_that("learner_ranger() is the default learner (identical results)", {
  m0 <- fit_lrn()
  m1 <- fit_lrn(learner = learner_ranger())
  expect_identical(m0$Indicators$Mean, m1$Indicators$Mean)
  expect_identical(m0$MSE_Estimates$Mean, m1$MSE_Estimates$Mean)
  expect_identical(m0$PopInfo$learner, "ranger")
  expect_s3_class(m0$MERFmodel$learner, "sae_learner")
})

test_that("a user-defined learner (linear model) runs through fit, bootstrap and prediction", {
  lrn_lm <- sae_learner(
    name = "lm",
    fit = function(x, y, args, threads) {
      d <- as.data.frame(x); d$.y <- y
      stats::lm(.y ~ ., data = d)
    },
    predict = function(model, x, threads) as.numeric(stats::predict(model, newdata = as.data.frame(x))),
    pack = function(model) { model$model <- NULL; model$qr$qr <- model$qr$qr[0, , drop = FALSE]; model$coefficients },
    unpack = function(packed) packed,
    info = function(model) list(type = "linear model")
  )
  # pack returns the coefficients only: predict must then handle a numeric vector
  lrn_lm$predict <- function(model, x, threads) {
    m <- if (is.matrix(x)) x else data.matrix(x)
    if (is.numeric(model)) return(as.numeric(model[1] + m[, names(model)[-1], drop = FALSE] %*% model[-1]))
    as.numeric(stats::predict(model, newdata = as.data.frame(x)))
  }
  m <- fit_lrn(learner = lrn_lm, oof_folds = 4)
  expect_identical(m$PopInfo$learner, "lm")
  expect_true(all(is.finite(m$Indicators$Mean)))
  expect_true(all(is.finite(m$MSE_Estimates$Mean) & m$MSE_Estimates$Mean > 0))
  expect_length(m$MERFmodel$OOBpred, nrow(m$MERFmodel$data))
  s <- summary(m)
  expect_s3_class(s, "summary.SAEforest")
})

test_that("point estimates of a learner are reproducible with a seed and the same across pools", {
  skip_if_not_installed("xgboost")
  m1 <- fit_lrn(learner = learner_xgboost(nrounds = 20), B = 0)
  m2 <- fit_lrn(learner = learner_xgboost(nrounds = 20), B = 0)
  expect_identical(m1$Indicators$Mean, m2$Indicators$Mean)
  expect_identical(m1$PopInfo$learner, "xgboost")
})

test_that("predefined learners fit and bootstrap on the example data", {
  for (l in list(if (requireNamespace("xgboost", quietly = TRUE)) learner_xgboost(nrounds = 20),
                 if (requireNamespace("lightgbm", quietly = TRUE)) learner_lightgbm(nrounds = 20),
                 if (requireNamespace("glmnet", quietly = TRUE)) learner_glmnet(),
                 if (requireNamespace("mgcv", quietly = TRUE)) learner_gam(k = 4))) {
    if (is.null(l)) next
    m <- fit_lrn(learner = l, B = 2, oof_folds = 3)
    expect_true(all(is.finite(m$Indicators$Mean)), info = l$name)
    expect_true(all(is.finite(m$MSE_Estimates$Mean)), info = l$name)
  }
})

test_that("an unbounded learner works with the Box-Cox transformation (widened tables)", {
  skip_if_not_installed("glmnet")
  data("eusilcA_smp", package = "SAEforestLearner", envir = environment())
  smp <- eusilcA_smp
  smp$ypos <- smp$eqIncome - min(smp$eqIncome) + 1e5
  m <- SAEforest_model(Y = smp$ypos, X = smp[, -c(1, 2, 16, 17, 18, ncol(smp))], dName = "district",
                       smp_data = smp, pop_data = eusilcA_pop, MSE = "nonparametric", B = 2, B_adj = 2,
                       seed = 11, cores = 2, verbose = FALSE, learner = learner_glmnet(), oof_folds = 3,
                       transform = "boxcox", bc_lambda = 0.5, bc_shift = 1)
  expect_false(m$MERFmodel$learner$bounded)
  expect_true(all(is.finite(m$Indicators$Mean)))
  expect_true(all(is.finite(m$MSE_Estimates$Mean)))
})

test_that("invalid learners are rejected", {
  expect_error(sae_learner("x", fit = 1, predict = function(m, x, t) 0), "functions")
  expect_error(fit_lrn(learner = list(name = "bad")), "sae_learner")
  expect_error(SAEforest_model(Y = eusilcA_smp$eqIncome, X = eusilcA_smp[, -c(1, 16, 17, 18)], dName = "district",
                               smp_data = eusilcA_smp, pop_data = eusilcA_pop, engine = "legacy",
                               learner = learner_ranger(), num.trees = 10), "learner")
})

test_that("learner_lm reproduces the nested-error linear mixed model (EBLUP)", {
  data("eusilcA_pop", package = "SAEforestLearner", envir = environment())
  data("eusilcA_smp", package = "SAEforestLearner", envir = environment())
  smp <- eusilcA_smp
  X <- smp[, -c(1, 2, 16, 17, 18)]
  m <- MERFranger(Y = smp$eqIncome, X = X, random = "(1|district)", data = smp, learner = learner_lm(),
                  ErrorTolerance = 1e-8, MaxIterations = 50)
  expect_identical(m$em_control$criterion, "robust")
  # exact leave-one-out predictions: compare with a refit without one unit
  i <- 7
  fit <- m$learner$fit(X, smp$eqIncome, list(), 1)
  f2 <- m$learner$fit(X[-i, ], smp$eqIncome[-i], list(), 1)
  expect_equal(fit$loo[i], m$learner$predict(f2, X[i, , drop = FALSE], 1), tolerance = 1e-8)
  # the fixed-point of the iterations is close to the ML fit of the nested-error model
  d <- cbind(y = smp$eqIncome, X, district = smp$district)
  lmm <- lme4::lmer(y ~ . - district + (1 | district), data = d, REML = FALSE)
  u_merf <- lme4::ranef(m$EffectModel)$district[, 1]
  u_lmm <- lme4::ranef(lmm)$district[, 1]
  expect_gt(cor(u_merf, u_lmm), 0.99)
  expect_equal(m$RanEffSD, as.data.frame(lme4::VarCorr(lmm))$sdcor[1], tolerance = 0.1)
  expect_equal(m$ErrorSD, stats::sigma(lmm), tolerance = 0.05)
})

test_that("the robust EM criterion keeps the best iteration and the original one is unchanged for ranger", {
  data("eusilcA_smp", package = "SAEforestLearner", envir = environment())
  smp <- eusilcA_smp
  X <- smp[, -c(1, 2, 16, 17, 18)]
  set.seed(1)
  m0 <- MERFranger(Y = smp$eqIncome, X = X, random = "(1|district)", data = smp, num.trees = 30, num.threads = 1)
  expect_identical(m0$em_control$criterion, "original")
  expect_identical(m0$BestIteration, m0$IterationsUsed)
  set.seed(1)
  m1 <- MERFranger(Y = smp$eqIncome, X = X, random = "(1|district)", data = smp, num.trees = 30, num.threads = 1,
                   em_criterion = "robust", em_control = list(patience = 2))
  expect_identical(m1$em_control$criterion, "robust")
  expect_true(m1$BestIteration <= m1$IterationsUsed)
  expect_equal(m1$LogLik[m1$BestIteration + 1], max(m1$LogLik[-1]))
})

test_that("spatial_folds balances folds and tune_learner runs a small grid (sequential and parallel)", {
  data("eusilcA_smp", package = "SAEforestLearner", envir = environment())
  smp <- eusilcA_smp
  X <- smp[, -c(1, 2, 16, 17, 18)]
  set.seed(2)
  cx <- runif(nrow(smp), 0, 100); cy <- runif(nrow(smp), 0, 100)
  f <- spatial_folds(cx, cy, block_size = 10, k = 4, seed = 1)
  expect_setequal(unique(f), 1:4)
  expect_lt(diff(range(table(f))), 0.25 * nrow(smp) / 4)
  tn <- tune_learner(Y = smp$eqIncome, X = X, dName = "district", smp_data = smp, learner_fun = learner_lm,
                     grid = NULL, folds = f, cores = 1, seed = 3)
  expect_equal(nrow(tn$results), 1)
  expect_true(is.finite(tn$best$RMSE))
  skip_if_not_installed("xgboost")
  g <- data.frame(nrounds = c(10, 20), eta = 0.3)
  tp <- tune_learner(Y = smp$eqIncome, X = X, dName = "district", smp_data = smp, learner_fun = learner_xgboost,
                     grid = g, folds = f, cores = 2, seed = 3, oof_folds = 3, MaxIterations = 3)
  expect_equal(nrow(tp$results), 2)
  expect_equal(nrow(tp$by_fold), 8)
  expect_true(all(is.finite(tp$results$RMSE)))
})
