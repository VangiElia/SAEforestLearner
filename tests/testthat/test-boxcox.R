# Tests of the Box-Cox transformation, the bootstrap error pool and the heteroscedastic bootstrap
# (scalable engine)

fit_eusilc <- function(...) {
  data("eusilcA_pop", package = "SAEforestLearner", envir = environment())
  data("eusilcA_smp", package = "SAEforestLearner", envir = environment())
  smp <- eusilcA_smp
  # income shifted far above 0, so that the truncation of the back-transformation never acts
  # (Y must be a column of the sample data)
  smp$ypos <- smp$eqIncome - min(smp$eqIncome) + 1e5
  SAEforest_model(Y = smp$ypos, X = smp[, -c(1, 16, 17, 18, ncol(smp))], dName = "district",
                  smp_data = smp, pop_data = eusilcA_pop, MSE = "nonparametric",
                  B = 3, B_adj = 3, num.trees = 30, seed = 11, cores = 2, verbose = FALSE, ...)
}

test_that("Box-Cox with lambda = 1, shift = 1 reproduces the untransformed model", {
  m0 <- fit_eusilc(boot_error_pool = "all")
  m1 <- fit_eusilc(transform = "boxcox", bc_lambda = 1, bc_shift = 1)
  expect_identical(m1$PopInfo$boot_error_pool, "all")
  expect_equal(m1$Indicators$Mean, m0$Indicators$Mean, tolerance = 1e-10)
  expect_equal(m1$MSE_Estimates$Mean, m0$MSE_Estimates$Mean, tolerance = 1e-8)
})

test_that("the default error pool without transformation is the original one", {
  m_def <- fit_eusilc()
  m_dom <- fit_eusilc(boot_error_pool = "domain")
  expect_identical(m_def$PopInfo$boot_error_pool, "domain")
  expect_identical(m_def$MSE_Estimates$Mean, m_dom$MSE_Estimates$Mean)
})

test_that("Box-Cox estimates are finite and do not depend on the error pool", {
  m_all <- fit_eusilc(transform = "boxcox", bc_lambda = 0.35)
  m_dom <- fit_eusilc(transform = "boxcox", bc_lambda = 0.35, boot_error_pool = "domain")
  expect_identical(m_all$PopInfo$transform$type, "boxcox")
  expect_true(all(is.finite(m_all$Indicators$Mean)))
  expect_true(all(is.finite(m_all$MSE_Estimates$Mean) & m_all$MSE_Estimates$Mean > 0))
  expect_identical(m_all$Indicators$Mean, m_dom$Indicators$Mean)
})

test_that("the heteroscedastic bootstrap leaves the point estimates unchanged", {
  m0 <- fit_eusilc()
  m_h <- fit_eusilc(boot_hetero = "unit")
  expect_identical(m_h$Indicators$Mean, m0$Indicators$Mean)
  expect_true(all(is.finite(m_h$MSE_Estimates$Mean) & m_h$MSE_Estimates$Mean > 0))
})

test_that("invalid Box-Cox settings are rejected", {
  expect_error(fit_eusilc(transform = "boxcox"), "bc_lambda")
  expect_error(fit_eusilc(transform = "boxcox", bc_lambda = 0.5, boot_hetero = "unit"), "boot_hetero")
})
