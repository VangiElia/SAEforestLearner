# Tests of the scalable extension (SAEforestLearner)

test_that("flattened forests predict exactly as ranger", {
  data("eusilcA_smp")
  X <- eusilcA_smp[, -c(1, 16, 17, 18)]
  rf <- ranger::ranger(x = X, y = eusilcA_smp$eqIncome, num.trees = 30, seed = 1, num.threads = 1)
  ff <- flatten_forest(rf)
  expect_s3_class(ff, "flat_forest")
  p_ranger <- predict(rf, X, num.threads = 1)$predictions
  p_flat <- predict_flat_forest(ff, X, num_threads = 2)
  expect_identical(p_flat, p_ranger)
})

test_that("scalable point estimates equal the original implementation", {
  data("eusilcA_pop")
  data("eusilcA_smp")
  X <- eusilcA_smp[, -c(1, 16, 17, 18)]
  set.seed(5)
  m_old <- SAEforest_model(Y = eusilcA_smp$eqIncome, X = X, dName = "district",
                           smp_data = eusilcA_smp, pop_data = eusilcA_pop,
                           engine = "legacy", num.trees = 30, num.threads = 1)
  m_new <- SAEforest_model(Y = eusilcA_smp$eqIncome, X = X, dName = "district",
                           smp_data = eusilcA_smp, pop_data = eusilcA_pop, seed = 5,
                           cores = 1, num.trees = 30, verbose = FALSE)
  expect_equal(m_new$Indicators$Mean, m_old$Indicators$Mean, tolerance = 1e-10)
})

test_that("bootstrap MSE, totals and SAEaggregate are consistent", {
  data("eusilcA_pop")
  data("eusilcA_smp")
  args <- list(Y = eusilcA_smp$eqIncome, X = eusilcA_smp[, -c(1, 16, 17, 18)],
               dName = "district", smp_data = eusilcA_smp, pop_data = eusilcA_pop,
               MSE = "nonparametric", B = 3, B_adj = 3, seed = 2, cores = 2,
               num.trees = 50, total_factor = 2, verbose = FALSE)
  m_d <- do.call(SAEforest_model, args)
  expect_true(all(m_d$MSE_Estimates$Mean > 0))
  N <- as.numeric(table(eusilcA_pop$district)[m_d$Indicators$district])
  expect_equal(m_d$Indicators$Total, m_d$Indicators$Mean * N * 2)
  expect_equal(m_d$MSE_Estimates$Total, m_d$MSE_Estimates$Mean * (N * 2)^2)

  m_s <- do.call(SAEforest_model, c(args, aggregate_to = "state"))
  agg <- SAEaggregate(m_d, unique(eusilcA_pop[, c("district", "state")]), name = "state")
  agg <- agg[match(m_s$Indicators$state, agg$state), ]
  expect_equal(agg$Mean, m_s$Indicators$Mean)
  expect_equal(agg$Mean_MSE, m_s$MSE_Estimates$Mean)
  expect_equal(agg$Total_MSE, m_s$MSE_Estimates$Total)
})

test_that("raster populations are cached and estimated", {
  skip_if_not_installed("terra")
  set.seed(3)
  r <- terra::rast(nrows = 40, ncols = 50, nlyrs = 3, xmin = 0, xmax = 50, ymin = 0, ymax = 40,
                   crs = "EPSG:3035")
  terra::values(r) <- cbind(runif(2000), rnorm(2000), rep(c(1, NA, 1, 1), 500))
  names(r) <- c("x1", "x2", "mask")
  f <- file.path(tempdir(), "pop_test.tif")
  terra::writeRaster(r, f, overwrite = TRUE)
  dom <- terra::as.polygons(terra::rast(nrows = 2, ncols = 2, xmin = 0, xmax = 50, ymin = 0,
                                        ymax = 40, crs = "EPSG:3035", vals = 1:4))
  names(dom) <- "area"
  pop <- SAEpopulation(f, domains = dom, domain_field = "area", mask = "mask",
                       cache_dir = file.path(tempdir(), "pop_cache"), cores = 2,
                       rows_per_task = 10, overwrite = TRUE, verbose = FALSE)
  expect_s3_class(pop, "SAEpop")
  expect_equal(pop$N, 1500)
  expect_equal(pop$predictors, c("x1", "x2"))

  xy <- terra::spatSample(r, 200, cells = TRUE, na.rm = TRUE)
  smp <- data.frame(xy[, c("x1", "x2")])
  smp$area <- terra::extract(terra::rasterize(dom, r, field = "area"), xy$cell)[, 1]
  smp$y <- 10 + 5 * smp$x1 + smp$area + rnorm(200)
  m <- SAEforest_model(Y = smp$y, X = smp[, c("x1", "x2")], dName = "area", smp_data = smp,
                       pop_data = pop, MSE = "nonparametric", B = 3, B_adj = 3, cores = 2,
                       num.trees = 50, seed = 1, verbose = FALSE, total_factor = 1)
  expect_equal(nrow(m$Indicators), 4)
  expect_equal(sum(m$PopInfo$N_total), 1500)
  expect_true(all(is.finite(m$MSE_Estimates$Total)))
})
