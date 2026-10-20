# Scalable point and MSE estimates for domain-level means (and totals) ------------------
#
# Same estimator and nonparametric bootstrap as point_mean() / MSE_SAEforest_mean()
# (Krennmair & Schmid, 2022), reorganised so that memory does not grow with
# population size x bootstrap replications:
#   * the population is streamed in chunks (data.frame or raster cache, see SAEpopulation);
#   * only domain-level sums of forest predictions are kept;
#   * random effects are looked up from ranef() (identical to predict.merMod on new data);
#   * bootstrap errors of all population units are drawn in C++ one replicate at a time
#     and summed per domain, only the errors of the bootstrap-sampled units are kept;
#   * the B bootstrap MERFs are fitted in parallel on the (small) bootstrap samples and
#     flattened to disk; one pass over the population then predicts all of them with a
#     fast multi-threaded C++ traversal.
# The bootstrap population, samples and estimates follow the same distribution as in
# the original implementation.

SAEforest_mean_scalable <- function(Y,
                                    X,
                                    dName,
                                    smp_data,
                                    pop,
                                    MSE = "none",
                                    B = 100,
                                    B_adj = 100,
                                    importance = "impurity",
                                    initialRandomEffects = 0,
                                    ErrorTolerance = 0.0001,
                                    MaxIterations = 25,
                                    aggregate_to = NULL,
                                    na.rm = TRUE,
                                    total_factor = NULL,
                                    cores = NULL,
                                    chunk_rows = 1e6,
                                    forest_batch = 50,
                                    boot_dir = NULL,
                                    seed = NULL,
                                    verbose = TRUE,
                                    boot_hetero = c("none", "unit", "unit_re"),
                                    hetero_min_level = 5,
                                    transform = c("none", "boxcox"),
                                    bc_lambda = NULL,
                                    bc_shift = 1,
                                    bt_nodes = 20,
                                    bt_grid = 4000,
                                    boot_error_pool = c("auto", "domain", "all"),
                                    learner = NULL,
                                    oof_folds = 5,
                                    em_criterion = "auto",
                                    em_control = list(),
                                    out_call,
                                    ...) {
  boot_hetero <- match.arg(boot_hetero)
  transform <- match.arg(transform)
  boxcox <- transform == "boxcox"
  boot_error_pool <- match.arg(boot_error_pool)
  if (boot_error_pool == "auto") boot_error_pool <- if (boxcox) "all" else "domain"
  if (boxcox) {
    if (is.null(bc_lambda) || !is.numeric(bc_lambda) || length(bc_lambda) != 1 || bc_lambda < 0) {
      stop("transform = 'boxcox' needs a non-negative bc_lambda.")
    }
    if (boot_hetero != "none") stop("boot_hetero cannot be combined with transform = 'boxcox'.")
  }
  t_start <- Sys.time()
  timings <- list()
  tick <- function(what, t0) {
    timings[[what]] <<- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  }
  say <- function(...) if (verbose) message(stamp(), ...)

  cores <- cores %||% max(1L, parallel::detectCores() - 1L)
  if (!is.null(seed)) set.seed(seed)

  # learner of the fixed part (ranger by default): the ... of the call are its parameters
  dots <- list(...)
  threads_main <- if (is.null(dots$num.threads)) cores else dots$num.threads
  dots$num.threads <- NULL
  if (is.null(learner)) learner <- learner_ranger()
  learner <- learner_prepare(learner, params = dots, importance = importance)

  if (!inherits(pop, "SAEpop")) {
    pop <- SAEpopulation(as.data.frame(pop), domain_field = dName,
                         extra_fields = aggregate_to,
                         predictors = intersect(colnames(X), names(pop)))
  }
  if (!identical(pop$domain_field, dName)) {
    stop("The domain field of the population (", pop$domain_field, ") must equal dName (", dName, ").")
  }
  miss <- setdiff(colnames(X), pop_predictors(pop))
  if (length(miss) > 0) {
    stop("Covariates of X not available in the population: ", paste(utils::head(miss, 10), collapse = ", "))
  }

  # sample -------------------------------------------------------------------------------
  smp_data <- as.data.frame(smp_data)
  X <- as.data.frame(X)
  if (na.rm) {
    comp <- stats::complete.cases(X) & !is.na(Y) & !is.na(smp_data[[dName]])
    smp_data <- smp_data[comp, , drop = FALSE]
    X <- X[comp, , drop = FALSE]
    Y <- Y[comp]
  }
  smp_key0 <- dom_key(smp_data[[dName]])
  smp_data[[dName]] <- factor(smp_key0, levels = unique(smp_key0))
  order_in <- order(smp_data[[dName]])
  smp_data <- smp_data[order_in, , drop = FALSE]
  X <- X[order_in, , drop = FALSE]
  Y <- Y[order_in]
  smp_key <- as.character(smp_data[[dName]])
  if (boxcox) {
    if (any(Y + bc_shift <= 0)) stop("Y + bc_shift must be positive for the Box-Cox transformation.")
    Y <- bc_fwd(Y, bc_lambda, bc_shift)   # the model is fitted on the transformed scale
  }

  # population domains / aggregation cells ----------------------------------------------
  t0 <- Sys.time()
  say("Reading population domain structure ...")
  pd <- pop_domain_codes(pop, aggregate_to)
  dom_pop <- pd$code
  dom_levels <- pd$levels
  D <- length(dom_levels)
  N_pop <- length(dom_pop)
  N_d <- tabulate(dom_pop, nbins = D)

  if (is.null(aggregate_to)) {
    cell_pop <- dom_pop
    cell_dom <- seq_len(D)
    cell_grp <- seq_len(D)
    out_levels <- dom_levels
    out_name <- dName
  } else {
    G0 <- length(pd$agg_levels)
    key <- (as.numeric(dom_pop) - 1) * G0 + pd$agg_code
    ukey <- sort(unique(key))
    cell_pop <- match(key, ukey)
    cell_dom <- as.integer((ukey - 1) %/% G0) + 1L
    cell_grp <- as.integer((ukey - 1) %% G0) + 1L
    rm(key)
    out_levels <- pd$agg_levels
    out_name <- aggregate_to
  }
  C <- length(cell_dom)
  G <- length(out_levels)
  N_c <- tabulate(cell_pop, nbins = C)
  N_g <- as.numeric(cpp_group_sum(as.numeric(N_c), cell_grp, G))

  smp_dom <- match(smp_key, dom_levels)  # pop domain of each sample unit (NA: not in pop)
  if (anyNA(smp_dom)) {
    say(sum(is.na(smp_dom)), " sample units lie in ", length(unique(smp_key[is.na(smp_dom)])),
        " domains without population units: used for model fitting only.")
  }
  n_d <- tabulate(smp_dom[!is.na(smp_dom)], nbins = D)
  tick("population_structure", t0)
  say("Population: ", format(N_pop, big.mark = ","), " units, ", D, " domains (",
      sum(n_d > 0), " in sample, ", sum(n_d == 0), " out of sample); sample: ", length(Y), " units.")

  # point estimate: MERF on the sample -----------------------------------------------------
  t0 <- Sys.time()
  say("Fitting the MERF (learner '", learner$name, "') on the sample ...")
  random <- paste0("(1|", dName, ")")
  em_ctrl <- em_control_prepare(em_criterion, em_control, learner)
  unit_model <- MERF_learner(Y = Y, X = X, random = random, data = smp_data, learner = learner,
                             threads = threads_main, initialRandomEffects = initialRandomEffects,
                             ErrorTolerance = ErrorTolerance, MaxIterations = MaxIterations,
                             na.rm = FALSE, oof_folds = oof_folds, em_control = em_ctrl)
  tick("merf_fit", t0)
  fit_info <- if (is.function(learner$info)) learner$info(unit_model$Forest) else NULL
  say("MERF converged after ", unit_model$IterationsUsed, " iterations",
      if (unit_model$BestIteration != unit_model$IterationsUsed) paste0(" (best iteration ", unit_model$BestIteration, " kept)"),
      if (!is.null(fit_info$r.squared)) paste0(" (OOB R2 of forest: ", round(fit_info$r.squared, 3), ")"), ".")

  # compact form of the fitted learner used for all population predictions
  ff0 <- learner$unpack(learner$pack(unit_model$Forest))
  var_names <- colnames(X)
  u_hat <- ranef_lookup(unit_model$EffectModel, dName, dom_levels)

  # bootstrap sample units (drawn before the population pass to collect their covariates)
  do_mse <- MSE != "none"
  if (do_mse) {
    t0 <- Sys.time()
    ord_pop <- order(dom_pop)
    start_d <- c(0L, cumsum(N_d)[-D])
    samp_doms <- which(n_d > 0)
    if (any(n_d > N_d)) {
      say("Some domains have more sample than population units: bootstrap samples drawn with replacement there.")
    }
    boot_rows <- vector("list", B)
    for (b in seq_len(B)) {
      boot_rows[[b]] <- unlist(lapply(samp_doms, function(d) {
        ord_pop[start_d[d] + sample.int(N_d[d], n_d[d], replace = n_d[d] > N_d[d])]
      }), use.names = FALSE)
    }
    rm(ord_pop)
    U <- sort(unique(unlist(boot_rows, use.names = FALSE)))
    tick("boot_sampling", t0)
  } else {
    U <- integer(0)
  }

  # pass 1 over the population: point forest sums (+ covariates of bootstrap units) -----
  t0 <- Sys.time()
  chunks <- pop_chunks(pop, chunk_rows)
  say("Predicting the population (", nrow(chunks), " chunks) ...")
  S0 <- numeric(C)
  if (do_mse) {
    XU <- matrix(NA_real_, nrow = length(U), ncol = length(var_names),
                 dimnames = list(NULL, var_names))
    f0U <- numeric(length(U))
  }
  # heteroscedastic bootstrap errors and the Box-Cox bootstrap need the forest prediction of every unit
  f_pop <- if (do_mse && (boot_hetero != "none" || boxcox)) numeric(N_pop) else numeric(0)
  if (boxcox) {
    # back-transformation tables: classes of domain sample size, predictive SD per class
    n_classes <- sort(unique(n_d))
    k_dom <- match(n_d, n_classes)
    su0 <- unit_model$RanEffSD
    se0 <- unit_model$ErrorSD
    tau0 <- bt_tau(n_classes, su0, se0)
    # grid of the linear predictor: the sample range (tree ensembles cannot predict outside it);
    # unbounded learners (linear models, splines) extrapolate, so the grid is widened by the
    # in-sample prediction range plus half of it on each side
    rng0 <- bt_range(Y, learner, if (learner$bounded) NULL else learner$predict(ff0, X, threads_main), su0)
    tb0 <- bt_tables(rng0, tau0, bc_lambda, bc_shift, n_grid = bt_grid, n_nodes = bt_nodes)
    n_clamped <- 0
  }
  for (k in seq_len(nrow(chunks))) {
    ch <- chunks[k, ]
    xk <- pop_read_chunk(pop, ch, var_names)
    pk <- learner$predict(ff0, xk, cores)
    cells_k <- cell_pop[ch$gfrom:ch$gto]
    if (boxcox) {
      bs <- cpp_bt_sum(pk, dom_pop[ch$gfrom:ch$gto], u_hat, k_dom, tb0$tab, tb0$t0, tb0$dt, cells_k, C)
      S0 <- S0 + bs$sum
      n_clamped <- n_clamped + bs$clamped
    } else {
      S0 <- S0 + cpp_group_sum(pk, cells_k, C)
    }
    if (length(f_pop) > 0) f_pop[ch$gfrom:ch$gto] <- pk
    if (do_mse) {
      pos <- which(U >= ch$gfrom & U <= ch$gto)
      if (length(pos) > 0) {
        loc <- U[pos] - ch$gfrom + 1
        XU[pos, ] <- xk[loc, , drop = FALSE]
        f0U[pos] <- pk[loc]
      }
    }
    rm(xk, pk)
    if (verbose && (k %% 20 == 0 || k == nrow(chunks))) say("  chunk ", k, "/", nrow(chunks))
  }
  tick("population_pass_point", t0)

  if (boxcox) {
    if (n_clamped > 0) warning(n_clamped, " population units outside the back-transformation grid (clamped).")
    mean_g <- cpp_group_sum(S0, cell_grp, G) / N_g   # random effects already in the back-transformed sums
  } else {
    mean_g <- cpp_group_sum(S0 + N_c * u_hat[cell_dom], cell_grp, G) / N_g
  }
  Indicators <- data.frame(out_levels, Mean = mean_g, stringsAsFactors = FALSE)
  names(Indicators)[1] <- out_name
  N_total_g <- total_units(pop, dom_levels, cell_dom, cell_grp, N_c, G, aggregate_to)
  if (!is.null(total_factor)) {
    Indicators$Total <- Indicators$Mean * N_total_g * total_factor
  }

  data_specs <- list(
    N_surv = length(Y),
    N_pop = N_pop,
    D_out = sum(n_d == 0),
    D_in = sum(n_d > 0),
    D_total = D,
    ni_smp = table(smp_data[[dName]]),
    ni_pop = as.table(stats::setNames(N_d, dom_levels)),
    dName = if (is.null(aggregate_to)) dName else aggregate_to
  )

  pop_info <- list(N_units = stats::setNames(N_g, out_levels),
                   N_total = stats::setNames(N_total_g, out_levels),
                   total_factor = total_factor, cores = cores, chunk_rows = chunk_rows,
                   learner = learner$name, learner_args = learner$args, oof_folds = oof_folds,
                   em_control = em_ctrl,
                   transform = if (boxcox) list(type = "boxcox", lambda = bc_lambda, shift = bc_shift,
                                                nodes = bt_nodes, grid = bt_grid) else NULL)

  if (!do_mse) {
    tick("total", t_start)
    result <- list(
      MERFmodel = c(unit_model, call = out_call, data_specs = list(data_specs), data = list(smp_data)),
      Indicators = sortAlpha(Indicators, dName = out_name),
      MSE_Estimates = NULL,
      AdjustedSD = NULL,
      NrCovar = NULL,
      PopInfo = c(pop_info, list(timings = timings))
    )
    class(result) <- c("SAEforest_mean", "SAEforest")
    say("Done (", format(round(difftime(Sys.time(), t_start, units = "mins"), 1)), ").")
    return(result)
  }

  # MSE: bias-corrected error SD -----------------------------------------------------------
  t0 <- Sys.time()
  say("Error SD bootstrap (B_adj = ", B_adj, ") ...")
  # in-sample predictions of the learner (note: predict.ranger drew a seed from the RNG here in
  # the original code, so bootstrap random streams differ from SAEforest 1.0.1 by two draws)
  f_in <- learner$predict(ff0, X, threads_main)
  adj_SD <- adjust_ErrorSD_learner(Y = Y, X = X, mod = unit_model, B = B_adj, learner = learner,
                                   threads = threads_main, pred_in = f_in, oof_folds = oof_folds)
  if (!is.finite(adj_SD)) {
    warning("The bias correction of the error SD (Mendez & Lohr) exceeds the residual variance ",
            "(few trees or B_adj too small?): the uncorrected error SD of the MERF is used.")
    adj_SD <- unit_model$ErrorSD
  }
  tick("adjust_errorSD", t0)
  say("Adjusted error SD: ", signif(adj_SD, 5), " (model: ", signif(unit_model$ErrorSD, 5), ")")

  # random components as in ran_comp() -----------------------------------------------------
  ran_obj <- ran_comp(Y = Y, smp_data = smp_data, mod = unit_model, ADJsd = adj_SD, dName = dName,
                      pred_in = f_in)
  ran_effs <- ran_obj$ran_effs
  forest_res <- ran_obj$forest_res
  g_u <- rep(1, D)          # domain scale of the random effects (1: homoscedastic)
  scale_pop <- numeric(0)   # unit scale of the errors (empty: homoscedastic)
  hetero_info <- NULL
  if (boot_hetero != "none") {
    hc <- ran_comp_hetero(ran_obj, Y = Y, smp_data = smp_data, mod = unit_model, ADJsd = adj_SD,
                          dName = dName, min_level = hetero_min_level)
    forest_res <- hc$pool
    scale_pop <- hc$sfun(f_pop)
    hetero_info <- hc[c("power", "const", "min_level", "sd_by_decile")]
    say(sprintf("Heteroscedastic bootstrap errors: SD proportional to level^%.2f.", hc$power))
    if (boot_hetero == "unit_re") {
      # random effects proportional to the synthetic domain mean, scaled to keep sigma_u
      lev_d <- pmax(as.numeric(cpp_group_sum(S0, cell_dom, D)) / pmax(N_d, 1), hetero_min_level)
      g_u <- lev_d / sqrt(mean(lev_d[n_d > 0]^2))
      rb <- hc$r_bar
      re_std <- rb / g_u[match(names(rb), dom_levels)]
      re_std <- re_std[is.finite(re_std)]
      re_std <- re_std / stats::sd(re_std) * unit_model$RanEffSD
      ran_effs <- as.numeric(re_std - mean(re_std))
      hetero_info$re_scale_range <- range(g_u)
      say(sprintf("Random effects proportional to the domain level (scale %.2f-%.2f).", min(g_u), max(g_u)))
    }
  }

  # error pools: own residuals for sampled domains, all residuals otherwise ("domain");
  # all residuals for every domain ("all")
  in_pop <- !is.na(smp_dom)
  o <- order(smp_dom[in_pop])
  pool <- c(forest_res[in_pop][o], forest_res)
  own <- n_d > 0 & boot_error_pool == "domain"
  K <- ifelse(own, n_d, length(forest_res))
  offset <- ifelse(own, c(0L, cumsum(n_d)[-D]), sum(in_pop))
  K <- as.integer(K)
  offset <- as.integer(offset)

  # bootstrap populations (domain sums) and samples ---------------------------------------
  t0 <- Sys.time()
  say("Generating ", B, " bootstrap populations and samples ...")
  tau_star <- matrix(NA_real_, nrow = G, ncol = B)
  boot_tasks <- vector("list", B)
  task_seeds <- sample.int(.Machine$integer.max, B)
  boot_dir <- boot_dir %||% file.path(tempdir(), paste0("SAEforest_boot_", Sys.getpid(),
                                                        "_", format(Sys.time(), "%H%M%S")))
  dir.create(boot_dir, recursive = TRUE, showWarnings = FALSE)
  xcols <- colnames(X)

  y_range_b <- matrix(NA_real_, nrow = B, ncol = 2)
  for (b in seq_len(B)) {
    u_b <- sample(ran_effs, size = D, replace = TRUE) * g_u
    rows_b <- boot_rows[[b]]
    o_b <- order(rows_b)
    if (boxcox) {
      # bootstrap population on the transformed scale, true domain sums on the original scale
      eb <- cpp_boot_errors(dom_pop, cell_pop, C, offset, K, pool, as.integer(rows_b[o_b]), numeric(0),
                            f_pop, u_b, bc_lambda, bc_shift)
      tau_c <- eb$E
    } else {
      eb <- cpp_boot_errors(dom_pop, cell_pop, C, offset, K, pool, as.integer(rows_b[o_b]), scale_pop)
      tau_c <- S0 + N_c * u_b[cell_dom] + eb$E
    }
    e_rows <- numeric(length(rows_b))
    e_rows[o_b] <- eb$e_sel
    tau_star[, b] <- cpp_group_sum(tau_c, cell_grp, G) / N_g

    iu <- match(rows_b, U)
    d_b <- dom_pop[rows_b]
    y_b <- f0U[iu] + u_b[d_b] + e_rows
    y_range_b[b, ] <- range(y_b)
    if (anyNA(y_b) || anyNA(tau_star[, b])) {
      stop("Missing values in bootstrap replicate ", b, " (random effects: ", sum(is.na(u_b)),
           ", unit errors: ", sum(is.na(e_rows)), ", domain means: ", sum(is.na(tau_star[, b])), ").")
    }
    x_b <- as.data.frame(XU[iu, xcols, drop = FALSE], optional = TRUE)
    boot_tasks[[b]] <- list(b = b, seed = task_seeds[b], y = y_b, X = x_b,
                            dom = dom_levels[d_b], file = file.path(boot_dir, sprintf("boot_%04d.rds", b)))
  }
  rm(XU, f0U, boot_rows, scale_pop, f_pop)
  tick("boot_populations", t0)

  # fit the bootstrap MERFs in parallel -----------------------------------------------------
  t0 <- Sys.time()
  n_workers <- max(1L, min(cores, B))
  threads_w <- max(1L, floor(cores / n_workers))
  say("Fitting ", B, " bootstrap MERFs on ", n_workers, " processes ...")
  learner_boot <- learner_detach(learner)   # self-contained functions for the workers
  learner_boot$args$importance <- "none"
  common <- list(dName = dName, random = random, learner = learner_boot, threads = threads_w,
                 initialRandomEffects = initialRandomEffects, oof_folds = oof_folds, em_control = em_ctrl,
                 ErrorTolerance = ErrorTolerance, MaxIterations = MaxIterations)

  cl <- parallel::makePSOCKcluster(n_workers)
  on.exit(parallel::stopCluster(cl), add = TRUE)
  parallel::clusterCall(cl, detach_env(boot_worker_init), detach_env(MERF_learner), common)
  boot_fits <- pbapply::pblapply(boot_tasks, detach_env(boot_worker_fit), cl = cl)
  parallel::stopCluster(cl)
  on.exit()
  rm(boot_tasks)
  tick("boot_merf_fits", t0)
  iters <- vapply(boot_fits, `[[`, numeric(1), "iter")
  say("Bootstrap MERFs done (iterations: median ", stats::median(iters), ", max ", max(iters), ").")

  # pass 2 over the population: bootstrap forest sums ---------------------------------------
  t0 <- Sys.time()
  S_b <- matrix(0, nrow = C, ncol = B)
  batches <- split(seq_len(B), ceiling(seq_len(B) / max(1, forest_batch)))
  total_steps <- B * nrow(chunks)
  done_steps <- 0
  for (ib in seq_along(batches)) {
    bb <- batches[[ib]]
    say("Predicting the population with bootstrap forests ", min(bb), "-", max(bb), " of ", B, " ...")
    t_batch <- Sys.time()
    forests <- lapply(boot_fits[bb], function(z) learner$unpack(readRDS(z$file)))
    if (boxcox) {
      # back-transformation of each bootstrap estimator: its own random effects and variance components
      bt_b <- lapply(bb, function(b) {
        re <- boot_fits[[b]]$ranef
        u_hat_b <- re[match(dom_levels, names(re))]
        u_hat_b[is.na(u_hat_b)] <- 0
        rng <- bt_range(y_range_b[b, ], learner, boot_fits[[b]]$pred_range, boot_fits[[b]]$sigma_u)
        list(u = u_hat_b, tb = bt_tables(rng, bt_tau(n_classes, boot_fits[[b]]$sigma_u, boot_fits[[b]]$sigma_e),
                                          bc_lambda, bc_shift, n_grid = bt_grid, n_nodes = bt_nodes))
      })
    }
    for (k in seq_len(nrow(chunks))) {
      ch <- chunks[k, ]
      xk <- pop_read_chunk(pop, ch, var_names)
      cells_k <- cell_pop[ch$gfrom:ch$gto]
      if (boxcox) dk <- dom_pop[ch$gfrom:ch$gto]
      for (j in seq_along(bb)) {
        pk <- learner$predict(forests[[j]], xk, cores)
        if (boxcox) {
          tb <- bt_b[[j]]$tb
          S_b[, bb[j]] <- S_b[, bb[j]] + cpp_bt_sum(pk, dk, bt_b[[j]]$u, k_dom, tb$tab, tb$t0, tb$dt, cells_k, C)$sum
        } else {
          S_b[, bb[j]] <- S_b[, bb[j]] + cpp_group_sum(pk, cells_k, C)
        }
      }
      rm(xk)
      if (verbose) {
        done_steps <- done_steps + length(bb)
        say(progress_bar(done_steps, total_steps, t0,
                         sprintf("batch %d/%d chunk %d/%d", ib, length(batches), k, nrow(chunks))))
      }
    }
    rm(forests)
    gc(verbose = FALSE)
  }
  tick("population_pass_boot", t0)
  unlink(boot_dir, recursive = TRUE)

  # bootstrap estimates and MSE --------------------------------------------------------------
  tau_b <- matrix(NA_real_, nrow = G, ncol = B)
  for (b in seq_len(B)) {
    if (boxcox) {   # random effects already in the back-transformed sums
      tau_b[, b] <- cpp_group_sum(S_b[, b], cell_grp, G) / N_g
      next
    }
    re <- boot_fits[[b]]$ranef
    u_hat_b <- re[match(dom_levels, names(re))]
    u_hat_b[is.na(u_hat_b)] <- 0
    tau_b[, b] <- cpp_group_sum(S_b[, b] + N_c * u_hat_b[cell_dom], cell_grp, G) / N_g
  }
  mse_g <- rowMeans((tau_star - tau_b)^2)

  MSE_Estimates <- data.frame(out_levels, Mean = mse_g, stringsAsFactors = FALSE)
  names(MSE_Estimates)[1] <- out_name
  if (!is.null(total_factor)) {
    MSE_Estimates$Total <- mse_g * (N_total_g * total_factor)^2
  }

  o_out <- order(as.character(out_levels))
  tick("total", t_start)
  result <- list(
    MERFmodel = c(unit_model, call = out_call, data_specs = list(data_specs), data = list(smp_data)),
    Indicators = sortAlpha(Indicators, dName = out_name),
    MSE_Estimates = sortAlpha(MSE_Estimates, dName = out_name),
    AdjustedSD = adj_SD,
    NrCovar = NULL,
    PopInfo = c(pop_info, list(timings = timings, boot_iterations = iters,
                               boot_hetero = boot_hetero, hetero = hetero_info,
                               boot_error_pool = boot_error_pool)),
    # bootstrap replicates of the domain means (true and estimated), rows in the order
    # of Indicators: allow exact bootstrap MSEs of any aggregation (see SAEaggregate)
    Bootstrap = list(tau_star = tau_star[o_out, , drop = FALSE], tau_b = tau_b[o_out, , drop = FALSE],
                     domains = out_levels[o_out])
  )
  class(result) <- c("SAEforest_mean", "SAEforest")
  say("Done (", format(round(difftime(Sys.time(), t_start, units = "mins"), 1)), ").")
  result
}


# Random intercepts of the domains in 'levels' (0 for domains not in the model) ---------
# Equivalent to predict(EffectModel, newdata, allow.new.levels = TRUE): the offset is not
# used for new data and the model has no fixed part.
ranef_lookup <- function(effect_model, dName, levels) {
  re <- lme4::ranef(effect_model)[[dName]]
  out <- re[match(levels, dom_key(rownames(re))), 1]
  out[is.na(out)] <- 0
  out
}


# Box-Cox transformation (transform = "boxcox") -------------------------------------------
# T(y) = ((y + shift)^lambda - 1) / lambda (log(y + shift) for lambda = 0), fitted on the
# transformed scale. Domain means on the original scale are EBP-type predictions
#   mean_j E[ Tinv(eta_ij + u_i + v) ],  v ~ N(0, tau_i^2),
# with u_i the predicted random effect (0 for out-of-sample domains) and
# tau_i^2 = sigma_e^2 + sigma_u^2 (1 - gamma_i), gamma_i = sigma_u^2 / (sigma_u^2 + sigma_e^2 / n_i)
# (gamma = 0 out of sample), i.e. the conditional predictive variance of a unit given the sample.
# The expectation is computed by Gauss-Hermite quadrature, tabulated on a fine grid of the
# linear predictor for each distinct tau (domains differ only through n_i) and interpolated
# linearly for every population unit.
# range of the linear predictor covered by the back-transformation tables: the target range plus
# 6 sigma_u + 1 (as in the original Box-Cox extension); for unbounded learners also the range of the
# predictions (pred_range: in-sample predictions, or their range from the workers) plus half of it
bt_range <- function(y_range, learner, pred_range, sigma_u) {
  r <- range(y_range) + c(-1, 1) * (6 * sigma_u + 1)
  if (!learner$bounded && !is.null(pred_range)) {
    pr <- range(pred_range)
    pr <- pr + c(-1, 1) * 0.5 * diff(pr)
    r <- range(c(r, pr))
  }
  r
}

bc_fwd <- function(y, lambda, shift) {
  if (lambda == 0) log(y + shift) else ((y + shift)^lambda - 1) / lambda
}
bc_inv <- function(t, lambda, shift) {
  if (lambda == 0) return(exp(t) - shift)
  pmax(lambda * t + 1, 0)^(1 / lambda) - shift
}
# nodes and weights of Gauss-Hermite quadrature for the standard normal (Golub-Welsch)
gh_normal <- function(n) {
  i <- seq_len(n - 1)
  J <- matrix(0, n, n)
  J[cbind(i, i + 1)] <- sqrt(i)
  J[cbind(i + 1, i)] <- sqrt(i)
  e <- eigen(J, symmetric = TRUE)
  list(nodes = e$values, weights = e$vectors[1, ]^2)
}
# tables h_k(t) = E[Tinv(t + tau_k Z)] on a regular grid of t, one column per tau
bt_tables <- function(range_t, tau, lambda, shift, n_grid = 4000, n_nodes = 20) {
  gh <- gh_normal(n_nodes)
  grid <- seq(range_t[1], range_t[2], length.out = n_grid)
  tab <- vapply(tau, function(s) {
    acc <- numeric(n_grid)
    for (q in seq_along(gh$nodes)) acc <- acc + gh$weights[q] * bc_inv(grid + s * gh$nodes[q], lambda, shift)
    acc
  }, numeric(n_grid))
  list(t0 = grid[1], dt = grid[2] - grid[1], n = n_grid, tab = matrix(tab, nrow = n_grid))
}
# interpolate h_k(t) for vectors t and k (class index); points outside the grid are clamped
bt_eval <- function(t, k, tb) {
  x <- (t - tb$t0) / tb$dt + 1
  x <- pmin(pmax(x, 1), tb$n - 1e-9)
  i <- floor(x)
  w <- x - i
  tb$tab[cbind(i, k)] * (1 - w) + tb$tab[cbind(i + 1L, k)] * w
}
# predictive SD of a unit given the sample, per class of domain sample size
bt_tau <- function(n_class, sigma_u, sigma_e) {
  gam <- ifelse(n_class > 0, sigma_u^2 / (sigma_u^2 + sigma_e^2 / pmax(n_class, 1)), 0)
  sqrt(sigma_e^2 + sigma_u^2 * (1 - gam))
}


# Heteroscedastic random components for the bootstrap (boot_hetero != "none") ------------
# The SD of the unit errors is modelled as a power of the level, SD = const * level^power:
# 'power' is estimated by regressing log(SD) on log(mean level) over deciles of the MERF
# fitted values (OOB forest + random effect), using the OOB residuals. The unit residuals
# of ran_comp() are standardised by level^power (level = in-bag forest prediction, as in
# ran_comp), centred and scaled to unit variance; 'const' makes the mean squared error
# over the sample equal to ADJsd^2 (Mendez & Lohr correction). Bootstrap errors are then
# const * z * level^power, with the level of each population unit given by its forest
# prediction. Domain residual means (r_bar) are returned for level-proportional random
# effects.
ran_comp_hetero <- function(ran_obj, Y, smp_data, mod, ADJsd, dName, min_level = 5) {
  sd_ <- ran_obj$smp_data
  lev_ib <- Y - sd_$forest_res
  u <- ranef_lookup(mod$EffectModel, dName, dom_key(smp_data[[dName]]))
  lev_oob <- mod$OOBpred + u
  r_oob <- Y - lev_oob
  br <- unique(stats::quantile(lev_oob, seq(0, 1, 0.1)))
  bin <- cut(lev_oob, br, include.lowest = TRUE)
  sdb <- data.frame(level = as.numeric(tapply(lev_oob, bin, mean)),
                    sd = as.numeric(tapply(r_oob, bin, stats::sd)))
  sdb <- sdb[is.finite(sdb$sd) & sdb$sd > 0 & sdb$level > 0, , drop = FALSE]
  if (nrow(sdb) < 3) stop("Too few distinct fitted values to estimate the variance function.")
  power <- unname(stats::coef(stats::lm(log(sd) ~ log(level), data = sdb))[2])
  sfun <- function(l) pmax(l, min_level)^power
  z <- sd_$forest_eij / sfun(lev_ib)
  z <- (z - mean(z)) / stats::sd(z)
  const <- ADJsd / sqrt(mean(sfun(lev_ib)^2))
  r_bar <- tapply(sd_$forest_res, as.character(sd_[[dName]]), mean)
  r_bar <- stats::setNames(as.numeric(r_bar), names(r_bar))   # plain vector, not a 1-d array
  list(pool = const * z, sfun = sfun, power = power, const = const, min_level = min_level,
       sd_by_decile = sdb, r_bar = r_bar)
}


# Number of units used for totals (mask units incl. those with missing covariates) ------
total_units <- function(pop, dom_levels, cell_dom, cell_grp, N_c, G, aggregate_to) {
  N_units <- as.numeric(cpp_group_sum(as.numeric(N_c), cell_grp, G))
  if (identical(pop$type, "memory") || is.null(pop$unit_counts) || sum(pop$unit_counts$N_na) == 0) {
    return(N_units)
  }
  # scale each domain's cells by N_mask / N so that units without covariates are counted
  uc <- pop$unit_counts
  keys <- dom_key(pop$lookup[[pop$domain_field]][match(uc$.code, pop$lookup$.code)])
  ratio <- (uc$N_mask / uc$N)[match(dom_levels, keys)]
  ratio[!is.finite(ratio)] <- 1
  as.numeric(cpp_group_sum(as.numeric(N_c) * ratio[cell_dom], cell_grp, G))
}


# Worker side of the bootstrap MERF fits (run in PSOCK workers) -------------------------

boot_worker_init <- function(merf_fun, common) {
  suppressPackageStartupMessages({
    requireNamespace("lme4")
    for (p in common$learner$packages) requireNamespace(p)
  })
  assign(".SAEboot", list(merf = merf_fun, common = common), envir = globalenv())
  invisible(TRUE)
}

boot_worker_fit <- function(task) {
  env <- get(".SAEboot", envir = globalenv())
  cm <- env$common
  set.seed(task$seed)
  d <- task$X
  d[[cm$dName]] <- factor(task$dom, levels = unique(task$dom))
  m <- env$merf(Y = task$y, X = task$X, random = cm$random, data = d, learner = cm$learner,
                threads = cm$threads, initialRandomEffects = cm$initialRandomEffects,
                ErrorTolerance = cm$ErrorTolerance, MaxIterations = cm$MaxIterations,
                na.rm = FALSE, oof_folds = cm$oof_folds, em_control = cm$em_control)
  re <- lme4::ranef(m$EffectModel)[[cm$dName]]
  saveRDS(cm$learner$pack(m$Forest), task$file, compress = FALSE)
  list(b = task$b, ranef = stats::setNames(re[, 1], rownames(re)), iter = m$IterationsUsed,
       sigma_u = m$RanEffSD, sigma_e = m$ErrorSD, file = task$file,
       pred_range = if (cm$learner$bounded) NULL else range(cm$learner$predict(m$Forest, task$X, cm$threads)))
}


# One-line text progress bar with elapsed time and ETA (log friendly) --------------------
progress_bar <- function(done, total, t0, label = "", width = 30) {
  frac <- min(1, done / total)
  el <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
  eta <- if (frac > 0) el / frac * (1 - frac) else NA
  n <- round(frac * width)
  sprintf("[%s%s] %5.1f%% | %s | elapsed %.1f min | ETA %.1f min",
          strrep("#", n), strrep("-", width - n), 100 * frac, label, el, eta)
}
