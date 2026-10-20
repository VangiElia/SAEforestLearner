# Population data for the scalable MERF estimators ---------------------------------------
#
# The population (unit-level covariates of every domain unit) can be a data.frame in
# memory or a raster. A raster is converted once into an on-disk cache of the valid
# units (mask cells with a domain and complete covariates), stored as 'fst' parts that
# are streamed chunk-wise by the estimators, so memory use does not grow with the
# population size.

`%||%` <- function(a, b) if (is.null(a)) b else a


#' Population data for SAEforest_model
#'
#' Builds an object of class \code{SAEpop} describing the unit-level population. For
#' rasters, all valid cells (inside \code{mask}, inside a domain polygon and with complete
#' covariates) are extracted in parallel, row block by row block, into an on-disk cache
#' (\code{cache_dir}) that is later streamed in chunks by \code{\link{SAEforest_model}}.
#' The extraction runs once; an existing cache is re-opened by passing its directory.
#'
#' @param x A \code{terra::SpatRaster} (file-backed) or raster file name whose layers are
#' the covariates; a data.frame of unit-level population data; or the directory of an
#' existing population cache.
#' @param domains Domain definition for rasters: polygons (\code{sf}, \code{SpatVector}
#' or vector file name) with the attribute \code{domain_field}; a single-layer raster
#' (or file name) with integer domain codes on the grid of \code{x}; or the name of a
#' layer of \code{x}. Polygons are rasterized on the grid of \code{x} by cell centre.
#' @param domain_field Name of the domain identifier. It becomes the domain column of
#' the population and must equal \code{dName} in \code{SAEforest_model}.
#' @param mask Optional name of a layer of \code{x}; only cells where it equals
#' \code{mask_value} belong to the population.
#' @param mask_value Value of \code{mask} identifying population cells.
#' @param predictors Covariate layers (default: all layers except \code{mask}, a domain
#' layer and \code{exclude}).
#' @param exclude Layers not used as covariates.
#' @param extra_fields Further polygon attributes (e.g. region codes) stored per domain;
#' they can be used as \code{aggregate_to} in \code{SAEforest_model}.
#' @param cache_dir Directory of the cache (created; must be empty unless
#' \code{overwrite = TRUE}).
#' @param cores Number of parallel processes for the extraction.
#' @param rows_per_task,rows_per_read Raster rows per parallel task and per read.
#' @param compress Compression level of the fst parts (0-100).
#' @param overwrite Logical, overwrite an existing cache.
#' @param verbose Logical, report progress.
#'
#' @return An object of class \code{SAEpop}.
#' @export
SAEpopulation <- function(x, domains = NULL, domain_field = NULL, mask = NULL, mask_value = 1,
                          predictors = NULL, exclude = NULL, extra_fields = NULL,
                          cache_dir = NULL, cores = NULL, rows_per_task = 64L,
                          rows_per_read = 16L, compress = 50, overwrite = FALSE,
                          verbose = TRUE) {

  # re-open an existing cache
  if (is.character(x) && length(x) == 1 && dir.exists(x)) {
    return(open_pop_cache(x))
  }

  # in-memory population
  if (is.data.frame(x)) {
    if (is.null(domain_field) || !domain_field %in% names(x)) {
      stop("For data.frame populations, domain_field must be a column of x.")
    }
    predictors <- predictors %||% setdiff(names(x), c(domain_field, extra_fields, exclude))
    miss <- setdiff(c(predictors, extra_fields), names(x))
    if (length(miss) > 0) stop("Columns not found in x: ", paste(miss, collapse = ", "))
    return(structure(list(type = "memory", data = x, domain_field = domain_field,
                          predictors = predictors, extra_fields = extra_fields),
                     class = "SAEpop"))
  }

  # raster population -> disk cache
  build_pop_cache(x = x, domains = domains, domain_field = domain_field, mask = mask,
                  mask_value = mask_value, predictors = predictors, exclude = exclude,
                  extra_fields = extra_fields, cache_dir = cache_dir, cores = cores,
                  rows_per_task = rows_per_task, rows_per_read = rows_per_read,
                  compress = compress, overwrite = overwrite, verbose = verbose)
}


# Build the raster cache -------------------------------------------------------------------

build_pop_cache <- function(x, domains, domain_field, mask, mask_value, predictors, exclude,
                            extra_fields, cache_dir, cores, rows_per_task, rows_per_read,
                            compress, overwrite, verbose, gdal_cache_mb = 512) {
  t_start <- Sys.time()
  terra::setGDALconfig("GDAL_CACHEMAX", as.character(gdal_cache_mb))
  if (is.null(cache_dir)) stop("cache_dir is required to build a population cache from a raster.")
  if (is.null(domain_field)) stop("domain_field is required.")

  r <- if (inherits(x, "SpatRaster")) x else terra::rast(x)
  src <- terra::sources(r)
  if (any(src == "")) {
    stop("The raster must be file-backed (write it to disk with terra::writeRaster first).")
  }
  lyr <- names(r)
  if (anyDuplicated(lyr)) stop("Raster layer names must be unique.")

  if (dir.exists(cache_dir) && length(list.files(cache_dir)) > 0) {
    if (!overwrite) {
      stop("cache_dir is not empty. Use overwrite = TRUE or open it with SAEpopulation(cache_dir).")
    }
    unlink(file.path(cache_dir, c("parts", "meta.rds", "domains.tif", "domains.tif.aux.xml")),
           recursive = TRUE)
  }
  dir.create(file.path(cache_dir, "parts"), recursive = TRUE, showWarnings = FALSE)

  # domains ------------------------------------------------------------------------------
  dom_layer <- NULL
  dom_file <- NULL
  lookup <- NULL

  if (is.character(domains) && length(domains) == 1 && domains %in% lyr) {
    dom_layer <- domains
  } else if (inherits(domains, c("sf", "sfc", "SpatVector")) ||
             (is.character(domains) && length(domains) == 1 &&
              grepl("\\.(shp|gpkg|geojson|json|sqlite|gdb)$", domains, ignore.case = TRUE))) {
    v <- if (inherits(domains, "SpatVector")) domains else terra::vect(domains)
    if (!domain_field %in% names(v)) stop("domain_field not found in the domain polygons.")
    miss <- setdiff(extra_fields, names(v))
    if (length(miss) > 0) stop("extra_fields not found in the domain polygons: ", paste(miss, collapse = ", "))
    v <- terra::project(v, terra::crs(r))
    att <- as.data.frame(v)[, c(domain_field, extra_fields), drop = FALSE]
    key <- att[[domain_field]]
    if (is.numeric(key) && all(key == round(key), na.rm = TRUE) &&
        all(abs(key) < .Machine$integer.max, na.rm = TRUE)) {
      code <- as.integer(key)
    } else {
      code <- match(as.character(key), unique(as.character(key)))
    }
    v$.code <- code
    lookup <- att
    lookup$.code <- code
    lookup <- lookup[!duplicated(lookup$.code) & !is.na(lookup$.code), c(".code", domain_field, extra_fields)]
    if (verbose) message(stamp(), "Rasterizing ", nrow(lookup), " domains on the population grid ...")
    dom_file <- file.path(cache_dir, "domains.tif")
    terra::rasterize(v, r[[1]], field = ".code", touches = FALSE, filename = dom_file,
                     datatype = "INT4S", overwrite = TRUE)
  } else if (inherits(domains, "SpatRaster") || (is.character(domains) && length(domains) == 1)) {
    dr <- if (inherits(domains, "SpatRaster")) domains else terra::rast(domains)
    if (terra::nlyr(dr) != 1) stop("The domain raster must have a single layer.")
    terra::compareGeom(r, dr)
    dom_file <- terra::sources(dr)
    if (dom_file == "") stop("The domain raster must be file-backed.")
  } else {
    stop("domains must be polygons, a domain raster or the name of a layer of x.")
  }

  predictors <- predictors %||% setdiff(lyr, c(mask, dom_layer, exclude))
  miss <- setdiff(c(predictors, mask), lyr)
  if (length(miss) > 0) stop("Layers not found in the raster: ", paste(miss, collapse = ", "))
  if (domain_field %in% predictors) stop("domain_field must not be a covariate.")

  # parallel extraction by row blocks ------------------------------------------------------
  nr <- terra::nrow(r)
  nc <- terra::ncol(r)
  starts <- seq(1L, nr, by = rows_per_task)
  cores <- cores %||% max(1L, min(16L, parallel::detectCores() - 1L))
  cores <- max(1L, min(cores, length(starts)))

  if (verbose) {
    message(stamp(), "Extracting ", length(predictors), " covariates from ", nr, " x ", nc,
            " cells in ", length(starts), " row blocks with ", cores, " processes ...")
  }

  job <- list(src = src, lyr = lyr, read_lyr = c(mask, predictors, dom_layer), mask = mask,
              mask_value = mask_value, predictors = predictors, dom_layer = dom_layer,
              dom_file = dom_file, nr = nr, nc = nc, rows_per_task = rows_per_task,
              rows_per_read = rows_per_read, out_dir = file.path(cache_dir, "parts"),
              compress = compress, gdal_cache_mb = gdal_cache_mb)

  cl <- parallel::makePSOCKcluster(cores)
  on.exit(parallel::stopCluster(cl), add = TRUE)
  parallel::clusterCall(cl, detach_env(pop_worker_init), job)
  res <- pbapply::pblapply(starts, detach_env(pop_worker_task), cl = cl)

  # assemble meta data ---------------------------------------------------------------------
  parts <- data.frame(
    row0 = vapply(res, `[[`, integer(1), "row0"),
    file = vapply(res, `[[`, character(1), "file"),
    n = vapply(res, `[[`, numeric(1), "n"),
    stringsAsFactors = FALSE
  )
  parts <- parts[parts$n > 0, ]
  parts <- parts[order(parts$row0), ]
  parts$offset <- c(0, cumsum(parts$n)[-nrow(parts)])
  parts$file <- basename(parts$file)
  rownames(parts) <- NULL

  sum_counts <- function(name) {
    tb <- do.call(rbind, lapply(res, `[[`, name))
    if (is.null(tb) || nrow(tb) == 0) return(data.frame(code = integer(0), n = numeric(0)))
    out <- stats::aggregate(n ~ code, data = tb, FUN = sum)
    out[order(out$code), ]
  }
  n_valid <- sum_counts("valid")  # mask cells with a domain
  n_na <- sum_counts("na")        # ... of which with incomplete covariates

  if (is.null(lookup)) {
    lookup <- data.frame(.code = n_valid$code)
    lookup[[domain_field]] <- n_valid$code
  }
  unit_counts <- data.frame(.code = lookup$.code)
  unit_counts$N_mask <- n_valid$n[match(lookup$.code, n_valid$code)]
  unit_counts$N_mask[is.na(unit_counts$N_mask)] <- 0
  unit_counts$N_na <- n_na$n[match(lookup$.code, n_na$code)]
  unit_counts$N_na[is.na(unit_counts$N_na)] <- 0
  unit_counts$N <- unit_counts$N_mask - unit_counts$N_na

  meta <- list(
    version = 1L,
    type = "cache",
    created = Sys.time(),
    raster = src,
    raster_layers = lyr,
    nrow = nr, ncol = nc,
    crs = terra::crs(r),
    extent = as.vector(terra::ext(r)),
    domain_field = domain_field,
    extra_fields = extra_fields,
    predictors = predictors,
    mask = mask, mask_value = mask_value,
    lookup = lookup,
    unit_counts = unit_counts,
    n_mask_nodomain = sum(vapply(res, `[[`, numeric(1), "nodom")),
    parts = parts,
    N = sum(parts$n)
  )
  saveRDS(meta, file.path(cache_dir, "meta.rds"))

  if (verbose) {
    message(stamp(), "Population cache ready: ", format(meta$N, big.mark = ","), " units in ",
            sum(unit_counts$N > 0), " domains; ", format(sum(unit_counts$N_na), big.mark = ","),
            " mask cells with missing covariates; ",
            format(meta$n_mask_nodomain, big.mark = ","), " mask cells outside all domains. (",
            format(round(difftime(Sys.time(), t_start, units = "mins"), 1)), ")")
  }
  open_pop_cache(cache_dir)
}


# Worker side of the raster extraction (run in PSOCK workers) ----------------------------

pop_worker_init <- function(job) {
  suppressPackageStartupMessages(requireNamespace("terra"))
  suppressPackageStartupMessages(requireNamespace("fst"))
  # GDAL keeps decompressed blocks in a cache of 5% of RAM per process by default:
  # with many workers this alone can exhaust memory.
  terra::setGDALconfig("GDAL_CACHEMAX", as.character(job$gdal_cache_mb))
  terra::terraOptions(memfrac = 0.1, progress = 0)
  r <- terra::rast(job$src)
  names(r) <- job$lyr
  r <- r[[job$read_lyr]]
  terra::readStart(r)
  dr <- NULL
  if (!is.null(job$dom_file)) {
    dr <- terra::rast(job$dom_file)
    terra::readStart(dr)
  }
  assign(".SAEpop_job", c(job, list(r = r, dr = dr)), envir = globalenv())
  invisible(TRUE)
}

pop_worker_task <- function(row0) {
  job <- get(".SAEpop_job", envir = globalenv())
  row_end <- min(job$nr, row0 + job$rows_per_task - 1L)
  out <- list()
  valid_tab <- list()
  na_tab <- list()
  nodom <- 0

  for (rr in seq(row0, row_end, by = job$rows_per_read)) {
    nrr <- min(job$rows_per_read, row_end - rr + 1L)
    v <- terra::readValues(job$r, row = rr, nrows = nrr, col = 1, ncols = job$nc, mat = TRUE)
    colnames(v) <- job$read_lyr
    if (!is.null(job$dr)) {
      dom <- terra::readValues(job$dr, row = rr, nrows = nrr, col = 1, ncols = job$nc, mat = FALSE)
    } else {
      dom <- v[, job$dom_layer]
    }
    inmask <- rep(TRUE, nrow(v))
    if (!is.null(job$mask)) {
      m <- v[, job$mask]
      inmask <- !is.na(m) & m == job$mask_value
    }
    nodom <- nodom + sum(inmask & is.na(dom))
    valid <- inmask & !is.na(dom)
    if (!any(valid)) next
    idx <- which(valid)
    vv <- v[idx, job$predictors, drop = FALSE]
    dd <- dom[idx]
    if (any(dd != round(dd))) stop("Domain codes must be integers.")
    dd <- as.integer(dd)
    cc <- stats::complete.cases(vv)
    valid_tab[[length(valid_tab) + 1]] <- as.data.frame(table(code = dd), stringsAsFactors = FALSE)
    if (any(!cc)) {
      na_tab[[length(na_tab) + 1]] <- as.data.frame(table(code = dd[!cc]), stringsAsFactors = FALSE)
    }
    if (any(cc)) {
      cells <- (rr - 1) * job$nc + idx  # terra cell numbers
      df <- data.frame(.cell = cells[cc], .dom = dd[cc], check.names = FALSE)
      df <- cbind(df, as.data.frame(vv[cc, , drop = FALSE], optional = TRUE))
      out[[length(out) + 1]] <- df
    }
    rm(v, vv)
  }

  tab_fix <- function(l) {
    if (length(l) == 0) return(data.frame(code = integer(0), n = numeric(0)))
    tb <- do.call(rbind, l)
    data.frame(code = as.integer(tb$code), n = as.numeric(tb$Freq))
  }

  file <- NA_character_
  n <- 0
  if (length(out) > 0) {
    df <- do.call(rbind, out)
    n <- nrow(df)
    file <- file.path(job$out_dir, sprintf("part_%06d.fst", row0))
    fst::write_fst(df, file, compress = job$compress)
  }
  list(row0 = as.integer(row0), file = file, n = as.numeric(n),
       valid = tab_fix(valid_tab), na = tab_fix(na_tab), nodom = as.numeric(nodom))
}


# Open an existing cache -------------------------------------------------------------------

open_pop_cache <- function(cache_dir) {
  meta_file <- file.path(cache_dir, "meta.rds")
  if (!file.exists(meta_file)) stop("No population cache (meta.rds) found in ", cache_dir)
  meta <- readRDS(meta_file)
  missing_parts <- !file.exists(file.path(cache_dir, "parts", meta$parts$file))
  if (any(missing_parts)) stop(sum(missing_parts), " cache parts are missing in ", cache_dir)
  structure(c(meta, list(cache_dir = normalizePath(cache_dir, winslash = "/"))), class = "SAEpop")
}


#' @export
print.SAEpop <- function(x, ...) {
  cat("SAEpop population data\n")
  if (identical(x$type, "memory")) {
    cat("  in memory:", format(nrow(x$data), big.mark = ","), "units,",
        length(x$predictors), "covariates, domain field:", x$domain_field, "\n")
  } else {
    cat("  cache:", x$cache_dir, "\n")
    cat("  raster:", paste(x$raster, collapse = ", "), "\n")
    cat("  units:", format(x$N, big.mark = ","), "in", sum(x$unit_counts$N > 0), "domains (",
        x$domain_field, ")\n")
    cat("  covariates:", length(x$predictors), "\n")
    if (!is.null(x$mask)) cat("  mask:", x$mask, "==", x$mask_value, "\n")
    cat("  mask cells with missing covariates:", format(sum(x$unit_counts$N_na), big.mark = ","), "\n")
    cat("  mask cells outside all domains:", format(x$n_mask_nodomain, big.mark = ","), "\n")
  }
  invisible(x)
}


# Internal interface used by the estimators ------------------------------------------------

# Domain (and optional aggregation) of every population unit, as integer codes.
pop_domain_codes <- function(pop, aggregate_to = NULL) {
  if (identical(pop$type, "memory")) {
    f <- factor(dom_key(pop$data[[pop$domain_field]]))
    out <- list(code = as.integer(f), levels = levels(f))
    if (!is.null(aggregate_to)) {
      fa <- factor(dom_key(pop$data[[aggregate_to]]))
      out$agg_code <- as.integer(fa)
      out$agg_levels <- levels(fa)
    }
    return(out)
  }

  files <- file.path(pop$cache_dir, "parts", pop$parts$file)
  raw <- unlist(lapply(files, function(f) fst::read_fst(f, columns = ".dom")[[1]]),
                use.names = FALSE)
  lk <- pop$lookup
  used <- sort(unique(raw))
  lk <- lk[match(used, lk$.code), , drop = FALSE]
  levels <- dom_key(lk[[pop$domain_field]])
  out <- list(code = match(raw, used), levels = levels)
  if (!is.null(aggregate_to)) {
    if (!aggregate_to %in% names(lk)) {
      stop("aggregate_to must be one of the extra_fields stored in the population cache.")
    }
    fa <- factor(dom_key(lk[[aggregate_to]]))
    out$agg_code <- as.integer(fa)[out$code]
    out$agg_levels <- levels(fa)
  }
  out
}

pop_size <- function(pop) {
  if (identical(pop$type, "memory")) nrow(pop$data) else pop$N
}

pop_predictors <- function(pop) pop$predictors

# Chunks of at most chunk_rows consecutive units (global row ranges; a chunk may span
# several cache parts, so that each prediction call gets enough rows to be efficient).
pop_chunks <- function(pop, chunk_rows) {
  chunk_rows <- max(1, floor(chunk_rows))
  n <- pop_size(pop)
  from <- seq(1, n, by = chunk_rows)
  data.frame(gfrom = from, gto = pmin(n, from + chunk_rows - 1))
}

# Covariate matrix of one chunk, columns in the order of var_names.
pop_read_chunk <- function(pop, chunk, var_names) {
  if (identical(pop$type, "memory")) {
    df <- pop$data[chunk$gfrom:chunk$gto, var_names, drop = FALSE]
    return(as_covariate_matrix(df, var_names))
  }
  parts <- pop$parts
  first <- parts$offset + 1
  last <- parts$offset + parts$n
  use <- which(last >= chunk$gfrom & first <= chunk$gto)
  x <- matrix(0, nrow = chunk$gto - chunk$gfrom + 1, ncol = length(var_names),
              dimnames = list(NULL, var_names))
  for (i in use) {
    lf <- max(chunk$gfrom, first[i]) - parts$offset[i]
    lt <- min(chunk$gto, last[i]) - parts$offset[i]
    df <- fst::read_fst(file.path(pop$cache_dir, "parts", parts$file[i]),
                        columns = var_names, from = lf, to = lt)
    rows <- (parts$offset[i] + lf - chunk$gfrom + 1):(parts$offset[i] + lt - chunk$gfrom + 1)
    for (j in seq_along(var_names)) x[rows, j] <- df[[j]]
    rm(df)
  }
  if (anyNA(x)) stop("Missing values in the covariates used for prediction.")
  x
}


stamp <- function() format(Sys.time(), "[%Y-%m-%d %H:%M:%S] ")

# Character keys for domain identifiers (integer-valued numbers without scientific notation)
dom_key <- function(x) {
  if (is.factor(x)) x <- as.character(x)
  if (is.numeric(x) && all(x == round(x), na.rm = TRUE) &&
      all(abs(x) < .Machine$integer.max, na.rm = TRUE)) {
    return(as.character(as.integer(x)))
  }
  as.character(x)
}

# Functions shipped to PSOCK workers must not reference the package namespace, so that
# workers also run when the package is loaded with devtools::load_all().
detach_env <- function(f) {
  environment(f) <- globalenv()
  f
}
