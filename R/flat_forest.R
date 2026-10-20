# Flattened ranger forests and fast multi-threaded prediction ----------------------------
#
# ranger::predict() sorts every covariate column on each call and re-converts the
# forest, which makes predicting millions of population units with many
# covariates very slow. A "flat_forest" stores the regression trees as plain
# vectors once; predict_flat_forest() traverses them in C++ (bit-identical to
# ranger for numeric covariates).

#' Flatten a ranger regression forest
#'
#' Converts the forest of a \code{ranger} regression model into a compact object
#' of class \code{flat_forest} that can be predicted with \code{\link{predict_flat_forest}}.
#' Forests using partitioned or re-ordered unordered factors are kept as
#' \code{ranger} objects and predicted with \code{ranger::predict}.
#'
#' @param object An object of class \code{ranger} (regression).
#' @return An object of class \code{flat_forest}.
#' @export
flatten_forest <- function(object) {
  if (!inherits(object, "ranger")) {
    stop("flatten_forest() needs a ranger object.")
  }
  f <- object$forest
  if (is.null(f)) {
    stop("The ranger object has no forest (write.forest = FALSE).")
  }
  if (!identical(f$treetype, "Regression")) {
    stop("Only regression forests are supported.")
  }

  fast <- all(f$is.ordered) &&
    (is.null(f$covariate.levels) || all(vapply(f$covariate.levels, is.null, logical(1))))

  if (!fast) {
    return(structure(list(ranger = object, var_names = f$independent.variable.names,
                          num_trees = f$num.trees), class = "flat_forest"))
  }

  nn <- lengths(f$split.varIDs)
  structure(list(
    var = as.integer(unlist(f$split.varIDs, use.names = FALSE)),
    val = as.numeric(unlist(f$split.values, use.names = FALSE)),
    left = as.integer(unlist(lapply(f$child.nodeIDs, `[[`, 1), use.names = FALSE)),
    right = as.integer(unlist(lapply(f$child.nodeIDs, `[[`, 2), use.names = FALSE)),
    tree_start = as.integer(c(0, cumsum(nn)[-length(nn)])),
    depth = NULL,
    var_names = f$independent.variable.names,
    num_trees = f$num.trees,
    ranger = NULL
  ), class = "flat_forest")
}


#' Predict with a flattened forest
#'
#' @param object An object of class \code{flat_forest} (or \code{ranger}).
#' @param newdata A numeric matrix or data.frame containing (at least) the covariates
#' used by the forest. Data.frames are converted with \code{data.matrix}, as in ranger.
#' @param num_threads Number of threads.
#' @param block_rows,group_mb Tuning parameters of the C++ traversal (rows per work
#' item and size in MB of the tree groups kept in cache).
#' @return Numeric vector of predictions.
#' @export
predict_flat_forest <- function(object, newdata, num_threads = 1L, block_rows = 128L,
                                group_mb = 3) {
  if (inherits(object, "ranger")) {
    object <- flatten_forest(object)
  }
  if (!is.null(object$ranger)) {
    return(stats::predict(object$ranger, newdata, num.threads = num_threads)$predictions)
  }

  x <- as_covariate_matrix(newdata, object$var_names)
  if (is.null(object$depth)) {
    object$depth <- cpp_tree_depth(object$left, object$right, object$tree_start)
  }
  cpp_predict_forest(x, object$var, object$val, object$left, object$right,
                     object$tree_start, object$depth, as.integer(num_threads),
                     as.integer(block_rows), as.numeric(group_mb))
}


# Numeric matrix with the covariates in forest order (no copy if already so) -----------
as_covariate_matrix <- function(newdata, var_names) {
  if (is.matrix(newdata) && is.double(newdata) && identical(colnames(newdata), var_names)) {
    x <- newdata
  } else {
    miss <- setdiff(var_names, colnames(newdata))
    if (length(miss) > 0) {
      stop("Covariates missing in the prediction data: ", paste(utils::head(miss, 10), collapse = ", "))
    }
    if (is.data.frame(newdata)) {
      x <- data.matrix(newdata[, var_names, drop = FALSE])
    } else {
      x <- newdata[, var_names, drop = FALSE]
    }
    if (!is.double(x)) {
      storage.mode(x) <- "double"
    }
  }
  if (anyNA(x)) {
    stop("Missing values in the covariates used for prediction.")
  }
  x
}


#' @export
print.flat_forest <- function(x, ...) {
  cat("Flattened ranger regression forest:", x$num_trees, "trees,",
      length(x$var_names), "covariates")
  if (!is.null(x$ranger)) cat(" (kept as ranger object)")
  cat("\n")
  invisible(x)
}
