#' Aggregate domain estimates to coarser levels with bootstrap MSE
#'
#' Aggregates the domain-level means and totals of a scalable \code{SAEforest_mean} model
#' (unit-level data, fitted with \code{MSE = "nonparametric"}) to groups of domains, e.g.
#' municipalities to provinces, regions or the whole country. Totals are sums of domain totals and
#' means are population-weighted means of domain means. MSEs are computed from the stored
#' bootstrap replicates of every domain, so that the correlation between domain estimates
#' induced by the common random forest is accounted for: the result is identical to
#' re-running \code{SAEforest_model} with \code{aggregate_to} (same seed).
#'
#' @param object An object of class \code{SAEforest_mean} with bootstrap replicates
#' (element \code{Bootstrap}), as returned by \code{\link{SAEforest_model}}.
#' @param groups Either a data.frame with two columns (domain identifier, group
#' identifier) or a vector of group identifiers named by domain. Use a single group
#' (e.g. \code{"country"}) for the overall estimate.
#' @param name Name of the group column in the output.
#'
#' @return A data.frame with, per group: number of domains, population units, \code{Mean},
#' \code{Mean_MSE}, \code{Mean_RMSE}, \code{Mean_RMSEp} (RMSE in \% of the estimate) and,
#' if the model was fitted with \code{total_factor}, the same for \code{Total}.
#' @export
SAEaggregate <- function(object, groups, name = "group") {
  if (!inherits(object, "SAEforest_mean")) {
    stop("object must be a SAEforest_mean object.")
  }
  if (is.null(object$PopInfo)) {
    stop("SAEaggregate() needs a model fitted with the scalable engine (engine = 'auto').")
  }
  dom_col <- colnames(object$Indicators)[1]
  doms <- dom_key(object$Indicators[[1]])

  if (is.data.frame(groups)) {
    if (ncol(groups) < 2) stop("groups must have a domain and a group column.")
    gmap <- stats::setNames(as.character(groups[[2]]), dom_key(groups[[1]]))
  } else {
    if (length(groups) == 1 && is.null(names(groups))) {
      gmap <- stats::setNames(rep(as.character(groups), length(doms)), doms)
    } else {
      gmap <- stats::setNames(as.character(groups), dom_key(names(groups)))
    }
  }
  gmap <- gmap[!duplicated(names(gmap))]
  g <- unname(gmap[doms])
  if (anyNA(g)) {
    stop(sum(is.na(g)), " domains have no group, e.g. ", paste(utils::head(doms[is.na(g)], 5), collapse = ", "))
  }

  f <- object$PopInfo$total_factor %||% 1
  N <- as.numeric(object$PopInfo$N_total[doms]) * f  # weight of each domain (e.g. forest area)
  mean_d <- object$Indicators$Mean

  gf <- factor(g)
  W <- as.numeric(tapply(N, gf, sum))
  tot <- as.numeric(tapply(mean_d * N, gf, sum))
  out <- data.frame(levels(gf), n_domains = as.integer(table(gf)),
                    N_units = as.numeric(tapply(N / f, gf, sum)),
                    Mean = tot / W, stringsAsFactors = FALSE)
  names(out)[1] <- name

  bt <- object$Bootstrap
  if (!is.null(bt)) {
    if (!identical(dom_key(bt$domains), doms)) {
      stop("Bootstrap replicates and indicators refer to different domains.")
    }
    # replicate totals per group: rows = groups, columns = bootstrap replicates
    tot_star <- rowsum(bt$tau_star * N, gf, reorder = TRUE)
    tot_b <- rowsum(bt$tau_b * N, gf, reorder = TRUE)
    mse_tot <- rowMeans((tot_star - tot_b)^2)
    out$Mean_MSE <- mse_tot / W^2
    out$Mean_RMSE <- sqrt(out$Mean_MSE)
    out$Mean_RMSEp <- 100 * out$Mean_RMSE / out$Mean
  }
  if (!is.null(object$PopInfo$total_factor)) {
    out$Total <- tot
    if (!is.null(bt)) {
      out$Total_MSE <- mse_tot
      out$Total_RMSE <- sqrt(mse_tot)
      out$Total_RMSEp <- 100 * out$Total_RMSE / out$Total
    }
  }
  rownames(out) <- NULL
  out
}
