# Convex-hull membership as a linear-programming feasibility test. lpSolve lives
# in Suggests, so callers gate access with rlang::check_installed().

#' Is a point inside the convex hull of a reference cloud?
#'
#' A query point `x` lies in the convex hull of the reference points
#' \eqn{r_1, \ldots, r_m} exactly when there exist weights
#' \eqn{\lambda_k \ge 0} with \eqn{\sum_k \lambda_k = 1} and
#' \eqn{\sum_k \lambda_k r_k = x}. This is a linear-programming feasibility
#' problem: the objective is irrelevant, so a zero objective is minimized subject
#' to the weight constraints and the solver status reports feasibility.
#'
#' @param query A numeric vector of length `p`, the point to test.
#' @param reference An `m` by `p` numeric matrix of reference points.
#'
#' @return A single logical value: `TRUE` when `query` is in the hull.
#' @keywords internal
#' @noRd
in_convex_hull <- function(query, reference) {
  m <- nrow(reference)
  const_mat <- rbind(t(reference), rep(1, m))
  const_rhs <- c(query, 1)
  solution <- lpSolve::lp(
    direction = "min",
    objective.in = rep(0, m),
    const.mat = const_mat,
    const.dir = rep("=", length(const_rhs)),
    const.rhs = const_rhs
  )
  # Status 0 is a feasible solution; a non-zero status reports infeasibility.
  solution$status == 0
}

#' Convex-hull membership of every unit against the opposite exposure group
#'
#' Tests each observation for membership in the convex hull of the opposite
#' exposure group, using the numeric covariates only. The hull test is defined
#' for continuous coordinates; categorical covariates carry no ordering and are
#' excluded, a choice documented in [check_extrapolation()].
#'
#' @param numeric_covariates An `n` by `p` numeric matrix of covariates.
#' @param group An integer or factor vector of length `n` giving the exposure
#'   group of each observation.
#'
#' @return A logical vector of length `n`.
#' @keywords internal
#' @noRd
hull_membership <- function(numeric_covariates, group) {
  n <- nrow(numeric_covariates)
  # Hull membership is invariant under a per-column affine rescaling, but the LP
  # feasibility tolerance is absolute, so columns are rescaled to [0, 1] to keep
  # large covariate offsets from misreporting membership.
  for (j in seq_len(ncol(numeric_covariates))) {
    column <- numeric_covariates[, j]
    column_min <- min(column)
    column_range <- max(column) - column_min
    numeric_covariates[, j] <- if (column_range == 0) {
      0
    } else {
      (column - column_min) / column_range
    }
  }
  in_hull <- logical(n)
  for (i in seq_len(n)) {
    opposite <- group != group[i]
    reference <- numeric_covariates[opposite, , drop = FALSE]
    in_hull[i] <- in_convex_hull(numeric_covariates[i, ], reference)
  }
  in_hull
}
