# Native Gower distance for mixed-type covariates and the geometric variability
# radius. A gower or cluster dependency is not warranted for these few lines.
# The geometric variability uses a closed form that never materialises the full
# distance matrix, and the per-unit statistics are computed in row chunks so that
# large samples stay within a bounded memory footprint.

#' The row-block size above which Gower distances are chunked
#'
#' Read from `options(positively.gower_chunk_threshold)`, defaulting to `5000`.
#' Samples larger than this are processed in row blocks of this size, keeping the
#' transient allocation to `chunk` by `n` rather than `n` by `n`. Lowering the
#' option in tests exercises the chunked path on a small sample.
#'
#' @return A single integer.
#' @keywords internal
#' @noRd
gower_chunk_threshold <- function() {
  getOption("positively.gower_chunk_threshold", default = 5000L)
}

#' The total Gower distance summed over every ordered pair of rows
#'
#' Returns \eqn{\sum_i \sum_j d(x_i, x_j)}, the sum of the full Gower distance
#' matrix, without ever forming that matrix. For a numeric covariate the sum of
#' absolute differences over ordered pairs is computed from the sorted values in
#' \eqn{O(n \log n)}; for a categorical covariate it is \eqn{n^2 - \sum_k m_k^2},
#' where \eqn{m_k} is the size of category \eqn{k}. A numeric covariate with zero
#' range (a constant column) contributes nothing.
#'
#' @param covariates A data frame of covariate columns, numeric or categorical.
#'
#' @return A single non-negative number, the total over all ordered pairs.
#' @keywords internal
#' @noRd
gower_total_distance <- function(covariates) {
  n <- nrow(covariates)
  p <- ncol(covariates)
  total <- 0
  for (j in seq_len(p)) {
    column <- covariates[[j]]
    if (is.numeric(column)) {
      column_range <- max(column) - min(column)
      if (column_range == 0) {
        next
      }
      sorted <- sort(column)
      # Sum over unordered pairs of |x_i - x_j| for sorted values; the ordered
      # sum, counting both directions, is twice this.
      unordered <- sum((2 * seq_len(n) - n - 1) * sorted)
      total <- total + 2 * unordered / column_range
    } else {
      counts <- table(column)
      total <- total + (n^2 - sum(counts^2))
    }
  }
  total / p
}

#' Geometric variability from the total Gower distance
#'
#' The geometric variability is one half the mean of the full Gower distance
#' matrix, the Cuadras-Arenas quantity used by the WhatIf diagnostics of King and
#' Zeng (2006). It sets the reference radius: an observation counts as supported
#' when an opposite-group unit lies within `nearby` geometric variabilities.
#'
#' @param total The total Gower distance over all ordered pairs.
#' @param n The number of observations.
#'
#' @return A single non-negative number.
#' @keywords internal
#' @noRd
geometric_variability <- function(total, n) {
  total / (2 * n^2)
}

#' The Gower distances from a block of rows to every row
#'
#' Builds the `length(rows)` by `n` block of the Gower distance matrix for the
#' requested rows, summing the per-covariate contributions. Constant numeric
#' covariates contribute zero.
#'
#' @param covariates A data frame of covariate columns.
#' @param rows An integer vector of row indices, the block to compute.
#'
#' @return A numeric matrix with `length(rows)` rows and `n` columns.
#' @keywords internal
#' @noRd
gower_block <- function(covariates, rows) {
  n <- nrow(covariates)
  p <- ncol(covariates)
  block <- matrix(0, nrow = length(rows), ncol = n)
  for (j in seq_len(p)) {
    column <- covariates[[j]]
    if (is.numeric(column)) {
      column_range <- max(column) - min(column)
      if (column_range == 0) {
        next
      }
      block <- block + abs(outer(column[rows], column, "-")) / column_range
    } else {
      column <- as.character(column)
      block <- block + outer(column[rows], column, "!=")
    }
  }
  block / p
}

#' Per-unit Gower support statistics against the opposite exposure group
#'
#' For each observation this computes the Gower distance to the nearest
#' opposite-group unit (`gower_min`), the mean Gower distance to the opposite
#' group (`gower_mean`), and the fraction of the opposite group lying within the
#' `radius` (`frac_nearby`). Rows are processed in blocks of `chunk_size` so that
#' at most a `chunk_size` by `n` distance block is held at once.
#'
#' @param covariates A data frame of covariate columns.
#' @param group An integer or factor vector of length `n` giving the exposure
#'   group of each observation.
#' @param radius The nearby radius, `nearby` times the geometric variability.
#' @param chunk_size The number of rows per block.
#'
#' @return A list with numeric vectors `gower_min`, `gower_mean`, and
#'   `frac_nearby`, each of length `n`.
#' @keywords internal
#' @noRd
gower_support_stats <- function(covariates, group, radius, chunk_size) {
  n <- nrow(covariates)
  gower_min <- numeric(n)
  gower_mean <- numeric(n)
  frac_nearby <- numeric(n)
  for (start in seq(1L, n, by = chunk_size)) {
    rows <- start:min(start + chunk_size - 1L, n)
    block <- gower_block(covariates, rows)
    for (r in seq_along(rows)) {
      i <- rows[r]
      distances <- block[r, group != group[i]]
      gower_min[i] <- min(distances)
      gower_mean[i] <- mean(distances)
      frac_nearby[i] <- mean(distances <= radius)
    }
  }
  list(
    gower_min = gower_min,
    gower_mean = gower_mean,
    frac_nearby = frac_nearby
  )
}
