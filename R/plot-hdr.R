# plot() calls autoplot() and prints, so both entry points draw the same view.

#' Plot an HDR non-overlap diagnostic
#'
#' Draws the non-overlap ratio \eqn{\hat{\tau}(a)} against the target exposure
#' value `a`, with a reference line at zero. For a sequential result the curves
#' are colored by time point.
#'
#' @param object An `hdr_result` from [check_hdr()] or [check_hdr_seq()].
#' @param ... Not used.
#'
#' @return A [ggplot2::ggplot] object.
#'
#' @examples
#' set.seed(1)
#' n <- 300
#' l <- rnorm(n)
#' dose <- rnorm(n, mean = l)
#' df <- data.frame(dose = dose, l = l)
#' result <- check_hdr(df, dose, l)
#'
#' ggplot2::autoplot(result)
#'
#' @name autoplot.hdr_result
#' @aliases plot.hdr_result
NULL

method(autoplot, hdr_result) <- function(object, ...) {
  results <- object@results
  if ("time" %in% names(results)) {
    autoplot_hdr_sequential(results)
  } else {
    autoplot_hdr_point(results)
  }
}

#' The point view of an HDR non-overlap diagnostic
#'
#' @param results The `@results` tibble of a point `hdr_result`.
#'
#' @return A [ggplot2::ggplot] object.
#' @keywords internal
#' @noRd
autoplot_hdr_point <- function(results) {
  ggplot2::ggplot(
    results,
    ggplot2::aes(x = .data$value, y = .data$nonoverlap)
  ) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed") +
    ggplot2::geom_line() +
    ggplot2::geom_point() +
    ggplot2::labs(
      x = "Target exposure value",
      y = "Non-overlap ratio",
      title = "HDR non-overlap across target exposure values"
    )
}

#' The sequential view of an HDR non-overlap diagnostic
#'
#' @param results The `@results` tibble of a sequential `hdr_result`.
#'
#' @return A [ggplot2::ggplot] object.
#' @keywords internal
#' @noRd
autoplot_hdr_sequential <- function(results) {
  plot_data <- results
  plot_data$time <- factor(plot_data$time)
  ggplot2::ggplot(
    plot_data,
    ggplot2::aes(
      x = .data$value,
      y = .data$nonoverlap,
      color = .data$time
    )
  ) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed") +
    ggplot2::geom_line() +
    ggplot2::geom_point() +
    ggplot2::labs(
      x = "Target exposure value",
      y = "Non-overlap ratio",
      color = "Time",
      title = "HDR non-overlap across target exposure values by time"
    )
}

method(plot, hdr_result) <- function(x, y, ...) {
  print(autoplot(x, ...))
  invisible(x)
}
