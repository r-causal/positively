# plot() calls autoplot() and prints, so both entry points draw the same view.

#' Plot an HDR non-overlap diagnostic
#'
#' Draws the non-overlap ratio \eqn{\hat{\tau}(a)} against the target exposure
#' value `a`, with a reference line at zero. For a sequential result the curves
#' are colored by time point. [pull_plot_data()] returns the tibble this view
#' draws.
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
#' autoplot(result)
#'
#' @name autoplot.hdr_result
#' @aliases plot.hdr_result
NULL

method(autoplot, hdr_result) <- function(object, ...) {
  plot_data <- pull_plot_data(object)
  if ("time" %in% names(plot_data)) {
    autoplot_hdr_sequential(plot_data)
  } else {
    autoplot_hdr_point(plot_data)
  }
}

#' The point view of an HDR non-overlap diagnostic
#'
#' @param plot_data The frame [pull_plot_data()] returns for a point
#'   `hdr_result`.
#'
#' @return A [ggplot2::ggplot] object.
#' @keywords internal
#' @noRd
autoplot_hdr_point <- function(plot_data) {
  ggplot2::ggplot(
    plot_data,
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
#' @param plot_data The frame [pull_plot_data()] returns for a sequential
#'   `hdr_result`.
#'
#' @return A [ggplot2::ggplot] object.
#' @keywords internal
#' @noRd
autoplot_hdr_sequential <- function(plot_data) {
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
