# plot() calls autoplot() and prints, so both entry points draw the same views.

#' Plot an extrapolation diagnostic
#'
#' Draws one of two views of a [check_extrapolation()] result.
#' `type = "distribution"` shows the distribution of `frac_nearby` as histograms
#' faceted by exposure group, the share of the opposite group that lies nearby
#' each observation. `type = "hull"` shows convex-hull membership as stacked bars
#' per exposure group, and is available only when the hull test ran.
#' [pull_plot_data()] returns the tibble either of these views draws.
#'
#' @param object An `extrapolation_result` from [check_extrapolation()].
#' @param type One of `"distribution"` or `"hull"`.
#' @param ... Not used.
#'
#' @return A [ggplot2::ggplot] object.
#'
#' @examples
#' set.seed(1)
#' n <- 100
#' x1 <- rnorm(n)
#' x2 <- rnorm(n)
#' exposure <- rbinom(n, 1, plogis(0.5 * x1 - 0.5 * x2))
#' df <- data.frame(exposure = exposure, x1 = x1, x2 = x2)
#' result <- check_extrapolation(df, exposure, c(x1, x2), hull = FALSE)
#'
#' autoplot(result, type = "distribution")
#'
#' @name autoplot.extrapolation_result
#' @aliases plot.extrapolation_result
NULL

method(autoplot, extrapolation_result) <- function(
  object,
  type = c("distribution", "hull"),
  ...
) {
  type <- resolve_arg_match(rlang::arg_match(type))
  if (type == "distribution") {
    autoplot_extrapolation_distribution(object)
  } else {
    autoplot_extrapolation_hull(object, call = rlang::current_call())
  }
}

#' The frac_nearby distribution view of an extrapolation diagnostic
#'
#' @param object An `extrapolation_result`.
#'
#' @return A [ggplot2::ggplot] object.
#' @keywords internal
#' @noRd
autoplot_extrapolation_distribution <- function(object) {
  plot_data <- pull_plot_data(object, type = "distribution")
  ggplot2::ggplot(
    plot_data,
    ggplot2::aes(x = .data$frac_nearby)
  ) +
    ggplot2::geom_histogram(
      bins = 30,
      boundary = 0,
      fill = "grey70",
      color = "white"
    ) +
    ggplot2::facet_wrap(
      ggplot2::vars(.data$exposure),
      ncol = 1,
      labeller = ggplot2::label_both
    ) +
    ggplot2::labs(
      x = "Fraction of the opposite group nearby",
      y = "Observations",
      title = "Nearby support by exposure group"
    )
}

#' The convex-hull membership view of an extrapolation diagnostic
#'
#' @param object An `extrapolation_result`.
#' @param call The calling environment, used to build the error's call.
#'
#' @return A [ggplot2::ggplot] object.
#' @keywords internal
#' @noRd
autoplot_extrapolation_hull <- function(object, call = rlang::caller_env()) {
  extrapolation_require_hull(object, call = call)
  plot_data <- pull_plot_data(object, type = "hull")
  ggplot2::ggplot(
    plot_data,
    ggplot2::aes(x = .data$exposure, fill = .data$membership)
  ) +
    ggplot2::geom_bar() +
    ggplot2::labs(
      x = "Exposure group",
      y = "Observations",
      fill = NULL,
      title = "Convex-hull membership against the opposite group"
    )
}

method(plot, extrapolation_result) <- function(x, y, ...) {
  print(autoplot(x, ...))
  invisible(x)
}
