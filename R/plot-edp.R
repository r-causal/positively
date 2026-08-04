# plot() calls autoplot() and prints, so both entry points draw the same views.

#' Plot an effective-data-points diagnostic
#'
#' Draws one view of a [check_edp()] result. `type = "histogram"` shows the
#' distribution of EDP faceted by intervention value, and `type = "ecdf"` shows
#' its empirical cumulative distribution colored by intervention value. For the
#' estimator variant, `type = "scatter"` plots `edp_outcome` against
#' `edp_treatment` colored by `ideal_weight`; it aborts for the data variant.
#'
#' @param object An `edp_result` from [check_edp()].
#' @param type One of `"histogram"`, `"ecdf"`, or `"scatter"`.
#' @param ... Not used.
#'
#' @return A [ggplot2::ggplot] object.
#'
#' @examples
#' set.seed(1)
#' n <- 100
#' x1 <- rnorm(n)
#' dose <- rnorm(n, mean = x1)
#' df <- data.frame(dose = dose, x1 = x1)
#' result <- check_edp(df, dose, x1, values = c(0, 1), exposure_type = "continuous")
#'
#' autoplot(result, type = "histogram")
#' autoplot(result, type = "ecdf")
#'
#' @name autoplot.edp_result
#' @aliases plot.edp_result
NULL

method(autoplot, edp_result) <- function(
  object,
  type = c("histogram", "ecdf", "scatter"),
  ...
) {
  type <- resolve_arg_match(rlang::arg_match(type))
  if (type == "histogram") {
    autoplot_edp_histogram(object)
  } else if (type == "ecdf") {
    autoplot_edp_ecdf(object)
  } else {
    autoplot_edp_scatter(object, call = rlang::current_call())
  }
}

#' The histogram view of an EDP diagnostic
#'
#' @param object An `edp_result`.
#'
#' @return A [ggplot2::ggplot] object.
#' @keywords internal
#' @noRd
autoplot_edp_histogram <- function(object) {
  measure <- edp_primary_column(object@results)
  plot_data <- object@results
  plot_data$value <- factor(plot_data$value)
  ggplot2::ggplot(plot_data, ggplot2::aes(x = .data[[measure]])) +
    ggplot2::geom_histogram(bins = 30, fill = "grey70", color = "white") +
    ggplot2::facet_wrap(
      ggplot2::vars(.data$value),
      labeller = ggplot2::label_both
    ) +
    ggplot2::labs(
      x = "Effective data points",
      y = "Observations",
      title = "Effective data points by intervention value"
    )
}

#' The empirical cumulative distribution view of an EDP diagnostic
#'
#' @param object An `edp_result`.
#'
#' @return A [ggplot2::ggplot] object.
#' @keywords internal
#' @noRd
autoplot_edp_ecdf <- function(object) {
  measure <- edp_primary_column(object@results)
  plot_data <- object@results
  plot_data$value <- factor(plot_data$value)
  ggplot2::ggplot(
    plot_data,
    ggplot2::aes(x = .data[[measure]], color = .data$value)
  ) +
    ggplot2::stat_ecdf() +
    ggplot2::labs(
      x = "Effective data points",
      y = "Cumulative proportion",
      color = "Intervention value",
      title = "Effective data points by intervention value"
    )
}

#' The estimator scatter view of an EDP diagnostic
#'
#' @param object An `edp_result`.
#' @param call The calling environment, used to build the error's call.
#'
#' @return A [ggplot2::ggplot] object.
#' @keywords internal
#' @noRd
autoplot_edp_scatter <- function(object, call = rlang::caller_env()) {
  if (object@variant != "estimator") {
    abort(
      c(
        "The scatter view needs the estimator variant.",
        i = "Rerun {.fn check_edp} with {.code variant = \"estimator\"}."
      ),
      error_class = "positively_variant_error",
      call = call
    )
  }
  results <- object@results
  finite_rows <- results[is.finite(results$ideal_weight), , drop = FALSE]
  infinite_rows <- results[is.infinite(results$ideal_weight), , drop = FALSE]

  plot <- ggplot2::ggplot(
    mapping = ggplot2::aes(x = .data$edp_treatment, y = .data$edp_outcome)
  )
  if (nrow(finite_rows) > 0) {
    plot <- plot +
      ggplot2::geom_point(
        data = finite_rows,
        mapping = ggplot2::aes(color = .data$ideal_weight),
        alpha = 0.6
      )
  }
  if (nrow(infinite_rows) > 0) {
    plot <- plot +
      ggplot2::geom_point(
        data = infinite_rows,
        shape = 4,
        color = "#B2182B",
        alpha = 0.6
      )
  }

  subtitle <- if (nrow(infinite_rows) > 0) {
    "Crosses mark infinite ideal weight, where outcome support is zero"
  }
  # The colour label belongs to the finite layer's mapping; without finite rows
  # there is no colour aesthetic, so naming it draws an unknown-label message.
  color_label <- if (nrow(finite_rows) > 0) {
    "Ideal weight"
  }
  plot +
    ggplot2::labs(
      x = "Treatment-model effective data points",
      y = "Outcome-model effective data points",
      color = color_label,
      title = "Outcome against treatment support",
      subtitle = subtitle
    )
}

method(plot, edp_result) <- function(x, y, ...) {
  print(autoplot(x, ...))
  invisible(x)
}
