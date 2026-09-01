# plot() calls autoplot() and prints, so both entry points draw the same views.

#' Plot an effective-data-points diagnostic
#'
#' Draws one view of a [check_edp()] result. `type = "boxplot"`, the default,
#' boxes EDP by intervention. Its whiskers reach the fifth and the ninety-fifth
#' percentile, following Ring and Schomaker (2026) rather than Tukey's fences,
#' and the values outside them are drawn as points. For the estimator variant
#' it draws the outcome-model and the treatment-model measure side by side in a
#' facet per intervention. `type = "histogram"` shows the distribution of EDP
#' faceted by intervention, `type = "ecdf"` shows its empirical cumulative
#' distribution colored by intervention, and `type = "density"` draws one
#' density ridgeline per intervention, which needs the \pkg{ggridges} package.
#' Each of those three draws one measure: `edp` for the data variant and
#' `edp_outcome` for the estimator variant, the measure carrying the exposure
#' dimension and so the one an intervention moves. `type = "scatter"` plots
#' `edp_outcome` against `edp_treatment` colored by `ideal_weight`; it needs
#' the estimator variant and aborts for the data variant.
#' [pull_plot_data()] returns the tibble any of these views draws.
#'
#' @param object An `edp_result` from [check_edp()].
#' @param type One of `"boxplot"`, `"histogram"`, `"ecdf"`, `"density"`, or
#'   `"scatter"`.
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
#' autoplot(result)
#' autoplot(result, type = "histogram")
#' autoplot(result, type = "ecdf")
#'
#' @name autoplot.edp_result
#' @aliases plot.edp_result
NULL

method(autoplot, edp_result) <- function(
  object,
  type = c("boxplot", "histogram", "ecdf", "density", "scatter"),
  ...
) {
  type <- resolve_arg_match(rlang::arg_match(type))
  switch(
    type,
    boxplot = autoplot_edp_boxplot(object),
    histogram = autoplot_edp_histogram(object),
    ecdf = autoplot_edp_ecdf(object),
    density = autoplot_edp_density(object, call = rlang::current_call()),
    scatter = autoplot_edp_scatter(object, call = rlang::current_call())
  )
}

#' The boxplot view of an EDP diagnostic
#'
#' The boxes are drawn from summaries computed ahead of the plot, because the
#' whiskers reach the fifth and the ninety-fifth percentile rather than Tukey's
#' fences, and the values outside them are a point layer of their own. That
#' layer draws the `outliers` list-column unchopped, so it holds no rows at all
#' when every value lies within the whiskers.
#'
#' @param object An `edp_result`.
#'
#' @return A [ggplot2::ggplot] object.
#' @keywords internal
#' @noRd
autoplot_edp_boxplot <- function(object) {
  estimator <- object@variant == "estimator"
  plot_data <- pull_plot_data(object, type = "boxplot")
  axis <- if (estimator) "measure" else "intervention"
  plot <- ggplot2::ggplot(
    plot_data,
    ggplot2::aes(
      x = .data[[axis]],
      ymin = .data$ymin,
      lower = .data$lower,
      middle = .data$middle,
      upper = .data$upper,
      ymax = .data$ymax
    )
  ) +
    ggplot2::geom_boxplot(stat = "identity") +
    ggplot2::geom_point(
      data = edp_boxplot_outliers(plot_data),
      mapping = ggplot2::aes(x = .data[[axis]], y = .data$outlier),
      inherit.aes = FALSE
    )
  if (!estimator) {
    return(
      plot +
        ggplot2::labs(
          x = "Intervention",
          y = "Effective data points",
          title = "Effective data points by intervention"
        )
    )
  }
  plot +
    ggplot2::facet_wrap(
      ggplot2::vars(.data$intervention),
      labeller = ggplot2::label_both
    ) +
    ggplot2::scale_x_discrete(
      labels = c(
        edp_outcome = "Outcome model",
        edp_treatment = "Treatment model"
      )
    ) +
    ggplot2::labs(
      x = "Model",
      y = "Effective data points",
      title = "Effective data points by model and intervention"
    )
}

#' The values outside the whiskers, one row each
#'
#' @param plot_data The box statistics from [pull_plot_data()].
#'
#' @return The keys of `plot_data` repeated once per outlier they carry, with
#'   the outliers themselves in `outlier`.
#' @keywords internal
#' @noRd
edp_boxplot_outliers <- function(plot_data) {
  statistics <- c("n", "ymin", "lower", "middle", "upper", "ymax", "outliers")
  keys <- plot_data[setdiff(names(plot_data), statistics)]
  outliers <- vctrs::vec_rep_each(keys, lengths(plot_data$outliers))
  # unlist() on a list of empty vectors returns NULL rather than an empty one.
  outliers$outlier <- as.double(unlist(plot_data$outliers, use.names = FALSE))
  outliers
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
  plot_data <- pull_plot_data(object, type = "histogram")
  ggplot2::ggplot(plot_data, ggplot2::aes(x = .data[[measure]])) +
    ggplot2::geom_histogram(bins = 30) +
    ggplot2::facet_wrap(
      ggplot2::vars(.data$intervention),
      labeller = ggplot2::label_both
    ) +
    ggplot2::labs(
      x = "Effective data points",
      y = "Observations",
      title = "Effective data points by intervention"
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
  plot_data <- pull_plot_data(object, type = "ecdf")
  ggplot2::ggplot(
    plot_data,
    ggplot2::aes(x = .data[[measure]], color = .data$intervention)
  ) +
    ggplot2::stat_ecdf() +
    ggplot2::labs(
      x = "Effective data points",
      y = "Cumulative proportion",
      color = "Intervention",
      title = "Effective data points by intervention"
    )
}

#' The ridgeline density view of an EDP diagnostic
#'
#' @param object An `edp_result`.
#' @param call The calling environment, used to build the error's call.
#'
#' @return A [ggplot2::ggplot] object.
#' @keywords internal
#' @noRd
autoplot_edp_density <- function(object, call = rlang::caller_env()) {
  if (!rlang::is_installed("ggridges")) {
    abort(
      c(
        "The density view requires the {.pkg ggridges} package.",
        i = "Install {.pkg ggridges}, or draw another view: {.code autoplot(x, type = \"histogram\")}."
      ),
      error_class = "positively_missing_package_error",
      call = call
    )
  }
  measure <- edp_primary_column(object@results)
  plot_data <- pull_plot_data(object, type = "density")
  ggplot2::ggplot(
    plot_data,
    ggplot2::aes(x = .data[[measure]], y = .data$intervention)
  ) +
    ggridges::geom_density_ridges(
      bandwidth = edp_joint_bandwidth(
        plot_data[[measure]],
        plot_data$intervention
      )
    ) +
    ggplot2::labs(
      x = "Effective data points",
      y = "Intervention",
      title = "Effective data points by intervention"
    )
}

#' The bandwidth the ridgelines share
#'
#' `ggridges::geom_density_ridges()` gives every ridgeline in a panel one
#' bandwidth, the mean of the per-group `stats::bw.nrd0()` rules, and reports
#' the value it settled on. Reproducing the rule here and passing the result in
#' keeps that report out of every transcript that draws the view without
#' changing the figure. Each intervention contributes one row per observation
#' and `check_edp()` requires at least two observations, so no group is ever too
#' small for `stats::bw.nrd0()`.
#'
#' @param values The measure the ridgelines smooth.
#' @param intervention A factor grouping `values`, one ridgeline per level.
#'
#' @return The joint bandwidth, a single positive number.
#' @keywords internal
#' @noRd
edp_joint_bandwidth <- function(values, intervention) {
  groups <- split(values, intervention)
  mean(vapply(groups, stats::bw.nrd0, numeric(1)))
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
  edp_require_estimator(object, call = call)
  results <- pull_plot_data(object, type = "scatter")
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
