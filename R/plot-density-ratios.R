# plot() calls autoplot() and prints, so both entry points draw the same views.

#' Plot a density-ratio diagnostic
#'
#' Draws one of two views of a [check_density_ratios()] result.
#' `type = "distribution"` shows the distribution of the raw ratios: a histogram
#' for a point treatment and a per-time boxplot for a time-varying treatment.
#' `type = "cumulative"` shows the cumulative-product quantiles across time
#' points, the sequential positivity signature, and is available for matrix
#' inputs.
#'
#' @param object A `density_ratios_result` from [check_density_ratios()].
#' @param type One of `"distribution"` or `"cumulative"`.
#' @param ... Not used.
#'
#' @return A [ggplot2::ggplot] object.
#'
#' @examples
#' m <- matrix(c(1, 1, 2, 1, 3, 1, 4, 1), nrow = 2)
#' result <- check_density_ratios(m)
#'
#' autoplot(result, type = "distribution")
#' autoplot(result, type = "cumulative")
#'
#' @name autoplot.density_ratios_result
#' @aliases plot.density_ratios_result
NULL

method(autoplot, density_ratios_result) <- function(
  object,
  type = c("distribution", "cumulative"),
  ...
) {
  type <- resolve_arg_match(rlang::arg_match(type))
  if (type == "cumulative") {
    autoplot_density_ratios_cumulative(object, call = rlang::current_call())
  } else if (length(object@ratios) == 1) {
    autoplot_density_ratios_histogram(object)
  } else {
    autoplot_density_ratios_boxplot(object)
  }
}

#' The point-treatment histogram view
#'
#' @param object A `density_ratios_result`.
#'
#' @return A [ggplot2::ggplot] object.
#' @keywords internal
#' @noRd
autoplot_density_ratios_histogram <- function(object) {
  plot_data <- tibble::tibble(ratio = object@ratios[[1]])
  ggplot2::ggplot(plot_data, ggplot2::aes(x = .data$ratio)) +
    ggplot2::geom_histogram(bins = 30, fill = "grey70", color = "white") +
    ggplot2::labs(
      x = "Density ratio",
      y = "Observations",
      title = "Distribution of density ratios"
    )
}

#' The time-varying boxplot view
#'
#' @param object A `density_ratios_result`.
#'
#' @return A [ggplot2::ggplot] object.
#' @keywords internal
#' @noRd
autoplot_density_ratios_boxplot <- function(object) {
  per_time <- purrr::imap(
    object@ratios,
    function(ratios, time) {
      tibble::tibble(time = as.integer(time), ratio = ratios)
    }
  )
  plot_data <- vctrs::vec_rbind(!!!per_time)
  ggplot2::ggplot(
    plot_data,
    ggplot2::aes(x = factor(.data$time), y = .data$ratio)
  ) +
    ggplot2::geom_boxplot() +
    ggplot2::labs(
      x = "Time point",
      y = "Density ratio",
      title = "Density ratios across time points"
    )
}

#' The cumulative-product quantile view
#'
#' @param object A `density_ratios_result`.
#'
#' @return A [ggplot2::ggplot] object.
#' @keywords internal
#' @noRd
autoplot_density_ratios_cumulative <- function(
  object,
  call = rlang::caller_env()
) {
  results <- object@results
  cumulative <- results[grepl("^cumulative_", results$statistic), ]
  if (nrow(cumulative) == 0) {
    abort(
      c(
        "{.code type = \"cumulative\"} requires a time-varying (matrix) input.",
        i = "This result summarizes a point treatment with a single time point."
      ),
      error_class = "positively_type_error",
      call = call
    )
  }
  is_quantile <- grepl("^cumulative_quantile_", cumulative$statistic)
  is_max <- cumulative$statistic == "cumulative_max"
  plot_data <- cumulative[is_quantile | is_max, ]
  plot_data$series <- cumulative_series(plot_data$statistic)
  ggplot2::ggplot(
    plot_data,
    ggplot2::aes(x = .data$time, y = .data$value, color = .data$series)
  ) +
    ggplot2::geom_line() +
    ggplot2::geom_point() +
    ggplot2::labs(
      x = "Time point",
      y = "Cumulative density ratio",
      color = "Summary",
      title = "Cumulative-product quantiles across time points"
    )
}

#' Readable, ordered series labels for the cumulative statistics
#'
#' Maps each cumulative statistic name to a human-readable label, ordered by the
#' underlying quantile probability with the maximum last, so the plot legend
#' reads as text rather than raw statistic keys.
#'
#' @param statistic A character vector of `cumulative_quantile_*` and
#'   `cumulative_max` names.
#'
#' @return An ordered factor of readable labels.
#' @keywords internal
#' @noRd
cumulative_series <- function(statistic) {
  is_max <- statistic == "cumulative_max"
  pct <- rep(100, length(statistic))
  pct[!is_max] <- as.numeric(
    sub("^cumulative_quantile_", "", statistic[!is_max])
  )
  label <- ifelse(
    is_max,
    "Maximum",
    paste0(vapply(pct, percentile_label, character(1)), " percentile")
  )
  factor(label, levels = unique(label[order(pct)]))
}

method(plot, density_ratios_result) <- function(x, y, ...) {
  print(autoplot(x, ...))
  invisible(x)
}
