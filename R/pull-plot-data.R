# S7 methods defined here are wired up at load time by S7::methods_register().
# The file sorts after every class it dispatches on.

#' Pull the data behind a plot
#'
#' `pull_plot_data()` returns the tibble one view of a diagnostic result draws.
#' A reader can rebuild the figure with [ggplot2::ggplot()] directly, take the
#' numbers somewhere else, or read what a view is showing without rendering it.
#'
#' @details
#' A view returns exactly one tibble, whatever the figure does with it. Where a
#' view splits its rows across layers, the whole frame comes back and the split
#' stays derivable from it.
#'
#' `type` names the view. Each method mirrors its class's [autoplot()] menu,
#' including the view that menu defaults to, so `pull_plot_data(x)` returns the
#' frame `autoplot(x)` draws. Every view is built from the frame this returns,
#' so the figure and the frame cannot disagree, and a view a result cannot draw
#' is refused here exactly as `autoplot()` refuses it.
#'
#' For a [check_edp()] result:
#'
#' * `"boxplot"` returns the box statistics themselves, one row per
#'   intervention, or one row per measure and intervention for the estimator
#'   variant. `ymin` and `ymax` are the fifth and the ninety-fifth percentile,
#'   `lower`, `middle`, and `upper` are the quartiles, `n` is the number of
#'   observations summarized, and `outliers` is a list-column of the values
#'   lying outside the whiskers, which the view draws as points. The estimator variant covers
#'   `edp_outcome` and `edp_treatment`; `ideal_weight` shares no scale with
#'   them and is left to the scatter view.
#' * `"histogram"`, `"ecdf"`, and `"density"` return the results at their own
#'   grain, one row per observation and intervention, with `intervention` as a
#'   factor whose levels follow the order the interventions were given in. The
#'   label is what separates the interventions rather than `value`, which a
#'   function intervention varies from observation to observation.
#' * `"scatter"` returns the whole estimator results tibble. The view draws the
#'   rows of finite `ideal_weight` in one layer and the infinite rows in
#'   another, and both are read off `ideal_weight`. It needs the estimator
#'   variant, and asking for it after a data-variant run is an error.
#'
#' @param x A [positivity_diagnostic].
#' @param ... Passed to methods. Every method takes `type`, the name of the view
#'   whose data to return, matched against that class's [autoplot()] menu.
#'
#' @return A [tibble][tibble::tibble] holding one view's data.
#'
#' @examples
#' set.seed(1)
#' n <- 100
#' x1 <- rnorm(n)
#' dose <- rnorm(n, mean = x1)
#' df <- data.frame(dose = dose, x1 = x1)
#' result <- check_edp(df, dose, x1, values = c(0, 1), exposure_type = "continuous")
#'
#' pull_plot_data(result)
#' pull_plot_data(result, type = "ecdf")
#' @export
pull_plot_data <- new_generic("pull_plot_data", "x", function(x, ...) {
  S7::S7_dispatch()
})

# ---- Effective data points -------------------------------------------------

#' The results tibble keyed on the intervention label
#'
#' Every view that separates the interventions keys on the label rather than on
#' `value`, which a function intervention varies from observation to
#' observation. The levels come from the normalized `values` rather than from
#' the rows, so a facet strip and an axis take the order the run was given
#' instead of an order the results happen to carry.
#'
#' Two interventions can carry the same label, either from repeated values or
#' from values that round to one label, so the levels are deduplicated. Such
#' interventions share a level and a view draws their rows together; every row
#' of the results is kept either way.
#'
#' @param x An `edp_result`.
#'
#' @return The results with `intervention` as a factor whose levels are the
#'   distinct labels in grid order.
#' @keywords internal
#' @noRd
edp_observation_frame <- function(x) {
  results <- x@results
  results$intervention <- factor(
    results$intervention,
    levels = unique(names(x@params$values))
  )
  results
}

#' Box statistics at the fifth and ninety-fifth percentiles
#'
#' Ring and Schomaker (2026) end the whiskers at the fifth and the ninety-fifth
#' percentile rather than at Tukey's fences, so the summaries are computed here
#' and the view draws them as they stand. Type 7 is the rule
#' `stats::quantile()` and [ggplot2::geom_boxplot()] both take by default, so
#' the quartiles are the ones an ordinary boxplot of the same values would
#' show.
#'
#' @param values The measure to summarize.
#' @param intervention A factor grouping `values`, whose levels fix the rows and
#'   their order.
#'
#' @return A tibble of one row per level of `intervention`, holding `n`, the
#'   five box statistics, and the outliers as a list-column.
#' @keywords internal
#' @noRd
edp_box_stats <- function(values, intervention) {
  groups <- unname(split(values, intervention))
  quantiles <- lapply(groups, function(group) {
    stats::quantile(
      group,
      c(0.05, 0.25, 0.5, 0.75, 0.95),
      type = 7,
      names = FALSE
    )
  })
  statistic <- function(position) {
    vapply(quantiles, `[[`, numeric(1), position)
  }
  tibble::tibble(
    n = lengths(groups),
    ymin = statistic(1L),
    lower = statistic(2L),
    middle = statistic(3L),
    upper = statistic(4L),
    ymax = statistic(5L),
    outliers = Map(
      function(group, quantile) {
        group[group < quantile[[1L]] | group > quantile[[5L]]]
      },
      groups,
      quantiles
    )
  )
}

#' The box statistics one row per box
#'
#' @param x An `edp_result`.
#'
#' @return A tibble keyed on the intervention, and on the measure as well for
#'   the estimator variant.
#' @keywords internal
#' @noRd
edp_boxplot_frame <- function(x) {
  results <- edp_observation_frame(x)
  interventions <- results$intervention
  keys <- factor(levels(interventions), levels = levels(interventions))
  if (x@variant != "estimator") {
    return(vctrs::vec_cbind(
      tibble::tibble(intervention = keys),
      edp_box_stats(results$edp, interventions)
    ))
  }
  # `ideal_weight` is a ratio rather than a count of effective data points, so
  # it shares no scale with the two EDP measures and no axis with them either.
  measures <- c("edp_outcome", "edp_treatment")
  boxes <- lapply(measures, function(measure) {
    vctrs::vec_cbind(
      tibble::tibble(
        intervention = keys,
        measure = factor(measure, levels = measures)
      ),
      edp_box_stats(results[[measure]], interventions)
    )
  })
  vctrs::vec_rbind(!!!boxes)
}

#' Require the estimator variant for the scatter view
#'
#' The scatter view reads `edp_outcome`, `edp_treatment`, and `ideal_weight`,
#' which only the estimator variant computes. Both the plot and the accessor
#' gate here, so a reader is told the same thing whichever they asked for and is
#' shown the call they made.
#'
#' @param x An `edp_result`.
#' @param call The call to blame, or the calling environment.
#'
#' @return `NULL`, invisibly, when the run is the estimator variant.
#' @keywords internal
#' @noRd
edp_require_estimator <- function(x, call = rlang::caller_env()) {
  if (x@variant != "estimator") {
    abort(
      c(
        "The scatter view needs the estimator variant.",
        i = "Rerun {.fn check_edp} with {.code variant = \"estimator\"}."
      ),
      error_class = "positively_variant_error",
      call = call
    )
  }
  invisible(NULL)
}

method(pull_plot_data, edp_result) <- function(
  x,
  type = c("boxplot", "histogram", "ecdf", "density", "scatter"),
  ...
) {
  # The reader's own call rather than the dispatched method, which names a
  # function nobody typed.
  type <- resolve_arg_match(
    rlang::arg_match(type),
    call = rlang::caller_call()
  )
  if (type == "boxplot") {
    edp_boxplot_frame(x)
  } else if (type == "scatter") {
    edp_require_estimator(x, call = rlang::caller_call())
    x@results
  } else {
    edp_observation_frame(x)
  }
}
