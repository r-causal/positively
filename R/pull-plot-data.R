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
#' * `"histogram"` and `"ecdf"` return the results at their own grain, one row
#'   per observation and intervention, with `intervention` as a factor whose
#'   levels follow the order the interventions were given in. The label is what
#'   separates the interventions rather than `value`, which a function
#'   intervention varies from observation to observation.
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
  type = c("histogram", "ecdf", "scatter"),
  ...
) {
  # The reader's own call rather than the dispatched method, which names a
  # function nobody typed.
  type <- resolve_arg_match(
    rlang::arg_match(type),
    call = rlang::caller_call()
  )
  if (type == "scatter") {
    edp_require_estimator(x, call = rlang::caller_call())
    x@results
  } else {
    edp_observation_frame(x)
  }
}
