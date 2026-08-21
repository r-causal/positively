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
#'   lying outside the whiskers, which the view draws as points. The estimator
#'   variant covers `edp_outcome` and `edp_treatment`; `ideal_weight` shares no
#'   scale with them and is left to the scatter view.
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
#' For a [check_density_ratios()] result:
#'
#' * `"distribution"` returns one row per ratio and time point, `time` and
#'   `ratio`. That is the grain the point-treatment histogram bins and the
#'   time-varying boxplot summarizes, and the raw ratios live nowhere else in
#'   the result.
#' * `"cumulative"` returns the cumulative-product summaries the series view
#'   draws, one row per time point and statistic, with `series` naming the
#'   summary in the words the key uses. It needs a time-varying input.
#'
#' For a [check_eta_bias()] result:
#'
#' * `"bootstrap"` returns one row per bootstrap draw, keyed on the estimand
#'   `term`, on the truncation `level`, and on the `label` its facet strip
#'   reads. The view bins the `estimate`. Each draw also carries the `truth` its
#'   term was aimed at, which the view draws as a reference line, one per term
#'   rather than one per draw, since a term's truth is the same in every draw.
#' * `"sweep"` returns one row per term and truncation level, holding the
#'   `truncation` the sweep is drawn against, the `bias`, and the `lower` and
#'   `upper` ends of the band two Monte Carlo standard errors either side of
#'   it. It needs a sweep of more than one level.
#'
#' For a [check_hat_values()] result:
#'
#' * `"null"` returns the null replicates in `phi`. The three numbers the
#'   figure marks beside them, `null_quantile`, `phi_hat`, and `conf_level`,
#'   ride as repeated columns.
#' * `"profile"` returns one row per exposure percentile the candidates were
#'   built at, `prob` and the `fraction` of that percentile's candidates
#'   reading as high leverage. The results hold one row per candidate, so the
#'   aggregation is what this view adds.
#'
#' For a [check_port()] or [check_port_seq()] result there is one view, so the
#' method takes `low_support_only` where a menu would take `type`. It returns
#' one row per bar: the reported subgroups, with a leaf that more than one tree
#' reported collapsed to a single row, the display `label`, the `width`
#' proportional to subgroup size, and the pair of thresholds the row's
#' prevalence was judged against in `beta_lower` and `beta_upper`.
#'
#' For a [check_hdr()] or [check_hdr_seq()] result there is one view, and the
#' frame is the results with `time` as a factor where a sequential run carries
#' one, because the wave separates the curves rather than measuring anything.
#'
#' For a [check_extrapolation()] result:
#'
#' * `"distribution"` returns the results with `exposure` as a factor, which is
#'   what panels the histograms.
#' * `"hull"` returns the same rows with `membership`, the reading of `in_hull`
#'   that the stacked bars and the key both use. It needs the hull test to have
#'   run.
#'
#' A [positivity_check] holds one frame per child rather than a frame of its
#' own, so name the diagnostic to pull: `pull_plot_data(check, "port")`. The
#' frame is that child's, unchanged, and a `type` passed alongside reaches it.
#'
#' @param x A [positivity_diagnostic] or a [positivity_check].
#' @param ... Passed to methods. A class drawing more than one view takes
#'   `type`, the name of the view whose data to return, matched against that
#'   class's [autoplot()] menu. A [check_port()] result takes
#'   `low_support_only` as [autoplot()] does, and a [positivity_check] takes the
#'   name of the diagnostic to pull.
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

# ---- Density ratios --------------------------------------------------------

#' The raw ratios, one row per ratio and time point
#'
#' `@ratios` is the only place the ratios themselves live, so this unchops them
#' into the grain both distribution figures read: a point treatment bins one
#' time point's ratios into a histogram, a time-varying treatment summarizes
#' each time point's own into a box. The `time` column stays either way, so one
#' view keeps one shape whichever figure it draws.
#'
#' @param x A `density_ratios_result`.
#'
#' @return A tibble of `time` and `ratio`.
#' @keywords internal
#' @noRd
density_ratios_distribution_frame <- function(x) {
  per_time <- purrr::imap(
    x@ratios,
    function(ratios, time) {
      tibble::tibble(time = as.integer(time), ratio = ratios)
    }
  )
  vctrs::vec_rbind(!!!per_time)
}

#' The cumulative-product series across time points
#'
#' Only the quantiles and the maximum, because a mean or a proportion shares no
#' scale with them.
#'
#' @param x A `density_ratios_result` whose results carry cumulative rows.
#'
#' @return The cumulative quantile and maximum rows with a `series` label.
#' @keywords internal
#' @noRd
density_ratios_cumulative_frame <- function(x) {
  results <- x@results
  cumulative <- results[grepl("^cumulative_", results$statistic), ]
  is_quantile <- grepl("^cumulative_quantile_", cumulative$statistic)
  is_max <- cumulative$statistic == "cumulative_max"
  plot_data <- cumulative[is_quantile | is_max, ]
  plot_data$series <- cumulative_series(plot_data$statistic)
  plot_data
}

#' Require cumulative rows for the cumulative view
#'
#' A point treatment has one time point and so no product across time points to
#' follow. Both the plot and the accessor gate here, so a reader is told the
#' same thing whichever they asked for and is shown the call they made.
#'
#' @param x A `density_ratios_result`.
#' @param call The call to blame, or the calling environment.
#'
#' @return `NULL`, invisibly, when the result carries cumulative rows.
#' @keywords internal
#' @noRd
density_ratios_require_cumulative <- function(x, call = rlang::caller_env()) {
  if (!any(grepl("^cumulative_", x@results$statistic))) {
    abort(
      c(
        "{.code type = \"cumulative\"} requires a time-varying (matrix) input.",
        i = "This result summarizes a point treatment with a single time point."
      ),
      error_class = "positively_type_error",
      call = call
    )
  }
  invisible(NULL)
}

method(pull_plot_data, density_ratios_result) <- function(
  x,
  type = c("distribution", "cumulative"),
  ...
) {
  call <- rlang::caller_call()
  type <- resolve_arg_match(rlang::arg_match(type), call = call)
  if (type == "cumulative") {
    density_ratios_require_cumulative(x, call = call)
    density_ratios_cumulative_frame(x)
  } else {
    density_ratios_distribution_frame(x)
  }
}

# ---- ETA.Bias --------------------------------------------------------------

#' The bootstrap draws, one row each
#'
#' `@boot_estimates` is the only place the draws live, so this unchops them and
#' repeats each reading's keys across its own draws. A draw is dropped per
#' reported reading when its estimate is not finite, so the readings can rest on
#' different numbers of draws and the keys are repeated by length rather than by
#' a shared count.
#'
#' The term, the level, and the level's label are all factors, so the panels
#' follow the sweep and the order the terms were built in. Left as an integer,
#' the level index would sort lexically and put a tenth level ahead of a second
#' one. The `label` is what the facet strips read, and the truth is a property
#' of the term alone, so every row of a term carries the one value the estimator
#' in it was aimed at.
#'
#' @param x An `eta_bias_result`.
#'
#' @return A tibble of `term`, `level`, `label`, `estimate`, and `truth`.
#' @keywords internal
#' @noRd
eta_bias_bootstrap_frame <- function(x) {
  results <- x@results
  n_terms <- length(x@truth)
  n_levels <- eta_bias_n_levels(results, n_terms)
  labels <- eta_bias_levels(x)$labels
  draws <- lengths(x@boot_estimates)
  level <- rep(eta_bias_level_index(n_levels, n_terms), times = draws)
  term_levels <- unique(results$term)
  frame <- tibble::tibble(
    term = factor(rep(results$term, times = draws), levels = term_levels),
    level = factor(level, levels = seq_len(n_levels)),
    label = factor(labels[level], levels = labels),
    estimate = unlist(x@boot_estimates, use.names = FALSE)
  )
  frame$truth <- unname(x@truth[as.character(frame$term)])
  frame
}

#' The truncation sweep with the band drawn around it
#'
#' @param x An `eta_bias_result` read over more than one truncation level.
#'
#' @return A tibble of `term`, `truncation`, `bias`, `lower`, and `upper`.
#' @keywords internal
#' @noRd
eta_bias_sweep_frame <- function(x) {
  results <- x@results
  n_terms <- length(x@truth)
  n_levels <- eta_bias_n_levels(results, n_terms)
  values <- eta_bias_levels(x)$values
  tibble::tibble(
    term = factor(results$term, levels = unique(results$term)),
    truncation = values[eta_bias_level_index(n_levels, n_terms)],
    bias = results$bias,
    lower = results$bias - 2 * results$mc_se,
    upper = results$bias + 2 * results$mc_se
  )
}

#' Require a sweep of more than one level for the sweep view
#'
#' An estimand of several terms fills several rows at one truncation level, so
#' the row count reads as a sweep that never ran and the level count is what
#' decides there is nothing to draw. Both the plot and the accessor gate here,
#' so a reader is told the same thing whichever they asked for and is shown the
#' call they made.
#'
#' @param x An `eta_bias_result`.
#' @param call The call to blame, or the calling environment.
#'
#' @return `NULL`, invisibly, when the result was read over several levels.
#' @keywords internal
#' @noRd
eta_bias_require_sweep <- function(x, call = rlang::caller_env()) {
  if (eta_bias_n_levels(x@results, length(x@truth)) == 1) {
    entry <- if (x@exposure_type == "continuous") {
      "quantile levels"
    } else {
      "lower bounds"
    }
    abort(
      c(
        "The sweep view needs a truncation sweep of more than one level.",
        i = "Rerun {.fn check_eta_bias} with a {.arg truncation_grid} of multiple {entry}."
      ),
      error_class = "positively_sweep_absent_error",
      call = call
    )
  }
  invisible(NULL)
}

method(pull_plot_data, eta_bias_result) <- function(
  x,
  type = c("bootstrap", "sweep"),
  ...
) {
  call <- rlang::caller_call()
  type <- resolve_arg_match(rlang::arg_match(type), call = call)
  if (type == "bootstrap") {
    eta_bias_bootstrap_frame(x)
  } else {
    eta_bias_require_sweep(x, call = call)
    eta_bias_sweep_frame(x)
  }
}

# ---- Hat values ------------------------------------------------------------

#' The null replicates and the readings marked beside them
#'
#' `@null_dist` is the only place the replicates live. The observed reading, the
#' null quantile it is judged against, and the confidence level that quantile
#' was taken at are single numbers the figure marks as lines, and they ride as
#' repeated columns because a view returns one tibble whatever the figure does
#' with it.
#'
#' @param x A `hat_values_result`.
#'
#' @return A tibble of `phi`, `null_quantile`, `phi_hat`, and `conf_level`.
#' @keywords internal
#' @noRd
hat_values_null_frame <- function(x) {
  tibble::tibble(
    phi = x@null_dist,
    null_quantile = x@null_quantile,
    phi_hat = x@phi_hat,
    conf_level = x@params$conf_level
  )
}

#' The share of high-leverage candidates at each exposure percentile
#'
#' The results hold one row per candidate and percentile, so the aggregation is
#' the whole of what the profile view adds. The percentiles are taken from the
#' results themselves rather than from the names the grouping carries, which
#' round-trip through text and back to a number that is close to the percentile
#' asked for rather than equal to it.
#'
#' @param x A `hat_values_result`.
#'
#' @return A tibble of `prob` and `fraction`, one row per percentile.
#' @keywords internal
#' @noRd
hat_values_profile_frame <- function(x) {
  results <- x@results
  # split() groups on the sorted distinct values, so the group means come back
  # in the order the percentiles are built in below.
  groups <- split(results$high_leverage, results$prob)
  tibble::tibble(
    prob = sort(unique(results$prob)),
    fraction = vapply(groups, mean, numeric(1), USE.NAMES = FALSE)
  )
}

method(pull_plot_data, hat_values_result) <- function(
  x,
  type = c("null", "profile"),
  ...
) {
  type <- resolve_arg_match(rlang::arg_match(type), call = rlang::caller_call())
  if (type == "null") {
    hat_values_null_frame(x)
  } else {
    hat_values_profile_frame(x)
  }
}

# ---- PoRT ------------------------------------------------------------------

#' The reported subgroups with the thresholds each was judged against
#'
#' The bars are the reported subgroups, a leaf reported by more than one tree
#' collapsed to one row; the pair of thresholds a row's prevalence was read
#' against comes last, because that pair is what shows why a subgroup reads as
#' low support. A sequential run resolves its thresholds per time point, and a
#' censoring run judges the exposure and censoring families on different risk
#' sets, so a row's thresholds follow its family as well as its wave.
#'
#' @param x A `port_result`.
#' @param low_support_only Keep only the low-support rows when `TRUE`.
#'
#' @return The plotting frame with `beta_lower` and `beta_upper` appended.
#' @keywords internal
#' @noRd
port_subgroup_frame <- function(x, low_support_only = FALSE) {
  frame <- port_plot_data(x, low_support_only = low_support_only)
  lower <- port_row_beta(x, frame)
  frame$beta_lower <- lower
  frame$beta_upper <- 1 - lower
  frame
}

#' The lower threshold each row was judged against
#'
#' @param x A `port_result`.
#' @param frame The rows to resolve a threshold for.
#'
#' @return A double vector, one entry per row of `frame`.
#' @keywords internal
#' @noRd
port_row_beta <- function(x, frame) {
  if (!("time" %in% names(frame))) {
    return(rep(x@beta[[1]], nrow(frame)))
  }
  if (!("type" %in% names(frame))) {
    return(as.double(x@beta[frame$time]))
  }
  as.double(ifelse(
    frame$type == "censoring",
    x@censoring_beta[frame$time],
    x@beta[frame$time]
  ))
}

method(pull_plot_data, port_result) <- function(
  x,
  low_support_only = FALSE,
  ...
) {
  validate_flag(
    low_support_only,
    arg_name = "low_support_only",
    call = rlang::caller_call()
  )
  port_subgroup_frame(x, low_support_only = low_support_only)
}

# ---- HDR -------------------------------------------------------------------

#' The non-overlap results keyed on the wave
#'
#' A sequential run draws one curve per time point, so the wave is a factor
#' here where the results carry it as the integer index. Left numeric it would
#' colour the curves off a continuous scale, which reads a wave as a quantity.
#'
#' @param x An `hdr_result`.
#'
#' @return The results, with `time` as a factor where there is one.
#' @keywords internal
#' @noRd
hdr_curve_frame <- function(x) {
  results <- x@results
  if ("time" %in% names(results)) {
    results$time <- factor(results$time)
  }
  results
}

method(pull_plot_data, hdr_result) <- function(x, ...) {
  hdr_curve_frame(x)
}

# ---- Extrapolation ---------------------------------------------------------

#' The results keyed on the exposure group
#'
#' The exposure panels the histograms rather than measuring anything, so it is a
#' factor here where the results carry the raw level.
#'
#' @param x An `extrapolation_result`.
#'
#' @return The results, with `exposure` as a factor.
#' @keywords internal
#' @noRd
extrapolation_distribution_frame <- function(x) {
  results <- x@results
  results$exposure <- factor(results$exposure)
  results
}

#' The results with hull membership read as the key reads it
#'
#' @param x An `extrapolation_result` whose hull test ran.
#'
#' @return The distribution frame with a `membership` factor appended.
#' @keywords internal
#' @noRd
extrapolation_hull_frame <- function(x) {
  results <- extrapolation_distribution_frame(x)
  results$membership <- factor(
    ifelse(results$in_hull, "Inside hull", "Outside hull"),
    levels = c("Inside hull", "Outside hull")
  )
  results
}

#' Require the hull test for the convex-hull view
#'
#' Both the plot and the accessor gate here, so a reader is told the same thing
#' whichever they asked for and is shown the call they made.
#'
#' @param x An `extrapolation_result`.
#' @param call The call to blame, or the calling environment.
#'
#' @return `NULL`, invisibly, when the hull test ran.
#' @keywords internal
#' @noRd
extrapolation_require_hull <- function(x, call = rlang::caller_env()) {
  if (!x@hull_run) {
    # A rerun with hull = TRUE only helps when numeric covariates were present;
    # without any, the missing precondition is what to report instead.
    advice <- if (isTRUE(x@params$n_numeric > 0)) {
      "Rerun {.fn check_extrapolation} with {.code hull = TRUE}."
    } else {
      "The hull test needs at least one numeric covariate."
    }
    abort(
      c(
        "The convex-hull view needs the hull test to have run.",
        i = advice
      ),
      error_class = "positively_hull_absent_error",
      call = call
    )
  }
  invisible(NULL)
}

method(pull_plot_data, extrapolation_result) <- function(
  x,
  type = c("distribution", "hull"),
  ...
) {
  call <- rlang::caller_call()
  type <- resolve_arg_match(rlang::arg_match(type), call = call)
  if (type == "distribution") {
    extrapolation_distribution_frame(x)
  } else {
    extrapolation_require_hull(x, call = call)
    extrapolation_hull_frame(x)
  }
}

# ---- Positivity check ------------------------------------------------------

method(pull_plot_data, positivity_check) <- function(x, diagnostic, ...) {
  call <- rlang::caller_call()
  # `autoplot()` reads a missing diagnostic as the whole panel and takes `NULL`
  # for it, so the two spellings arrive here alike. There is no panel of frames,
  # so both are refused rather than one of them guessed at.
  if (missing(diagnostic) || is.null(diagnostic)) {
    abort(
      c(
        "{.arg diagnostic} must name a diagnostic in this container.",
        i = "Available diagnostics are {.val {names(x)}}."
      ),
      error_class = "positively_diagnostic_error",
      call = call
    )
  }
  # extract_named_check() is what `$` and `[[` refuse an unknown name with, so
  # every way of naming a child spells the refusal the same.
  pull_plot_data(extract_named_check(x, diagnostic, call = call), ...)
}
