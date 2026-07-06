# S7 methods defined here are wired up at load time by S7::methods_register().

#' The extrapolation diagnostic result class
#'
#' `extrapolation_result` is the S7 class returned by [check_extrapolation()]. It
#' extends [positivity_diagnostic] with the geometric variability radius and a
#' flag recording whether the convex-hull test ran, and it exposes a per-group
#' `@summary` tibble through a read-only getter. It is created internally and is
#' not constructed directly by users.
#'
#' @keywords internal
#' @noRd
extrapolation_result <- new_class(
  "extrapolation_result",
  parent = positivity_diagnostic,
  properties = list(
    geometric_variability = class_double,
    hull_run = class_logical,
    summary = new_property(
      getter = function(self) extrapolation_summary(self)
    )
  ),
  validator = function(self) {
    if (length(self@geometric_variability) != 1) {
      "@geometric_variability must be a single value"
    } else if (length(self@hull_run) != 1) {
      "@hull_run must be a single value"
    }
  }
)

#' Diagnose extrapolation for binary exposures with Gower distances and the
#' convex hull
#'
#' `check_extrapolation()` implements the extrapolation diagnostics of King and
#' Zeng (2006) for binary point exposures. It measures, for every observation,
#' how well the opposite exposure group covers its position in covariate space,
#' using Gower distances that handle mixed numeric and categorical covariates,
#' and it optionally tests convex-hull membership against the opposite group. It
#' reports a diagnostic only and assigns no severity grade.
#'
#' @details
#' Estimating a causal contrast requires comparing treated and untreated units
#' with similar covariates. Where one exposure group has no nearby members of the
#' other, the contrast rests on extrapolation rather than data, the concern
#' framed by Petersen and colleagues (2012) and by King and Zeng (2006).
#'
#' The Gower distance between two rows is the mean over covariates of a
#' per-covariate distance: for a numeric covariate, the absolute difference
#' divided by the pooled sample range; for a factor or character covariate, one
#' when the categories differ and zero when they agree. Distances therefore lie
#' in `[0, 1]` and mixed covariate types combine on a common scale. Because the
#' denominator is the sample range, a single extreme covariate value rescales
#' every distance on that axis.
#'
#' The geometric variability is one half the mean of the full Gower distance
#' matrix, the Cuadras-Arenas quantity used by the WhatIf software. It sets the
#' reference radius: an observation is counted as nearby-supported when an
#' opposite-group unit lies within `nearby` geometric variabilities. Three
#' per-unit statistics are reported: `gower_min` (distance to the nearest
#' opposite-group unit), `gower_mean` (mean distance to the opposite group), and
#' `frac_nearby` (the fraction of the opposite group within `nearby` geometric
#' variabilities). A planted support gap raises `gower_min` and lowers
#' `frac_nearby` in both directions. The per-group `@summary` column
#' `prop_supported`, and the `prop_supported` statistic in [glance()] and the
#' print method, share one fixed definition: the fraction of units whose nearest
#' opposite-group unit lies within one geometric variability
#' (`gower_min <= geometric_variability`), independent of `nearby`.
#'
#' When every covariate is constant, all rows coincide, the geometric
#' variability is zero, and each observation is at distance zero from the
#' opposite group. Such an observation is fully supported: `gower_min` is zero
#' and `frac_nearby` is one. Constant numeric covariates likewise contribute
#' zero distance on their own axis without affecting the others.
#'
#' The convex-hull test asks whether an observation lies inside the convex hull
#' of the opposite group, solved as a linear-programming feasibility problem
#' (weights that are non-negative, sum to one, and reproduce the point). The hull
#' test uses the numeric covariates only, since categorical covariates carry no
#' ordering to interpolate over. In-hull membership collapses toward zero as the
#' covariate dimension grows even under perfect overlap, so the test is
#' degenerate in high dimensions (D'Amour and colleagues, 2021). It therefore
#' runs automatically only when there are at most ten numeric covariates and
#' \pkg{lpSolve} is installed; above that it is skipped and `in_hull` is `NA`.
#' Setting `hull = TRUE` forces the test to run and warns when the dimension is
#' high; `hull = FALSE` skips it.
#'
#' @param .data A data frame.
#' @param .exposure The binary exposure column, selected with data-masking.
#'   `check_extrapolation()` aborts for continuous or categorical exposures.
#' @param .covariates The covariate columns, selected with tidyselect. Numeric
#'   and categorical covariates are both supported by the Gower distance.
#' @param nearby The reference radius as a multiple of the geometric variability.
#'   An observation counts an opposite-group unit as nearby when their Gower
#'   distance is at most `nearby` times the geometric variability. A single
#'   positive number, defaulting to `1`.
#' @param hull Whether to run the convex-hull membership test. `NULL` (the
#'   default) runs it automatically when there is at least one and at most ten
#'   numeric covariates and \pkg{lpSolve} is installed, and skips it otherwise.
#'   `TRUE` forces it to run, aborting when \pkg{lpSolve} is not installed and
#'   warning when there are no numeric covariates or the dimension is high.
#'   `FALSE` skips it.
#'
#' @return An `extrapolation_result` object, an S7 subclass of
#'   [positivity_diagnostic]. Its `@results` tibble has one row per observation
#'   with columns `.id` (the row index), `exposure` (the exposure value),
#'   `frac_nearby`, `gower_min`, `gower_mean`, and `in_hull` (the hull membership
#'   flag, or `NA` when the hull test did not run). It also carries the scalar
#'   properties `@geometric_variability` and `@hull_run`, and a per-group
#'   `@summary` tibble with one row per exposure level and columns `exposure`,
#'   `n`, `mean_frac_nearby`, `mean_gower_min`, `prop_supported` (the one
#'   geometric variability fraction described in Details), and `prop_in_hull`
#'   (the fraction inside the opposite group's hull, `NA` when the hull test did
#'   not run).
#'
#' @references
#' King G, Zeng L (2006). The Dangers of Extreme Counterfactuals.
#' *Political Analysis*, 14(2):131-159. \doi{10.1093/pan/mpj004}
#'
#' Petersen ML, Porter KE, Gruber S, Wang Y, van der Laan MJ (2012). Diagnosing
#' and Responding to Violations in the Positivity Assumption.
#' *Statistical Methods in Medical Research*, 21(1):31-54.
#' \doi{10.1177/0962280210386207}
#'
#' D'Amour A, Ding P, Feller A, Lei L, Sekhon J (2021). Overlap in Observational
#' Studies with High-Dimensional Covariates. *Journal of Econometrics*,
#' 221(2):644-654. \doi{10.1016/j.jeconom.2019.10.014}
#'
#' @examples
#' set.seed(1)
#' n <- 100
#' x1 <- rnorm(n)
#' x2 <- rnorm(n)
#' exposure <- rbinom(n, 1, plogis(0.5 * x1 - 0.5 * x2))
#' df <- data.frame(exposure = exposure, x1 = x1, x2 = x2)
#'
#' # hull = FALSE keeps the example free of the optional lpSolve dependency.
#' result <- check_extrapolation(df, exposure, c(x1, x2), hull = FALSE)
#' result
#'
#' @export
check_extrapolation <- function(
  .data,
  .exposure,
  .covariates,
  nearby = 1,
  hull = NULL
) {
  validate_data_frame(.data)
  validate_positive_number(nearby, arg_name = "nearby")
  if (!is.null(hull) && !rlang::is_bool(hull)) {
    abort(
      "{.arg hull} must be `NULL`, `TRUE`, or `FALSE`.",
      error_class = "positively_type_error"
    )
  }

  exposure_pos <- tidyselect::eval_select(rlang::enquo(.exposure), .data)
  if (length(exposure_pos) != 1) {
    abort(
      "{.arg .exposure} must select exactly one column, not {length(exposure_pos)}.",
      error_class = "positively_selection_error"
    )
  }
  exposure_name <- names(exposure_pos)
  exposure_vec <- .data[[exposure_pos]]
  if (anyNA(exposure_vec)) {
    abort(
      "{.arg .exposure} must not contain missing values.",
      error_class = "positively_missing_error"
    )
  }

  n <- nrow(.data)
  if (n < 2) {
    abort(
      "{.arg .data} must have at least two observations, not {n}.",
      error_class = "positively_size_error"
    )
  }

  exposure_type <- detect_exposure_type(exposure_vec)
  if (exposure_type != "binary") {
    abort(
      c(
        "{.fn check_extrapolation} supports binary exposures only.",
        i = "{.arg .exposure} was detected as {.val {exposure_type}}."
      ),
      error_class = "positively_exposure_type_error"
    )
  }

  covariate_pos <- tidyselect::eval_select(rlang::enquo(.covariates), .data)
  validate_column_selection(covariate_pos, ".covariates")
  covariate_names <- names(covariate_pos)
  missing_covariates <- covariate_names[vapply(
    covariate_names,
    function(column) anyNA(.data[[column]]),
    logical(1)
  )]
  if (length(missing_covariates) > 0) {
    abort(
      c(
        "{.arg .covariates} must not contain missing values.",
        x = "Missing values in {.val {missing_covariates}}."
      ),
      error_class = "positively_missing_error"
    )
  }

  covariates <- .data[covariate_pos]
  gv <- geometric_variability(gower_total_distance(covariates), n)
  radius <- nearby * gv

  group <- match(exposure_vec, sort(unique(exposure_vec)))
  chunk_size <- gower_chunk_threshold()
  if (n > chunk_size) {
    alert_info(
      "Computing Gower distances for {n} observations in row chunks of {chunk_size}."
    )
  }
  stats <- gower_support_stats(covariates, group, radius, chunk_size)

  is_numeric_column <- vapply(covariates, is.numeric, logical(1))
  n_numeric <- sum(is_numeric_column)
  hull_run <- resolve_hull(hull, n_numeric)
  if (hull_run) {
    numeric_covariates <- as.matrix(covariates[is_numeric_column])
    in_hull <- hull_membership(numeric_covariates, group)
  } else {
    in_hull <- rep(NA, n)
  }

  results <- tibble::tibble(
    .id = seq_len(n),
    exposure = exposure_vec,
    frac_nearby = stats$frac_nearby,
    gower_min = stats$gower_min,
    gower_mean = stats$gower_mean,
    in_hull = in_hull
  )

  extrapolation_result(
    results = results,
    exposure = exposure_name,
    exposure_type = exposure_type,
    n = as.integer(n),
    params = list(nearby = nearby, hull = hull),
    call = rlang::current_call(),
    geometric_variability = gv,
    hull_run = hull_run
  )
}

# The largest numeric-covariate dimension at which the convex-hull test runs
# automatically; above it the finite-sample hull is degenerate.
hull_max_p <- 10L

#' Decide whether the convex-hull test should run
#'
#' Resolves the `hull` argument against the numeric-covariate dimension and the
#' availability of \pkg{lpSolve}. `NULL` runs the test automatically for at most
#' `hull_max_p` numeric covariates when \pkg{lpSolve} is installed; `TRUE` forces
#' it, aborting without \pkg{lpSolve} and warning at high dimension; `FALSE`
#' skips it.
#'
#' @param hull The user-supplied `hull` argument: `NULL`, `TRUE`, or `FALSE`.
#' @param n_numeric The number of numeric covariates.
#'
#' @return A single logical value: whether to run the hull test.
#' @keywords internal
#' @noRd
resolve_hull <- function(hull, n_numeric, call = rlang::caller_env()) {
  if (isFALSE(hull)) {
    return(FALSE)
  }
  if (n_numeric < 1) {
    if (isTRUE(hull)) {
      warn(
        "The convex-hull test needs at least one numeric covariate; skipping it.",
        warning_class = "positively_hull_covariate_warning",
        call = call
      )
    }
    return(FALSE)
  }
  if (isTRUE(hull)) {
    if (!rlang::is_installed("lpSolve")) {
      abort(
        c(
          "The convex-hull test requires the {.pkg lpSolve} package.",
          i = "Install {.pkg lpSolve} or set {.code hull = FALSE}."
        ),
        error_class = "positively_missing_package_error",
        call = call
      )
    }
    if (n_numeric > hull_max_p) {
      warn(
        c(
          "The convex-hull test is degenerate above {hull_max_p} numeric covariates.",
          i = "{n_numeric} numeric covariates were supplied; membership will be near zero."
        ),
        warning_class = "positively_hull_dimension_warning",
        call = call
      )
    }
    return(TRUE)
  }
  # hull is NULL: run automatically only when the test is defensible.
  if (n_numeric > hull_max_p) {
    alert_info(
      "Skipping the convex-hull test: {n_numeric} numeric covariates exceed the maximum of {hull_max_p}."
    )
    return(FALSE)
  }
  if (!rlang::is_installed("lpSolve")) {
    alert_info(
      "Skipping the convex-hull test: the {.pkg lpSolve} package is not installed."
    )
    return(FALSE)
  }
  TRUE
}

#' The per-group summary of an extrapolation diagnostic
#'
#' Aggregates the per-unit results to one row per exposure level, reporting the
#' group size, the mean support statistics, the fraction of units with an
#' opposite-group neighbour within one geometric variability, and the fraction
#' inside the opposite group's convex hull.
#'
#' @param self An `extrapolation_result`.
#'
#' @return A tibble with one row per exposure level.
#' @keywords internal
#' @noRd
extrapolation_summary <- function(self) {
  results <- self@results
  gv <- self@geometric_variability
  levels_exposure <- sort(unique(results$exposure))
  rows <- lapply(levels_exposure, function(level) {
    in_group <- results$exposure == level
    tibble::tibble(
      exposure = level,
      n = sum(in_group),
      mean_frac_nearby = mean(results$frac_nearby[in_group]),
      mean_gower_min = mean(results$gower_min[in_group]),
      prop_supported = mean(results$gower_min[in_group] <= gv),
      prop_in_hull = mean(results$in_hull[in_group])
    )
  })
  vctrs::vec_rbind(!!!rows)
}

# ---- Methods --------------------------------------------------------------

method(glance, extrapolation_result) <- function(x, ...) {
  tibble::tibble(
    exposure = paste(x@exposure, collapse = ", "),
    exposure_type = x@exposure_type,
    n = x@n,
    nearby = x@params$nearby,
    geometric_variability = x@geometric_variability,
    prop_supported = mean(x@results$gower_min <= x@geometric_variability),
    hull_run = x@hull_run
  )
}

method(print, extrapolation_result) <- function(x, ...) {
  results <- x@results
  gv <- x@geometric_variability
  n_supported <- sum(results$gower_min <= gv)
  cat_cli({
    cli::cli_h1("{S7::S7_class(x)@name}")
    cli::cli_text("Exposure: {.val {x@exposure}} ({x@exposure_type})")
    cli::cli_text("Observations: {x@n}")
    cli::cli_text("Geometric variability: {round(gv, 3)}")
    cli::cli_text(
      "Nearby radius ({x@params$nearby} x gv): {round(x@params$nearby * gv, 3)}"
    )
    cli::cli_text(
      "Mean fraction nearby: {round(mean(results$frac_nearby), 3)}"
    )
    cli::cli_text(
      "Nearest opposite within one geometric variability: {n_supported} of {x@n}"
    )
    if (x@hull_run) {
      n_in_hull <- sum(results$in_hull)
      cli::cli_text("In opposite-group hull: {n_in_hull} of {x@n}")
    } else {
      cli::cli_text("Convex-hull test: not run")
    }
  })
  invisible(x)
}
