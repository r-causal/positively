# S7 methods defined here are wired up at load time by S7::methods_register().

#' The HDR non-overlap diagnostic result class
#'
#' `hdr_result` is the S7 class returned by [check_hdr()] and [check_hdr_seq()].
#' It extends [positivity_diagnostic] with the HDR probability mass and a label
#' for the conditional-density estimator. It is created internally and is not
#' constructed directly by users.
#'
#' @usage NULL
#' @keywords internal
#' @noRd
hdr_result <- new_class(
  "hdr_result",
  parent = positivity_diagnostic,
  properties = list(
    mass = class_double,
    density_estimator = class_character
  ),
  validator = function(self) {
    if (length(self@mass) != 1) {
      "@mass must be a single value"
    } else if (length(self@density_estimator) != 1) {
      "@density_estimator must be a single label"
    }
  }
)

#' Diagnose positivity for continuous exposures with HDR non-overlap
#'
#' `check_hdr()` computes the highest-density-region (HDR) non-overlap ratio of
#' Bao and Schomaker (2025) for continuous exposures. For a common target
#' exposure value, it reports the fraction of covariate profiles whose
#' highest-density region excludes that value, that is, the share of the
#' population for which the target dose is not supported. It reports a diagnostic
#' only and assigns no severity grade.
#'
#' @details
#' For a continuous exposure \eqn{A} with covariates \eqn{L}, fix a probability
#' mass `mass` (the paper's \eqn{\alpha}). At covariate profile \eqn{l} the HDR
#' is the smallest set of exposure values capturing that mass of the conditional
#' density,
#' \deqn{A_\alpha(l) = \{ a : f(a \mid l) \ge f_\alpha(l) \},
#' \qquad P(A \in A_\alpha(l) \mid L = l) = \mathrm{mass}.}
#' The non-overlap ratio at a target `a` is the fraction of profiles whose HDR
#' excludes it,
#' \deqn{\hat{\tau}(a) = \frac{1}{n} \sum_{j=1}^{n}
#' \mathbf{1}\{ a \notin A_\alpha(l_j) \},}
#' a value in \eqn{[0, 1]}: zero when `a` is supported everywhere, one when it is
#' supported nowhere.
#'
#' The conditional density is supplied by `density_estimator`. The default
#' [hdr_density_normal()] fits `lm(exposure ~ covariates)`, treats the density
#' as Gaussian with the residual standard deviation \eqn{\hat{\sigma}}, and uses
#' the closed-form cutoff \eqn{\mathrm{dnorm}(z) / \hat{\sigma}} with
#' \eqn{z = \Phi^{-1}((1 + \mathrm{mass}) / 2)}. Membership then reduces to the
#' interval test \eqn{|a - \hat{\mu}(l)| \le z\,\hat{\sigma}}, so
#' \eqn{\hat{\tau}(a)} is the fraction of fitted means more than
#' \eqn{z\,\hat{\sigma}} from `a`.
#'
#' The non-overlap ratio is a population-level, common-target quantity. It is not
#' small merely because the density model fits well: when the exposure is
#' strongly covariate-driven, even the best central target has a nonzero floor,
#' because setting the whole population to one common dose is infeasible. The default normal estimator detects mean-shift support gaps,
#' where a stratum's supported dose moves away from a target. It does not detect
#' multimodal gaps: a hole between two modes of the true conditional density is
#' filled by the fitted normal and reported as supported. Supply a flexible
#' estimator through [new_hdr_density()] when multimodal structure is expected.
#'
#' @param .data A data frame.
#' @param .exposure The continuous exposure column, selected with data-masking.
#'   `check_hdr()` aborts for binary or categorical exposures.
#' @param .covariates The covariate columns, selected with tidyselect.
#' @param mass The HDR probability mass, a single value strictly between 0 and 1.
#'   Defaults to `0.95`.
#' @param values A numeric vector of target exposure values at which to evaluate
#'   the non-overlap ratio. `NULL` (the default) uses a 100-point grid spanning
#'   the observed exposure range.
#' @param density_estimator An `hdr_density` conditional-density estimator, built
#'   with [hdr_density_normal()] (the default) or [new_hdr_density()].
#'
#' @return An `hdr_result` object, an S7 subclass of [positivity_diagnostic].
#'   Its `@results` tibble has one row per target value with columns `value` (the
#'   target `a`) and `nonoverlap` (the ratio \eqn{\hat{\tau}(a)}). It also
#'   carries the properties `@mass` and `@density_estimator`.
#'
#' @references
#' Bao Y, Schomaker M (2025). Feasible Dose-Response Curves for Continuous
#' Treatments Under Positivity Violations.
#'
#' @examples
#' set.seed(1)
#' n <- 300
#' l <- rnorm(n)
#' dose <- rnorm(n, mean = l)
#' df <- data.frame(dose = dose, l = l)
#'
#' result <- check_hdr(df, dose, l, values = c(-2, 0, 2))
#' result
#'
#' @export
check_hdr <- function(
  .data,
  .exposure,
  .covariates,
  mass = 0.95,
  values = NULL,
  density_estimator = hdr_density_normal()
) {
  validate_data_frame(.data)
  validate_probability(mass, arg_name = "mass")
  validate_hdr_estimator(density_estimator)

  exposure_name <- select_single_exposure(rlang::enquo(.exposure), .data)
  exposure_vec <- .data[[exposure_name]]

  exposure_type <- detect_exposure_type(exposure_vec)
  if (exposure_type != "continuous") {
    abort(
      c(
        "{.fn check_hdr} supports continuous exposures only.",
        i = "{.arg .exposure} was detected as {.val {exposure_type}}."
      ),
      error_class = "positively_exposure_type_error"
    )
  }

  covariate_pos <- tidyselect::eval_select(rlang::enquo(.covariates), .data)
  validate_column_selection(covariate_pos, ".covariates")
  covariate_names <- names(covariate_pos)

  validate_numeric_columns(.data, exposure_name, ".exposure")
  validate_numeric_columns(.data, covariate_names, ".covariates")

  n <- nrow(.data)
  if (n < 2) {
    abort(
      "{.arg .data} must have at least two observations, not {n}.",
      error_class = "positively_size_error"
    )
  }

  targets <- resolve_targets(values, exposure_vec)

  analysis <- .data[c(exposure_name, covariate_names)]
  nonoverlap <- hdr_nonoverlap(
    density_estimator,
    exposure_name = exposure_name,
    covariate_names = covariate_names,
    analysis = analysis,
    exposure = as.double(exposure_vec),
    mass = mass,
    targets = targets
  )

  results <- tibble::tibble(value = targets, nonoverlap = nonoverlap)

  hdr_result(
    results = results,
    exposure = exposure_name,
    exposure_type = exposure_type,
    n = as.integer(n),
    params = list(mass = mass, values = values),
    call = rlang::current_call(),
    mass = mass,
    density_estimator = density_estimator@label
  )
}

# ---- Shared internals -----------------------------------------------------

#' Validate that a density estimator is an `hdr_density`
#'
#' @param density_estimator The object supplied to `density_estimator`.
#' @param call The calling environment, used to build the error's call.
#'
#' @return `density_estimator`, invisibly, when it is an `hdr_density`.
#' @keywords internal
#' @noRd
validate_hdr_estimator <- function(
  density_estimator,
  call = rlang::caller_env()
) {
  if (!S7::S7_inherits(density_estimator, hdr_density)) {
    supplied <- class(density_estimator)[1]
    abort(
      c(
        "{.arg density_estimator} must be an {.cls hdr_density} object.",
        i = "Build one with {.fn hdr_density_normal} or {.fn new_hdr_density}.",
        x = "Received a {.cls {supplied}}."
      ),
      error_class = "positively_type_error",
      call = call
    )
  }
  invisible(density_estimator)
}

#' Validate that every named exposure column is continuous
#'
#' Runs the exposure-type detector on each column and aborts when any is binary
#' or categorical, so the sequential variant holds the same continuous-only line
#' as [check_hdr()].
#'
#' @param .data The data frame.
#' @param exposure_names The exposure column names.
#' @param call The calling environment, used to build the error's call.
#'
#' @return `.data`, invisibly, when every exposure is continuous.
#' @keywords internal
#' @noRd
validate_continuous_exposures <- function(
  .data,
  exposure_names,
  call = rlang::caller_env()
) {
  types <- vapply(
    exposure_names,
    function(name) detect_exposure_type(.data[[name]], announce = FALSE),
    character(1)
  )
  non_continuous <- exposure_names[types != "continuous"]
  if (length(non_continuous) > 0) {
    abort(
      c(
        "{.fn check_hdr_seq} supports continuous exposures only.",
        x = "{.val {non_continuous}} {?is/are} not continuous."
      ),
      error_class = "positively_exposure_type_error",
      call = call
    )
  }
  invisible(.data)
}

#' Resolve a data-masked selection to a single exposure column name
#'
#' @param exposure_quo A quosure capturing the `.exposure` selection.
#' @param .data The data frame.
#' @param call The calling environment, used to build the error's call.
#'
#' @return The name of the single selected exposure column.
#' @keywords internal
#' @noRd
select_single_exposure <- function(
  exposure_quo,
  .data,
  call = rlang::caller_env()
) {
  exposure_pos <- tidyselect::eval_select(exposure_quo, .data)
  if (length(exposure_pos) != 1) {
    abort(
      "{.arg .exposure} must select exactly one column, not {length(exposure_pos)}.",
      error_class = "positively_selection_error",
      call = call
    )
  }
  names(exposure_pos)
}

#' Resolve the target grid of exposure values
#'
#' Returns a validated user grid, or a default 100-point grid spanning the
#' observed exposure range when `values` is `NULL`.
#'
#' @param values The user-supplied `values`, or `NULL`.
#' @param exposure The observed exposure vector.
#' @param call The calling environment, used to build the error's call.
#'
#' @return A numeric vector of target values.
#' @keywords internal
#' @noRd
resolve_targets <- function(values, exposure, call = rlang::caller_env()) {
  if (is.null(values)) {
    return(unique(seq(min(exposure), max(exposure), length.out = 100)))
  }
  if (!is.numeric(values)) {
    values_class <- class(values)[1]
    abort(
      "{.arg values} must be numeric, not a {.cls {values_class}}.",
      error_class = "positively_type_error",
      call = call
    )
  }
  if (length(values) == 0) {
    abort(
      "{.arg values} must contain at least one value.",
      error_class = "positively_empty_error",
      call = call
    )
  }
  if (anyNA(values)) {
    abort(
      "{.arg values} must not contain missing values.",
      error_class = "positively_missing_error",
      call = call
    )
  }
  if (!all(is.finite(values))) {
    abort(
      "{.arg values} must be finite.",
      error_class = "positively_range_error",
      call = call
    )
  }
  as.double(values)
}

#' The HDR non-overlap ratio at each target value
#'
#' Fits the conditional-density estimator, resolves the per-row HDR density
#' cutoff, then for every target value returns the fraction of rows whose
#' density at that value falls below the cutoff, that is, whose HDR excludes it.
#'
#' @param density_estimator An `hdr_density` estimator.
#' @param exposure_name The exposure column name.
#' @param covariate_names The covariate column names.
#' @param analysis A data frame holding the exposure and covariate columns.
#' @param exposure The observed exposure vector as a double.
#' @param mass The HDR probability mass.
#' @param targets The target exposure values.
#' @param call The calling environment, used to build the error's call.
#'
#' @return A numeric vector of non-overlap ratios, one per target value.
#' @keywords internal
#' @noRd
hdr_nonoverlap <- function(
  density_estimator,
  exposure_name,
  covariate_names,
  analysis,
  exposure,
  mass,
  targets,
  call = rlang::caller_env()
) {
  formula <- stats::reformulate(covariate_names, response = exposure_name)
  state <- density_estimator@fit(formula, analysis)
  cutoff <- hdr_thresholds(
    density_estimator,
    state = state,
    newdata = analysis,
    mass = mass,
    exposure = exposure,
    call = call
  )
  vapply(
    targets,
    function(a) {
      density <- density_estimator@density(state, a, analysis)
      mean(density < cutoff)
    },
    numeric(1)
  )
}

#' Diagnose sequential positivity for continuous exposures with HDR non-overlap
#'
#' `check_hdr_seq()` applies the HDR non-overlap ratio of Bao and Schomaker
#' (2025) to a time-varying continuous exposure, one time point at a time. At
#' each time point it fits the conditional density of that exposure on the
#' history entering its conditioning set, following the sequential-support logic
#' of sPoRT (Chatton et al.), and reports the non-overlap ratio at the target
#' values. It reports a diagnostic only and assigns no severity grade.
#'
#' @details
#' Data are supplied in wide form, one row per subject with one exposure column
#' per time point. For time point \eqn{t} the diagnostic fits the conditional
#' density of the exposure at \eqn{t} on its conditioning set, then computes the
#' point-in-time non-overlap ratio \eqn{\hat{\tau}_t(a)} exactly as
#' [check_hdr()] does. The conditioning set is the baseline covariates, the
#' time-varying covariates for each time point within the `lag` history window,
#' and the exposures at earlier time points within that window. With the default
#' `lag = Inf` the full history enters.
#'
#' Because each time point is treated separately, a support gap that moves a
#' stratum's supported dose away from a target is isolated to the time point
#' where it occurs. The default normal estimator sees mean-shift gaps only; a
#' multimodal gap at a single time point is invisible to it, as for
#' [check_hdr()].
#'
#' @param .data A data frame in wide form, one row per subject.
#' @param .exposures An ordered tidyselect of exposure columns, one per time
#'   point.
#' @param .covariates A list of tidyselect expressions, one per time point, of
#'   the time-varying covariates. A length-one list is recycled across time
#'   points. Write it as a literal `list()` call, for example
#'   `list(l0, l1, l2)`, rather than a pre-built list held in a variable.
#' @param .baseline A tidyselect of baseline covariates always included in every
#'   conditioning set. Defaults to `NULL`.
#' @param mass The HDR probability mass, a single value strictly between 0 and 1.
#'   Defaults to `0.95`.
#' @param values A numeric vector of target exposure values. `NULL` (the
#'   default) uses a 100-point grid spanning the pooled exposure range.
#' @param density_estimator An `hdr_density` conditional-density estimator, built
#'   with [hdr_density_normal()] (the default) or [new_hdr_density()].
#' @param lag The history window: the number of earlier time points whose
#'   covariates and exposures enter each conditioning set. Defaults to `Inf`, the
#'   full history.
#'
#' @return An `hdr_result` object, an S7 subclass of [positivity_diagnostic].
#'   Its `@results` tibble has one row per time point and target value with
#'   columns `time`, `value`, and `nonoverlap`.
#'
#' @references
#' Bao Y, Schomaker M (2025). Feasible Dose-Response Curves for Continuous
#' Treatments Under Positivity Violations.
#'
#' Chatton A, Schomaker M, Luque-Fernandez MA, Platt RW, Schnitzer ME (2025).
#' Is checking for sequential positivity violations getting you down? Try
#' sPoRT!
#'
#' @examples
#' set.seed(1)
#' n <- 200
#' df <- data.frame(
#'   l0 = rnorm(n),
#'   a1 = rnorm(n),
#'   l1 = rnorm(n),
#'   a2 = rnorm(n),
#'   l2 = rnorm(n),
#'   a3 = rnorm(n)
#' )
#'
#' result <- check_hdr_seq(
#'   df,
#'   c(a1, a2, a3),
#'   list(l0, l1, l2),
#'   values = c(-2, 0, 2)
#' )
#' result
#'
#' @export
check_hdr_seq <- function(
  .data,
  .exposures,
  .covariates,
  .baseline = NULL,
  mass = 0.95,
  values = NULL,
  density_estimator = hdr_density_normal(),
  lag = Inf
) {
  validate_data_frame(.data)
  validate_probability(mass, arg_name = "mass")
  validate_hdr_estimator(density_estimator)
  validate_lag(lag)

  exposure_pos <- tidyselect::eval_select(rlang::enquo(.exposures), .data)
  validate_column_selection(exposure_pos, ".exposures")
  exposure_names <- names(exposure_pos)
  n_times <- length(exposure_names)

  validate_continuous_exposures(.data, exposure_names)

  covariate_sets <- parse_covariate_list(
    rlang::enquo(.covariates),
    .data,
    n_times = n_times
  )
  baseline_names <- eval_optional_selection(rlang::enquo(.baseline), .data)

  validate_numeric_columns(.data, exposure_names, ".exposures")
  validate_numeric_columns(
    .data,
    unique(unlist(covariate_sets)),
    ".covariates"
  )
  if (length(baseline_names) > 0) {
    validate_numeric_columns(.data, baseline_names, ".baseline")
  }

  n <- nrow(.data)
  if (n < 2) {
    abort(
      "{.arg .data} must have at least two observations, not {n}.",
      error_class = "positively_size_error"
    )
  }

  pooled_exposure <- as.double(unlist(.data[exposure_names], use.names = FALSE))
  targets <- resolve_targets(values, pooled_exposure)

  seq_call <- rlang::current_env()
  per_time <- purrr::map(seq_len(n_times), function(t) {
    conditioning <- conditioning_set(
      time = t,
      covariate_sets = covariate_sets,
      exposure_names = exposure_names,
      baseline_names = baseline_names,
      lag = lag
    )
    exposure_name <- exposure_names[[t]]
    analysis <- .data[c(exposure_name, conditioning)]
    nonoverlap <- hdr_nonoverlap(
      density_estimator,
      exposure_name = exposure_name,
      covariate_names = conditioning,
      analysis = analysis,
      exposure = as.double(.data[[exposure_name]]),
      mass = mass,
      targets = targets,
      call = seq_call
    )
    tibble::tibble(
      time = rep(t, length(targets)),
      value = targets,
      nonoverlap = nonoverlap
    )
  })
  results <- vctrs::vec_rbind(!!!per_time)

  hdr_result(
    results = results,
    exposure = exposure_names,
    exposure_type = "continuous",
    n = as.integer(n),
    params = list(mass = mass, values = values, lag = lag),
    call = rlang::current_call(),
    mass = mass,
    density_estimator = density_estimator@label
  )
}

#' Parse a list of tidyselect expressions into per-time covariate names
#'
#' Accepts the `list(...)` form of `.covariates`, one selection per time point,
#' and returns a list of character vectors of resolved column names. A single
#' selection (length-one list or a bare selection) is recycled across all time
#' points.
#'
#' @param covariates_quo A quosure capturing `.covariates`.
#' @param .data The data frame.
#' @param n_times The number of time points.
#' @param call The calling environment, used to build the error's call.
#'
#' @return A list of length `n_times` of character vectors of column names.
#' @keywords internal
#' @noRd
parse_covariate_list <- function(
  covariates_quo,
  .data,
  n_times,
  call = rlang::caller_env()
) {
  expr <- rlang::quo_get_expr(covariates_quo)
  env <- rlang::quo_get_env(covariates_quo)
  is_list_call <- rlang::is_call(expr, "list")
  elements <- if (is_list_call) rlang::call_args(expr) else list(expr)

  sets <- lapply(elements, function(element) {
    element_quo <- rlang::new_quosure(element, env)
    selection <- tidyselect::eval_select(element_quo, .data)
    names(selection)
  })

  if (length(sets) == 1 && n_times > 1) {
    sets <- rep(sets, n_times)
  }
  if (length(sets) != n_times) {
    abort(
      c(
        "{.arg .covariates} must supply one selection per time point.",
        i = "Found {length(sets)} selection{?s} for {n_times} exposure{?s}."
      ),
      error_class = "positively_selection_error",
      call = call
    )
  }
  sets
}

#' Evaluate an optional data-masked selection to column names
#'
#' @param selection_quo A quosure capturing the selection.
#' @param .data The data frame.
#'
#' @return A character vector of column names, empty when the selection is
#'   `NULL`.
#' @keywords internal
#' @noRd
eval_optional_selection <- function(selection_quo, .data) {
  if (rlang::quo_is_null(selection_quo)) {
    return(character(0))
  }
  names(tidyselect::eval_select(selection_quo, .data))
}

#' The conditioning set for one time point
#'
#' Assembles the baseline covariates, the time-varying covariates for each time
#' point within the `lag` window, and the exposures at earlier time points
#' within that window.
#'
#' @param time The current time point index.
#' @param covariate_sets The per-time covariate name list.
#' @param exposure_names The exposure column names, time-ordered.
#' @param baseline_names The baseline covariate names.
#' @param lag The history window.
#'
#' @return A character vector of unique conditioning-set column names.
#' @keywords internal
#' @noRd
conditioning_set <- function(
  time,
  covariate_sets,
  exposure_names,
  baseline_names,
  lag
) {
  earliest <- if (is.infinite(lag)) 1L else max(1L, time - lag)
  covariate_names <- unlist(
    covariate_sets[earliest:time],
    use.names = FALSE
  )
  prior_exposures <- if (time > earliest) {
    exposure_names[earliest:(time - 1)]
  } else {
    character(0)
  }
  unique(c(baseline_names, covariate_names, prior_exposures))
}

# ---- Methods --------------------------------------------------------------

method(print, hdr_result) <- function(x, ...) {
  nonoverlap <- x@results$nonoverlap
  n_target <- length(unique(x@results$value))
  cat_cli({
    cli::cli_h1("{S7::S7_class(x)@name}")
    cli::cli_text("Exposure: {.val {x@exposure}} ({x@exposure_type})")
    cli::cli_text("Observations: {x@n}")
    cli::cli_text("HDR mass: {x@mass}")
    cli::cli_text("Density estimator: {x@density_estimator}")
    if ("time" %in% names(x@results)) {
      cli::cli_text("Time points: {length(unique(x@results$time))}")
    }
    cli::cli_text(
      "Non-overlap over {n_target} target{?s}: {round(min(nonoverlap), 3)} to {round(max(nonoverlap), 3)}"
    )
  })
  invisible(x)
}
