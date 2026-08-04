# S7 methods defined here are wired up at load time by S7::methods_register().

#' The effective-data-points diagnostic result class
#'
#' `edp_result` is the S7 class returned by [check_edp()]. It extends
#' [positivity_diagnostic] with the variant that produced it and the resolved
#' kernel half-distances. It is created internally and is not constructed
#' directly by users.
#'
#' @keywords internal
#' @noRd
edp_result <- new_class(
  "edp_result",
  parent = positivity_diagnostic,
  properties = list(
    variant = class_character,
    bandwidths = class_list
  ),
  validator = function(self) {
    if (
      length(self@variant) != 1 || !self@variant %in% c("data", "estimator")
    ) {
      "@variant must be a single value: \"data\" or \"estimator\""
    }
  }
)

#' Diagnose positivity with effective data points
#'
#' `check_edp()` computes the effective-data-points (EDP) positivity diagnostic
#' of Ring and Schomaker (2026) for binary, categorical, and continuous
#' exposures. For every observation moved to a candidate intervention value, it
#' sums a product kernel over the sample to measure how much observed support
#' surrounds that intervened-on point.
#'
#' @details
#' For each observation \eqn{i} the diagnostic forms the intervened-on point
#' \eqn{o^*_i = (l_i, a^*)}, holding the covariates at their observed values and
#' setting the exposure to an intervention value \eqn{a^*}. The effective data
#' points at that point are
#' \deqn{\mathrm{EDP}_i(a^*) = \sum_{j=1}^{n} k(o_j, o^*_i),}
#' a sum over the whole sample of a product kernel that multiplies one factor per
#' dimension. Each factor lies in \eqn{[0, 1]}, so \eqn{0 \le \mathrm{EDP} \le
#' n}. A continuous dimension with half-distance \eqn{h} contributes
#' \eqn{0.5^{(\delta / h)^2}}, which equals `1` at \eqn{\delta = 0} and exactly
#' `0.5` at \eqn{|\delta| = h}; this is a reparameterized Gaussian kernel. A
#' categorical dimension contributes `1` on a match and `categorical_similarity`
#' otherwise. Numeric covariates and continuous exposures are continuous
#' dimensions; factor and character covariates, together with binary and
#' categorical exposures, are categorical dimensions.
#'
#' The half-distances follow the paper's rule of thumb: `1 * sd` for each
#' numeric covariate and `0.5 * sd` for a continuous exposure. These are marginal
#' scales, so a mid-support gap in a strongly multimodal exposure can be masked
#' by an inflated standard deviation; supply `bw_exposure` or `bw_covariates`
#' directly when a robust scale is preferable. As the half-distances grow without
#' bound every kernel factor approaches `1` and EDP approaches `n`, which renders
#' the diagnostic uninformative; at a half-distance of `0` a continuous dimension
#' counts only exact matches.
#'
#' EDP is a product kernel, so adding dimensions can only lower it. It is
#' meaningful relative across observations and strata for a fixed covariate set,
#' never against a universal threshold; compare a stratum against another stratum
#' or a candidate value against another candidate value rather than against a
#' fixed cutoff.
#'
#' The estimator variant reframes EDP around the two nuisance models of a causal
#' estimator. `edp_outcome` conditions on the exposure and the outcome-model
#' covariates \eqn{(A, L)}; `edp_treatment` conditions on the treatment-model
#' covariates \eqn{(L)} alone. Because the outcome measure carries the extra
#' exposure dimension, `edp_outcome <= edp_treatment` whenever the two covariate
#' sets are equal. The reported `ideal_weight` is `edp_treatment / edp_outcome`,
#' proportional to the inverse of the estimated treatment density at the
#' intervention, so it grows where treatment support is thin. Where a target has
#' no outcome-model support, `edp_outcome` is zero and `ideal_weight` is
#' infinite, the honest value of one over an estimated density of zero.
#'
#' @param .data A data frame.
#' @param .exposure The exposure column, selected with data-masking.
#' @param .covariates The covariate columns, selected with tidyselect.
#' @param .outcome_covariates The outcome-model covariate columns for the
#'   estimator variant, selected with tidyselect. Defaults to `.covariates`.
#' @param .treatment_covariates The treatment-model covariate columns for the
#'   estimator variant, selected with tidyselect. Defaults to `.covariates`.
#' @param values The intervention values \eqn{a^*}. `NULL` (the default) uses
#'   every level of a factor exposure and every observed value of a binary or
#'   character exposure, and the deciles of the observed exposure for continuous
#'   exposures.
#' @param variant One of `"data"` (report a single `edp` per point, the default)
#'   or `"estimator"` (report `edp_outcome`, `edp_treatment`, and
#'   `ideal_weight`).
#' @param kernel The continuous-dimension kernel. Only `"gaussian"` is currently
#'   available; categorical dimensions always use the match kernel.
#' @param bw_exposure The continuous-exposure half-distance. `NULL` (the default)
#'   uses `0.5 * sd(exposure)`. A half-distance of `0` counts exact matches.
#' @param bw_covariates The numeric-covariate half-distance, applied to every
#'   numeric covariate. `NULL` (the default) uses `1 * sd` per covariate.
#' @param categorical_similarity The kernel value in `[0, 1]` for non-matching
#'   categories, shared by categorical exposures and factor or character
#'   covariates. Defaults to `0`. A positive value masks categorical violations.
#' @param exposure_type One of `"auto"` (detect from the data, the default),
#'   `"binary"`, `"categorical"`, or `"continuous"`. A supplied type is
#'   authoritative and detection is not consulted, so it is rejected only when
#'   the exposure column cannot carry it: `"continuous"` needs a numeric column
#'   and `"binary"` needs exactly two distinct values, while `"categorical"` asks
#'   nothing of the column.
#'
#' @return An `edp_result` object, an S7 subclass of [positivity_diagnostic].
#'   Its `@results` tibble has one row per observation and intervention value.
#'   The data variant carries columns `.id` (the observation row index), `value`
#'   (the intervention value \eqn{a^*}), and `edp`. The estimator variant
#'   replaces `edp` with `edp_outcome`, `edp_treatment`, and `ideal_weight`. It
#'   also carries the `@variant` and `@bandwidths` properties.
#'
#'   [generics::glance()] returns a one-row tibble with `n` (the sample size),
#'   `variant`, `n_values` (the number of intervention values), and the observed
#'   range of every measure the variant reports. That is `edp_min` and `edp_max`
#'   for the data variant, and `edp_outcome_min`, `edp_outcome_max`,
#'   `edp_treatment_min`, `edp_treatment_max`, `ideal_weight_min`, and
#'   `ideal_weight_max` for the estimator variant.
#'
#' @references
#' Ring C, Schomaker M (2026). A Diagnostic to Find and Help Combat Stochastic
#' Positivity Issues, with a Focus on Continuous Treatments.
#'
#' @examples
#' set.seed(1)
#' n <- 100
#' x1 <- rnorm(n)
#' dose <- rnorm(n, mean = x1)
#' df <- data.frame(dose = dose, x1 = x1)
#'
#' result <- check_edp(df, dose, x1, values = c(0, 1), exposure_type = "continuous")
#' result
#'
#' @export
check_edp <- function(
  .data,
  .exposure,
  .covariates,
  .outcome_covariates = NULL,
  .treatment_covariates = NULL,
  values = NULL,
  variant = c("data", "estimator"),
  kernel = c("gaussian"),
  bw_exposure = NULL,
  bw_covariates = NULL,
  categorical_similarity = 0,
  exposure_type = c("auto", "binary", "categorical", "continuous")
) {
  validate_data_frame(.data)
  variant <- resolve_arg_match(rlang::arg_match(variant))
  kernel <- resolve_arg_match(rlang::arg_match(kernel))
  validate_unit_scalar(categorical_similarity, "categorical_similarity")
  if (!is.null(bw_exposure)) {
    validate_nonnegative_scalar(bw_exposure, "bw_exposure")
  }
  if (!is.null(bw_covariates)) {
    validate_nonnegative_scalar(bw_covariates, "bw_covariates")
  }

  exposure_pos <- eval_select_column(
    rlang::enquo(.exposure),
    .data,
    ".exposure"
  )
  exposure_name <- names(exposure_pos)
  exposure_vec <- .data[[exposure_pos]]

  covariate_pos <- eval_select_columns(
    rlang::enquo(.covariates),
    .data,
    ".covariates"
  )
  validate_column_selection(covariate_pos, ".covariates")
  covariate_names <- names(covariate_pos)

  outcome_covariates_quo <- rlang::enquo(.outcome_covariates)
  treatment_covariates_quo <- rlang::enquo(.treatment_covariates)
  outcome_names <- resolve_covariate_names(
    outcome_covariates_quo,
    .data,
    covariate_names,
    ".outcome_covariates"
  )
  treatment_names <- resolve_covariate_names(
    treatment_covariates_quo,
    .data,
    covariate_names,
    ".treatment_covariates"
  )
  if (variant == "data") {
    warn_unused_estimator_covariates(
      outcome_covariates_quo,
      treatment_covariates_quo
    )
  }

  validate_complete_columns(.data, exposure_name, ".exposure")
  validate_finite_columns(.data, exposure_name, ".exposure")
  validate_complete_columns(.data, covariate_names, ".covariates")
  validate_finite_columns(.data, covariate_names, ".covariates")
  if (variant == "estimator") {
    validate_complete_columns(.data, outcome_names, ".outcome_covariates")
    validate_finite_columns(.data, outcome_names, ".outcome_covariates")
    validate_complete_columns(.data, treatment_names, ".treatment_covariates")
    validate_finite_columns(.data, treatment_names, ".treatment_covariates")
  }

  n <- nrow(.data)
  if (n < 2) {
    abort(
      "{.arg .data} must have at least two rows, not {n}.",
      error_class = "positively_size_error"
    )
  }

  exposure_type <- resolve_exposure_type(
    exposure_type,
    exposure_vec,
    supported = diagnostic_supported_types()[["edp"]],
    fn = "check_edp"
  )

  if (!is.null(values)) {
    validate_values(values, exposure_type)
  }

  compute_names <- if (variant == "estimator") {
    unique(c(outcome_names, treatment_names))
  } else {
    covariate_names
  }
  warn_unused_bandwidths(
    bw_exposure,
    bw_covariates,
    exposure_type,
    .data,
    compute_names
  )

  bw_exposure_resolved <- if (exposure_type == "continuous") {
    bw_exposure %||% (0.5 * stats::sd(as.double(exposure_vec)))
  } else {
    NA_real_
  }
  values_grid <- resolve_values(values, exposure_vec, exposure_type)

  exposure_kernels <- lapply(
    values_grid,
    function(a) {
      exposure_kernel(
        exposure_vec,
        a,
        exposure_type,
        bw_exposure_resolved,
        categorical_similarity
      )
    }
  )

  ids <- rep(seq_len(n), times = length(values_grid))
  value_column <- rep(values_grid, each = n)

  if (variant == "data") {
    cov_bw <- resolve_covariate_bandwidths(
      .data,
      covariate_names,
      bw_covariates
    )
    similarity <- covariate_similarity(
      .data,
      covariate_names,
      cov_bw,
      categorical_similarity
    )
    edp <- unlist(
      lapply(exposure_kernels, function(e) as.double(similarity %*% e)),
      use.names = FALSE
    )
    results <- tibble::tibble(.id = ids, value = value_column, edp = edp)
    bandwidths <- list(exposure = bw_exposure_resolved, covariates = cov_bw)
  } else {
    outcome_bw <- resolve_covariate_bandwidths(
      .data,
      outcome_names,
      bw_covariates
    )
    treatment_bw <- resolve_covariate_bandwidths(
      .data,
      treatment_names,
      bw_covariates
    )
    outcome_similarity <- covariate_similarity(
      .data,
      outcome_names,
      outcome_bw,
      categorical_similarity
    )
    treatment_similarity <- covariate_similarity(
      .data,
      treatment_names,
      treatment_bw,
      categorical_similarity
    )
    edp_outcome <- unlist(
      lapply(exposure_kernels, function(e) as.double(outcome_similarity %*% e)),
      use.names = FALSE
    )
    edp_treatment <- rep(
      rowSums(treatment_similarity),
      times = length(values_grid)
    )
    results <- tibble::tibble(
      .id = ids,
      value = value_column,
      edp_outcome = edp_outcome,
      edp_treatment = edp_treatment,
      ideal_weight = edp_treatment / edp_outcome
    )
    bandwidths <- list(
      exposure = bw_exposure_resolved,
      outcome_covariates = outcome_bw,
      treatment_covariates = treatment_bw
    )
  }

  edp_result(
    results = results,
    exposure = exposure_name,
    exposure_type = exposure_type,
    n = as.integer(n),
    params = list(
      variant = variant,
      kernel = kernel,
      values = values_grid,
      categorical_similarity = categorical_similarity,
      bw_exposure = bw_exposure_resolved,
      bw_covariates = bw_covariates,
      covariates = covariate_names,
      outcome_covariates = outcome_names,
      treatment_covariates = treatment_names
    ),
    call = rlang::current_call(),
    variant = variant,
    bandwidths = bandwidths
  )
}

# ---- Argument helpers -----------------------------------------------------

#' Resolve an optional covariate selection, defaulting to the base covariates
#'
#' @param quo A quosure from an optional data-masking covariate argument.
#' @param .data The data frame.
#' @param default The base covariate names used when the argument is `NULL`.
#' @param arg_name The argument name used in error messages.
#' @param call The calling environment, used to build the error's call.
#'
#' @return A character vector of covariate column names.
#' @keywords internal
#' @noRd
resolve_covariate_names <- function(
  quo,
  .data,
  default,
  arg_name,
  call = rlang::caller_env()
) {
  if (rlang::quo_is_null(quo)) {
    return(default)
  }
  selection <- eval_select_columns(quo, .data, arg_name, call = call)
  validate_column_selection(selection, arg_name, call = call)
  names(selection)
}

#' Resolve the intervention-value grid
#'
#' @param values The user-supplied values, or `NULL` for the default grid.
#' @param exposure_vec The observed exposure vector.
#' @param exposure_type The resolved exposure type.
#'
#' @return A vector of intervention values.
#' @keywords internal
#' @noRd
resolve_values <- function(values, exposure_vec, exposure_type) {
  if (!is.null(values)) {
    return(values)
  }
  if (exposure_type == "continuous") {
    unique(stats::quantile(
      as.double(exposure_vec),
      probs = seq(0.1, 0.9, by = 0.1),
      names = FALSE,
      type = 7
    ))
  } else if (is.factor(exposure_vec)) {
    levels(exposure_vec)
  } else {
    sort(unique(exposure_vec))
  }
}

#' Validate a user-supplied intervention-value grid
#'
#' @param values The user-supplied intervention values.
#' @param exposure_type The resolved exposure type.
#' @param arg_name The argument name used in error messages.
#' @param call The calling environment, used to build the error's call.
#'
#' @return `values`, invisibly, when it is a non-empty, complete grid that is
#'   numeric for a continuous exposure.
#' @keywords internal
#' @noRd
validate_values <- function(
  values,
  exposure_type,
  arg_name = "values",
  call = rlang::caller_env()
) {
  if (length(values) == 0) {
    abort(
      "{.arg {arg_name}} must contain at least one value.",
      error_class = "positively_empty_error",
      call = call
    )
  }
  if (anyNA(values)) {
    abort(
      "{.arg {arg_name}} must not contain missing values.",
      error_class = "positively_missing_error",
      call = call
    )
  }
  if (exposure_type == "continuous" && !is.numeric(values)) {
    values_class <- class(values)[1]
    abort(
      "{.arg {arg_name}} must be numeric for a continuous exposure, not a {.cls {values_class}}.",
      error_class = "positively_type_error",
      call = call
    )
  }
  invisible(values)
}

#' Warn when a supplied bandwidth reaches no continuous dimension
#'
#' Raises a classed warning when `bw_exposure` is supplied but the exposure is
#' not continuous, or when `bw_covariates` is supplied but none of the covariates
#' that enter the computation is numeric. In either case the half-distance would
#' be silently discarded.
#'
#' @param bw_exposure The user-supplied exposure half-distance, or `NULL`.
#' @param bw_covariates The user-supplied covariate half-distance, or `NULL`.
#' @param exposure_type The resolved exposure type.
#' @param .data The data frame.
#' @param compute_names The covariate columns that enter the computation.
#' @param call The calling environment, used to build the warning's call.
#'
#' @return Invisibly `NULL`; raises warnings as a side effect.
#' @keywords internal
#' @noRd
warn_unused_bandwidths <- function(
  bw_exposure,
  bw_covariates,
  exposure_type,
  .data,
  compute_names,
  call = rlang::caller_env()
) {
  if (!is.null(bw_exposure) && exposure_type != "continuous") {
    warn(
      c(
        "{.arg bw_exposure} was supplied but the exposure is not continuous.",
        i = "The half-distance is ignored for a {exposure_type} exposure."
      ),
      warning_class = "positively_unused_arg_warning",
      call = call
    )
  }
  has_numeric <- any(vapply(
    compute_names,
    function(column) is.numeric(.data[[column]]),
    logical(1)
  ))
  if (!is.null(bw_covariates) && !has_numeric) {
    warn(
      c(
        "{.arg bw_covariates} was supplied but no covariate is numeric.",
        i = "The half-distance is ignored when every covariate is categorical."
      ),
      warning_class = "positively_unused_arg_warning",
      call = call
    )
  }
  invisible()
}

#' Warn when estimator covariates are supplied to the data variant
#'
#' The `.outcome_covariates` and `.treatment_covariates` selections shape the
#' estimator variant alone. Supplying either under the data variant has no
#' effect, so a classed warning flags the discarded selection.
#'
#' @param outcome_quo The captured `.outcome_covariates` selection.
#' @param treatment_quo The captured `.treatment_covariates` selection.
#' @param call The calling environment, used to build the warning's call.
#'
#' @return Invisibly `NULL`; raises a warning as a side effect.
#' @keywords internal
#' @noRd
warn_unused_estimator_covariates <- function(
  outcome_quo,
  treatment_quo,
  call = rlang::caller_env()
) {
  supplied <- character(0)
  if (!rlang::quo_is_null(outcome_quo)) {
    supplied <- c(supplied, ".outcome_covariates")
  }
  if (!rlang::quo_is_null(treatment_quo)) {
    supplied <- c(supplied, ".treatment_covariates")
  }
  if (length(supplied) > 0) {
    warn(
      c(
        "{.arg {supplied}} {?applies/apply} only to the estimator variant.",
        i = "The {cli::qty(supplied)}selection{?s} {?is/are} ignored when {.arg variant} is {.val {\"data\"}}."
      ),
      warning_class = "positively_unused_arg_warning",
      call = call
    )
  }
  invisible()
}

#' Validate a single number in the closed unit interval
#'
#' @param x A candidate value.
#' @param arg_name The argument name used in error messages.
#' @param call The calling environment, used to build the error's call.
#'
#' @return `x`, invisibly, when it is a single value in `[0, 1]`.
#' @keywords internal
#' @noRd
validate_unit_scalar <- function(x, arg_name, call = rlang::caller_env()) {
  if (!is.numeric(x) || length(x) != 1 || is.na(x) || !is.finite(x)) {
    abort(
      "{.arg {arg_name}} must be a single number.",
      error_class = "positively_type_error",
      call = call
    )
  }
  if (x < 0 || x > 1) {
    abort(
      c(
        "{.arg {arg_name}} must be between {.val {0}} and {.val {1}}.",
        x = "Found {.val {x}}."
      ),
      error_class = "positively_range_error",
      call = call
    )
  }
  invisible(x)
}

#' Validate a single non-negative number
#'
#' @param x A candidate value.
#' @param arg_name The argument name used in error messages.
#' @param call The calling environment, used to build the error's call.
#'
#' @return `x`, invisibly, when it is a single finite value of at least `0`.
#' @keywords internal
#' @noRd
validate_nonnegative_scalar <- function(
  x,
  arg_name,
  call = rlang::caller_env()
) {
  if (!is.numeric(x) || length(x) != 1 || is.na(x) || !is.finite(x)) {
    abort(
      "{.arg {arg_name}} must be a single finite number.",
      error_class = "positively_type_error",
      call = call
    )
  }
  if (x < 0) {
    abort(
      c(
        "{.arg {arg_name}} must not be negative.",
        x = "Found {.val {x}}."
      ),
      error_class = "positively_range_error",
      call = call
    )
  }
  invisible(x)
}

# ---- Kernel machinery -----------------------------------------------------

#' The exposure-dimension kernel for one intervention value
#'
#' @param exposure_vec The observed exposure vector.
#' @param a A single intervention value.
#' @param exposure_type The resolved exposure type.
#' @param bw The continuous-exposure half-distance (ignored when discrete).
#' @param categorical_similarity The non-matching kernel value for discrete
#'   exposures.
#'
#' @return A numeric vector, one kernel value per observation.
#' @keywords internal
#' @noRd
exposure_kernel <- function(
  exposure_vec,
  a,
  exposure_type,
  bw,
  categorical_similarity
) {
  if (exposure_type == "continuous") {
    continuous_kernel(as.double(exposure_vec) - a, bw)
  } else {
    match_kernel(exposure_vec, a, categorical_similarity)
  }
}

#' The reparameterized Gaussian kernel over a difference vector
#'
#' Returns `0.5 ^ ((delta / bw) ^ 2)`, which equals `1` at `delta == 0` and
#' `0.5` at `abs(delta) == bw`. A half-distance of `0` counts exact matches.
#'
#' @param delta A numeric vector of differences.
#' @param bw The half-distance.
#'
#' @return A numeric vector in `[0, 1]`.
#' @keywords internal
#' @noRd
continuous_kernel <- function(delta, bw) {
  if (bw == 0) {
    return(as.double(delta == 0))
  }
  0.5^((delta / bw)^2)
}

#' Blend match indicators with the non-matching kernel value
#'
#' Maps a logical match indicator to `1` on a match and `categorical_similarity`
#' otherwise. This is the single source of the categorical kernel formula, shared
#' by the one-value match kernel and the pairwise similarity matrix.
#'
#' @param matches A logical vector or matrix of match indicators.
#' @param categorical_similarity The kernel value for non-matching categories.
#'
#' @return A numeric object of the same shape as `matches`.
#' @keywords internal
#' @noRd
blend_match <- function(matches, categorical_similarity) {
  categorical_similarity + (1 - categorical_similarity) * matches
}

#' The match kernel comparing a vector against one value
#'
#' @param x A vector of observed categories.
#' @param value A single comparison value.
#' @param categorical_similarity The kernel value for non-matching categories.
#'
#' @return A numeric vector: `1` on a match, `categorical_similarity` otherwise.
#' @keywords internal
#' @noRd
match_kernel <- function(x, value, categorical_similarity) {
  blend_match(as.character(x) == as.character(value), categorical_similarity)
}

#' Resolve per-covariate half-distances for the numeric covariates
#'
#' @param .data The data frame.
#' @param columns The covariate column names.
#' @param bw_covariates The supplied half-distance, or `NULL` for `1 * sd`.
#'
#' @return A named numeric vector of half-distances for the numeric covariates.
#' @keywords internal
#' @noRd
resolve_covariate_bandwidths <- function(.data, columns, bw_covariates) {
  numeric_columns <- columns[vapply(
    columns,
    function(column) is.numeric(.data[[column]]),
    logical(1)
  )]
  vapply(
    numeric_columns,
    function(column) {
      if (is.null(bw_covariates)) {
        stats::sd(.data[[column]])
      } else {
        bw_covariates
      }
    },
    numeric(1)
  )
}

#' The pairwise covariate similarity matrix
#'
#' Multiplies one per-dimension kernel matrix per covariate: numeric covariates
#' use the reparameterized Gaussian kernel, factor and character covariates use
#' the match kernel. Entry `[i, j]` is the product-kernel similarity between the
#' covariates of observations `i` and `j`.
#'
#' @param .data The data frame.
#' @param columns The covariate column names.
#' @param cov_bw A named numeric vector of half-distances for numeric covariates.
#' @param categorical_similarity The non-matching kernel value.
#'
#' @return An `n` by `n` numeric matrix.
#' @keywords internal
#' @noRd
covariate_similarity <- function(
  .data,
  columns,
  cov_bw,
  categorical_similarity
) {
  n <- nrow(.data)
  similarity <- matrix(1, nrow = n, ncol = n)
  for (column in columns) {
    values <- .data[[column]]
    if (is.numeric(values)) {
      dimension <- continuous_kernel(
        outer(values, values, "-"),
        cov_bw[[column]]
      )
    } else {
      characters <- as.character(values)
      dimension <- blend_match(
        outer(characters, characters, "=="),
        categorical_similarity
      )
    }
    similarity <- similarity * dimension
  }
  similarity
}

#' The measures an EDP variant reports
#'
#' @param variant Either `"data"` or `"estimator"`.
#'
#' @return A character vector of `@results` column names.
#' @keywords internal
#' @noRd
edp_measures <- function(variant) {
  if (variant == "estimator") {
    c("edp_outcome", "edp_treatment", "ideal_weight")
  } else {
    "edp"
  }
}

# ---- Methods --------------------------------------------------------------

method(glance, edp_result) <- function(x, ...) {
  ranges <- do.call(
    c,
    lapply(edp_measures(x@variant), function(measure) {
      values <- x@results[[measure]]
      bounds <- list(min(values), max(values))
      names(bounds) <- paste0(measure, c("_min", "_max"))
      bounds
    })
  )
  tibble::as_tibble(c(
    list(
      n = x@n,
      variant = x@variant,
      n_values = length(x@params$values)
    ),
    ranges
  ))
}

method(diagnostic_label, edp_result) <- function(x) {
  "Effective data points"
}

method(diagnostic_headline, edp_result) <- function(x) {
  measure <- edp_primary_column(x@results)
  values <- x@results[[measure]]
  variant <- if (x@variant == "data") "Data" else "Estimator"
  n_values <- length(x@params$values)
  cli::format_inline(
    "{variant} variant over {n_values} intervention value{?s}; {measure} {round(min(values), 3)} to {round(max(values), 3)}"
  )
}

method(print, edp_result) <- function(x, ...) {
  results <- x@results
  columns <- edp_measures(x@variant)
  cat_cli({
    cli::cli_h1("{diagnostic_label(x)}")
    cli::cli_text("Exposure: {.val {x@exposure}} ({x@exposure_type})")
    cli::cli_text("Observations: {x@n}")
    cli::cli_text("Variant: {x@variant}")
    cli::cli_text("Intervention values: {length(x@params$values)}")
    for (column in columns) {
      values <- results[[column]]
      cli::cli_text(
        "{column} range: {round(min(values), 3)} to {round(max(values), 3)}"
      )
    }
  })
  invisible(x)
}

#' The primary EDP measure column of a results tibble
#'
#' The data variant reports `edp`; the estimator variant reports `edp_outcome`.
#' This is the single source of that rule, shared by the print and plot methods.
#'
#' @param results An `edp_result` results tibble.
#'
#' @return The name of the primary EDP column.
#' @keywords internal
#' @noRd
edp_primary_column <- function(results) {
  if ("edp" %in% names(results)) "edp" else "edp_outcome"
}
