# S7 methods defined here are wired up at load time by S7::methods_register().

#' The ETA.Bias diagnostic result class
#'
#' `eta_bias_result` is the S7 class returned by [check_eta_bias()]. It extends
#' [positivity_diagnostic] with the estimator label, the `truth` against which
#' the bootstrap estimates are compared (one value per estimand term, named by
#' term), and the list of bootstrap-estimate vectors (one per row of `@results`,
#' that is per term and truncation level). It is created internally and is not
#' constructed directly by users.
#'
#' @keywords internal
#' @noRd
eta_bias_result <- new_class(
  "eta_bias_result",
  parent = positivity_diagnostic,
  properties = list(
    estimator = class_character,
    truth = class_double,
    boot_estimates = class_list
  ),
  validator = function(self) {
    n_terms <- length(unique(self@results$term))
    if (length(self@estimator) != 1) {
      "@estimator must be a single value"
    } else if (length(self@truth) != n_terms) {
      "@truth must have one value per term in @results"
    } else if (length(self@boot_estimates) != nrow(self@results)) {
      "@boot_estimates must have one vector per row of @results"
    }
  }
)

#' Diagnose positivity bias with the parametric bootstrap (ETA.Bias)
#'
#' `check_eta_bias()` implements the ETA.Bias parametric-bootstrap diagnostic of
#' Petersen et al. (2012) for point exposures. It estimates the bias that
#' a chosen causal estimator would incur under the data's own fitted treatment
#' and outcome mechanisms, isolating the component of bias that stems from
#' positivity violations, weight truncation, and finite-sample sparsity.
#'
#' @details
#' The diagnostic (the paper's Section 4.1) treats the fitted data-generating
#' mechanism as a known truth and asks how a candidate estimator behaves under
#' it. The algorithm has three steps.
#'
#' First, it fits the treatment mechanism \eqn{g_n(1 \mid W)}, a main-effects
#' logistic model of the exposure on the covariates, and the outcome regression
#' \eqn{\bar{Q}_n(A, W)}, a main-effects linear or logistic model of the outcome
#' on the exposure and covariates. Factor and character covariates enter these
#' models as indicator terms, exactly as [stats::glm()] expands them. An
#' exposure of more than two levels fits the treatment mechanism as a
#' multinomial logit and so requires the \pkg{nnet} package. The target
#' of inference is a single scalar,
#' the G-computation estimate on the observed data,
#' \eqn{\psi = \frac{1}{n}\sum_i [\bar{Q}_n(1, W_i) - \bar{Q}_n(0, W_i)]}. This
#' `truth` depends on the outcome model alone. It is identical across estimators
#' and constant across truncation levels.
#'
#' Second, it draws `n_boot` bootstrap datasets. Each resamples the covariate
#' rows with replacement, draws a bootstrap exposure
#' \eqn{A^* \sim \mathrm{Bernoulli}(g_n(1 \mid W^*))} from the untruncated fitted
#' propensity, and draws a bootstrap outcome from \eqn{\bar{Q}_n(A^*, W^*)}: for
#' continuous outcomes by adding normal or empirical-residual error, and for
#' binary outcomes by a Bernoulli draw. Truncation is a property of the
#' estimator, not of the data-generating mechanism, so the bootstrap exposure is
#' always drawn from the untruncated propensity.
#'
#' Third, it refits both nuisance models on each bootstrap dataset and applies
#' the candidate estimator, bounding the fitted propensity into the truncation
#' interval where the estimator uses it. `ETA.Bias` is the mean of the bootstrap
#' estimates minus `truth`, and its Monte Carlo standard error is the standard
#' deviation of the bootstrap estimates divided by the square root of the number
#' of retained draws. A bootstrap draw whose estimate is non-finite, which can
#' happen when a resampled exposure lands in a single arm or a refit propensity
#' reaches exactly 0 or 1, is dropped with a warning before the summaries are
#' formed.
#'
#' The three estimators behave differently under positivity violations.
#' G-computation ignores the propensity entirely, so its ETA.Bias is Monte Carlo
#' noise around zero and is flat across a truncation grid. Inverse probability
#' weighting divides by the fitted propensity, so extreme scores inflate its
#' ETA.Bias, and tightening the truncation bound trades bias for variance. The
#' augmented (doubly robust) estimator stays near zero in bias while carrying
#' more variance than G-computation.
#'
#' When `truncation_grid` is supplied, every grid point reuses one shared set of
#' bootstrap draws and nuisance refits, so the sweep isolates the effect of the
#' truncation bound alone: the inverse-probability-weighting bias rises and its
#' bootstrap spread shrinks as the bound tightens.
#'
#' ETA.Bias captures only the positivity, truncation, and sparsity component of
#' bias and excludes model-misspecification bias by construction. It is a red
#' flag for a fitted mechanism, not a bias correction.
#'
#' @param .data A data frame.
#' @param .exposure The exposure column, selected with data-masking. A binary
#'   exposure has exactly two distinct values; a categorical exposure has as
#'   many levels as it has distinct values, and every non-reference level is
#'   contrasted against the reference; a continuous exposure has no levels, and
#'   is summarized by the working model `msm_formula` sets.
#' @param .outcome The outcome column, selected with data-masking.
#' @param .covariates The covariate columns, selected with tidyselect. Numeric,
#'   logical, factor, and character covariates are all supported. The fitted
#'   treatment and outcome models expand factor and character covariates into
#'   indicator terms, and the bootstrap resamples covariate rows, so factor
#'   levels are preserved throughout.
#' @param estimator The causal estimator to diagnose, one of `"ipw"` (inverse
#'   probability weighting, matching [propensity::ipw()]), `"gcomp"`
#'   (G-computation), or `"aipw"` (augmented, doubly robust). The augmented
#'   estimator's correction is a sum over the exposure levels, so it is an error
#'   for a continuous exposure.
#' @param exposure_formula An optional model formula for the treatment mechanism.
#'   `NULL` (the default) fits a main-effects logistic model of the exposure on
#'   the covariates.
#' @param outcome_formula An optional model formula for the outcome regression.
#'   `NULL` (the default) fits a main-effects model of the outcome on the
#'   exposure and covariates. Terms that recompute a data-dependent basis from
#'   the exposure column (for example `poly()` of the exposure) are not
#'   supported, because the counterfactual predictions rebuild the design with
#'   the exposure set to each of its levels in turn.
#' @param outcome_type One of `"auto"` (detect from the outcome, the default),
#'   `"continuous"`, or `"binary"`.
#' @param truncation An optional length-two numeric vector `c(lower, upper)`
#'   bounding the fitted propensity for a single run. `NULL` (the default)
#'   applies no truncation. A pair of bounds describes the one-dimensional
#'   simplex of a two-level exposure, so it applies to binary exposures only and
#'   is an error past two levels or on a continuous exposure.
#' @param truncation_grid An optional numeric vector, each entry in `[0, 1/k)`
#'   for an exposure with `k` levels, which is `[0, 0.5)` for a binary exposure,
#'   and in `[0, 1)` for a continuous one. All levels share one set of bootstrap
#'   draws, and the grid overrides `truncation`. What an entry means depends on
#'   the exposure type: at two levels it is a lower bound and becomes the swept
#'   truncation level `c(lower, 1 - lower)`; past two levels it raises every
#'   fitted probability below it and renormalizes each row, leaving no upper
#'   bound to apply or report; and for a continuous exposure it is a quantile
#'   level `g` that caps the stabilized weight at its `1 - g` quantile, read once
#'   from the observed fit and held fixed across the bootstrap draws.
#' @param n_boot The number of bootstrap datasets. Defaults to `500`. Must be at
#'   least 2, since the Monte Carlo standard error needs two draws.
#' @param error_dist The bootstrap error model for continuous outcomes, one of
#'   `"normal"` (add Gaussian noise matched to the residual standard deviation,
#'   the default) or `"empirical"` (resample residuals). Ignored for binary
#'   outcomes.
#' @param exposure_type One of `"auto"` (detect from the data, the default),
#'   `"binary"`, `"categorical"`, or `"continuous"`. A supplied type is
#'   authoritative and detection is not consulted, so it is rejected only when
#'   the exposure column cannot carry it: `"binary"` needs exactly two distinct
#'   values, `"continuous"` needs a numeric column, and `"categorical"` needs
#'   nothing beyond the two levels every discrete run needs. Within the discrete
#'   types the run is decided by the number of levels, so declaring
#'   `"categorical"` on a two-level column returns the binary run unchanged.
#' @param reference_level The exposure level every contrast is taken against.
#'   `NULL` (the default) takes the first level, which is the lower of the two
#'   values for a binary exposure and the first factor level otherwise. A
#'   continuous exposure has no levels to contrast, so supplying it there is an
#'   error.
#' @param msm_formula An optional one-sided formula for the marginal structural
#'   model a continuous exposure is summarized by, for example `~ a + I(a^2)`.
#'   `NULL` (the default) fits a model linear in the exposure, whose one
#'   coefficient is named after the exposure column. A discrete estimand is a set
#'   of contrasts and imposes no working model, so supplying it there is an
#'   error.
#'
#' @return An `eta_bias_result` object, an S7 subclass of
#'   [positivity_diagnostic]. Its `@results` tibble has one row per estimand term
#'   and truncation level, with the columns `term`, `truncation_lower`,
#'   `truncation_upper`, `bias` (ETA.Bias), `mc_se` (the Monte Carlo standard
#'   error), and `boot_mean` (the mean bootstrap estimate), in that order. `term`
#'   names the estimand contrast as `<level> - <reference level>`, so an exposure
#'   coded 0/1 has the single term `"1 - 0"`; for a continuous exposure it names
#'   a working-model coefficient instead. A continuous run has no lower cap to
#'   place on a density ratio, so `truncation_lower` is `NA_real_` and
#'   `truncation_upper` holds the weight cap, which is `Inf` where no cap was
#'   applied. The object also carries the properties `@estimator`, `@truth` (a
#'   named numeric vector holding one truth per term), and `@boot_estimates` (a
#'   list of bootstrap-estimate vectors, one per row of `@results`).
#'
#'   [generics::glance()] returns a one-row tibble with `n` (the sample size),
#'   `estimator`, and `n_boot`. A single-term estimand adds `truth`; an estimand
#'   of more than one term has no single truth, so it reports `n_terms` in its
#'   place. A one-row result then adds `bias` and `mc_se`, while a result of more
#'   rows replaces those two with `n_levels`, `bias_min`, and `bias_max`. The row
#'   count is what decides this, not the argument used, so a length-one
#'   `truncation_grid` takes the one-row form.
#'
#' @references
#' Petersen ML, Porter KE, Gruber S, Wang Y, van der Laan MJ (2012). Diagnosing
#' and responding to violations in the positivity assumption. *Statistical
#' Methods in Medical Research*, 21(1):31--54. \doi{10.1177/0962280210386207}
#'
#' @examples
#' set.seed(1)
#' n <- 300
#' x1 <- rnorm(n)
#' x2 <- rnorm(n)
#' a <- rbinom(n, 1, plogis(x1 + x2))
#' y <- a + x1 + x2 + rnorm(n)
#' df <- data.frame(a = a, y = y, x1 = x1, x2 = x2)
#'
#' # n_boot is small here to keep the example fast.
#' check_eta_bias(df, a, y, c(x1, x2), estimator = "ipw", n_boot = 25)
#'
#' @export
check_eta_bias <- function(
  .data,
  .exposure,
  .outcome,
  .covariates,
  estimator = c("ipw", "gcomp", "aipw"),
  exposure_formula = NULL,
  outcome_formula = NULL,
  outcome_type = c("auto", "continuous", "binary"),
  truncation = NULL,
  truncation_grid = NULL,
  n_boot = 500,
  error_dist = c("normal", "empirical"),
  exposure_type = c("auto", "binary", "categorical", "continuous"),
  reference_level = NULL,
  msm_formula = NULL
) {
  validate_data_frame(.data)
  estimator <- resolve_arg_match(rlang::arg_match(estimator))
  outcome_type <- resolve_arg_match(rlang::arg_match(outcome_type))
  error_dist <- resolve_arg_match(rlang::arg_match(error_dist))
  validate_count(n_boot, arg_name = "n_boot", min = 2)

  exposure_pos <- eval_select_column(
    rlang::enquo(.exposure),
    .data,
    ".exposure"
  )
  outcome_pos <- eval_select_column(rlang::enquo(.outcome), .data, ".outcome")
  covariate_pos <- eval_select_columns(
    rlang::enquo(.covariates),
    .data,
    ".covariates"
  )
  validate_column_selection(covariate_pos, ".covariates")

  exposure_name <- names(exposure_pos)
  outcome_name <- names(outcome_pos)
  covariate_names <- names(covariate_pos)

  # The treatment model would separate perfectly on the exposure and the outcome
  # model would absorb the whole treatment effect, so neither may appear among
  # the covariates.
  overlap <- intersect(covariate_pos, c(exposure_pos, outcome_pos))
  if (length(overlap) > 0) {
    overlapping <- names(.data)[overlap]
    abort(
      c(
        "{.arg .covariates} must not include the exposure or outcome column.",
        x = "{.val {overlapping}} {?is/are} also selected by {.arg .exposure} or {.arg .outcome}."
      ),
      error_class = "positively_selection_error"
    )
  }

  validate_eta_covariates(.data, covariate_names)
  validate_numeric_columns(.data, outcome_name, ".outcome")

  # The exposure is not required to be numeric: a factor, character, or logical
  # column is a valid discrete exposure, matching check_extrapolation(). Only
  # missingness and the resolved type are enforced here.
  exposure_vec <- .data[[exposure_pos]]
  if (anyNA(exposure_vec)) {
    abort(
      "{.arg .exposure} must not contain missing values.",
      error_class = "positively_missing_error"
    )
  }

  # The resolved type is checked against the column on both paths, which is what
  # rules out a one-level factor: detection reads it as binary, and eta_spec()
  # would otherwise code its single level to 0, leaving a treatment mechanism
  # with no second arm and a finite, unwarned, meaningless bias.
  exposure_type <- resolve_exposure_type(
    exposure_type,
    exposure_vec,
    supported = diagnostic_supported_types()[["eta_bias"]],
    fn = "check_eta_bias"
  )

  validate_eta_estimand_args(
    exposure_type,
    estimator = estimator,
    reference_level = reference_level,
    msm_formula = msm_formula
  )

  spec <- eta_spec(
    exposure_vec,
    exposure_type,
    reference_level = reference_level,
    msm_formula = msm_formula,
    exposure_name = exposure_name
  )

  # The multinomial treatment mechanism is the only part of the run that needs
  # nnet, and it is refit once per bootstrap draw, so availability is settled
  # here rather than inside the loop. Two levels take a logistic fit and never
  # reach it, and a continuous exposure has no level to model.
  needs_multinom <- spec$type != "continuous" && spec$k > 2
  if (needs_multinom && !rlang::is_installed("nnet")) {
    abort(
      c(
        "A categorical exposure of more than two levels requires the {.pkg nnet} package.",
        i = "Install {.pkg nnet} to run {.fn check_eta_bias} on this exposure."
      ),
      error_class = "positively_missing_package_error"
    )
  }

  # The truncation bounds are resolved after the exposure type, because what a
  # grid point bounds, and how large it may be, both depend on the type.
  levels <- resolve_truncation_levels(
    truncation,
    truncation_grid,
    spec,
    call = rlang::current_env()
  )

  a <- spec$a
  y <- as.double(.data[[outcome_pos]])
  covariates <- .data[covariate_pos]

  outcome_type <- resolve_outcome_type(outcome_type, y)
  if (outcome_type == "binary" && !all(y %in% c(0, 1))) {
    bad <- utils::head(unique(y[!y %in% c(0, 1)]), 5)
    abort(
      c(
        "{.arg .outcome} must contain only {.val {0}} and {.val {1}} when {.arg outcome_type} is {.val {\"binary\"}}.",
        x = "Found other values, including {.val {bad}}."
      ),
      error_class = "positively_range_error"
    )
  }

  formulas <- resolve_eta_formulas(
    exposure_formula,
    outcome_formula,
    exposure_name,
    outcome_name,
    covariate_names
  )

  run <- eta_bias_bootstrap(
    a = a,
    y = y,
    covariates = covariates,
    exposure_name = exposure_name,
    outcome_name = outcome_name,
    estimator = estimator,
    outcome_type = outcome_type,
    error_dist = error_dist,
    levels = levels,
    n_boot = n_boot,
    formulas = formulas,
    spec = spec,
    call = rlang::current_env()
  )

  # A continuous run resolves its weight cap from the observed fit, so the
  # bounds the results report come back from the bootstrap rather than from the
  # arguments alone.
  levels <- run$levels
  results <- tibble::tibble(
    term = run$term,
    truncation_lower = levels$lower[run$level],
    truncation_upper = levels$upper[run$level],
    bias = run$boot_mean - unname(run$truth[run$term]),
    mc_se = run$boot_sd / sqrt(run$n_valid),
    boot_mean = run$boot_mean
  )

  eta_bias_result(
    results = results,
    exposure = exposure_name,
    exposure_type = exposure_type,
    n = nrow(.data),
    params = list(
      estimator = estimator,
      outcome_type = outcome_type,
      error_dist = error_dist,
      n_boot = n_boot,
      truncation = truncation,
      truncation_grid = truncation_grid,
      reference_level = reference_level,
      msm_formula = msm_formula
    ),
    call = rlang::current_call(),
    estimator = estimator,
    truth = run$truth,
    boot_estimates = run$boot_estimates
  )
}

# ---- Argument resolution --------------------------------------------------

#' Reject the estimand arguments the resolved exposure type leaves nothing for
#'
#' A discrete estimand is a set of contrasts against a reference level; a
#' continuous one is the coefficient set of a working model. Each type therefore
#' has one argument that shapes the estimand and one that says nothing on it, and
#' the augmented estimator is defined for the discrete estimand alone. Each
#' unusable argument is refused rather than ignored, because silently dropping any
#' of them would report a number the caller did not ask for.
#'
#' @param exposure_type The resolved exposure type.
#' @param estimator One of `"ipw"`, `"gcomp"`, or `"aipw"`.
#' @param reference_level The supplied reference level, or `NULL`.
#' @param msm_formula The supplied working-model formula, or `NULL`.
#' @param call The calling environment, used to build the error's call.
#'
#' @return `NULL`, invisibly, when every supplied argument applies.
#' @keywords internal
#' @noRd
validate_eta_estimand_args <- function(
  exposure_type,
  estimator,
  reference_level,
  msm_formula,
  call = rlang::caller_env()
) {
  if (exposure_type != "continuous") {
    if (!is.null(msm_formula)) {
      abort(
        c(
          "{.arg msm_formula} sets the working model a continuous exposure is
           summarized by.",
          x = "{.arg .exposure} is {exposure_type}.",
          i = "A discrete estimand contrasts each level against
               {.arg reference_level}, which imposes no working model."
        ),
        error_class = "positively_args_error",
        call = call
      )
    }
    return(invisible(NULL))
  }

  if (estimator == "aipw") {
    abort(
      c(
        "{.arg estimator} {.val aipw} is not available for a continuous
         exposure.",
        i = "The augmentation term integrates the outcome residual over the
             exposure grid, which this path does not compute.",
        i = "Use {.val ipw} or {.val gcomp}."
      ),
      error_class = "positively_args_error",
      call = call
    )
  }
  if (!is.null(reference_level)) {
    abort(
      c(
        "{.arg reference_level} names the level every contrast is taken against,
         which a continuous exposure does not have.",
        i = "A continuous estimand is the coefficient set of the working model
             set by {.arg msm_formula}."
      ),
      error_class = "positively_args_error",
      call = call
    )
  }
  invisible(NULL)
}

#' Resolve the truncation levels for a run or a sweep
#'
#' Returns the per-level lower and upper probability bounds. A `truncation_grid`
#' takes precedence, a single `truncation` bound is honoured verbatim, and when
#' both are `NULL` the run applies no truncation.
#'
#' The upper bound belongs to two levels alone. There the simplex is
#' one-dimensional, so a lower bound on one level is an upper bound on the other
#' and a grid point expands to `c(lower, 1 - lower)`. Past two levels truncation
#' raises every level from below and renormalizes, which leaves no upper bound to
#' apply and none to report, so the upper bound is `NA_real_` at every level.
#'
#' A continuous exposure bounds a stabilized density ratio rather than a
#' probability, so a grid point is the quantile level at which that weight is
#' capped and the pair of probability bounds has nothing to bound. The cap itself
#' depends on the observed fit, so it is filled in by the bootstrap; only the
#' quantile levels are settled here.
#'
#' @param truncation A length-two bound or `NULL`.
#' @param truncation_grid A vector of lower bounds, or of quantile levels for a
#'   continuous exposure, or `NULL`.
#' @param spec The exposure specification from `eta_spec()`, whose type and level
#'   count decide what a grid point bounds and how large it may be.
#'
#' @return A list with numeric vectors `lower` and `upper`, one entry per level,
#'   and for a continuous exposure the `quantile` levels the caps come from.
#' @keywords internal
#' @noRd
resolve_truncation_levels <- function(
  truncation,
  truncation_grid,
  spec,
  call = rlang::caller_env()
) {
  if (spec$type == "continuous") {
    return(resolve_weight_truncation_levels(
      truncation,
      truncation_grid,
      call = call
    ))
  }
  k <- spec$k
  two_sided <- k == 2
  if (!is.null(truncation_grid)) {
    validate_truncation_grid(truncation_grid, k = k, call = call)
    upper <- if (two_sided) {
      1 - truncation_grid
    } else {
      rep(NA_real_, length(truncation_grid))
    }
    list(lower = truncation_grid, upper = upper)
  } else if (!is.null(truncation)) {
    if (!two_sided) {
      abort(
        c(
          "{.arg truncation} bounds a fitted probability from both sides, which
           describes a two-level exposure only.",
          x = "{.arg .exposure} has {k} levels.",
          i = "Use {.arg truncation_grid} to raise every level from below."
        ),
        error_class = "positively_args_error",
        call = call
      )
    }
    validate_truncation(truncation, call = call)
    list(lower = truncation[[1]], upper = truncation[[2]])
  } else {
    list(lower = 0, upper = if (two_sided) 1 else NA_real_)
  }
}

#' Resolve the quantile levels a continuous run caps its weights at
#'
#' The estimator's weight is the stabilized density ratio, so there is no lower
#' cap to place on it and the reported lower bound is missing at every level. The
#' upper bound is the cap, which the bootstrap resolves from the observed fit and
#' writes back over the placeholder left here.
#'
#' @inheritParams resolve_truncation_levels
#'
#' @return A list of `lower`, `upper`, and `quantile`, one entry per level.
#' @keywords internal
#' @noRd
resolve_weight_truncation_levels <- function(
  truncation,
  truncation_grid,
  call = rlang::caller_env()
) {
  if (!is.null(truncation)) {
    abort(
      c(
        "{.arg truncation} is a pair of probability bounds, and a continuous
         exposure weights by a density ratio rather than by a probability.",
        i = "Use {.arg truncation_grid} to cap the stabilized weight at a
             quantile of its own distribution."
      ),
      error_class = "positively_args_error",
      call = call
    )
  }
  quantile_levels <- if (is.null(truncation_grid)) {
    0
  } else {
    validate_truncation_grid(
      truncation_grid,
      exposure_type = "continuous",
      call = call
    )
  }
  n_levels <- length(quantile_levels)
  list(
    lower = rep(NA_real_, n_levels),
    upper = rep(NA_real_, n_levels),
    quantile = quantile_levels
  )
}

#' Validate the covariate columns for the ETA.Bias bootstrap
#'
#' The treatment and outcome models accept any covariate type that
#' [stats::glm()] can expand into model terms: numeric, logical, factor, and
#' character. Other column types, such as list or date columns, are rejected.
#' Missing values are rejected regardless of type, since the bootstrap resamples
#' complete covariate rows.
#'
#' @param .data The data frame.
#' @param columns A character vector of covariate column names.
#' @param call The calling environment, used to build the error's call.
#'
#' @return `.data`, invisibly, when every covariate is a supported type and free
#'   of missing values.
#' @keywords internal
#' @noRd
validate_eta_covariates <- function(
  .data,
  columns,
  call = rlang::caller_env()
) {
  is_supported <- function(column) {
    value <- .data[[column]]
    is.numeric(value) ||
      is.logical(value) ||
      is.factor(value) ||
      is.character(value)
  }
  unsupported <- columns[!vapply(columns, is_supported, logical(1))]
  if (length(unsupported) > 0) {
    abort(
      c(
        "{.arg .covariates} must select numeric, logical, factor, or character columns.",
        x = "{.val {unsupported}} {?is/are} of an unsupported type."
      ),
      error_class = "positively_type_error",
      call = call
    )
  }
  missing_columns <- columns[vapply(
    columns,
    function(column) anyNA(.data[[column]]),
    logical(1)
  )]
  if (length(missing_columns) > 0) {
    abort(
      c(
        "{.arg .covariates} must not contain missing values.",
        x = "Missing values in {.val {missing_columns}}."
      ),
      error_class = "positively_missing_error",
      call = call
    )
  }
  invisible(.data)
}

#' Resolve the outcome type, detecting it when requested
#'
#' @param outcome_type One of `"auto"`, `"continuous"`, or `"binary"`.
#' @param y The outcome vector.
#'
#' @return `"continuous"` or `"binary"`.
#' @keywords internal
#' @noRd
resolve_outcome_type <- function(outcome_type, y) {
  if (outcome_type != "auto") {
    return(outcome_type)
  }
  values <- unique(y)
  if (length(values) == 2 && all(values %in% c(0, 1))) {
    "binary"
  } else {
    "continuous"
  }
}

#' Resolve the treatment and outcome model formulas
#'
#' Builds the default main-effects formulas when the user supplies none, and
#' records whether both formulas are the defaults so the bootstrap can take its
#' fast matrix path.
#'
#' @param exposure_formula A user formula for the treatment model or `NULL`.
#' @param outcome_formula A user formula for the outcome model or `NULL`.
#' @param exposure_name The exposure column name.
#' @param outcome_name The outcome column name.
#' @param covariate_names The covariate column names.
#'
#' @return A list with the resolved `exposure` and `outcome` formulas and a
#'   logical `default` flag.
#' @keywords internal
#' @noRd
resolve_eta_formulas <- function(
  exposure_formula,
  outcome_formula,
  exposure_name,
  outcome_name,
  covariate_names
) {
  default <- is.null(exposure_formula) && is.null(outcome_formula)
  # Backticks keep non-syntactic names intact through parsing.
  backticked <- function(name) paste0("`", name, "`")
  covariate_terms <- paste(backticked(covariate_names), collapse = " + ")
  exposure_formula <- exposure_formula %||%
    stats::as.formula(paste(backticked(exposure_name), "~", covariate_terms))
  outcome_formula <- outcome_formula %||%
    stats::as.formula(paste(
      backticked(outcome_name),
      "~",
      backticked(exposure_name),
      "+",
      covariate_terms
    ))
  list(
    exposure = exposure_formula,
    outcome = outcome_formula,
    default = default
  )
}

# ---- Exposure specification -----------------------------------------------

#' Describe the exposure grid, coding, and estimand for one run
#'
#' The three exposure types differ only in the grid of exposure values at which
#' the outcome regression is evaluated and in the marginal structural model that
#' summarizes the resulting counterfactual means. Collecting both in one object
#' keeps the bootstrap engine and the estimators free of exposure-type branches.
#'
#' For a discrete exposure the marginal structural model is the saturated factor
#' model, so summarizing its counterfactual means is a contrast against the
#' reference level and imposes nothing. A binary exposure is the two-level case.
#' A continuous exposure has no levels to saturate, so the working model is a
#' genuine restriction and the grid is the exposure's own deciles.
#'
#' @param exposure_vec The exposure column, before coding.
#' @param exposure_type `"binary"`, `"categorical"`, or `"continuous"`.
#' @param reference_level The level contrasts are taken against, or `NULL` for
#'   the first level.
#' @param msm_formula The one-sided working-model formula for a continuous
#'   exposure, or `NULL` for a model linear in the exposure.
#' @param exposure_name The exposure column name, which names the variable the
#'   working model is written in.
#' @param call The calling environment, used to build the error's call.
#'
#' @return A list with the coded exposure `a`, the `grid` of exposure values, and
#'   the coefficient `terms`. A discrete specification adds the level `labels`,
#'   their count `k`, and the `reference` index; a continuous one adds the
#'   working-model `design()` closure, the `grid_design` it builds over the grid,
#'   and the `estimand` columns of that design.
#' @keywords internal
#' @noRd
eta_spec <- function(
  exposure_vec,
  exposure_type,
  reference_level = NULL,
  msm_formula = NULL,
  exposure_name = NULL,
  call = rlang::caller_env()
) {
  if (exposure_type == "continuous") {
    return(eta_spec_continuous(exposure_vec, msm_formula, exposure_name))
  }
  exposure_factor <- factor(exposure_vec)
  labels <- levels(exposure_factor)
  # A declared categorical type imposes no structure of its own, so a one-level
  # column reaches here whenever it was declared rather than detected. There is
  # no contrast to take against the only level there is.
  if (length(labels) < 2) {
    abort(
      c(
        "{.arg .exposure} must have at least two levels.",
        x = "It has {length(labels)} level{?s}."
      ),
      error_class = "positively_exposure_type_error",
      call = call
    )
  }
  reference <- resolve_reference_level(reference_level, labels, call = call)

  # Codes run 0 to K-1 in level order, so a binary exposure keeps the 0/1 coding
  # the logistic treatment model expects and the two-level case reduces exactly.
  list(
    type = exposure_type,
    a = as.double(as.integer(exposure_factor) - 1L),
    grid = as.double(seq_along(labels) - 1L),
    k = length(labels),
    labels = labels,
    reference = reference,
    terms = paste0(labels[-reference], " - ", labels[[reference]])
  )
}

#' Describe the grid and working model of a continuous exposure
#'
#' The grid is the exposure's own deciles, deduplicated, so equal weights over it
#' weight the working model by the marginal density of the exposure, the standard
#' choice. The estimand is the working model's coefficients other than the
#' intercept, which is what the counterfactual surface over that grid is projected
#' onto. `k`, `labels`, and `reference` are the vocabulary of levels, which a
#' continuous exposure has none of, so `k` is missing rather than one and any
#' branch that reads it without checking the type fails loudly.
#'
#' @inheritParams eta_spec
#'
#' @return The continuous exposure specification.
#' @keywords internal
#' @noRd
eta_spec_continuous <- function(exposure_vec, msm_formula, exposure_name) {
  a <- as.double(exposure_vec)
  grid <- unique(stats::quantile(
    a,
    probs = seq(0.1, 0.9, by = 0.1),
    names = FALSE,
    type = 7
  ))
  msm_formula <- msm_formula %||%
    stats::as.formula(paste0("~ `", exposure_name, "`"))

  design <- function(values) {
    frame <- stats::model.frame(
      msm_formula,
      data = stats::setNames(list(values), exposure_name)
    )
    stats::model.matrix(attr(frame, "terms"), frame)
  }
  grid_design <- design(grid)
  # A working model written without an intercept contributes every column to the
  # estimand, so the columns are chosen by name rather than by dropping the first.
  estimand <- which(colnames(grid_design) != "(Intercept)")

  list(
    type = "continuous",
    a = a,
    grid = grid,
    k = NA_integer_,
    labels = NULL,
    reference = NA_integer_,
    terms = colnames(grid_design)[estimand],
    design = design,
    grid_design = grid_design,
    estimand = estimand
  )
}

#' Resolve which level contrasts are taken against
#'
#' @param reference_level The requested level, or `NULL` for the first.
#' @param labels The exposure level labels, in order.
#' @param call The calling environment, used to build the error's call.
#'
#' @return The index of the reference level.
#' @keywords internal
#' @noRd
resolve_reference_level <- function(
  reference_level,
  labels,
  call = rlang::caller_env()
) {
  if (is.null(reference_level)) {
    return(1L)
  }
  if (length(reference_level) != 1) {
    abort(
      "{.arg reference_level} must be a single value.",
      error_class = "positively_args_error",
      call = call
    )
  }
  index <- match(as.character(reference_level), labels)
  if (is.na(index)) {
    abort(
      c(
        "{.arg reference_level} must be one of the exposure's levels.",
        x = "{.val {as.character(reference_level)}} is not among {.val {labels}}."
      ),
      error_class = "positively_args_error",
      call = call
    )
  }
  index
}

#' Code an exposure as level indicators
#'
#' The `n` by `K` matrix of indicators for the levels of `spec$grid`, in that
#' order. The treatment mechanism takes the whole matrix as its multinomial
#' response, and `eta_exposure_design()` takes every column but the first as the
#' saturated coding of the exposure in the outcome regression.
#'
#' Which column the design drops is a choice of parameterization and not of
#' estimand: any full-rank coding of a saturated factor gives the same fit and
#' the same counterfactual predictions, so `reference_level` bears on the
#' contrasts alone.
#'
#' @param a The exposure codes.
#' @param spec The exposure specification from `eta_spec()`.
#'
#' @return An `n` by `K` numeric matrix of zeros and ones.
#' @keywords internal
#' @noRd
eta_exposure_indicator <- function(a, spec) {
  indicators <- vapply(
    spec$grid,
    function(value) as.double(a == value),
    numeric(length(a))
  )
  matrix(indicators, nrow = length(a), ncol = spec$k)
}

#' Code an exposure for the outcome regression's design matrix
#'
#' Drops the first level's indicator, leaving the `K - 1` columns a saturated
#' factor contributes to a design that already carries an intercept. A two-level
#' exposure is coded 0 and 1, so its one remaining indicator is the exposure
#' itself and the design is the one the binary path has always built. A
#' continuous exposure enters the outcome regression as itself, one column.
#'
#' @inheritParams eta_exposure_indicator
#'
#' @return An `n` by `K - 1` numeric matrix, or an `n` by one matrix for a
#'   continuous exposure.
#' @keywords internal
#' @noRd
eta_exposure_design <- function(a, spec) {
  if (spec$type == "continuous") {
    return(matrix(a, ncol = 1L))
  }
  eta_exposure_indicator(a, spec)[, -1L, drop = FALSE]
}

#' The stabilized density ratio a continuous exposure weights by
#'
#' `f(A) / f(A | W)` under normal treatment mechanisms, mirroring
#' `propensity:::ate_continuous()` exactly: both densities carry the
#' maximum-likelihood variance rather than the `n - 1` estimate, the numerator's
#' marginal fit as much as the denominator's conditional one.
#'
#' @param a The exposure vector.
#' @param fitted The fitted conditional mean of the exposure.
#' @param sigma The pooled residual standard deviation of the treatment model.
#'
#' @return A numeric vector of stabilized weights.
#' @keywords internal
#' @noRd
eta_stabilized_weights <- function(a, fitted, sigma) {
  centre <- mean(a)
  stats::dnorm(a, centre, sqrt(mean((a - centre)^2))) /
    stats::dnorm(a, fitted, sigma)
}

#' Resolve the weight caps a continuous truncation sweep applies
#'
#' A quantile level `g` caps the stabilized weight at its `1 - g` quantile. The
#' cap is read once from the observed fit and then held fixed across bootstrap
#' draws, which is what leaves one number per level to report and what keeps a
#' sweep isolating the bound rather than mixing per-draw quantile noise into it.
#' The binary path's probability bounds are fixed across draws for the same
#' reason.
#'
#' A level of zero is no truncation, so its cap is infinite rather than the
#' observed maximum: the observed maximum is a no-op on the observed weights but
#' would still clip a bootstrap refit, whose heaviest weights run an order of
#' magnitude above it.
#'
#' @param weights The stabilized weights under the observed fit.
#' @param quantile_levels The quantile levels swept over.
#'
#' @return A numeric vector of caps, one per level.
#' @keywords internal
#' @noRd
eta_resolve_weight_caps <- function(weights, quantile_levels) {
  vapply(
    quantile_levels,
    function(level) {
      if (level == 0) {
        Inf
      } else {
        unname(stats::quantile(weights, 1 - level, type = 7))
      }
    },
    numeric(1)
  )
}

#' Bound a probability matrix away from zero
#'
#' With more than two levels this mirrors the matrix method of
#' [propensity::ps_trunc()] under `method = "ps"`: every entry below `lower` is
#' raised to `lower`, each row is renormalized to sum to one, and no upper bound
#' is applied or recorded. A lower pin on every level already bounds every weight
#' denominator, which is the failure mode truncation exists for.
#'
#' Two levels differ because the simplex is one-dimensional there and the binary
#' convention bounds both sides. Pinning the second column into `[lower, upper]`
#' and taking the first as its complement leaves the row summing to one, so there
#' is nothing to renormalize, and renormalizing would distort an asymmetric bound.
#'
#' @param gmat An `n` by `K` matrix of fitted probabilities.
#' @param lower The lower bound, applied at every number of levels.
#' @param upper The upper bound. Used only when `K` is two, and ignored
#'   otherwise, since a run over more levels carries no upper bound.
#'
#' @return The bounded matrix, whose rows sum to one.
#' @keywords internal
#' @noRd
truncate_probabilities <- function(gmat, lower, upper) {
  if (ncol(gmat) == 2) {
    bounded <- pmin(pmax(gmat[, 2], lower), upper)
    return(cbind(1 - bounded, bounded))
  }
  bounded <- pmax(gmat, lower)
  bounded / rowSums(bounded)
}

# ---- Bootstrap engine -----------------------------------------------------

#' Run the ETA.Bias parametric bootstrap
#'
#' Fits the nuisance models once on the observed data to define the target
#' `truth` and the data-generating mechanism, then draws `n_boot` bootstrap
#' datasets and refits both models on each. Every truncation level reuses the
#' one shared set of draws and refits, so a sweep isolates the truncation bound.
#'
#' @param a The coded exposure vector, taking values from `spec$grid`.
#' @param y The outcome vector.
#' @param covariates A data frame of covariates.
#' @param exposure_name,outcome_name The exposure and outcome column names.
#' @param estimator One of `"ipw"`, `"gcomp"`, or `"aipw"`.
#' @param outcome_type `"continuous"` or `"binary"`.
#' @param error_dist `"normal"` or `"empirical"`; used for continuous outcomes.
#' @param levels A list of `lower` and `upper` truncation vectors.
#' @param n_boot The number of bootstrap datasets.
#' @param formulas The resolved model formulas and default flag.
#' @param spec The exposure specification from `eta_spec()`.
#'
#' @return A list describing one cell per truncation level and estimand term.
#'   `truth` is a numeric vector named by term, one entry per term. `term` and
#'   `level` label the cells, with the level varying fastest and the term
#'   slowest. `boot_estimates` holds one estimate vector per cell, with the
#'   non-finite draws removed, and `boot_mean`, `boot_sd`, and `n_valid` hold
#'   that vector's mean, standard deviation, and length, one entry per cell.
#'   `levels` is the truncation specification the run applied, which for a
#'   continuous exposure carries the weight caps read off the observed fit.
#' @keywords internal
#' @noRd
eta_bias_bootstrap <- function(
  a,
  y,
  covariates,
  exposure_name,
  outcome_name,
  estimator,
  outcome_type,
  error_dist,
  levels,
  n_boot,
  formulas,
  spec,
  call = rlang::caller_env()
) {
  q_family <- if (outcome_type == "binary") {
    stats::binomial()
  } else {
    stats::gaussian()
  }
  fitters <- eta_fitters(
    covariates = covariates,
    exposure_name = exposure_name,
    outcome_name = outcome_name,
    q_family = q_family,
    formulas = formulas,
    spec = spec
  )
  continuous <- spec$type == "continuous"
  observed <- fitters$observed(a, y)
  # The truth is G-computation under the observed fit, which is exact because
  # that fit is the truth of the bootstrap data-generating mechanism.
  truth <- if (continuous) {
    eta_estimate_continuous(
      estimator = "gcomp",
      y = y,
      qmat = observed$qmat,
      design = NULL,
      weights = NULL,
      spec = spec
    )
  } else {
    eta_estimate(
      estimator = "gcomp",
      a = a,
      y = y,
      gmat = observed$gmat,
      qmat = observed$qmat,
      spec = spec
    )
  }
  # The cap is read from the observed fit's weights once, so every draw is
  # truncated by the same rule and the sweep isolates the bound.
  if (continuous) {
    levels$upper <- eta_resolve_weight_caps(
      eta_stabilized_weights(a, observed$ps_fitted, observed$ps_sigma),
      levels$quantile
    )
  }

  n <- length(a)
  n_levels <- length(levels$lower)
  n_terms <- length(spec$terms)
  # One row per bootstrap draw, then one slice per truncation level and term.
  estimates <- array(NA_real_, dim = c(n_boot, n_levels, n_terms))
  # All covariate-row resamples are drawn together, one column per draw.
  index_matrix <- matrix(
    sample.int(n, n * n_boot, replace = TRUE),
    nrow = n,
    ncol = n_boot
  )

  for (b in seq_len(n_boot)) {
    idx <- index_matrix[, b]
    # A discrete counterfactual mean is a lookup on the grid the outcome model
    # was already evaluated over. A continuous draw lands between grid points, so
    # its mean is evaluated at the drawn value itself.
    if (continuous) {
      a_star <- stats::rnorm(n, observed$ps_fitted[idx], observed$ps_sigma)
      mu <- observed$predict_at(idx, a_star)
    } else {
      a_star <- draw_bootstrap_exposure(
        observed$gmat[idx, , drop = FALSE],
        spec
      )
      mu <- observed$qmat[cbind(idx, match(a_star, spec$grid))]
    }
    y_star <- draw_bootstrap_outcome(
      mu = mu,
      outcome_type = outcome_type,
      error_dist = error_dist,
      residuals = observed$residuals,
      sigma = observed$sigma
    )
    refit <- fitters$refit(idx, a_star, y_star)
    if (continuous) {
      weights <- eta_stabilized_weights(
        a_star,
        refit$ps_fitted,
        refit$ps_sigma
      )
      design_star <- spec$design(a_star)
    }
    for (level in seq_len(n_levels)) {
      estimates[b, level, ] <- if (continuous) {
        eta_estimate_continuous(
          estimator = estimator,
          y = y_star,
          qmat = refit$qmat,
          design = design_star,
          weights = pmin(weights, levels$upper[[level]]),
          spec = spec
        )
      } else {
        eta_estimate(
          estimator = estimator,
          a = a_star,
          y = y_star,
          gmat = truncate_probabilities(
            refit$gmat,
            levels$lower[[level]],
            levels$upper[[level]]
          ),
          qmat = refit$qmat,
          spec = spec
        )
      }
    }
  }

  # A bootstrap exposure that lands in a single arm, or a refit propensity of
  # exactly 0 or 1, yields a non-finite estimate; such draws carry no
  # information and are dropped per truncation level before summarizing.
  # Terms vary slowest so a sweep reads as one block of levels per term.
  cells <- expand.grid(level = seq_len(n_levels), term = seq_len(n_terms))
  boot_estimates <- lapply(seq_len(nrow(cells)), function(cell) {
    column <- estimates[, cells$level[[cell]], cells$term[[cell]]]
    column[is.finite(column)]
  })
  n_valid <- lengths(boot_estimates)
  n_dropped <- max(n_boot - n_valid)
  if (n_dropped > 0) {
    warn(
      c(
        "Dropped {n_dropped} of {n_boot} bootstrap draws with a non-finite estimate.",
        i = "Such draws arise when a bootstrap exposure lands in a single arm or a refit propensity reaches exactly 0 or 1."
      ),
      warning_class = "positively_degenerate_boot_warning",
      call = call
    )
  }
  list(
    truth = truth,
    term = spec$terms[cells$term],
    level = cells$level,
    levels = levels,
    boot_mean = vapply(boot_estimates, mean, numeric(1)),
    boot_sd = vapply(boot_estimates, stats::sd, numeric(1)),
    boot_estimates = boot_estimates,
    n_valid = n_valid
  )
}

#' Draw a bootstrap exposure from the fitted treatment mechanism
#'
#' A two-level exposure keeps the Bernoulli draw, so a binary run consumes the
#' random stream exactly as it did before the grid was generalized. More levels
#' take an inverse-CDF draw from each row of the fitted probability matrix.
#'
#' @param gmat_rows The fitted probabilities at the resampled rows, `n` by `K`.
#' @param spec The exposure specification from `eta_spec()`.
#'
#' @return A numeric vector of exposure codes drawn from `spec$grid`.
#' @keywords internal
#' @noRd
draw_bootstrap_exposure <- function(gmat_rows, spec) {
  n <- nrow(gmat_rows)
  if (spec$k == 2) {
    return(as.double(stats::rbinom(n, 1L, gmat_rows[, 2])))
  }
  # Accumulating column by column performs the same additions in the same
  # left-to-right order as a per-row cumsum, so the result is bit-identical
  # while avoiding n R-level cumsum calls and a transpose.
  cumulative <- gmat_rows
  for (j in seq_len(spec$k)[-1L]) {
    cumulative[, j] <- cumulative[, j - 1L] + gmat_rows[, j]
  }
  draws <- rowSums(matrix(stats::runif(n), n, spec$k) > cumulative)
  # Rounding can leave the last cumulative probability a hair below one, which
  # would otherwise index one past the final level.
  spec$grid[pmin(draws, spec$k - 1L) + 1L]
}

#' Draw a bootstrap outcome from the fitted outcome mechanism
#'
#' @param mu The fitted outcome mean at the bootstrap exposure and covariates.
#' @param outcome_type `"continuous"` or `"binary"`.
#' @param error_dist `"normal"` or `"empirical"`.
#' @param residuals The observed-model residuals, resampled under `"empirical"`.
#' @param sigma The residual standard deviation, used under `"normal"`.
#'
#' @return A numeric vector of bootstrap outcomes.
#' @keywords internal
#' @noRd
draw_bootstrap_outcome <- function(
  mu,
  outcome_type,
  error_dist,
  residuals,
  sigma
) {
  n <- length(mu)
  if (outcome_type == "binary") {
    return(as.double(stats::rbinom(n, 1L, mu)))
  }
  error <- if (error_dist == "empirical") {
    sample(residuals, size = n, replace = TRUE)
  } else {
    stats::rnorm(n, mean = 0, sd = sigma)
  }
  mu + error
}

#' Build the observed-fit and bootstrap-refit closures
#'
#' Selects the fast matrix path when both formulas are the main-effects defaults
#' and every covariate is numeric, and the general formula path otherwise. A
#' factor, character, or logical covariate, or a user-supplied formula, takes the
#' formula path, where [stats::glm()] expands non-numeric covariates into model
#' terms. Both closures return `gmat`, the `n` by `K` matrix of fitted exposure
#' probabilities, and `qmat`, the `n` by `K` matrix of outcome predictions over
#' `spec$grid`, which are the only quantities the estimators consume. `observed()`
#' adds the residuals and the residual standard deviation, which the bootstrap
#' outcome draw needs.
#'
#' @inheritParams eta_bias_bootstrap
#'
#' @return A list of two functions: `observed(a, y)` and
#'   `refit(idx, a_star, y_star)`.
#' @keywords internal
#' @noRd
eta_fitters <- function(
  covariates,
  exposure_name,
  outcome_name,
  q_family,
  formulas,
  spec
) {
  numeric_covariates <- all(vapply(covariates, is.numeric, logical(1)))
  if (formulas$default && numeric_covariates) {
    eta_fitters_matrix(as.matrix(covariates), q_family, spec)
  } else {
    eta_fitters_formula(
      covariates = covariates,
      exposure_name = exposure_name,
      outcome_name = outcome_name,
      q_family = q_family,
      formulas = formulas,
      spec = spec
    )
  }
}

#' Fit a multinomial logit treatment mechanism
#'
#' The treatment mechanism past two levels, wrapped so that it returns the same
#' `n` by `K` matrix of fitted probabilities the two-level path returns, in
#' `spec$grid` order. [nnet::multinom()] writes an iteration trace to the console
#' unless `trace` is off, and this runs once per bootstrap draw.
#'
#' The response is an indicator matrix rather than a factor because
#' [nnet::multinom()] drops a level no row holds, which a bootstrap resample of a
#' sparse arm can produce, and returns one column too few for the caller to read.
#' An indicator response keeps every column, and a level no row holds takes a
#' fitted probability near zero, which is the reading that draw has earned.
#'
#' @param formula The treatment-model formula, whose response is the indicator
#'   matrix of exposure levels.
#' @param frame A list or data frame holding the response and the predictors.
#' @param n The number of rows.
#' @param k The number of exposure levels.
#'
#' @return An `n` by `k` matrix of fitted probabilities.
#' @keywords internal
#' @noRd
eta_multinom_probs <- function(formula, frame, n, k) {
  fit <- nnet::multinom(formula, data = frame, trace = FALSE)
  fitted <- stats::fitted(fit)
  # Reshaping to the expected dimensions would recycle silently if a future
  # nnet returned a different shape, and a recycled probability matrix is a
  # plausible-looking wrong answer rather than a failure.
  if (!is.matrix(fitted) || nrow(fitted) != n || ncol(fitted) != k) {
    shape <- paste(dim(as.matrix(fitted)), collapse = " by ")
    abort(
      c(
        "The fitted multinomial probabilities must be {n} by {k}.",
        x = "{.fn nnet::multinom} returned {shape}."
      ),
      error_class = "positively_size_error"
    )
  }
  matrix(fitted, nrow = n, ncol = k)
}

#' Fast matrix fitters for the default main-effects models
#'
#' Fits the treatment model and the outcome model on prebuilt design matrices,
#' avoiding the per-draw cost of building a data frame and parsing a formula.
#' This is the path a default run with numeric covariates takes. The outcome
#' model is a [stats::glm.fit()] throughout; the treatment model is a two-level
#' [stats::glm.fit()], the multinomial logit of `eta_multinom_probs()`, or a
#' linear model with a pooled residual standard deviation.
#'
#' @param x_cov The observed covariate matrix, `n` by `q`.
#' @param q_family The outcome-model family.
#' @param spec The exposure specification from `eta_spec()`, whose `grid` fixes
#'   the exposure values the outcome predictions are made at.
#'
#' @return A list of `observed(a, y)` and `refit(idx, a_star, y_star)`.
#' @keywords internal
#' @noRd
eta_fitters_matrix <- function(x_cov, q_family, spec) {
  n <- nrow(x_cov)
  intercept <- rep.int(1, n)
  binomial_family <- stats::binomial()
  gaussian_family <- stats::gaussian()

  # glm.fit warns on separation and near-degenerate fits, which are expected
  # under positivity violations, the very condition the diagnostic measures.
  # glm.fit returns NA coefficients for aliased (constant or collinear) columns.
  # Zeroing them reproduces the aliased-column predictions that stats::glm()
  # makes on the formula path, so the two paths agree on rank-deficient designs.
  fit_ps <- if (spec$type == "continuous") {
    function(design, exposure) {
      coef <- suppressWarnings(
        stats::glm.fit(design, exposure, family = gaussian_family)$coefficients
      )
      coef[is.na(coef)] <- 0
      fitted <- drop(design %*% coef)
      eta_normal_mechanism(exposure, fitted)
    }
  } else if (spec$k == 2) {
    function(design, exposure) {
      coef <- suppressWarnings(
        stats::glm.fit(design, exposure, family = binomial_family)$coefficients
      )
      coef[is.na(coef)] <- 0
      positive <- binomial_family$linkinv(drop(design %*% coef))
      list(gmat = cbind(1 - positive, positive))
    }
  } else {
    function(design, exposure) {
      # The response and the covariates enter as whole matrices, so the model
      # frame holds two variables and no covariate name can collide with them.
      frame <- list(.a = eta_exposure_indicator(exposure, spec), .x = design)
      list(
        gmat = eta_multinom_probs(.a ~ .x - 1, frame, nrow(design), spec$k)
      )
    }
  }
  fit_q <- function(x_cov_rows, exposure, response) {
    design <- cbind(intercept, eta_exposure_design(exposure, spec), x_cov_rows)
    fit <- suppressWarnings(
      stats::glm.fit(design, response, family = q_family)
    )
    coef <- fit$coefficients
    coef[is.na(coef)] <- 0
    counterfactual <- vapply(
      spec$grid,
      function(value) {
        at_value <- cbind(
          intercept,
          eta_exposure_design(rep_len(value, n), spec),
          x_cov_rows
        )
        q_family$linkinv(drop(at_value %*% coef))
      },
      numeric(nrow(x_cov_rows))
    )
    list(coef = coef, rank = fit$rank, qmat = counterfactual)
  }

  observed <- function(a, y) {
    fit <- fit_ps(cbind(intercept, x_cov), a)
    q <- fit_q(x_cov, a, y)
    design <- cbind(intercept, eta_exposure_design(a, spec), x_cov)
    fitted_y <- q_family$linkinv(drop(design %*% q$coef))
    residuals <- y - fitted_y
    # Residual degrees of freedom count estimated parameters only. An aliased
    # column contributes no estimated parameter, so it does not deflate sigma,
    # matching the residual standard deviation on the formula path.
    df_residual <- n - q$rank
    fit$qmat <- q$qmat
    fit$residuals <- residuals
    fit$sigma <- sqrt(sum(residuals^2) / df_residual)
    fit$predict_at <- function(idx, values) {
      at_values <- cbind(
        intercept,
        eta_exposure_design(values, spec),
        x_cov[idx, , drop = FALSE]
      )
      q_family$linkinv(drop(at_values %*% q$coef))
    }
    fit
  }

  refit <- function(idx, a_star, y_star) {
    x_rows <- x_cov[idx, , drop = FALSE]
    fit <- fit_ps(cbind(intercept, x_rows), a_star)
    fit$qmat <- fit_q(x_rows, a_star, y_star)$qmat
    fit
  }

  list(observed = observed, refit = refit)
}

#' General formula fitters for user-supplied models
#'
#' Fits the treatment and outcome models from arbitrary formulas by building the
#' design matrix with [stats::glm.fit()], rebuilding the model frame on each
#' bootstrap dataset. Fitting the design directly, rather than through
#' [stats::glm()], avoids dropping unused factor levels, so a resample missing a
#' rare level keeps its design column instead of aborting. [nnet::multinom()],
#' which the treatment model takes past two levels, leaves them undropped
#' already. Used only when the user overrides a model formula or supplies a
#' non-numeric covariate.
#'
#' @inheritParams eta_bias_bootstrap
#'
#' @return A list of `observed(a, y)` and `refit(idx, a_star, y_star)`.
#' @keywords internal
#' @noRd
eta_fitters_formula <- function(
  covariates,
  exposure_name,
  outcome_name,
  q_family,
  formulas,
  spec
) {
  binomial_family <- stats::binomial()
  gaussian_family <- stats::gaussian()
  # A bootstrap resample can drop every row that holds a rare character level,
  # re-factoring the covariate to a single level and aborting the contrasts.
  # Fixing character covariates to factors once on the observed data keeps their
  # levels constant across resamples, so the design columns never change.
  covariates[] <- lapply(covariates, function(column) {
    if (is.character(column)) factor(column) else column
  })

  # Past two levels the outcome model needs the exposure expanded into level
  # indicators, which stats::model.matrix() does for a factor and not for the
  # numeric codes. Two levels keep the numeric coding, whose one indicator is
  # the code itself, and a continuous exposure enters as itself.
  exposure_column <- function(values) {
    if (spec$type == "continuous" || spec$k == 2) {
      values
    } else {
      factor(values, levels = spec$grid)
    }
  }

  assemble <- function(exposure, response) {
    frame <- covariates
    frame[[exposure_name]] <- exposure_column(exposure)
    frame[[outcome_name]] <- response
    frame
  }

  # stats::glm() forces drop.unused.levels = TRUE in its model.frame call, so a
  # resample missing a factor level collapses that factor to one level and
  # `contrasts<-` aborts. Building the model frame directly leaves
  # drop.unused.levels at its FALSE default, retaining an absent level as an
  # all-zero design column whose aliased coefficient is zeroed, matching the
  # matrix path on rank-deficient designs.
  fit_glm <- function(formula, frame, family) {
    frame <- stats::model.frame(formula, data = frame)
    design <- stats::model.matrix(attr(frame, "terms"), frame)
    fit <- suppressWarnings(stats::glm.fit(
      design,
      stats::model.response(frame),
      family = family,
      offset = stats::model.offset(frame)
    ))
    coef <- fit$coefficients
    coef[is.na(coef)] <- 0
    list(
      formula = formula,
      family = family,
      coef = coef,
      df_residual = fit$df.residual
    )
  }
  predict_fit <- function(fit, frame) {
    frame <- stats::model.frame(fit$formula, data = frame)
    design <- stats::model.matrix(attr(frame, "terms"), frame)
    eta <- drop(design %*% fit$coef)
    offset <- stats::model.offset(frame)
    if (!is.null(offset)) {
      eta <- eta + offset
    }
    as.double(fit$family$linkinv(eta))
  }

  # The treatment model's response is the exposure itself, so the multinomial
  # fit reads the level indicators from a frame of its own rather than the shared
  # one, whose exposure column is coded for the outcome model.
  fit_ps <- if (spec$type == "continuous") {
    function(frame, exposure) {
      ps_fit <- fit_glm(formulas$exposure, frame, gaussian_family)
      eta_normal_mechanism(exposure, predict_fit(ps_fit, frame))
    }
  } else if (spec$k == 2) {
    function(frame, exposure) {
      ps_fit <- fit_glm(formulas$exposure, frame, binomial_family)
      positive <- predict_fit(ps_fit, frame)
      list(gmat = cbind(1 - positive, positive))
    }
  } else {
    function(frame, exposure) {
      ps_frame <- as.list(frame)
      ps_frame[[exposure_name]] <- eta_exposure_indicator(exposure, spec)
      list(
        gmat = eta_multinom_probs(
          formulas$exposure,
          ps_frame,
          nrow(frame),
          spec$k
        )
      )
    }
  }

  fit_pair <- function(frame, exposure) {
    fit <- fit_ps(frame, exposure)
    q_fit <- fit_glm(formulas$outcome, frame, q_family)
    fit$qmat <- vapply(
      spec$grid,
      function(value) {
        counterfactual_frame <- frame
        counterfactual_frame[[exposure_name]] <- exposure_column(
          rep_len(value, nrow(frame))
        )
        predict_fit(q_fit, counterfactual_frame)
      },
      numeric(nrow(frame))
    )
    fit$q_fit <- q_fit
    fit
  }

  observed <- function(a, y) {
    frame <- assemble(a, y)
    fit <- fit_pair(frame, a)
    # Bound outside `fit` so the prediction closure keeps the outcome fit after
    # it is dropped from the returned list.
    q_fit <- fit$q_fit
    fit$q_fit <- NULL
    residuals <- y - predict_fit(q_fit, frame)
    fit$residuals <- residuals
    fit$sigma <- sqrt(sum(residuals^2) / q_fit$df_residual)
    fit$predict_at <- function(idx, values) {
      at_values <- assemble(values, y)
      at_values[seq_along(covariates)] <- covariates[idx, , drop = FALSE]
      predict_fit(q_fit, at_values)
    }
    fit
  }

  refit <- function(idx, a_star, y_star) {
    frame <- assemble(a_star, y_star)
    frame[seq_along(covariates)] <- covariates[idx, , drop = FALSE]
    fit <- fit_pair(frame, a_star)
    fit$q_fit <- NULL
    fit
  }

  list(observed = observed, refit = refit)
}

#' Summarize a normal treatment mechanism for a continuous exposure
#'
#' The mechanism is the fitted conditional mean paired with the pooled
#' maximum-likelihood residual standard deviation, which is the spread
#' `propensity` derives for itself when it is handed no sigma.
#'
#' @param exposure The exposure vector.
#' @param fitted The fitted conditional mean.
#'
#' @return A list of `ps_fitted` and `ps_sigma`.
#' @keywords internal
#' @noRd
eta_normal_mechanism <- function(exposure, fitted) {
  list(ps_fitted = fitted, ps_sigma = sqrt(mean((exposure - fitted)^2)))
}

#' Apply a causal estimator to a bootstrap dataset
#'
#' Returns one estimate per coefficient of the marginal structural model, which
#' for a discrete exposure is the saturated factor model and so yields one
#' contrast per non-reference level. G-computation ignores the fitted
#' probabilities; inverse probability weighting is the Hajek difference in
#' weighted means, matching [propensity::ipw()]; the augmented estimator adds the
#' doubly robust correction. A binary exposure is the two-level case of each.
#'
#' @param estimator One of `"ipw"`, `"gcomp"`, or `"aipw"`.
#' @param a The bootstrap exposure, coded as in `spec$grid`.
#' @param y The bootstrap outcome.
#' @param gmat The truncated fitted probabilities, `n` by `K`.
#' @param qmat The outcome predictions at each grid value, `n` by `K`.
#' @param spec The exposure specification from `eta_spec()`.
#'
#' @return A named numeric vector, one entry per coefficient.
#' @keywords internal
#' @noRd
eta_estimate <- function(estimator, a, y, gmat, qmat, spec) {
  reference <- spec$reference
  others <- seq_len(spec$k)[-reference]

  if (estimator == "ipw") {
    hajek <- function(level) {
      weight <- (a == spec$grid[[level]]) / gmat[, level]
      sum(weight * y) / sum(weight)
    }
    mean_reference <- hajek(reference)
    estimates <- vapply(
      others,
      function(level) hajek(level) - mean_reference,
      numeric(1)
    )
  } else if (estimator == "gcomp") {
    estimates <- vapply(
      others,
      function(level) mean(qmat[, level] - qmat[, reference]),
      numeric(1)
    )
  } else {
    residual <- function(level) {
      (a == spec$grid[[level]]) / gmat[, level] * (y - qmat[, level])
    }
    residual_reference <- residual(reference)
    estimates <- vapply(
      others,
      function(level) {
        mean(
          qmat[, level] -
            qmat[, reference] +
            residual(level) -
            residual_reference
        )
      },
      numeric(1)
    )
  }

  names(estimates) <- spec$terms
  estimates
}

#' Apply a causal estimator to a bootstrap dataset with a continuous exposure
#'
#' A continuous estimand is the working model's coefficients rather than a set of
#' contrasts, so both estimators are a least-squares projection onto that model.
#' G-computation projects the counterfactual surface, the mean outcome prediction
#' at each grid value, with equal weight on every grid point; because the grid is
#' a quantile grid, that weights by the marginal density of the exposure.
#' Inverse probability weighting projects the observed outcomes at the exposure
#' values actually drawn, weighted by the truncated stabilized density ratio.
#'
#' @param estimator One of `"ipw"` or `"gcomp"`.
#' @param y The bootstrap outcome.
#' @param qmat The outcome predictions at each grid value, `n` by the grid size.
#' @param design The working-model design at the drawn exposure, used by inverse
#'   probability weighting alone.
#' @param weights The truncated stabilized weights, used by inverse probability
#'   weighting alone.
#' @param spec The exposure specification from `eta_spec()`.
#'
#' @return A named numeric vector, one entry per coefficient.
#' @keywords internal
#' @noRd
eta_estimate_continuous <- function(estimator, y, qmat, design, weights, spec) {
  coefficients <- if (estimator == "gcomp") {
    eta_msm_coefficients(spec$grid_design, colMeans(qmat))
  } else {
    eta_msm_coefficients(design, y, weights)
  }
  estimates <- coefficients[spec$estimand]
  names(estimates) <- spec$terms
  estimates
}

#' Fit a working model by weighted least squares
#'
#' @param design The working-model design matrix.
#' @param response The values the model is fitted to.
#' @param weights The observation weights. Defaults to equal weights.
#'
#' @return The coefficient vector, including the intercept where the design
#'   carries one.
#' @keywords internal
#' @noRd
eta_msm_coefficients <- function(design, response, weights = 1) {
  weighted <- design * weights
  drop(solve(crossprod(design, weighted), crossprod(weighted, response)))
}

# ---- Methods --------------------------------------------------------------

#' Count the truncation levels a result was read over
#'
#' The results tibble holds one row per term and truncation level, so the level
#' count is the row count divided by the term count. Counting the distinct lower
#' bounds instead would read one level for a continuous run, whose lower bound is
#' missing at every level because there is no lower cap on a density ratio.
#'
#' @param results The results tibble.
#' @param n_terms The number of estimand terms.
#'
#' @return A single integer.
#' @keywords internal
#' @noRd
eta_bias_n_levels <- function(results, n_terms) {
  nrow(results) %/% n_terms
}

method(glance, eta_bias_result) <- function(x, ...) {
  results <- x@results
  n_terms <- length(x@truth)
  columns <- list(
    n = x@n,
    estimator = x@estimator,
    # @params holds n_boot as supplied, which is a double.
    n_boot = as.integer(x@params$n_boot)
  )
  # An estimand of more than one term has one truth per term, so no single truth
  # describes the run and the term count stands in its place.
  if (n_terms == 1) {
    columns$truth <- unname(x@truth)
  } else {
    columns$n_terms <- n_terms
  }
  if (nrow(results) == 1) {
    columns$bias <- results$bias
    columns$mc_se <- results$mc_se
  } else {
    # A sweep has one bias per truncation level, so no single bias and no single
    # Monte Carlo error beside it describe the run.
    columns$n_levels <- eta_bias_n_levels(results, n_terms)
    columns$bias_min <- min(results$bias)
    columns$bias_max <- max(results$bias)
  }
  tibble::as_tibble(columns)
}

# The heading names the diagnostic in prose; the literature's `ETA.Bias` spelling
# stays on the line that reports the statistic itself.
method(diagnostic_label, eta_bias_result) <- function(x) {
  "ETA bias"
}

method(diagnostic_headline, eta_bias_result) <- function(x) {
  results <- x@results
  n_terms <- length(x@truth)
  n_levels <- eta_bias_n_levels(results, n_terms)
  estimand <- if (n_terms == 1) {
    cli::format_inline("against a truth of {round(unname(x@truth), 3)}")
  } else {
    cli::format_inline("across {n_terms} estimand term{?s}")
  }
  span <- round(range(results$bias), 3)
  # The bias belongs to an estimator under a truncation rule rather than to the
  # data, so the estimator is named alongside it.
  if (nrow(results) == 1) {
    cli::format_inline(
      "{x@estimator} estimator {estimand}; ETA.Bias {round(results$bias, 3)} (MC SE {round(results$mc_se, 3)}) from {x@params$n_boot} bootstrap draw{?s}"
    )
  } else if (n_levels == 1) {
    # More than one row over one truncation level is a range across terms, and
    # naming the one level it was read at would report a sweep that never ran.
    cli::format_inline(
      "{x@estimator} estimator {estimand}; ETA.Bias {span[[1]]} to {span[[2]]}"
    )
  } else {
    cli::format_inline(
      "{x@estimator} estimator {estimand}; ETA.Bias {span[[1]]} to {span[[2]]} over {n_levels} truncation level{?s}"
    )
  }
}

method(print, eta_bias_result) <- function(x, ...) {
  results <- x@results
  n_terms <- length(x@truth)
  n_levels <- eta_bias_n_levels(results, n_terms)

  # A sweep reports the span of the bias over its levels; a single level reports
  # the one bias it read, beside the Monte Carlo error of that reading.
  reading <- function(rows) {
    if (n_levels == 1) {
      cli::format_inline(
        "{round(rows$bias, 3)} (MC SE {round(rows$mc_se, 3)})"
      )
    } else {
      cli::format_inline(
        "{round(min(rows$bias), 3)} to {round(max(rows$bias), 3)}"
      )
    }
  }

  cat_cli({
    cli::cli_h1("{diagnostic_label(x)}")
    cli::cli_text("Exposure: {.val {x@exposure}} ({x@exposure_type})")
    cli::cli_text("Observations: {x@n}")
    cli::cli_text("Estimator: {x@estimator}")
    cli::cli_text("Bootstrap draws: {x@params$n_boot}")
    if (n_terms == 1) {
      cli::cli_text("Truth: {round(unname(x@truth), 3)}")
    } else {
      cli::cli_text("Estimand terms: {n_terms}")
    }
    if (n_levels > 1) {
      cli::cli_text("Truncation levels: {n_levels}")
    }
    # Each term is its own estimand with its own truth, so a run of several
    # states each one rather than a span over quantities that are not comparable.
    if (n_terms == 1) {
      cli::cli_text("ETA.Bias: {reading(results)}")
    } else {
      for (term in unique(results$term)) {
        term_reading <- reading(results[results$term == term, ])
        cli::cli_text("ETA.Bias ({term}): {term_reading}")
      }
    }
  })
  invisible(x)
}
