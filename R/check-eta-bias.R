# S7 methods defined here are wired up at load time by S7::methods_register().

#' The ETA.Bias diagnostic result class
#'
#' `eta_bias_result` is the S7 class returned by [check_eta_bias()]. It extends
#' [positivity_diagnostic] with the estimator label, the single scalar `truth`
#' against which the bootstrap estimates are compared, and the list of
#' bootstrap-estimate vectors (one per truncation level). It is created
#' internally and is not constructed directly by users.
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
    if (length(self@estimator) != 1) {
      "@estimator must be a single value"
    } else if (length(self@truth) != 1) {
      "@truth must be a single value"
    } else if (length(self@boot_estimates) != nrow(self@results)) {
      "@boot_estimates must have one vector per row of @results"
    }
  }
)

#' Diagnose positivity bias with the parametric bootstrap (ETA.Bias)
#'
#' `check_eta_bias()` implements the ETA.Bias parametric-bootstrap diagnostic of
#' Petersen et al. (2012) for binary point exposures. It estimates the bias that
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
#' models as indicator terms, exactly as [stats::glm()] expands them. The target
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
#' @param .exposure The binary exposure column, selected with data-masking.
#'   `check_eta_bias()` aborts unless the exposure has exactly two distinct
#'   values, whether that type is detected or declared through `exposure_type`.
#' @param .outcome The outcome column, selected with data-masking.
#' @param .covariates The covariate columns, selected with tidyselect. Numeric,
#'   logical, factor, and character covariates are all supported. The fitted
#'   treatment and outcome models expand factor and character covariates into
#'   indicator terms, and the bootstrap resamples covariate rows, so factor
#'   levels are preserved throughout.
#' @param estimator The causal estimator to diagnose, one of `"ipw"` (inverse
#'   probability weighting, matching [propensity::ipw()]), `"gcomp"`
#'   (G-computation), or `"aipw"` (augmented, doubly robust).
#' @param exposure_formula An optional model formula for the treatment mechanism.
#'   `NULL` (the default) fits a main-effects logistic model of the exposure on
#'   the covariates.
#' @param outcome_formula An optional model formula for the outcome regression.
#'   `NULL` (the default) fits a main-effects model of the outcome on the
#'   exposure and covariates. Terms that recompute a data-dependent basis from
#'   the exposure column (for example `poly()` of the exposure) are not
#'   supported, because the counterfactual predictions rebuild the design with
#'   the exposure set to 0 or 1.
#' @param outcome_type One of `"auto"` (detect from the outcome, the default),
#'   `"continuous"`, or `"binary"`.
#' @param truncation An optional length-two numeric vector `c(lower, upper)`
#'   bounding the fitted propensity for a single run. `NULL` (the default)
#'   applies no truncation.
#' @param truncation_grid An optional numeric vector of lower bounds, each in
#'   `[0, 0.5)`. Each becomes a swept truncation level `c(lower, 1 - lower)`, and
#'   all levels share one set of bootstrap draws. Overrides `truncation`.
#' @param n_boot The number of bootstrap datasets. Defaults to `500`. Must be at
#'   least 2, since the Monte Carlo standard error needs two draws.
#' @param error_dist The bootstrap error model for continuous outcomes, one of
#'   `"normal"` (add Gaussian noise matched to the residual standard deviation,
#'   the default) or `"empirical"` (resample residuals). Ignored for binary
#'   outcomes.
#' @param exposure_type One of `"auto"` (detect from the data, the default) or
#'   `"binary"`. A supplied type is authoritative and detection is not consulted,
#'   so `"binary"` is rejected only when the exposure column does not have
#'   exactly two distinct values. Declaring it unlocks no call that `"auto"`
#'   would refuse, since detection reads any two-valued column as binary; it is
#'   accepted so that every diagnostic takes the same argument.
#'
#' @return An `eta_bias_result` object, an S7 subclass of
#'   [positivity_diagnostic]. Its `@results` tibble has one row per truncation
#'   level with columns `truncation_lower`, `truncation_upper`, `bias`
#'   (ETA.Bias), `mc_se` (the Monte Carlo standard error), and `boot_mean` (the
#'   mean bootstrap estimate). It also carries the properties `@estimator`,
#'   `@truth`, and `@boot_estimates` (a list of bootstrap-estimate vectors, one
#'   per row of `@results`).
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
  exposure_type = c("auto", "binary")
) {
  validate_data_frame(.data)
  estimator <- resolve_arg_match(rlang::arg_match(estimator))
  outcome_type <- resolve_arg_match(rlang::arg_match(outcome_type))
  error_dist <- resolve_arg_match(rlang::arg_match(error_dist))
  validate_count(n_boot, arg_name = "n_boot", min = 2)

  levels <- resolve_truncation_levels(
    truncation,
    truncation_grid,
    call = rlang::current_env()
  )

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

  # The exposure is not required to be numeric: a two-level factor, character,
  # or logical column is a valid binary exposure, matching check_extrapolation().
  # Only missingness and a binary type are enforced here.
  exposure_vec <- .data[[exposure_pos]]
  if (anyNA(exposure_vec)) {
    abort(
      "{.arg .exposure} must not contain missing values.",
      error_class = "positively_missing_error"
    )
  }

  # exposure_type grants no new capability here: detect_exposure_type() reads any
  # two-valued column as binary, so declaring "binary" can never unlock a call
  # that "auto" rejects. It is accepted so that every diagnostic takes the same
  # argument. Resolving through resolve_exposure_type() is not inert, though: the
  # resolved type is checked against the column on both paths, and a one-level
  # factor, which detection also calls binary, would otherwise pass through
  # binary_to_01() as all zeros and yield a finite, unwarned, meaningless bias.
  exposure_type <- resolve_exposure_type(
    exposure_type,
    exposure_vec,
    supported = diagnostic_supported_types()[["eta_bias"]],
    fn = "check_eta_bias"
  )

  a <- binary_to_01(exposure_vec)
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
    call = rlang::current_env()
  )

  results <- tibble::tibble(
    truncation_lower = levels$lower,
    truncation_upper = levels$upper,
    bias = run$boot_mean - run$truth,
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
      truncation_grid = truncation_grid
    ),
    call = rlang::current_call(),
    estimator = estimator,
    truth = run$truth,
    boot_estimates = run$boot_estimates
  )
}

# ---- Argument resolution --------------------------------------------------

#' Resolve the truncation levels for a run or a sweep
#'
#' Returns the per-level lower and upper propensity bounds. A `truncation_grid`
#' takes precedence and expands each lower bound to `c(lower, 1 - lower)`. A
#' single `truncation` bound is honoured verbatim. When both are `NULL` the run
#' spans the whole unit interval, `c(0, 1)`.
#'
#' @param truncation A length-two bound or `NULL`.
#' @param truncation_grid A vector of lower bounds or `NULL`.
#'
#' @return A list with numeric vectors `lower` and `upper`, one entry per level.
#' @keywords internal
#' @noRd
resolve_truncation_levels <- function(
  truncation,
  truncation_grid,
  call = rlang::caller_env()
) {
  if (!is.null(truncation_grid)) {
    validate_truncation_grid(truncation_grid, call = call)
    list(lower = truncation_grid, upper = 1 - truncation_grid)
  } else if (!is.null(truncation)) {
    validate_truncation(truncation, call = call)
    list(lower = truncation[[1]], upper = truncation[[2]])
  } else {
    list(lower = 0, upper = 1)
  }
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

#' Map a binary exposure to a 0/1 double
#'
#' The lower of the two distinct levels becomes `0` and the higher becomes `1`,
#' matching the coding a logistic model expects for its response.
#'
#' @param exposure A two-level exposure vector.
#'
#' @return A double vector of `0` and `1`.
#' @keywords internal
#' @noRd
binary_to_01 <- function(exposure) {
  as.double(as.integer(factor(exposure)) - 1L)
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

# ---- Bootstrap engine -----------------------------------------------------

#' Run the ETA.Bias parametric bootstrap
#'
#' Fits the nuisance models once on the observed data to define the target
#' `truth` and the data-generating mechanism, then draws `n_boot` bootstrap
#' datasets and refits both models on each. Every truncation level reuses the
#' one shared set of draws and refits, so a sweep isolates the truncation bound.
#'
#' @param a The 0/1 exposure vector.
#' @param y The outcome vector.
#' @param covariates A data frame of covariates.
#' @param exposure_name,outcome_name The exposure and outcome column names.
#' @param estimator One of `"ipw"`, `"gcomp"`, or `"aipw"`.
#' @param outcome_type `"continuous"` or `"binary"`.
#' @param error_dist `"normal"` or `"empirical"`; used for continuous outcomes.
#' @param levels A list of `lower` and `upper` truncation vectors.
#' @param n_boot The number of bootstrap datasets.
#' @param formulas The resolved model formulas and default flag.
#'
#' @return A list with `truth` (a scalar), `boot_mean` and `boot_sd` (one entry
#'   per level), `boot_estimates` (a list of per-level estimate vectors with
#'   non-finite draws removed), and `n_valid` (the number of retained draws per
#'   level).
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
    formulas = formulas
  )
  observed <- fitters$observed(a, y)

  n <- length(a)
  n_levels <- length(levels$lower)
  # One row per bootstrap draw, one column per truncation level.
  estimates <- matrix(NA_real_, nrow = n_boot, ncol = n_levels)
  # All covariate-row resamples are drawn together, one column per draw.
  index_matrix <- matrix(
    sample.int(n, n * n_boot, replace = TRUE),
    nrow = n,
    ncol = n_boot
  )

  for (b in seq_len(n_boot)) {
    idx <- index_matrix[, b]
    a_star <- stats::rbinom(n, 1L, observed$ps[idx])
    mu <- ifelse(a_star == 1, observed$q1[idx], observed$q0[idx])
    y_star <- draw_bootstrap_outcome(
      mu = mu,
      outcome_type = outcome_type,
      error_dist = error_dist,
      residuals = observed$residuals,
      sigma = observed$sigma
    )
    refit <- fitters$refit(idx, a_star, y_star)
    for (level in seq_len(n_levels)) {
      ps_level <- pmin(
        pmax(refit$ps, levels$lower[[level]]),
        levels$upper[[level]]
      )
      estimates[b, level] <- eta_estimate(
        estimator = estimator,
        a = a_star,
        y = y_star,
        ps = ps_level,
        q1 = refit$q1,
        q0 = refit$q0
      )
    }
  }

  # A bootstrap exposure that lands in a single arm, or a refit propensity of
  # exactly 0 or 1, yields a non-finite estimate; such draws carry no
  # information and are dropped per truncation level before summarizing.
  boot_estimates <- lapply(seq_len(n_levels), function(level) {
    column <- estimates[, level]
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
    truth = observed$truth,
    boot_mean = vapply(boot_estimates, mean, numeric(1)),
    boot_sd = vapply(boot_estimates, stats::sd, numeric(1)),
    boot_estimates = boot_estimates,
    n_valid = n_valid
  )
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
#' terms. Both closures return the fitted propensity and the outcome predictions
#' at the two exposure values, the only quantities the estimators consume.
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
  formulas
) {
  numeric_covariates <- all(vapply(covariates, is.numeric, logical(1)))
  if (formulas$default && numeric_covariates) {
    eta_fitters_matrix(as.matrix(covariates), q_family)
  } else {
    eta_fitters_formula(
      covariates = covariates,
      exposure_name = exposure_name,
      outcome_name = outcome_name,
      q_family = q_family,
      formulas = formulas
    )
  }
}

#' Fast matrix fitters for the default main-effects models
#'
#' Fits the logistic treatment model and the outcome model by [stats::glm.fit()]
#' on prebuilt design matrices, avoiding the per-draw cost of building a data
#' frame and parsing a formula. This is the path a default run with numeric
#' covariates takes.
#'
#' @param x_cov The observed covariate matrix, `n` by `q`.
#' @param q_family The outcome-model family.
#'
#' @return A list of `observed(a, y)` and `refit(idx, a_star, y_star)`.
#' @keywords internal
#' @noRd
eta_fitters_matrix <- function(x_cov, q_family) {
  n <- nrow(x_cov)
  intercept <- rep.int(1, n)
  binomial_family <- stats::binomial()

  # glm.fit warns on separation and near-degenerate fits, which are expected
  # under positivity violations, the very condition the diagnostic measures.
  # glm.fit returns NA coefficients for aliased (constant or collinear) columns.
  # Zeroing them reproduces the aliased-column predictions that stats::glm()
  # makes on the formula path, so the two paths agree on rank-deficient designs.
  fit_ps <- function(design, response) {
    coef <- suppressWarnings(
      stats::glm.fit(design, response, family = binomial_family)$coefficients
    )
    coef[is.na(coef)] <- 0
    binomial_family$linkinv(drop(design %*% coef))
  }
  fit_q <- function(x_cov_rows, exposure, response) {
    design <- cbind(intercept, exposure, x_cov_rows)
    fit <- suppressWarnings(
      stats::glm.fit(design, response, family = q_family)
    )
    coef <- fit$coefficients
    coef[is.na(coef)] <- 0
    design1 <- cbind(intercept, 1, x_cov_rows)
    design0 <- cbind(intercept, 0, x_cov_rows)
    list(
      coef = coef,
      rank = fit$rank,
      q1 = q_family$linkinv(drop(design1 %*% coef)),
      q0 = q_family$linkinv(drop(design0 %*% coef))
    )
  }

  observed <- function(a, y) {
    ps <- fit_ps(cbind(intercept, x_cov), a)
    q <- fit_q(x_cov, a, y)
    fitted_y <- q_family$linkinv(drop(cbind(intercept, a, x_cov) %*% q$coef))
    residuals <- y - fitted_y
    # Residual degrees of freedom count estimated parameters only. An aliased
    # column contributes no estimated parameter, so it does not deflate sigma,
    # matching the residual standard deviation on the formula path.
    df_residual <- n - q$rank
    list(
      ps = ps,
      q1 = q$q1,
      q0 = q$q0,
      truth = mean(q$q1 - q$q0),
      residuals = residuals,
      sigma = sqrt(sum(residuals^2) / df_residual)
    )
  }

  refit <- function(idx, a_star, y_star) {
    x_rows <- x_cov[idx, , drop = FALSE]
    ps <- fit_ps(cbind(intercept, x_rows), a_star)
    q <- fit_q(x_rows, a_star, y_star)
    list(ps = ps, q1 = q$q1, q0 = q$q0)
  }

  list(observed = observed, refit = refit)
}

#' General formula fitters for user-supplied models
#'
#' Fits the treatment and outcome models from arbitrary formulas by building the
#' design matrix with [stats::glm.fit()], rebuilding the model frame on each
#' bootstrap dataset. Fitting the design directly, rather than through
#' [stats::glm()], avoids dropping unused factor levels, so a resample missing a
#' rare level keeps its design column instead of aborting. Used only when the
#' user overrides a model formula or supplies a non-numeric covariate.
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
  formulas
) {
  binomial_family <- stats::binomial()
  # A bootstrap resample can drop every row that holds a rare character level,
  # re-factoring the covariate to a single level and aborting the contrasts.
  # Fixing character covariates to factors once on the observed data keeps their
  # levels constant across resamples, so the design columns never change.
  covariates[] <- lapply(covariates, function(column) {
    if (is.character(column)) factor(column) else column
  })

  assemble <- function(exposure, response) {
    frame <- covariates
    frame[[exposure_name]] <- exposure
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

  fit_pair <- function(frame) {
    ps_fit <- fit_glm(formulas$exposure, frame, binomial_family)
    q_fit <- fit_glm(formulas$outcome, frame, q_family)
    frame1 <- frame
    frame0 <- frame
    frame1[[exposure_name]] <- 1
    frame0[[exposure_name]] <- 0
    list(
      ps = predict_fit(ps_fit, frame),
      q1 = predict_fit(q_fit, frame1),
      q0 = predict_fit(q_fit, frame0),
      q_fit = q_fit
    )
  }

  observed <- function(a, y) {
    frame <- assemble(a, y)
    fit <- fit_pair(frame)
    residuals <- y - predict_fit(fit$q_fit, frame)
    list(
      ps = fit$ps,
      q1 = fit$q1,
      q0 = fit$q0,
      truth = mean(fit$q1 - fit$q0),
      residuals = residuals,
      sigma = sqrt(sum(residuals^2) / fit$q_fit$df_residual)
    )
  }

  refit <- function(idx, a_star, y_star) {
    frame <- assemble(a_star, y_star)
    frame[seq_along(covariates)] <- covariates[idx, , drop = FALSE]
    fit <- fit_pair(frame)
    list(ps = fit$ps, q1 = fit$q1, q0 = fit$q0)
  }

  list(observed = observed, refit = refit)
}

#' Apply a causal estimator to a bootstrap dataset
#'
#' Computes the average treatment effect under the chosen estimator, given the
#' bootstrap exposure and outcome, the truncated propensity, and the outcome
#' predictions at both exposure values. G-computation ignores the propensity;
#' inverse probability weighting is the Hajek difference in weighted means,
#' matching [propensity::ipw()]; the augmented estimator adds the doubly robust
#' correction.
#'
#' @param estimator One of `"ipw"`, `"gcomp"`, or `"aipw"`.
#' @param a The bootstrap 0/1 exposure.
#' @param y The bootstrap outcome.
#' @param ps The truncated fitted propensity.
#' @param q1,q0 The outcome predictions at exposure `1` and `0`.
#'
#' @return A single numeric estimate.
#' @keywords internal
#' @noRd
eta_estimate <- function(estimator, a, y, ps, q1, q0) {
  if (estimator == "gcomp") {
    return(mean(q1 - q0))
  }
  if (estimator == "ipw") {
    weight_treated <- a / ps
    weight_control <- (1 - a) / (1 - ps)
    mean_treated <- sum(weight_treated * y) / sum(weight_treated)
    mean_control <- sum(weight_control * y) / sum(weight_control)
    return(mean_treated - mean_control)
  }
  mean(
    q1 - q0 + a / ps * (y - q1) - (1 - a) / (1 - ps) * (y - q0)
  )
}

# ---- Methods --------------------------------------------------------------

method(print, eta_bias_result) <- function(x, ...) {
  results <- x@results
  n_levels <- nrow(results)
  cat_cli({
    cli::cli_h1("{S7::S7_class(x)@name}")
    cli::cli_text("Exposure: {.val {x@exposure}} ({x@exposure_type})")
    cli::cli_text("Observations: {x@n}")
    cli::cli_text("Estimator: {x@estimator}")
    cli::cli_text("Bootstrap draws: {x@params$n_boot}")
    cli::cli_text("Truth: {round(x@truth, 3)}")
    if (n_levels == 1) {
      cli::cli_text(
        "ETA.Bias: {round(results$bias, 3)} (MC SE {round(results$mc_se, 3)})"
      )
    } else {
      cli::cli_text("Truncation levels: {n_levels}")
      cli::cli_text(
        "ETA.Bias: {round(min(results$bias), 3)} to {round(max(results$bias), 3)}"
      )
    }
  })
  invisible(x)
}
