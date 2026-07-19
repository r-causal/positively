# S7 methods defined here are wired up at load time by S7::methods_register().

#' The hat-value diagnostic result class
#'
#' `hat_values_result` is the S7 class returned by [check_hat_values()]. It
#' extends [positivity_diagnostic] with the scalar summaries of the null
#' comparison. It is created internally and is not constructed directly by
#' users.
#'
#' @noRd
hat_values_result <- new_class(
  "hat_values_result",
  parent = positivity_diagnostic,
  properties = list(
    phi_hat = class_double,
    null_dist = class_double,
    null_quantile = class_double,
    exceeds_null = class_logical,
    p = class_integer
  ),
  validator = function(self) {
    if (length(self@phi_hat) != 1) {
      "@phi_hat must be a single value"
    } else if (length(self@null_quantile) != 1) {
      "@null_quantile must be a single value"
    } else if (length(self@exceeds_null) != 1) {
      "@exceeds_null must be a single value"
    } else if (length(self@p) != 1 || is.na(self@p)) {
      "@p must be a single integer"
    }
  }
)

#' Diagnose positivity for continuous exposures with hat values
#'
#' `check_hat_values()` is the hat-value (leverage) positivity diagnostic of
#' Moodie and Schulz (2025) for continuous exposures. It measures how far a set
#' of candidate exposure-covariate combinations sits from the observed data by
#' the leverage those combinations would exert on a linear model, then compares
#' the observed leverage profile against a null in which the exposure is drawn
#' independently of the covariates. It reports a diagnostic only and assigns no
#' severity grade.
#'
#' @details
#' The method (the paper's Box 1) fits the linear design matrix
#' \eqn{M = [1, D, X]}, where `D` is the continuous exposure and `X` the
#' covariates, so the number of parameters is \eqn{p = q + 2} for `q`
#' covariates. For each requested exposure percentile \eqn{d_q} and each
#' observed covariate row \eqn{x_i}, it forms the candidate point
#' \eqn{x_* = (1, d_q, x_i)} and computes its leverage (hat value)
#' \eqn{h_* = x_*^\top (M^\top M)^{-1} x_*}. A candidate is flagged as high
#' leverage when \eqn{h_* > \mathrm{threshold} \times p / n}, with the default
#' `threshold` of 2 reproducing the paper's \eqn{2p/n} rule. The summary
#' statistic \eqn{\hat{\phi}} is the proportion of candidates that are high
#' leverage.
#'
#' Because a high value of \eqn{\hat{\phi}} can arise from the shape of the
#' exposure distribution alone, the diagnostic calibrates it against a null in
#' which the exposure is drawn independently of the covariates. Each of
#' `null_reps` replicates redraws the exposure, rebuilds the design matrix, and
#' recomputes \eqn{\hat{\phi}} at the same candidate points. The observed
#' \eqn{\hat{\phi}} is compared against the `conf_level` quantile of this null
#' distribution; when it exceeds that quantile, `exceeds_null` is `TRUE`.
#'
#' Three null-resampling schemes are available through `null_method`.
#' `"permutation"` (the default) permutes the observed exposures, `"bootstrap"`
#' resamples them with replacement, and `"gaussian"` draws from a normal
#' distribution matched to the observed exposure mean and standard deviation.
#' The permutation and bootstrap schemes are empirical: they preserve the
#' observed exposure distribution exactly. A uniform null is deliberately not
#' offered, because a uniform draw combined with a skewed exposure fabricates
#' violations that are not present, whereas the empirical and Gaussian nulls
#' stay quiet in that case.
#'
#' The per-candidate `high_leverage` flag is a conservative leverage indicator
#' that feeds \eqn{\hat{\phi}}, not a calibrated per-unit accept or reject rule.
#' The absolute value of \eqn{\hat{\phi}} depends steeply on the candidate
#' percentile grid, so the null comparison, not the raw magnitude, carries the
#' diagnostic signal.
#'
#' @param .data A data frame.
#' @param .exposure The continuous exposure column, selected with data-masking.
#'   `check_hat_values()` aborts for binary or categorical exposures.
#' @param .covariates The covariate columns, selected with tidyselect.
#' @param probs A numeric vector of exposure percentiles, each in `[0, 1]`, at
#'   which candidate points are formed. Defaults to `seq(0.05, 0.95, by =
#'   0.05)`.
#' @param threshold The leverage cutoff multiplier. A candidate is flagged when
#'   its hat value exceeds `threshold * p / n`. Defaults to `2`.
#' @param null_reps The number of null replicates used to calibrate
#'   \eqn{\hat{\phi}}. Defaults to `500`.
#' @param null_method The null-resampling scheme, one of `"permutation"`
#'   (permute observed exposures, the default), `"bootstrap"` (resample observed
#'   exposures with replacement), or `"gaussian"` (draw from a matched normal).
#'   A uniform null is deliberately not offered; see Details.
#' @param conf_level The null quantile, strictly between 0 and 1, against which the observed
#'   \eqn{\hat{\phi}} is compared. Defaults to `0.95`.
#'
#' @return A `hat_values_result` object, an S7 subclass of
#'   [positivity_diagnostic]. Its `@results` tibble has one row per candidate
#'   point with columns `.id` (the covariate row index), `prob` (the exposure
#'   percentile), `value` (the candidate exposure value \eqn{d_q}), `hat_value`
#'   (the leverage), and `high_leverage` (the logical flag). It also carries the
#'   scalar properties `@phi_hat`, `@null_dist`, `@null_quantile`,
#'   `@exceeds_null`, and `@p`, the number of model parameters
#'   \eqn{p = q + 2} that sets the \eqn{2p/n} leverage cutoff.
#'
#' @references
#' Moodie EEM, Schulz J (2025). A Simple Diagnostic for the Positivity
#' Assumption for Continuous Exposures. *Statistics in Medicine*, 44:e70194.
#' \doi{10.1002/sim.70194}
#'
#' @examples
#' set.seed(1)
#' n <- 200
#' x1 <- rnorm(n)
#' # The exposure depends on x1, so extreme doses are implausible for some
#' # covariate values.
#' dose <- rnorm(n, mean = x1)
#' df <- data.frame(dose = dose, x1 = x1)
#'
#' # null_reps is small here to keep the example fast.
#' result <- check_hat_values(df, dose, x1, null_reps = 50)
#' result
#'
#' @export
check_hat_values <- function(
  .data,
  .exposure,
  .covariates,
  probs = seq(0.05, 0.95, by = 0.05),
  threshold = 2,
  null_reps = 500,
  null_method = c("permutation", "bootstrap", "gaussian"),
  conf_level = 0.95
) {
  validate_data_frame(.data)
  validate_prob(probs, arg_name = "probs")
  validate_probability(conf_level, arg_name = "conf_level")
  validate_positive_number(threshold, arg_name = "threshold")
  validate_count(null_reps, arg_name = "null_reps")
  null_method <- rlang::arg_match(null_method)

  exposure_pos <- tidyselect::eval_select(rlang::enquo(.exposure), .data)
  if (length(exposure_pos) != 1) {
    abort(
      "{.arg .exposure} must select exactly one column, not {length(exposure_pos)}.",
      error_class = "positively_selection_error"
    )
  }
  exposure_name <- names(exposure_pos)
  exposure_vec <- .data[[exposure_pos]]

  exposure_type <- detect_exposure_type(exposure_vec)
  if (exposure_type != "continuous") {
    abort(
      c(
        "{.fn check_hat_values} supports continuous exposures only.",
        i = "{.arg .exposure} was detected as {.val {exposure_type}}."
      ),
      error_class = "positively_exposure_type_error"
    )
  }

  covariate_pos <- tidyselect::eval_select(rlang::enquo(.covariates), .data)
  validate_column_selection(covariate_pos, ".covariates")

  covariate_names <- names(covariate_pos)
  validate_numeric_columns(.data, covariate_names, ".covariates")
  validate_numeric_columns(.data, exposure_name, ".exposure")

  dose <- as.double(exposure_vec)
  covariates <- as.matrix(.data[covariate_pos])
  n <- nrow(.data)

  design <- cbind(1, dose, covariates)
  p <- ncol(design)
  design_qr <- qr(design)
  if (design_qr$rank < p) {
    abort(
      c(
        "The design matrix formed from {.arg .exposure} and {.arg .covariates} is not full rank.",
        i = "Check for constant or collinear columns in {.arg .covariates}."
      ),
      error_class = "positively_rank_error"
    )
  }
  gram_inverse <- gram_inverse_from_qr(design_qr, p)
  cutoff <- threshold * p / n

  candidate_values <- stats::quantile(
    dose,
    probs = probs,
    names = FALSE,
    type = 7
  )
  hat_values <- candidate_hat_values(gram_inverse, candidate_values, covariates)

  high_leverage <- hat_values > cutoff
  results <- tibble::tibble(
    .id = rep(seq_len(n), times = length(probs)),
    prob = rep(probs, each = n),
    value = rep(candidate_values, each = n),
    hat_value = hat_values,
    high_leverage = high_leverage
  )
  phi_hat <- mean(high_leverage)

  null_dist <- null_phi_distribution(
    dose = dose,
    covariates = covariates,
    candidate_values = candidate_values,
    cutoff = cutoff,
    null_reps = null_reps,
    null_method = null_method
  )
  null_quantile <- stats::quantile(
    null_dist,
    probs = conf_level,
    names = FALSE,
    type = 7
  )

  hat_values_result(
    results = results,
    exposure = exposure_name,
    exposure_type = exposure_type,
    n = as.integer(n),
    params = list(
      probs = probs,
      threshold = threshold,
      null_reps = null_reps,
      null_method = null_method,
      conf_level = conf_level
    ),
    call = rlang::current_call(),
    phi_hat = phi_hat,
    null_dist = null_dist,
    null_quantile = null_quantile,
    exceeds_null = phi_hat > null_quantile,
    p = as.integer(p)
  )
}

#' Invert a design cross-product from its QR decomposition
#'
#' Forms \eqn{(M^\top M)^{-1}} from the QR decomposition of \eqn{M} as
#' \eqn{(R^\top R)^{-1}}, undoing the column pivoting so the result is in the
#' original column order. Working from \eqn{R} keeps the effective condition
#' number at that of \eqn{M} rather than its square.
#'
#' @param design_qr The QR decomposition of the design matrix.
#' @param p The number of design columns.
#'
#' @return The `p` by `p` inverse cross-product matrix.
#' @noRd
gram_inverse_from_qr <- function(design_qr, p) {
  pivoted_inverse <- chol2inv(qr.R(design_qr))
  gram_inverse <- matrix(0, p, p)
  gram_inverse[design_qr$pivot, design_qr$pivot] <- pivoted_inverse
  gram_inverse
}

#' Leverage of every candidate point against a fitted design
#'
#' Computes the hat value \eqn{x_*^\top A x_*} for every candidate point
#' \eqn{x_* = (1, d_q, x_i)}, formed by crossing each candidate exposure value
#' with every observed covariate row, given a precomputed
#' \eqn{A = (M^\top M)^{-1}}.
#'
#' @param gram_inverse The inverse of the design cross-product, a `p` by `p`
#'   matrix.
#' @param candidate_values A numeric vector of candidate exposure values.
#' @param covariates The observed covariate matrix, `n` by `q`.
#'
#' @return A numeric vector of length `length(candidate_values) * n`, ordered so
#'   that all `n` covariate rows for the first candidate value come first.
#' @noRd
candidate_hat_values <- function(gram_inverse, candidate_values, covariates) {
  n <- nrow(covariates)
  intercept <- rep(1, n)
  per_value <- lapply(candidate_values, function(value) {
    candidates <- cbind(intercept, rep(value, n), covariates)
    rowSums((candidates %*% gram_inverse) * candidates)
  })
  unlist(per_value, use.names = FALSE)
}

#' The null distribution of the leverage summary statistic
#'
#' Draws `null_reps` replicates in which the exposure is redrawn independently
#' of the covariates according to `null_method`, then recomputes
#' \eqn{\hat{\phi}} at the same candidate points.
#'
#' @param dose The observed exposure vector.
#' @param covariates The observed covariate matrix.
#' @param candidate_values The candidate exposure values, held fixed from the
#'   observed exposure.
#' @param cutoff The leverage cutoff `threshold * p / n`.
#' @param null_reps The number of replicates.
#' @param null_method One of `"permutation"`, `"bootstrap"`, or `"gaussian"`.
#'
#' @return A numeric vector of length `null_reps` of null \eqn{\hat{\phi}}
#'   values.
#' @noRd
null_phi_distribution <- function(
  dose,
  covariates,
  candidate_values,
  cutoff,
  null_reps,
  null_method
) {
  n <- length(dose)
  dose_mean <- mean(dose)
  dose_sd <- stats::sd(dose)
  draw_null_dose <- switch(
    null_method,
    permutation = function() sample(dose),
    bootstrap = function() sample(dose, size = n, replace = TRUE),
    gaussian = function() stats::rnorm(n, mean = dose_mean, sd = dose_sd)
  )
  vapply(
    seq_len(null_reps),
    function(rep) {
      null_design <- cbind(1, draw_null_dose(), covariates)
      null_gram_inverse <- gram_inverse_from_qr(
        qr(null_design),
        ncol(null_design)
      )
      null_hat_values <- candidate_hat_values(
        null_gram_inverse,
        candidate_values,
        covariates
      )
      mean(null_hat_values > cutoff)
    },
    numeric(1)
  )
}

# ---- Methods --------------------------------------------------------------

method(print, hat_values_result) <- function(x, ...) {
  n_flagged <- sum(x@results$high_leverage)
  n_candidates <- nrow(x@results)
  cat_cli({
    cli::cli_h1("{S7::S7_class(x)@name}")
    cli::cli_text("Exposure: {.val {x@exposure}} ({x@exposure_type})")
    cli::cli_text("Observations: {x@n}")
    cli::cli_text(
      "Null: {x@params$null_method} ({x@params$null_reps} reps), cutoff {x@params$threshold}p/n"
    )
    cli::cli_text("phi-hat: {round(x@phi_hat, 3)}")
    cli::cli_text(
      "Null {x@params$conf_level} quantile: {round(x@null_quantile, 3)}"
    )
    cli::cli_text("Exceeds null: {x@exceeds_null}")
    cli::cli_text(
      "High-leverage candidates: {n_flagged} of {n_candidates}"
    )
  })
  invisible(x)
}
