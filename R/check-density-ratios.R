# S7 methods defined here are wired up at load time by S7::methods_register().

#' The density-ratio diagnostic result class
#'
#' `density_ratios_result` is the S7 class returned by
#' [check_density_ratios()]. It extends [positivity_diagnostic] with the raw
#' per-time ratios that produced the summaries. It is created internally and is
#' not constructed directly by users.
#'
#' @keywords internal
#' @noRd
density_ratios_result <- new_class(
  "density_ratios_result",
  parent = positivity_diagnostic,
  properties = list(
    ratios = class_list
  ),
  validator = function(self) {
    if (!all(vapply(self@ratios, is.numeric, logical(1)))) {
      "@ratios must contain only numeric vectors"
    } else if (length(self@ratios) != length(unique(self@results$time))) {
      "@ratios must have one element per time point in @results"
    }
  }
)

#' Summarize density ratios for positivity diagnosis
#'
#' `check_density_ratios()` summarizes user-supplied density ratios, the
#' Radon-Nikodym derivatives that arise as weights in modified treatment policy
#' and longitudinal causal analyses. It is a bring-your-own-numbers diagnostic:
#' a numeric vector is treated as a point treatment, an n-by-T matrix as a
#' time-varying treatment whose cumulative products across time points are
#' summarized alongside the per-time ratios, and an \pkg{lmtp} fit is read for
#' its density-ratio component. It reports a diagnostic only and assigns no
#' severity grade.
#'
#' @details
#' A density ratio is \eqn{r = g^*(a \mid h) / g(a \mid h)}, the intervention
#' treatment law relative to the observed treatment law. It is the weight that
#' reweights the observed data to the interventional world, and it is large
#' exactly where the observed conditional treatment density is small, that is,
#' where positivity is threatened. Petersen and colleagues (2012) recommend
#' inspecting the distribution of these weights, because near violations of
#' positivity surface as a heavy upper tail, while values of the observed
#' density close to zero for some treatment level signal a violation directly.
#'
#' Because a density ratio is a Radon-Nikodym derivative, its mean is one by
#' construction whenever the intervention law is absolutely continuous with
#' respect to the observed law. The mean is therefore a specification check, not
#' a violation detector: well-behaved weights alone do not guarantee positivity.
#' The signal instead lives in the upper tail (the quantiles, the maximum, and
#' the proportion of ratios above the thresholds) and in the Kish effective
#' sample size, which collapses as a few large weights come to dominate. A mean
#' visibly below one is itself informative: it marks a structural violation in
#' which intervention probability mass leaks outside the observed support, so
#' that some ratios are exactly zero.
#'
#' For each time point the diagnostic reports the mean, the quantiles at `probs`
#' (named `quantile_50`, `quantile_90`, and so on, with the `probs = 1` quantile
#' reported as `max`), the proportion of ratios exceeding each of `thresholds`
#' (named `prop_gt_10`, `prop_gt_50`), the proportion of exact-zero ratios
#' (`prop_zero`), the Kish effective sample size `ess`, and its fraction of the
#' sample size `ess_fraction`. The Kish effective sample size is
#' \eqn{(\sum_i r_i)^2 / \sum_i r_i^2}; it equals the sample size when the ratios
#' are equal and is invariant to their overall scale.
#'
#' For a matrix input the same summaries are also applied to the cumulative
#' product of the ratios across time points, under the `cumulative_` prefix.
#' This is the sequential positivity signature: mild per-step non-overlap
#' compounds multiplicatively, so the cumulative effective sample size can
#' collapse even while each per-time fraction stays healthy. For a single-column
#' matrix the cumulative summaries coincide with the per-time summaries. The
#' default thresholds of 10 and 50 are conventional weight-magnitude heuristics
#' rather
#' than values sanctioned by either source paper.
#'
#' @param ratios The density ratios. A numeric vector for a point treatment (a
#'   [propensity::psw] object is accepted and read through
#'   [vctrs::vec_data()]), an n-by-T numeric matrix whose columns are time
#'   points for a time-varying treatment, or a fitted \pkg{lmtp} object whose
#'   density-ratio component is read out. Ratios must be non-negative, finite,
#'   and free of missing values.
#' @param probs A numeric vector of quantile probabilities in `[0, 1]`. Defaults
#'   to `c(0.5, 0.9, 0.95, 0.99, 1)`; the `1` quantile is reported as the
#'   maximum.
#' @param thresholds A numeric vector of positive exceedance thresholds.
#'   Defaults to `c(10, 50)`.
#' @param ... Passed to methods.
#'
#' @return A `density_ratios_result` object, an S7 subclass of
#'   [positivity_diagnostic]. Its `@results` tibble is long, with columns `time`
#'   (the time point, `1` for a point treatment), `statistic` (the summary
#'   name), and `value`. It also carries the raw per-time ratios as the list
#'   property `@ratios`, with one element per time point.
#'
#' @references
#' Petersen ML, Porter KE, Gruber S, Wang Y, van der Laan MJ (2012). Diagnosing
#' and responding to violations in the positivity assumption.
#' *Statistical Methods in Medical Research*, 21(1):31-54.
#' \doi{10.1177/0962280210386207}
#'
#' Ring C, Schomaker M (2026). A Diagnostic to Find and Help Combat Stochastic
#' Positivity Issues, with a Focus on Continuous Treatments.
#'
#' @examples
#' # Point treatment: a numeric vector of density ratios.
#' ratios <- c(0.5, 1, 1, 2, 8)
#' check_density_ratios(ratios)
#'
#' # Time-varying treatment: one column per time point. The cumulative products
#' # across columns are summarized alongside the per-time ratios.
#' m <- matrix(c(1, 1, 2, 1, 3, 1), nrow = 2)
#' check_density_ratios(m)
#'
#' @export
check_density_ratios <- function(ratios, ...) {
  UseMethod("check_density_ratios")
}

#' @rdname check_density_ratios
#' @export
check_density_ratios.numeric <- function(
  ratios,
  probs = c(0.5, 0.9, 0.95, 0.99, 1),
  thresholds = c(10, 50),
  ...
) {
  validate_ratios(ratios)
  validate_summary_probs(probs)
  validate_thresholds(thresholds)

  ratios <- as.double(ratios)
  results <- density_ratio_rows(ratios, time = 1L, probs, thresholds)

  build_density_ratios_result(
    results = results,
    ratios = list(ratios),
    n = length(ratios),
    probs = probs,
    thresholds = thresholds,
    call = rlang::current_call()
  )
}

#' @export
check_density_ratios.psw <- function(ratios, ...) {
  check_density_ratios(vctrs::vec_data(ratios), ...)
}

#' @rdname check_density_ratios
#' @export
check_density_ratios.matrix <- function(
  ratios,
  probs = c(0.5, 0.9, 0.95, 0.99, 1),
  thresholds = c(10, 50),
  ...
) {
  validate_ratios(ratios)
  validate_summary_probs(probs)
  validate_thresholds(thresholds)

  ratios <- matrix(as.double(ratios), nrow = nrow(ratios), ncol = ncol(ratios))
  n_times <- ncol(ratios)
  cumulative <- cumulative_products(ratios)

  per_time <- purrr::map(
    seq_len(n_times),
    function(time) density_ratio_rows(ratios[, time], time, probs, thresholds)
  )
  cumulative_time <- purrr::map(
    seq_len(n_times),
    function(time) {
      density_ratio_rows(
        cumulative[, time],
        time,
        probs,
        thresholds,
        prefix = "cumulative_"
      )
    }
  )
  results <- vctrs::vec_rbind(!!!per_time, !!!cumulative_time)

  build_density_ratios_result(
    results = results,
    ratios = purrr::map(seq_len(n_times), function(time) ratios[, time]),
    n = nrow(ratios),
    probs = probs,
    thresholds = thresholds,
    call = rlang::current_call()
  )
}

#' @export
check_density_ratios.lmtp <- function(ratios, ...) {
  density_ratios <- ratios[["density_ratios"]]
  if (is.null(density_ratios)) {
    abort(
      c(
        "The {.pkg lmtp} fit has no density-ratio component.",
        i = "Expected a {.field density_ratios} element in the fitted object.",
        i = "Pass the density ratios directly as a numeric vector or matrix."
      ),
      error_class = "positively_lmtp_error"
    )
  }
  check_density_ratios(density_ratios, ...)
}

#' @export
check_density_ratios.default <- function(ratios, ...) {
  ratios_class <- class(ratios)[1]
  abort(
    c(
      "{.arg ratios} must be a numeric vector, a numeric matrix, or an {.pkg lmtp} fit.",
      x = "A {.cls {ratios_class}} is not supported."
    ),
    error_class = "positively_type_error"
  )
}

# ---- Summary machinery ----------------------------------------------------

#' Construct a density-ratio result object
#'
#' @param results The long results tibble.
#' @param ratios A list of the raw per-time ratio vectors.
#' @param n The number of observations.
#' @param probs,thresholds The resolved summary parameters.
#' @param call The originating call.
#'
#' @return A `density_ratios_result`.
#' @keywords internal
#' @noRd
build_density_ratios_result <- function(
  results,
  ratios,
  n,
  probs,
  thresholds,
  call
) {
  density_ratios_result(
    results = results,
    exposure = "density ratios",
    exposure_type = "continuous",
    n = as.integer(n),
    params = list(probs = probs, thresholds = thresholds),
    call = call,
    ratios = ratios
  )
}

#' The cumulative product of ratios across time points
#'
#' @param ratios An n-by-T numeric matrix.
#'
#' @return An n-by-T numeric matrix whose column `t` is the product of the first
#'   `t` columns of `ratios`.
#' @keywords internal
#' @noRd
cumulative_products <- function(ratios) {
  cumulative <- ratios
  for (time in seq_len(ncol(ratios))[-1]) {
    cumulative[, time] <- cumulative[, time - 1] * ratios[, time]
  }
  cumulative
}

#' The long summary rows for one vector of ratios
#'
#' @param ratios A numeric vector of density ratios.
#' @param time The time-point label.
#' @param probs,thresholds The summary parameters.
#' @param prefix A statistic-name prefix, `"cumulative_"` for cumulative rows.
#'
#' @return A tibble with columns `time`, `statistic`, and `value`.
#' @keywords internal
#' @noRd
density_ratio_rows <- function(
  ratios,
  time,
  probs,
  thresholds,
  prefix = ""
) {
  statistics <- density_ratio_statistics(ratios, probs, thresholds)
  tibble::tibble(
    time = as.integer(time),
    statistic = paste0(prefix, names(statistics)),
    value = unname(statistics)
  )
}

#' The named density-ratio summary statistics for one vector
#'
#' @param ratios A numeric vector of density ratios.
#' @param probs,thresholds The summary parameters.
#'
#' @return A named numeric vector of summaries.
#' @keywords internal
#' @noRd
density_ratio_statistics <- function(ratios, probs, thresholds) {
  n <- length(ratios)
  quantiles <- stats::quantile(ratios, probs = probs, names = FALSE, type = 7)
  exceedance <- vapply(
    thresholds,
    function(threshold) mean(ratios > threshold),
    numeric(1)
  )
  ess <- kish_ess(ratios)
  c(
    mean = mean(ratios),
    stats::setNames(quantiles, quantile_statistic_names(probs)),
    stats::setNames(exceedance, paste0("prop_gt_", number_label(thresholds))),
    ess = ess,
    ess_fraction = ess / n,
    prop_zero = mean(ratios == 0)
  )
}

#' A minimal, precision-preserving label for a number
#'
#' Formats a number to at most eight decimal places and drops trailing zeros, so
#' that whole numbers stay whole (`50` becomes `"50"`) while distinct fractional
#' values remain distinct (`99.5` becomes `"99.5"`). This keeps statistic names
#' unique for arbitrary `probs` and `thresholds`.
#'
#' @param x A numeric vector.
#'
#' @return A character vector of labels.
#' @keywords internal
#' @noRd
number_label <- function(x) {
  labels <- formatC(round(x, 8), format = "f", digits = 8)
  labels <- sub("0+$", "", labels)
  sub("\\.$", "", labels)
}

#' The English ordinal suffix for a whole number
#'
#' Returns `"st"`, `"nd"`, `"rd"`, or `"th"`, with the 11, 12, and 13 exception
#' that always takes `"th"`.
#'
#' @param n A whole number.
#'
#' @return A single-character suffix.
#' @keywords internal
#' @noRd
ordinal_suffix <- function(n) {
  n <- abs(n)
  if ((n %% 100) %in% c(11, 12, 13)) {
    return("th")
  }
  switch(
    as.character(n %% 10),
    "1" = "st",
    "2" = "nd",
    "3" = "rd",
    "th"
  )
}

#' An ordinal percentile label
#'
#' Formats a percentile for display: a whole-number percentile takes its English
#' ordinal suffix (`"1st"`, `"21st"`), while a fractional percentile takes a
#' plain `"th"` (`"99.5th"`).
#'
#' @param pct A percentile, the probability times 100.
#'
#' @return A single string, such as `"95th"`.
#' @keywords internal
#' @noRd
percentile_label <- function(pct) {
  suffix <- if (pct == round(pct)) ordinal_suffix(round(pct)) else "th"
  paste0(number_label(pct), suffix)
}

#' The statistic names for a set of quantile probabilities
#'
#' Maps each probability to `quantile_<percent>`, preserving the documented
#' default names (`quantile_50`, `quantile_90`, `quantile_95`, `quantile_99`)
#' while staying unique for arbitrary probabilities. The `probs = 1` quantile is
#' the maximum and is named `max`.
#'
#' @param probs A numeric vector of probabilities.
#'
#' @return A character vector of statistic names.
#' @keywords internal
#' @noRd
quantile_statistic_names <- function(probs) {
  ifelse(probs == 1, "max", paste0("quantile_", number_label(probs * 100)))
}

#' The Kish effective sample size of a set of weights
#'
#' Computes \eqn{(\sum_i w_i)^2 / \sum_i w_i^2}, which equals the sample size
#' for equal weights and is invariant to the overall scale. Returns zero when
#' every weight is zero.
#'
#' @param weights A non-negative numeric vector.
#'
#' @return A single numeric value.
#' @keywords internal
#' @noRd
kish_ess <- function(weights) {
  if (length(weights) == 0) {
    return(0)
  }
  peak <- max(weights)
  if (peak == 0) {
    return(0)
  }
  # Scaling by the maximum keeps the sums representable when the squared weights
  # would otherwise overflow or underflow double precision; the ratio itself is
  # scale invariant.
  scaled <- weights / peak
  sum(scaled)^2 / sum(scaled^2)
}

#' Validate a set of density ratios
#'
#' @param ratios A numeric vector or matrix of density ratios.
#' @param call The calling environment, used to build the error's call.
#'
#' @return `ratios`, invisibly, when it is numeric, non-empty, complete,
#'   non-negative, and finite.
#' @keywords internal
#' @noRd
validate_ratios <- function(ratios, call = rlang::caller_env()) {
  if (!is.numeric(ratios)) {
    ratios_class <- class(ratios)[1]
    abort(
      "{.arg ratios} must be numeric, not a {.cls {ratios_class}}.",
      error_class = "positively_type_error",
      call = call
    )
  }
  values <- as.vector(ratios)
  if (length(values) == 0) {
    abort(
      "{.arg ratios} must contain at least one value.",
      error_class = "positively_empty_error",
      call = call
    )
  }
  if (anyNA(values)) {
    abort(
      "{.arg ratios} must not contain missing values.",
      error_class = "positively_missing_error",
      call = call
    )
  }
  if (any(values < 0)) {
    abort(
      c(
        "{.arg ratios} must be non-negative.",
        x = "Density ratios are Radon-Nikodym derivatives and cannot be negative."
      ),
      error_class = "positively_range_error",
      call = call
    )
  }
  if (any(is.infinite(values))) {
    abort(
      c(
        "{.arg ratios} must be finite.",
        x = "An infinite ratio means the observed treatment density was estimated as zero.",
        i = "Inspect the density estimates before summarizing the ratios."
      ),
      error_class = "positively_range_error",
      call = call
    )
  }
  invisible(ratios)
}

#' Validate a set of exceedance thresholds
#'
#' @param thresholds A numeric vector of thresholds.
#' @param call The calling environment, used to build the error's call.
#'
#' @return `thresholds`, invisibly, when it is numeric, non-empty, complete, and
#'   strictly positive.
#' @keywords internal
#' @noRd
validate_thresholds <- function(thresholds, call = rlang::caller_env()) {
  if (!is.numeric(thresholds)) {
    thresholds_class <- class(thresholds)[1]
    abort(
      "{.arg thresholds} must be numeric, not a {.cls {thresholds_class}}.",
      error_class = "positively_type_error",
      call = call
    )
  }
  if (length(thresholds) == 0) {
    abort(
      "{.arg thresholds} must contain at least one value.",
      error_class = "positively_empty_error",
      call = call
    )
  }
  if (anyNA(thresholds)) {
    abort(
      "{.arg thresholds} must not contain missing values.",
      error_class = "positively_missing_error",
      call = call
    )
  }
  if (any(thresholds <= 0)) {
    bad <- thresholds[thresholds <= 0]
    abort(
      c(
        "{.arg thresholds} must be positive.",
        x = "Found {.val {bad}}."
      ),
      error_class = "positively_range_error",
      call = call
    )
  }
  if (anyDuplicated(thresholds)) {
    dupes <- unique(thresholds[duplicated(thresholds)])
    abort(
      c(
        "{.arg thresholds} must not contain duplicate values.",
        x = "Found {.val {dupes}} more than once."
      ),
      error_class = "positively_duplicate_error",
      call = call
    )
  }
  labels <- number_label(thresholds)
  if (anyDuplicated(labels)) {
    dupes <- unique(labels[duplicated(labels)])
    abort(
      c(
        "{.arg thresholds} must be distinct to eight decimal places.",
        x = "{.val {paste0('prop_gt_', dupes)}} would name more than one statistic."
      ),
      error_class = "positively_duplicate_error",
      call = call
    )
  }
  invisible(thresholds)
}

#' Validate the quantile probabilities for a density-ratio summary
#'
#' Reuses `validate_prob()` for the unit-interval checks and additionally
#' rejects duplicates, since each probability becomes a primary key in the tidy
#' output.
#'
#' @param probs A numeric vector of probabilities.
#' @param call The calling environment, used to build the error's call.
#'
#' @return `probs`, invisibly, when it is valid and free of duplicates.
#' @keywords internal
#' @noRd
validate_summary_probs <- function(probs, call = rlang::caller_env()) {
  validate_prob(probs, arg_name = "probs", call = call)
  if (anyDuplicated(probs)) {
    dupes <- unique(probs[duplicated(probs)])
    abort(
      c(
        "{.arg probs} must not contain duplicate values.",
        x = "Found {.val {dupes}} more than once."
      ),
      error_class = "positively_duplicate_error",
      call = call
    )
  }
  statistic_names <- quantile_statistic_names(probs)
  if (anyDuplicated(statistic_names)) {
    dupes <- unique(statistic_names[duplicated(statistic_names)])
    abort(
      c(
        "{.arg probs} must be distinct to eight decimal places of the percent scale.",
        x = "{.val {dupes}} would name more than one statistic."
      ),
      error_class = "positively_duplicate_error",
      call = call
    )
  }
  invisible(probs)
}

# ---- Methods --------------------------------------------------------------

method(glance, density_ratios_result) <- function(x, ...) {
  n_times <- length(x@ratios)
  # The headline effective sample size and maximum are the cumulative-product
  # values at the final time point, which for a point treatment is the single
  # per-time vector. Both are computed directly so the tibble is one row for any
  # probs.
  final <- final_ratios(x@ratios)
  ess <- kish_ess(final)
  tibble::tibble(
    exposure = x@exposure,
    exposure_type = x@exposure_type,
    n = x@n,
    n_times = n_times,
    ess_fraction = ess / length(final),
    max = max(final)
  )
}

method(print, density_ratios_result) <- function(x, ...) {
  n_times <- length(x@ratios)
  multi_time <- n_times > 1
  probs <- x@params$probs
  thresholds <- x@params$thresholds
  # The per-time summaries at the final time point are computed directly from
  # the raw ratios so the displayed lines never depend on which statistics the
  # requested probs happened to name.
  displayed <- x@ratios[[n_times]]
  quantile_probs <- probs[probs != 1]
  quantile_values <- if (length(quantile_probs) > 0) {
    stats::quantile(displayed, quantile_probs, names = FALSE, type = 7)
  } else {
    numeric(0)
  }
  ess <- kish_ess(displayed)
  cat_cli({
    cli::cli_h1("{S7::S7_class(x)@name}")
    cli::cli_text(
      "Density ratios: {x@n} observation{?s} across {n_times} time point{?s}"
    )
    if (multi_time) {
      cli::cli_text("Summaries shown for time point {n_times}")
    }
    cli::cli_text("Mean ratio: {round(mean(displayed), 3)}")
    for (i in seq_along(quantile_probs)) {
      cli::cli_text(
        "{percentile_label(quantile_probs[i] * 100)} percentile: {round(quantile_values[i], 3)}"
      )
    }
    cli::cli_text("Maximum: {round(max(displayed), 3)}")
    for (threshold in thresholds) {
      cli::cli_text(
        "Proportion > {threshold}: {round(mean(displayed > threshold), 4)}"
      )
    }
    cli::cli_text(
      "Proportion exactly zero: {round(mean(displayed == 0), 4)}"
    )
    cli::cli_text("Kish ESS fraction: {round(ess / length(displayed), 3)}")
    if (multi_time) {
      final <- final_ratios(x@ratios)
      cli::cli_text(
        "Cumulative ESS fraction: {round(kish_ess(final) / length(final), 3)}"
      )
    }
  })
  invisible(x)
}

#' The cumulative-product ratios at the final time point
#'
#' Multiplies the per-time ratio vectors elementwise, giving the density ratio
#' of the whole treatment history. For a point treatment this is the single
#' per-time vector.
#'
#' @param ratios A list of per-time ratio vectors.
#'
#' @return A numeric vector.
#' @keywords internal
#' @noRd
final_ratios <- function(ratios) {
  Reduce(`*`, ratios)
}
