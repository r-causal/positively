# S7 methods for port_result are defined in check-port.R. This file adds the
# sequential entry point and its per-time machinery.

#' Locate sequential positivity violations with the sPoRT algorithm
#'
#' `check_port_seq()` applies the PoRT reading rule one time point at a time
#' along a treatment regime, the sequential Positivity Regression Tree (sPoRT)
#' algorithm of Chatton et al. Under the stratified strategy it fits a tree at
#' each time point among the subjects still following the regime, so that a
#' violation masked by already-treated subjects is exposed. It reports a
#' diagnostic only and assigns no severity grade.
#'
#' @details
#' Data are supplied in wide form, one row per subject with one exposure column
#' per time point. Under `strategy = "stratified"` the risk set at time \eqn{t}
#' is the set of subjects still following the rule. For a monotone binary
#' treatment this is the subjects not yet initiated, that is, those whose
#' immediately preceding treatment is the untreated level; every subject is in
#' the risk set at the first time point. Restricting to these followers is what
#' reveals a violation that pooling over already-treated subjects would mask,
#' because the already-treated carry their treatment forward and inflate the
#' subgroup prevalence.
#'
#' At each time point the diagnostic runs [check_port()] with the time-\eqn{t}
#' exposure as the response and a conditioning set of the baseline covariates,
#' the time-varying covariates within the `lag` history window, and the earlier
#' exposures within that window. The reading rule and the combination search
#' behave exactly as in the point diagnostic.
#'
#' The untreated level is resolved once from the pooled exposure across all time
#' points, so the follower set is defined against the regime's untreated level
#' rather than whichever levels happen to survive at a given wave. A non-binary
#' exposure carries no follower restriction, so every subject stays in every
#' risk set. When a risk set at some time point is empty or too small to support
#' the reading rule, for example a fully treated preceding wave or a single
#' follower under the Gruber bound, that time point is skipped with a warning,
#' contributes no rows, and records `NA` in `@beta`.
#'
#' When `beta = "gruber"` the prevalence threshold is resolved per time point
#' from the risk-set size \eqn{n_t}, as \eqn{5 / (\sqrt{n_t} \, \ln n_t)}. Under
#' a monotone rule the risk set shrinks over time, so the per-time threshold
#' grows: smaller follower sets are held to a looser bound. The resolved
#' thresholds are stored in `@beta`, ordered by time point.
#'
#' @param .data A data frame in wide form, one row per subject.
#' @param .exposures An ordered tidyselect of exposure columns, one per time
#'   point.
#' @param .covariates A list of tidyselect expressions, one per time point, of
#'   the time-varying covariates. A length-one list is recycled across time
#'   points.
#' @param .baseline A tidyselect of baseline covariates always included in every
#'   conditioning set. Defaults to `NULL`.
#' @param .censoring An ordered tidyselect of censoring indicators. Reserved for
#'   a later release and currently unused. Defaults to `NULL`.
#' @param strategy The sequential strategy. Only `"stratified"` (per-time trees
#'   among the regime followers, the default) is implemented; `"pooled"` is
#'   reserved for a later release.
#' @param lag The history window: the number of earlier time points whose
#'   covariates and exposures enter each conditioning set. Defaults to `Inf`,
#'   the full history.
#' @param alpha,beta,gamma,n_bins,breaks As in [check_port()]. `beta` is
#'   resolved per time point when it is `"gruber"`.
#' @param ... Passed to [rpart::rpart.control()].
#'
#' @return A `port_result` object, an S7 subclass of [positivity_diagnostic].
#'   Its `@results` tibble carries the point columns plus a leading `time`
#'   column. `@beta` holds the per-time thresholds ordered by time point.
#'
#' @references
#' Danelian G, Foucher Y, Léger M, Le Borgne F, Chatton A (2023). Identification
#' of in-sample positivity violations using regression trees: the PoRT
#' algorithm.
#'
#' Chatton A, Schomaker M, Luque-Fernandez MA, Platt RW, Schnitzer ME (2025).
#' Is checking for sequential positivity violations getting you down? Try
#' sPoRT!
#'
#' @examples
#' set.seed(1)
#' n <- 400
#' c1 <- rnorm(n)
#' a1 <- rbinom(n, 1, 0.3)
#' # Monotone treatment: once treated, stays treated.
#' a2 <- a1
#' followers <- a1 == 0
#' # Subjects with c1 > 1 never initiate at time 2.
#' a2[followers] <- rbinom(sum(followers), 1, ifelse(c1[followers] > 1, 0, 0.4))
#' df <- data.frame(c1 = c1, a1 = a1, a2 = a2)
#'
#' result <- check_port_seq(df, c(a1, a2), list(c1))
#' result
#'
#' @export
check_port_seq <- function(
  .data,
  .exposures,
  .covariates,
  .baseline = NULL,
  .censoring = NULL,
  strategy = c("stratified", "pooled"),
  lag = Inf,
  alpha = 0.05,
  beta = 0.05,
  gamma = 2,
  n_bins = 3,
  breaks = NULL,
  ...
) {
  validate_data_frame(.data)
  validate_probability(alpha, "alpha")
  beta_spec <- validate_beta(beta)
  validate_count(gamma, "gamma")
  validate_count(n_bins, "n_bins")
  validate_breaks(breaks)
  validate_lag(lag)
  strategy <- rlang::arg_match(strategy)

  exposure_pos <- tidyselect::eval_select(rlang::enquo(.exposures), .data)
  validate_column_selection(exposure_pos, ".exposures")
  exposure_names <- names(exposure_pos)
  n_times <- length(exposure_names)

  covariate_sets <- parse_covariate_list(
    rlang::enquo(.covariates),
    .data,
    n_times = n_times
  )
  baseline_names <- eval_optional_selection(rlang::enquo(.baseline), .data)

  if (!rlang::quo_is_null(rlang::enquo(.censoring))) {
    abort(
      c(
        "{.arg .censoring} is not yet implemented.",
        i = "Censoring-indicator trees are planned for a later release."
      ),
      error_class = "positively_censoring_error"
    )
  }

  if (strategy == "pooled") {
    abort(
      c(
        "{.arg strategy} {.val pooled} is not yet implemented.",
        i = "Use {.val stratified}, the default."
      ),
      error_class = "positively_strategy_error"
    )
  }

  pooled_exposure <- unlist(.data[exposure_names], use.names = FALSE)
  resolved_type <- detect_exposure_type(pooled_exposure)
  untreated_level <- resolve_untreated_level(pooled_exposure, resolved_type)

  control <- port_control(...)
  per_time <- purrr::map(seq_len(n_times), function(t) {
    port_seq_time(
      .data = .data,
      time = t,
      exposure_names = exposure_names,
      covariate_sets = covariate_sets,
      baseline_names = baseline_names,
      exposure_type = resolved_type,
      untreated_level = untreated_level,
      lag = lag,
      alpha = alpha,
      beta_spec = beta_spec,
      gamma = gamma,
      n_bins = n_bins,
      breaks = breaks,
      control = control
    )
  })

  rows <- do.call(c, purrr::map(per_time, "rows"))
  trees <- do.call(c, purrr::map(per_time, "trees"))
  beta_by_time <- vapply(per_time, function(x) x$beta, numeric(1))

  results <- assemble_port_results(rows, sequential = TRUE)

  port_result(
    results = results,
    exposure = exposure_names,
    exposure_type = resolved_type,
    n = as.integer(nrow(.data)),
    params = list(
      alpha = alpha,
      beta = beta,
      gamma = gamma,
      n_bins = n_bins,
      breaks = breaks,
      strategy = strategy,
      lag = lag
    ),
    call = rlang::current_call(),
    trees = trees,
    alpha = as.double(alpha),
    beta = beta_by_time,
    gamma = as.double(gamma)
  )
}

#' Run sPoRT for one time point
#'
#' Subsets to the regime followers at this time point, builds the conditioning
#' set, resolves the per-time prevalence threshold, and runs the PoRT search on
#' the restricted risk set.
#'
#' @param .data The full wide data frame.
#' @param time The current time point index.
#' @param exposure_names The exposure column names, time-ordered.
#' @param covariate_sets The per-time covariate name list.
#' @param baseline_names The baseline covariate names.
#' @param exposure_type The resolved exposure type.
#' @param untreated_level The regime's untreated level, resolved once globally.
#' @param lag The history window.
#' @param alpha,gamma,n_bins,breaks,control As in [check_port()].
#' @param beta_spec Either `"gruber"` or a numeric threshold.
#'
#' @return A list with `rows` (a list of per-subgroup tibbles carrying the
#'   `time` column), `trees` (fitted `rpart` objects), and `beta` (the resolved
#'   per-time threshold, or `NA` for a skipped time point).
#' @keywords internal
#' @noRd
port_seq_time <- function(
  .data,
  time,
  exposure_names,
  covariate_sets,
  baseline_names,
  exposure_type,
  untreated_level,
  lag,
  alpha,
  beta_spec,
  gamma,
  n_bins,
  breaks,
  control
) {
  indices <- risk_set_indices(
    .data,
    exposure_names,
    time,
    exposure_type,
    untreated_level
  )
  risk_set <- .data[indices, , drop = FALSE]
  n_t <- nrow(risk_set)
  beta_value <- resolve_beta_scalar(beta_spec, n_t)

  if (n_t < 2 || !is.finite(beta_value) || beta_value >= 0.5) {
    warn(
      c(
        "Time point {time} was skipped: its risk set cannot support the reading rule.",
        i = "The risk set holds {n_t} subject{?s} and the prevalence threshold resolved to {.val {beta_value}}."
      ),
      warning_class = "positively_risk_set_warning"
    )
    return(list(
      rows = list(),
      trees = list(),
      beta = NA_real_
    ))
  }

  conditioning <- conditioning_set(
    time = time,
    covariate_sets = covariate_sets,
    exposure_names = exposure_names,
    baseline_names = baseline_names,
    lag = lag
  )

  exposure_vec <- risk_set[[exposure_names[[time]]]]
  level_pairs <- port_level_responses(
    exposure_vec,
    exposure_type,
    n_bins = n_bins,
    breaks = breaks
  )

  search <- run_port(
    covariate_data = risk_set[conditioning],
    predictors = conditioning,
    level_pairs = level_pairs,
    alpha = alpha,
    beta_value = beta_value,
    gamma = gamma,
    n_total = n_t,
    control = control
  )

  rows <- lapply(search$rows, function(subgroup_rows) {
    vctrs::vec_cbind(
      tibble::tibble(time = rep(time, nrow(subgroup_rows))),
      subgroup_rows
    )
  })

  list(rows = rows, trees = search$trees, beta = beta_value)
}

#' The regime's untreated level
#'
#' Resolves the untreated (reference) exposure level once, from the pooled
#' exposure across every time point, so that the follower restriction uses a
#' single regime-wide level rather than whichever levels happen to survive at a
#' given wave. Only binary regimes have an untreated level; other exposure types
#' return `NA`.
#'
#' @param pooled_exposure The exposure values pooled across time points.
#' @param exposure_type The resolved exposure type.
#'
#' @return The untreated level, or `NA` for a non-binary exposure.
#' @keywords internal
#' @noRd
resolve_untreated_level <- function(pooled_exposure, exposure_type) {
  if (exposure_type != "binary") {
    return(NA)
  }
  min(unique(pooled_exposure[!is.na(pooled_exposure)]))
}

#' The risk-set row indices for one sPoRT time point
#'
#' At the first time point every subject is in the risk set. For a monotone
#' binary treatment at later time points the risk set is the subjects not yet
#' initiated, whose preceding treatment is the regime's untreated level. For a
#' non-binary exposure no follower restriction applies, so every subject is
#' retained.
#'
#' @param .data The full wide data frame.
#' @param exposure_names The exposure column names, time-ordered.
#' @param time The current time point index.
#' @param exposure_type The resolved exposure type.
#' @param untreated_level The regime's untreated level.
#'
#' @return An integer vector of row indices.
#' @keywords internal
#' @noRd
risk_set_indices <- function(
  .data,
  exposure_names,
  time,
  exposure_type,
  untreated_level
) {
  n <- nrow(.data)
  if (time == 1L || exposure_type != "binary") {
    return(seq_len(n))
  }
  prior <- .data[[exposure_names[[time - 1L]]]]
  which(prior == untreated_level)
}
