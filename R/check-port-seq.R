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
#' Supplying `.censoring` adds a second family of per-time trees for the
#' censoring process (Chatton et al. 2025). Each censoring indicator is analyzed
#' as a binary exposure with the same two-sided reading rule, which flags a
#' subgroup whose censoring prevalence is below `beta` or above `1 - beta`. The
#' upper tail flags a subgroup nearly certain to be censored, just as the
#' exposure trees flag a subgroup nearly certain to be treated; the lower tail
#' flags a subgroup nearly certain to remain uncensored, the never-censored
#' subgroups on which the counterfactual under no censoring rests (Chatton et
#' al. 2025, page 11). Because real censoring is light in most subgroups, the
#' lower tail usually accounts for most of the flagged censoring rows. The
#' censoring tree at time \eqn{t} runs on the subjects uncensored through the
#' previous time point, the never-censored risk set, which is the full uncensored
#' cohort rather than the treatment-regime followers, since the treatment rule
#' plays no role in the censoring process. Its conditioning set is the same
#' baseline and `lag`-window history the exposure trees use, with the current
#' treatment added as a covariate. Its prevalence threshold is resolved from the
#' size of the censoring risk set, which shrinks as earlier censoring accrues.
#' Censoring rows are marked in the results by a `type` column with the value
#' `"censoring"`, against `"exposure"` for the exposure trees, and this column
#' appears only when `.censoring` is supplied.
#'
#' Supplying `.censoring` also restricts the exposure risk sets. At time \eqn{t}
#' the exposure risk set becomes the regime followers who are also uncensored
#' through time \eqn{t - 1}, since a subject censored earlier has no observed
#' exposure at time \eqn{t} and cannot follow the regime (Chatton et al. 2025,
#' Figure 2). The exposure risk sets therefore shrink faster than under the
#' follower restriction alone, and the per-time thresholds in `@beta` grow
#' accordingly. A censored subject may carry an `NA` exposure at the time points
#' after it leaves the study; the exposure type is resolved from the non-missing
#' pooled values, and such subjects fall outside the restricted exposure risk
#' sets, so their missing exposures never enter a tree.
#'
#' @param .data A data frame in wide form, one row per subject.
#' @param .exposures An ordered tidyselect of exposure columns, one per time
#'   point.
#' @param .covariates A list of tidyselect expressions, one per time point, of
#'   the time-varying covariates. A length-one list is recycled across time
#'   points.
#' @param .baseline A tidyselect of baseline covariates always included in every
#'   conditioning set. Defaults to `NULL`.
#' @param .censoring An ordered tidyselect of censoring indicators, one per time
#'   point, each coded `1` when the subject is censored and `0` otherwise. When
#'   supplied, each indicator is analyzed as its own binary-exposure tree, in
#'   addition to the exposure trees. The selection must match `.exposures` in
#'   length and every indicator must be binary; monotonicity is not required of
#'   the data, since a subject censored at a time point is dropped from later
#'   censoring risk sets. Supplying it also restricts each exposure risk set to
#'   the followers uncensored so far, and the exposure type is detected from the
#'   non-missing pooled values so that `NA` exposures for censored subjects do
#'   not distort it. Only the default `strategy = "stratified"` supports
#'   censoring. Defaults to `NULL`, which runs the exposure trees alone.
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
#'   column, and a `type` column distinguishing exposure from censoring rows when
#'   `.censoring` is supplied. `@beta` holds the per-time exposure thresholds
#'   ordered by time point, and `@censoring_beta` holds the per-time censoring
#'   thresholds when `.censoring` is supplied.
#'
#' @references
#' Danelian G, Foucher Y, Léger M, Le Borgne F, Chatton A (2023). Identification
#' of in-sample positivity violations using regression trees: the PoRT
#' algorithm. \doi{10.1515/jci-2022-0032}
#'
#' Chatton A, Schomaker M, Luque-Fernandez MA, Platt RW, Schnitzer ME (2025).
#' Is checking for sequential positivity violations getting you down? Try
#' sPoRT!
#'
#' @examples
#' set.seed(1)
#' n <- 1000
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
#' # Add censoring indicators, coded 1 when censored. Censoring is near-certain
#' # at time 2 where c1 > 1, which the censoring trees flag.
#' cens1 <- rbinom(n, 1, 0.08)
#' cens2 <- rbinom(n, 1, ifelse(c1 > 1, 0.98, 0.1))
#' df$cens1 <- cens1
#' df$cens2 <- cens2
#'
#' check_port_seq(df, c(a1, a2), list(c1), .censoring = c(cens1, cens2))
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
  validate_count(n_bins, "n_bins", min = 2)
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

  if (strategy == "pooled") {
    abort(
      c(
        "{.arg strategy} {.val pooled} is not yet implemented.",
        i = "Use {.val stratified}, the default."
      ),
      error_class = "positively_strategy_error"
    )
  }

  censoring_quo <- rlang::enquo(.censoring)
  has_censoring <- !rlang::quo_is_null(censoring_quo)
  censoring_names <- character(0)
  if (has_censoring) {
    censoring_pos <- tidyselect::eval_select(censoring_quo, .data)
    censoring_names <- names(censoring_pos)
    validate_censoring_indicators(.data, censoring_names, n_times)
  }

  pooled_exposure <- unlist(.data[exposure_names], use.names = FALSE)
  pooled_exposure <- pooled_exposure[!is.na(pooled_exposure)]
  resolved_type <- detect_exposure_type(pooled_exposure)
  untreated_level <- resolve_untreated_level(pooled_exposure, resolved_type)

  # Supplying censoring restricts the exposure risk sets to the followers who
  # are also uncensored through the previous time point; without it the exposure
  # pass sees no censoring information and keeps every follower.
  exposure_censoring_names <- if (has_censoring) censoring_names else NULL

  control <- port_control(...)
  per_time <- purrr::map(seq_len(n_times), function(t) {
    port_seq_time(
      .data = .data,
      time = t,
      exposure_names = exposure_names,
      censoring_names = exposure_censoring_names,
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

  exposure_rows <- do.call(c, purrr::map(per_time, "rows"))
  trees <- do.call(c, purrr::map(per_time, "trees"))
  beta_by_time <- vapply(per_time, function(x) x$beta, numeric(1))

  if (has_censoring) {
    per_time_censoring <- purrr::map(seq_len(n_times), function(t) {
      port_seq_censoring_time(
        .data = .data,
        time = t,
        exposure_names = exposure_names,
        censoring_names = censoring_names,
        covariate_sets = covariate_sets,
        baseline_names = baseline_names,
        lag = lag,
        alpha = alpha,
        beta_spec = beta_spec,
        gamma = gamma,
        n_bins = n_bins,
        breaks = breaks,
        control = control
      )
    })
    censoring_rows <- do.call(c, purrr::map(per_time_censoring, "rows"))
    trees <- c(trees, do.call(c, purrr::map(per_time_censoring, "trees")))
    censoring_beta_by_time <- vapply(
      per_time_censoring,
      function(x) x$beta,
      numeric(1)
    )
    results <- assemble_port_seq_results(exposure_rows, censoring_rows)
  } else {
    censoring_beta_by_time <- numeric(0)
    results <- assemble_port_results(exposure_rows, sequential = TRUE)
  }

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
    censoring_beta = censoring_beta_by_time,
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
#' @param censoring_names The censoring indicator column names, time-ordered, or
#'   `NULL` when no censoring is supplied. When present the risk set is further
#'   restricted to the subjects uncensored through the previous time point.
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
  censoring_names,
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
    untreated_level,
    censoring_names
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

#' Run the censoring sPoRT pass for one time point
#'
#' Subsets to the subjects uncensored through the previous time point, builds the
#' conditioning set with the current treatment added, resolves the per-time
#' prevalence threshold from the censoring risk-set size, and runs the PoRT
#' search with the time-point censoring indicator as a binary response.
#'
#' @param .data The full wide data frame.
#' @param time The current time point index.
#' @param exposure_names The exposure column names, time-ordered.
#' @param censoring_names The censoring indicator column names, time-ordered.
#' @param covariate_sets The per-time covariate name list.
#' @param baseline_names The baseline covariate names.
#' @param lag The history window.
#' @param alpha,gamma,n_bins,breaks,control As in [check_port()].
#' @param beta_spec Either `"gruber"` or a numeric threshold.
#'
#' @return A list with `rows` (a list of per-subgroup tibbles carrying the `time`
#'   column), `trees` (fitted `rpart` objects), and `beta` (the resolved per-time
#'   threshold, or `NA` for a skipped time point).
#' @keywords internal
#' @noRd
port_seq_censoring_time <- function(
  .data,
  time,
  exposure_names,
  censoring_names,
  covariate_sets,
  baseline_names,
  lag,
  alpha,
  beta_spec,
  gamma,
  n_bins,
  breaks,
  control
) {
  indices <- censoring_risk_set_indices(.data, censoring_names, time)
  risk_set <- .data[indices, , drop = FALSE]
  n_t <- nrow(risk_set)
  beta_value <- resolve_beta_scalar(beta_spec, n_t)

  if (n_t < 2 || !is.finite(beta_value) || beta_value >= 0.5) {
    warn(
      c(
        "Censoring at time point {time} was skipped: its risk set cannot support the reading rule.",
        i = "The risk set holds {n_t} subject{?s} and the prevalence threshold resolved to {.val {beta_value}}."
      ),
      warning_class = "positively_risk_set_warning"
    )
    return(list(rows = list(), trees = list(), beta = NA_real_))
  }

  response_vec <- risk_set[[censoring_names[[time]]]]
  if (length(unique(response_vec[!is.na(response_vec)])) < 2) {
    return(list(rows = list(), trees = list(), beta = beta_value))
  }

  conditioning <- unique(c(
    conditioning_set(
      time = time,
      covariate_sets = covariate_sets,
      exposure_names = exposure_names,
      baseline_names = baseline_names,
      lag = lag
    ),
    exposure_names[[time]]
  ))

  level_pairs <- port_level_responses(
    response_vec,
    "binary",
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

#' The censoring risk-set row indices for one sPoRT time point
#'
#' At the first time point every subject is in the censoring risk set. At later
#' time points the risk set is the subjects uncensored through the previous time
#' point, that is, those whose earlier censoring indicators are all zero. This is
#' the never-censored cohort and carries no treatment-follower restriction.
#'
#' @param .data The full wide data frame.
#' @param censoring_names The censoring indicator column names, time-ordered.
#' @param time The current time point index.
#'
#' @return An integer vector of row indices.
#' @keywords internal
#' @noRd
censoring_risk_set_indices <- function(.data, censoring_names, time) {
  n <- nrow(.data)
  if (time == 1L) {
    return(seq_len(n))
  }
  prior <- censoring_names[seq_len(time - 1L)]
  uncensored <- rep(TRUE, n)
  for (name in prior) {
    uncensored <- uncensored & .data[[name]] == 0
  }
  which(uncensored)
}

#' Validate the censoring indicator selection
#'
#' The selection must supply one indicator per time point and every indicator
#' must be binary, coded 0 or 1. Monotonicity is not checked, since a subject
#' censored at a time point is dropped from later censoring risk sets.
#'
#' @param .data The data frame.
#' @param censoring_names The resolved censoring indicator column names.
#' @param n_times The number of exposure time points.
#' @param call The calling environment, used to build the error's call.
#'
#' @return `censoring_names`, invisibly, when the selection is valid.
#' @keywords internal
#' @noRd
validate_censoring_indicators <- function(
  .data,
  censoring_names,
  n_times,
  call = rlang::caller_env()
) {
  if (length(censoring_names) != n_times) {
    abort(
      c(
        "{.arg .censoring} must select one indicator per time point.",
        i = "Found {length(censoring_names)} indicator{?s} for {n_times} exposure{?s}."
      ),
      error_class = "positively_selection_error",
      call = call
    )
  }
  non_binary <- censoring_names[
    !vapply(
      censoring_names,
      function(name) {
        values <- .data[[name]]
        all(values[!is.na(values)] %in% c(0, 1))
      },
      logical(1)
    )
  ]
  if (length(non_binary) > 0) {
    abort(
      c(
        "{.arg .censoring} must select binary 0/1 indicators.",
        x = "{.val {non_binary}} {?is/are} not coded 0/1.",
        i = "Code each censoring indicator as 1 when censored and 0 otherwise."
      ),
      error_class = "positively_type_error",
      call = call
    )
  }
  invisible(censoring_names)
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
#' When censoring indicators are supplied the risk set is the intersection of
#' the regime followers with the subjects uncensored through the previous time
#' point, since a subject censored earlier has no observed exposure at this time
#' point and is not a follower of the regime.
#'
#' @param .data The full wide data frame.
#' @param exposure_names The exposure column names, time-ordered.
#' @param time The current time point index.
#' @param exposure_type The resolved exposure type.
#' @param untreated_level The regime's untreated level.
#' @param censoring_names The censoring indicator column names, time-ordered, or
#'   `NULL` when no censoring restriction applies.
#'
#' @return An integer vector of row indices.
#' @keywords internal
#' @noRd
risk_set_indices <- function(
  .data,
  exposure_names,
  time,
  exposure_type,
  untreated_level,
  censoring_names = NULL
) {
  n <- nrow(.data)
  followers <- if (time == 1L || exposure_type != "binary") {
    seq_len(n)
  } else {
    prior <- .data[[exposure_names[[time - 1L]]]]
    which(prior == untreated_level)
  }
  if (is.null(censoring_names)) {
    return(followers)
  }
  uncensored <- censoring_risk_set_indices(.data, censoring_names, time)
  intersect(followers, uncensored)
}

#' Assemble exposure and censoring rows into the sequential results tibble
#'
#' Tags each row with a `type` column, `"exposure"` or `"censoring"`, and orders
#' the columns with `time` and `type` leading. Used only when `.censoring` is
#' supplied, so the exposure-only path keeps its original columns.
#'
#' @param exposure_rows A list of exposure subgroup tibbles.
#' @param censoring_rows A list of censoring subgroup tibbles.
#'
#' @return A tibble with the fixed sequential PoRT columns and a leading `type`
#'   column.
#' @keywords internal
#' @noRd
assemble_port_seq_results <- function(exposure_rows, censoring_rows) {
  cols <- c(
    "time",
    "type",
    "subgroup",
    "description",
    "exposure_level",
    "n",
    "proportion",
    "prevalence",
    "flagged"
  )
  tagged <- c(
    tag_port_rows(exposure_rows, "exposure"),
    tag_port_rows(censoring_rows, "censoring")
  )
  populated <- tagged[vapply(tagged, nrow, integer(1)) > 0]
  if (length(populated) == 0) {
    base <- empty_port_tibble(sequential = TRUE)
    return(vctrs::vec_cbind(
      base["time"],
      tibble::tibble(type = character(0)),
      base[setdiff(names(base), "time")]
    )[cols])
  }
  combined <- vctrs::vec_rbind(!!!populated)
  combined[cols]
}

#' Tag a list of subgroup tibbles with a row type
#'
#' @param rows A list of subgroup tibbles.
#' @param type The type label to attach, `"exposure"` or `"censoring"`.
#'
#' @return The list with a `type` column added to every non-empty tibble.
#' @keywords internal
#' @noRd
tag_port_rows <- function(rows, type) {
  lapply(rows, function(row) {
    row$type <- rep(type, nrow(row))
    row
  })
}
