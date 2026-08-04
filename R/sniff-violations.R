# S7 methods defined here are wired up at load time by S7::methods_register().
# The file sorts after every class it dispatches on.

#' Report what a set of diagnostics found
#'
#' `sniff_violations()` answers what a run turned up, one row per finding,
#' across every diagnostic it holds. It is the row-level companion to
#' [summary()], which reports a statistic per diagnostic whether or not anything
#' was found.
#'
#' @details
#' A diagnostic reports a finding by declaring one. Each participating class
#' states which of its statistics it reports and under which label; a class that
#' declares nothing reports nothing, and no column is read as a flag on its
#' behalf. A diagnostic defined outside positively therefore reports no rows
#' until it defines a method, whatever its results hold.
#'
#' `scope` names the kind of thing a row is about, and rows of different kinds
#' are not comparable:
#'
#' * `"subgroup"`, a covariate region, which is something a reader can act on.
#' * `"overall"`, a reading taken across the whole run.
#' * `"unit"`, a single row of a distribution whose aggregate is the reading.
#'
#' Unit rows are reachable but are never the default, because a diagnostic that
#' reports per unit produces far more rows than one that reports per subgroup,
#' and stacking the two buries the smaller set.
#'
#' A reading contributes a row only when it fires, meaning at least one unit
#' falls on the wrong side of it. A run in which every unit is supported reports
#' nothing rather than reporting how well supported it was.
#'
#' `threshold` carries the number behind a row so that a reader can see the cut
#' came from an argument they passed rather than from the package. It is `NA`
#' where a reading has no such number: [check_extrapolation()] compares each
#' unit against a geometric variability computed from the data, which is a
#' property of the sample rather than a choice.
#'
#' @param x A [positivity_diagnostic] or a [positivity_check].
#' @param scope The kinds of row to report, any of `"subgroup"`, `"overall"`,
#'   and `"unit"`. Defaults to the first two.
#' @param ... Passed to methods.
#'
#' @return A [tibble][tibble::tibble] with one row per finding and the columns
#'   `diagnostic`, `scope`, `label`, `n`, `statistic`, `value`, and `threshold`.
#'   `value` and `threshold` are numeric, and `n` is an integer count that is
#'   `NA` where a row is not about a set of rows.
#'
#' @examples
#' set.seed(1)
#' n <- 300
#' x1 <- rnorm(n)
#' df <- data.frame(exposure = rbinom(n, 1, plogis(3 * x1)), x1 = x1)
#'
#' res <- check_positivity(
#'   df,
#'   exposure,
#'   x1,
#'   diagnostics = c("port", "extrapolation")
#' )
#'
#' sniff_violations(res)
#'
#' # Unit rows are available but are never the default.
#' sniff_violations(res, scope = "unit")
#' @export
sniff_violations <- new_generic(
  "sniff_violations",
  "x",
  function(x, scope = c("subgroup", "overall"), ...) {
    # Matched here rather than in each method, so an unrecognized value is
    # refused once and the refusal names the call the reader made rather than
    # the method it dispatched to.
    scope <- resolve_scope(scope)
    S7::S7_dispatch()
  }
)

#' Match a requested scope against the kinds a row can have
#'
#' Every method resolves `scope` through here, so an unrecognized value is
#' refused the same way wherever it was supplied.
#'
#' @param scope The requested scopes.
#' @param call The calling environment, used to build the error's call.
#'
#' @return The matched scopes.
#' @keywords internal
#' @noRd
resolve_scope <- function(scope, call = rlang::caller_env()) {
  resolve_arg_match(
    rlang::arg_match(
      scope,
      values = c("subgroup", "overall", "unit"),
      multiple = TRUE,
      error_arg = "scope",
      error_call = call
    ),
    call = call
  )
}

#' Build findings in the shared row schema
#'
#' The single place the column set and its types are decided, so that the rows
#' of one diagnostic stack against another's without a type being widened by one
#' method's choice.
#'
#' @param diagnostic The diagnostic the row came from.
#' @param scope The kind of thing the row is about.
#' @param label What the row names, in the diagnostic's own terms.
#' @param n How many rows of `@results` the finding covers, or `NA` where the
#'   row is not about a set of rows.
#' @param statistic The statistic the row reports.
#' @param value The statistic's value.
#' @param threshold The number the value was compared against, or `NA`.
#'
#' @return A [tibble][tibble::tibble] in the shared schema.
#' @keywords internal
#' @noRd
sniff_rows <- function(
  diagnostic,
  scope,
  label,
  n,
  statistic,
  value,
  threshold
) {
  tibble::tibble(
    diagnostic = as.character(diagnostic),
    scope = as.character(scope),
    label = as.character(label),
    n = as.integer(n),
    statistic = as.character(statistic),
    value = as.double(value),
    threshold = as.double(threshold)
  )
}

#' The shared row schema with no rows in it
#'
#' Reporting nothing is the common answer, so it comes back in the same shape as
#' a populated one rather than needing a guard at every call site.
#'
#' @return A zero-row [tibble][tibble::tibble] in the shared schema.
#' @keywords internal
#' @noRd
empty_sniff_rows <- function() {
  sniff_rows(
    diagnostic = character(0),
    scope = character(0),
    label = character(0),
    n = integer(0),
    statistic = character(0),
    value = double(0),
    threshold = double(0)
  )
}

#' Stack findings, keeping the schema when there are none
#'
#' @param rows A list of tibbles in the shared schema.
#'
#' @return One tibble in the shared schema.
#' @keywords internal
#' @noRd
bind_sniff_rows <- function(rows) {
  populated <- rows[vapply(rows, nrow, integer(1)) > 0]
  if (length(populated) == 0) {
    return(empty_sniff_rows())
  }
  vctrs::vec_rbind(!!!populated)
}

# ---- Methods ---------------------------------------------------------------

method(sniff_violations, positivity_diagnostic) <- function(
  x,
  scope = c("subgroup", "overall"),
  ...
) {
  empty_sniff_rows()
}

#' The prevalence cut behind each reported PoRT subgroup
#'
#' A sequential run resolves one threshold per wave, and the censoring side
#' resolves its own from the uncensored-so-far risk sets, so a row's cut is
#' looked up by the wave and the side it came from rather than taken from the
#' run as a whole.
#'
#' @param x A `port_result`.
#' @param found The reported subgroups.
#'
#' @return A numeric vector, one threshold per row of `found`.
#' @keywords internal
#' @noRd
port_row_thresholds <- function(x, found) {
  if (!"time" %in% names(found)) {
    return(rep(common_beta(x@beta), nrow(found)))
  }
  is_censoring <- if ("type" %in% names(found)) {
    found$type == "censoring"
  } else {
    logical(nrow(found))
  }
  thresholds <- x@beta[found$time]
  thresholds[is_censoring] <- x@censoring_beta[found$time[is_censoring]]
  thresholds
}

method(sniff_violations, port_result) <- function(
  x,
  scope = c("subgroup", "overall"),
  ...
) {
  if (!"subgroup" %in% scope) {
    return(empty_sniff_rows())
  }
  found <- x@results[x@results$low_support, , drop = FALSE]
  if (nrow(found) == 0) {
    return(empty_sniff_rows())
  }
  sniff_rows(
    diagnostic = "port",
    scope = "subgroup",
    label = found$description,
    n = found$n,
    statistic = "prevalence",
    value = found$prevalence,
    threshold = port_row_thresholds(x, found)
  )
}

method(sniff_violations, hat_values_result) <- function(
  x,
  scope = c("subgroup", "overall"),
  ...
) {
  rows <- list()
  if ("overall" %in% scope && x@exceeds_null) {
    rows$overall <- sniff_rows(
      diagnostic = "hat_values",
      scope = "overall",
      label = paste("phi-hat exceeds the", x@params$null_method, "null"),
      n = NA_integer_,
      statistic = "phi_hat",
      value = x@phi_hat,
      threshold = x@null_quantile
    )
  }
  if ("unit" %in% scope) {
    high_leverage <- x@results[x@results$high_leverage, , drop = FALSE]
    if (nrow(high_leverage) > 0) {
      rows$unit <- sniff_rows(
        diagnostic = "hat_values",
        scope = "unit",
        # The column is a conservative indicator feeding phi-hat rather than a
        # calibrated accept or reject rule for one unit, so the rows are named
        # for the column instead of for a violation.
        label = "high_leverage",
        n = NA_integer_,
        statistic = "hat_value",
        value = high_leverage$hat_value,
        threshold = x@params$threshold * x@p / x@n
      )
    }
  }
  bind_sniff_rows(rows)
}

method(sniff_violations, extrapolation_result) <- function(
  x,
  scope = c("subgroup", "overall"),
  ...
) {
  results <- x@results
  beyond <- results$low_support
  # in_hull is missing throughout when the hull test did not run, so the reading
  # fires for nobody rather than for everybody.
  outside_hull <- if (x@hull_run) !results$in_hull else logical(nrow(results))

  rows <- list()
  if ("overall" %in% scope) {
    if (any(beyond)) {
      rows$supported <- sniff_rows(
        diagnostic = "extrapolation",
        scope = "overall",
        label = "beyond one geometric variability of an opposite unit",
        n = sum(beyond),
        statistic = "prop_supported",
        value = mean(!beyond),
        threshold = NA_real_
      )
    }
    if (any(outside_hull)) {
      rows$hull <- sniff_rows(
        diagnostic = "extrapolation",
        scope = "overall",
        label = "outside the opposite-exposure hull",
        n = sum(outside_hull),
        statistic = "prop_in_hull",
        value = mean(results$in_hull),
        threshold = NA_real_
      )
    }
  }
  # A unit contributes a row for the reading that fires on it and no row for the
  # one that does not, so the two readings cost what they found rather than one
  # row per unit each.
  if ("unit" %in% scope) {
    if (any(beyond)) {
      rows$beyond_units <- sniff_rows(
        diagnostic = "extrapolation",
        scope = "unit",
        label = "beyond one geometric variability of an opposite unit",
        n = NA_integer_,
        statistic = "gower_min",
        value = results$gower_min[beyond],
        threshold = NA_real_
      )
    }
    if (any(outside_hull)) {
      rows$hull_units <- sniff_rows(
        diagnostic = "extrapolation",
        scope = "unit",
        label = "outside the opposite-exposure hull",
        n = NA_integer_,
        statistic = "in_hull",
        # Hull membership has no per-unit quantity behind it, and these rows
        # exist only where it is FALSE, so a zero would be determined by the row
        # existing rather than measured. It would also read as the opposite of
        # what it says: a Gower distance of zero in the same column is a unit
        # sitting on top of an opposite-exposure unit. The missing value is
        # repeated rather than left to recycle, because it is the only column
        # here that would otherwise carry the number of rows.
        value = rep(NA_real_, sum(outside_hull)),
        threshold = NA_real_
      )
    }
  }
  bind_sniff_rows(rows)
}

#' @rdname sniff_violations
#' @usage NULL
#' @export
method(sniff_violations, positivity_check) <- function(
  x,
  scope = c("subgroup", "overall"),
  ...
) {
  rows <- purrr::map2(
    names(x@checks),
    x@checks,
    function(name, check) {
      found <- sniff_violations(check, scope = scope)
      # The container names a child by the name it holds it under, which is what
      # `names()` lists and `$` extracts by, rather than by whatever the child
      # calls itself.
      found$diagnostic <- rep(name, nrow(found))
      found
    }
  )
  bind_sniff_rows(rows)
}
