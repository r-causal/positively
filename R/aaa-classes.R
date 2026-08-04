# aaa-classes.R sorts first alphabetically so the abstract parent exists when
# subclass definitions execute at load. S7 methods are registered here and wired
# up at load time by S7::methods_register() in .onLoad().

#' The abstract parent class for positivity diagnostics
#'
#' `positivity_diagnostic` is the abstract S7 parent that every positively
#' diagnostic result inherits from. It cannot be instantiated directly. It fixes
#' the shared property set (a tidy results tibble, the exposure column names and
#' type, the sample size, the resolved method parameters, and the originating
#' call) and supplies the [generics::tidy()], [generics::glance()],
#' [summary()], and [print()] behavior that every diagnostic reuses. Package
#' developers extend it when adding a new diagnostic.
#'
#' @details
#' [generics::tidy()] returns `@results` unchanged. The inherited
#' [generics::glance()] returns a one-row [tibble][tibble::tibble] holding the
#' single column `n`, the sample size, because a subclass that adds no
#' statistics of its own has nothing else to report. Each diagnostic overrides
#' it to state its own statistics beside `n`, typed as they are computed.
#'
#' [summary()] reports those statistics in long form, with the columns
#' `statistic`, `value`, and `threshold`. Not every [generics::glance()] column
#' earns a row. `value` is numeric throughout, so the character and logical
#' statistics stay behind in the wide output, as do the sample size, the
#' resolved method parameters, and any column that another statistic reports as
#' its threshold, which would otherwise be stated twice.
#'
#' `threshold` is the one number behind the row, stated in the units of the
#' quantity it cuts, which are not always the units of `value`. Where the
#' statistic is itself a reading, the cut applies to that reading, so `phi_hat`
#' sits beside the null quantile it was compared against. Where the statistic
#' counts how many rows crossed a cut, the cut applies to the per-row quantity
#' rather than to the count, so a count of low-support subgroups sits beside a
#' subgroup prevalence and a count of high-leverage candidates beside a hat
#' value.
#'
#' `threshold` is `NA` wherever there is no one number to state. That covers a
#' statistic that is a raw magnitude, a statistic governed by several
#' parameters with no single cut, a run that resolved a different cut at each
#' wave rather than one throughout, and a statistic whose cut is real but is
#' not one of the pairings the package tracks. An `NA` therefore reports that
#' the summary has no cut to show, not that the diagnostic compared nothing
#' against anything. The pairings are fixed inside positively; a diagnostic
#' defined outside the package reports `NA` throughout.
#'
#' @usage NULL
#' @param results A [tibble][tibble::tibble] of tidy diagnostic output.
#' @param exposure The exposure column name or names, time-ordered when the
#'   diagnostic is sequential.
#' @param exposure_type One of `"binary"`, `"categorical"`, or `"continuous"`.
#' @param n The number of observations, an integer.
#' @param params A list of method parameters as supplied and resolved.
#' @param call The originating call.
#'
#' @return An abstract class object. Construction of a subclass returns a
#'   `positivity_diagnostic`.
#' @export
positivity_diagnostic <- new_class(
  "positivity_diagnostic",
  abstract = TRUE,
  properties = list(
    results = class_data.frame,
    exposure = class_character,
    exposure_type = class_character,
    n = class_integer,
    params = class_list,
    call = class_call
  ),
  validator = function(self) {
    if (!tibble::is_tibble(self@results)) {
      "@results must be a tibble"
    } else if (length(self@exposure) == 0) {
      "@exposure must contain at least one column name"
    } else if (
      length(self@exposure_type) != 1 ||
        !self@exposure_type %in% c("binary", "categorical", "continuous")
    ) {
      "@exposure_type must be a single value: binary, categorical, or continuous"
    } else if (length(self@n) != 1 || is.na(self@n)) {
      "@n must be a single integer"
    }
  }
)

#' A container for a set of positivity diagnostics
#'
#' `positivity_check` is the S7 container returned by [check_positivity()]. It
#' holds one [positivity_diagnostic] child per diagnostic, named for the
#' diagnostic that produced it, together with the exposure, covariates, exposure
#' type, sample size, and call the run resolved. Its print method shows each
#' child's summary in its own section.
#'
#' @details
#' [summary()] gives the overview of the run: one row per statistic per child,
#' with the columns `diagnostic`, `statistic`, `value`, and `threshold`. Each
#' child contributes the statistics its [generics::glance()] computed, so the
#' overview reads the same way whichever diagnostics were run.
#'
#' Extract a child diagnostic with `x[["port"]]` or `x$port`, and list the
#' diagnostics the container holds with `names(x)`. Both extractors reject a
#' name the container does not hold rather than returning `NULL`. A child
#' carries its own [generics::tidy()] and [generics::glance()] methods, which
#' are where a diagnostic's full results and its wide, typed statistics are
#' read.
#'
#' @param checks A named list of [positivity_diagnostic] objects, in the order
#'   the diagnostics were requested. Each name is the diagnostic that produced
#'   the child.
#' @param exposure The exposure column name or names, time-ordered when the
#'   diagnostics are sequential.
#' @param exposure_type One of `"binary"`, `"categorical"`, or `"continuous"`.
#' @param covariates The covariate column names.
#' @param n The number of observations, an integer.
#' @param call The originating call.
#'
#' @return A `positivity_check` object. Its [summary()] method returns a
#'   [tibble][tibble::tibble] with one row per statistic per diagnostic.
#'
#' @examples
#' set.seed(1)
#' n <- 150
#' x1 <- rnorm(n)
#' df <- data.frame(exposure = rbinom(n, 1, plogis(0.5 * x1)), x1 = x1)
#'
#' res <- check_positivity(df, exposure, x1, diagnostics = "port")
#'
#' # The overview across every child.
#' summary(res)
#'
#' # List the diagnostics and extract one child.
#' names(res)
#' res[["port"]]
#' res$port
#'
#' # A child carries its own results and its own statistics.
#' tidy(res$port)
#' glance(res$port)
#' @order 1
#' @export
positivity_check <- new_class(
  "positivity_check",
  properties = list(
    checks = class_list,
    exposure = class_character,
    exposure_type = class_character,
    covariates = class_character,
    n = class_integer,
    call = class_call
  ),
  validator = function(self) {
    is_diagnostic <- vapply(
      self@checks,
      function(check) S7::S7_inherits(check, positivity_diagnostic),
      logical(1)
    )
    # names(list()) is NULL, so an empty container cannot be asked to be named.
    check_names <- names(self@checks)
    unnamed <- length(self@checks) > 0 &&
      (is.null(check_names) || any(is.na(check_names) | check_names == ""))

    if (!all(is_diagnostic)) {
      "@checks must all be positivity_diagnostic objects"
    } else if (unnamed) {
      "@checks must be a fully named list, one name per diagnostic"
    } else if (anyDuplicated(check_names) > 0) {
      "@checks must not repeat a diagnostic name"
    } else if (length(self@exposure) == 0) {
      "@exposure must contain at least one column name"
    } else if (
      length(self@exposure_type) != 1 ||
        !self@exposure_type %in% c("binary", "categorical", "continuous")
    ) {
      "@exposure_type must be a single value: binary, categorical, or continuous"
    } else if (length(self@covariates) == 0) {
      "@covariates must contain at least one column name"
    } else if (length(self@n) != 1 || is.na(self@n)) {
      "@n must be a single integer"
    }
  }
)

# ---- Shared methods -------------------------------------------------------

method(tidy, positivity_diagnostic) <- function(x, ...) {
  x@results
}

method(glance, positivity_diagnostic) <- function(x, ...) {
  tibble::tibble(n = x@n)
}

# The glance() columns that get no summary() row. `n` is sample-size metadata,
# which a container states once rather than once per child. alpha through n_boot
# are settings the caller chose rather than findings the run produced.
# null_quantile is neither: it is dropped because phi_hat already reports it as
# a threshold, and a column that appears in the threshold column would otherwise
# be stated twice. Being a cut is not on its own grounds for dropping a column:
# geometric_variability is one, is paired with nothing here, and keeps its row.
non_statistic_columns <- function() {
  c("n", "alpha", "beta", "gamma", "nearby", "mass", "n_boot", "null_quantile")
}

#' The cut behind each of a diagnostic's statistics
#'
#' Names the one number behind a statistic, for the statistics that have one,
#' stated in the units of the quantity it cuts. For a reported reading that
#' quantity is the reading itself, as `phi_hat` is compared against the null
#' quantile. For a count it is the per-row quantity being counted, as
#' `n_low_support` counts the subgroups whose prevalence fell below `beta` or
#' rose above `1 - beta`, so a threshold need not share the units of the value
#' beside it.
#'
#' A statistic absent from the returned vector reports `NA`, and absence is not
#' a claim that no cut exists. It covers a raw magnitude with nothing behind it,
#' a statistic governed by several parameters with no single cut, and a
#' statistic whose cut is real but unpaired here: `prop_supported` is the
#' fraction of units within one `@geometric_variability` and is not paired. A
#' named entry can also resolve to `NA` on its own: `common_beta()` returns `NA`
#' when a sequential run settled on a different prevalence threshold at each
#' wave, so there is no one number to state.
#'
#' The generic is deliberately unexported, so the pairings are a
#' package-internal facility. A diagnostic subclassed outside positively can
#' override [generics::glance()] and [summary()] but reports `NA` thresholds.
#'
#' @param x A [positivity_diagnostic].
#'
#' @return A named numeric vector, keyed by statistic name.
#' @keywords internal
#' @noRd
statistic_thresholds <- new_generic("statistic_thresholds", "x")

method(statistic_thresholds, positivity_diagnostic) <- function(x) {
  rlang::set_names(numeric(0), character(0))
}

#' A diagnostic's statistics in long form
#'
#' Pivots the numeric statistics of [generics::glance()] into one row each and
#' pairs every row with its threshold. Shared by the default [summary()] method
#' and by the container, which builds its overview from its children's
#' statistics rather than from their [summary()] output, since a child may
#' summarize something else entirely.
#'
#' @param x A [positivity_diagnostic].
#'
#' @return A tibble with columns `statistic`, `value`, and `threshold`.
#' @keywords internal
#' @noRd
glance_statistics <- function(x) {
  glanced <- generics::glance(x)
  is_statistic <- vapply(glanced, is.numeric, logical(1)) &
    !names(glanced) %in% non_statistic_columns()
  statistics <- glanced[is_statistic]
  thresholds <- statistic_thresholds(x)
  tibble::tibble(
    statistic = names(statistics),
    value = unname(vapply(
      statistics,
      function(column) as.numeric(column[[1]]),
      numeric(1)
    )),
    threshold = unname(thresholds[names(statistics)])
  )
}

method(summary, positivity_diagnostic) <- function(object, ...) {
  glance_statistics(object)
}

# Render a cli block to a string and write it to stdout. Print methods need
# their output on stdout so that both `print()` at the console and testthat's
# output capture see it; cli would otherwise divert to the message stream
# whenever stdout is redirected.
cat_cli <- function(expr) {
  cat(cli::cli_fmt(expr), sep = "\n")
}

method(print, positivity_diagnostic) <- function(x, ...) {
  cat_cli({
    cli::cli_h1("{S7::S7_class(x)@name}")
    cli::cli_text("Exposure: {.val {x@exposure}} ({x@exposure_type})")
    cli::cli_text("Observations: {x@n}")
    cli::cli_text(
      "Results: {nrow(x@results)} row{?s}, {ncol(x@results)} column{?s}"
    )
  })
  invisible(x)
}

method(print, positivity_check) <- function(x, ...) {
  cat_cli(cli::cli_h1("Positivity check"))
  check_names <- names(x@checks)
  for (i in seq_along(x@checks)) {
    cat_cli(cli::cli_h2("{check_names[[i]]}"))
    print(x@checks[[i]])
  }
  invisible(x)
}

#' @rdname positivity_check
#' @usage NULL
#' @order 3
method(`[[`, positivity_check) <- function(x, i, ...) {
  if (length(i) != 1) {
    abort(
      c(
        "{.arg i} must select a single diagnostic.",
        x = "You supplied {length(i)} value{?s}."
      ),
      error_class = "positively_diagnostic_error"
    )
  }
  if (is.character(i)) {
    return(extract_named_check(x, i))
  }
  if (!is.numeric(i) || is.na(i) || i != trunc(i)) {
    abort(
      "{.arg i} must be a diagnostic name or a whole-number position.",
      error_class = "positively_diagnostic_error"
    )
  }
  if (i < 1 || i > length(x@checks)) {
    abort(
      c(
        "Index {i} is out of bounds.",
        i = "The container holds {length(x@checks)} diagnostic{?s}."
      ),
      error_class = "positively_bounds_error"
    )
  }
  x@checks[[i]]
}

#' @rdname positivity_check
#' @usage NULL
#' @order 4
method(`$`, positivity_check) <- function(x, name) {
  extract_named_check(x, name)
}

# `method<-` is a replacement function, so the registration above also leaves
# `$` bound in this namespace. codetools reads that binding as a redefinition of
# base `$`, drops the handler that keeps `x$field` from being treated as a
# variable reference, and then reports every field name used anywhere in the
# package as an undefined global under R CMD check. Dropping the leftover
# binding restores the handler; the S3 method itself is registered in base's
# method table and is unaffected.
rm(`$`)

#' @rdname positivity_check
#' @usage NULL
#' @order 5
method(names, positivity_check) <- function(x) {
  names(x@checks) %||% character(0)
}

#' Extract a child diagnostic by name
#'
#' Shared by `[[` and `$` so that both spell the same refusal for a name the
#' container does not hold. Neither extractor partially matches, which keeps the
#' two in agreement on a typo.
#'
#' @param x A [positivity_check] object.
#' @param name A single diagnostic name.
#' @param call The execution environment used to build the error's call, which
#'   is the extractor's frame so that the error is blamed on `x[["port"]]` or
#'   `x$port` rather than on this helper.
#'
#' @return The named [positivity_diagnostic] child.
#' @keywords internal
#' @noRd
extract_named_check <- function(x, name, call = rlang::caller_env()) {
  idx <- match(name, names(x@checks))
  if (is.na(idx)) {
    abort(
      c(
        "{.val {name}} is not a diagnostic in this container.",
        i = "Available diagnostics are {.val {names(x@checks)}}."
      ),
      error_class = "positively_diagnostic_error",
      call = call
    )
  }
  x@checks[[idx]]
}
