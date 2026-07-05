# aaa-classes.R sorts first alphabetically so the abstract parent exists when
# subclass definitions execute at load. S7 methods are registered here and wired
# up at load time by S7::methods_register() in .onLoad().

#' The abstract parent class for positivity diagnostics
#'
#' `positivity_diagnostic` is the abstract S7 parent that every positively
#' diagnostic result inherits from. It cannot be instantiated directly. It fixes
#' the shared property set (a tidy results tibble, the exposure column names and
#' type, the sample size, the resolved method parameters, and the originating
#' call) and supplies the [generics::tidy()], [generics::glance()], and
#' [print()] behaviour that every diagnostic reuses. Package developers extend
#' it when adding a new diagnostic.
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
#' `positivity_check` is the S7 container returned by `check_positivity()`. It
#' bundles one or more [positivity_diagnostic] children with the names of the
#' diagnostics that produced them. Its print method shows each child's summary
#' in its own section. It renders no verdicts.
#'
#' @param checks A list of [positivity_diagnostic] objects.
#' @param diagnostics A character vector naming the diagnostics, aligned with
#'   `checks`.
#'
#' @return A `positivity_check` object.
#' @export
positivity_check <- new_class(
  "positivity_check",
  properties = list(
    checks = class_list,
    diagnostics = class_character
  ),
  validator = function(self) {
    is_diagnostic <- vapply(
      self@checks,
      function(check) S7::S7_inherits(check, positivity_diagnostic),
      logical(1)
    )
    if (!all(is_diagnostic)) {
      "@checks must all be positivity_diagnostic objects"
    } else if (length(self@checks) != length(self@diagnostics)) {
      "@diagnostics must have one name per element of @checks"
    }
  }
)

# ---- Shared methods -------------------------------------------------------

method(tidy, positivity_diagnostic) <- function(x, ...) {
  x@results
}

method(glance, positivity_diagnostic) <- function(x, ...) {
  tibble::tibble(
    exposure = paste(x@exposure, collapse = ", "),
    exposure_type = x@exposure_type,
    n = x@n,
    n_rows = nrow(x@results)
  )
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
  for (i in seq_along(x@checks)) {
    cat_cli(cli::cli_h2("{x@diagnostics[[i]]}"))
    print(x@checks[[i]])
  }
  invisible(x)
}
