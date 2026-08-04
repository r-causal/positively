# plot() calls autoplot() and prints, so both entry points draw the same views.

#' Plot the diagnostics a positivity check holds
#'
#' Draws each child's default view. With no `diagnostic`, the views are composed
#' into one panel; naming a diagnostic draws that child alone.
#'
#' @details
#' The panel needs the \pkg{patchwork} package, which composes the children into
#' a single object that is still a [ggplot2::ggplot] and so can be built on
#' further. Without it, naming a diagnostic still works, because that path draws
#' one child's plot and composes nothing.
#'
#' A diagnostic is named the way [names()] lists it and `$` extracts it, which is
#' also the name the printed report heads its sections with, so a reader can go
#' from a section straight to its plot.
#'
#' @param object A [positivity_check] from [check_positivity()].
#' @param diagnostic Optionally, the name of one diagnostic to draw. Defaults to
#'   `NULL`, which panels every child.
#' @param ... Passed to a child's own [autoplot()] method. The panel
#'   passes them to every child, so an argument only one diagnostic accepts
#'   belongs with a named `diagnostic` rather than with the panel.
#'
#' @return A [ggplot2::ggplot] object.
#'
#' @examplesIf rlang::is_installed("patchwork")
#' set.seed(1)
#' n <- 300
#' x1 <- rnorm(n)
#' df <- data.frame(exposure = rbinom(n, 1, plogis(3 * x1)), x1 = x1)
#'
#' check <- check_positivity(df, exposure, x1, diagnostics = c("port", "edp"))
#'
#' autoplot(check)
#' autoplot(check, "port")
#'
#' @name autoplot.positivity_check
#' @aliases plot.positivity_check
NULL

method(autoplot, positivity_check) <- function(
  object,
  diagnostic = NULL,
  ...
) {
  if (!is.null(diagnostic)) {
    # extract_named_check() is what `$` and `[[` refuse an unknown name with, so
    # every way of naming a child spells the refusal the same.
    return(autoplot(extract_named_check(object, diagnostic), ...))
  }
  if (length(object@checks) == 0) {
    abort(
      "There are no diagnostics to plot.",
      error_class = "positively_empty_error"
    )
  }
  if (!rlang::is_installed("patchwork")) {
    single <- paste0("autoplot(check, \"", names(object@checks)[[1]], "\")")
    abort(
      c(
        "Plotting a whole positivity check requires the {.pkg patchwork} package.",
        i = "Install {.pkg patchwork}, or plot one diagnostic: {.code {single}}."
      ),
      error_class = "positively_missing_package_error"
    )
  }
  patchwork::wrap_plots(lapply(
    object@checks,
    function(check) autoplot(check, ...)
  ))
}

method(plot, positivity_check) <- function(x, y, ...) {
  # The second positional argument names a diagnostic, as it does for
  # autoplot(), rather than carrying the base generic's sense of `y`.
  if (missing(y)) {
    print(autoplot(x, ...))
  } else {
    print(autoplot(x, y, ...))
  }
  invisible(x)
}
