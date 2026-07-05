#' Should informational messages be suppressed?
#'
#' Reads `options(positively.quiet)`, defaulting to `FALSE`. When `TRUE`,
#' `alert_info()` stays silent so that diagnostics can run without console
#' chatter.
#'
#' @return A single logical value.
#' @keywords internal
#' @noRd
be_quiet <- function() {
  getOption("positively.quiet", default = FALSE)
}

#' Signal a classed positively error
#'
#' A thin wrapper around [cli::cli_abort()] that attaches the
#' `positively_error` condition class to every error the package raises. A
#' caller may layer an additional, more specific subclass through
#' `error_class` so that individual failures can be caught or snapshot tested.
#'
#' @param ... Passed to [cli::cli_abort()]. The first element is the message,
#'   which may use cli inline markup.
#' @param error_class Optional character vector of extra condition subclasses
#'   layered ahead of `positively_error`.
#' @param call The execution environment used to build the error's call for
#'   display.
#' @param .envir The environment in which to evaluate cli glue expressions.
#'
#' @return Never returns; always raises an error.
#' @keywords internal
#' @noRd
abort <- function(
  ...,
  error_class = NULL,
  call = rlang::caller_env(),
  .envir = parent.frame()
) {
  cli::cli_abort(
    ...,
    class = c(error_class, "positively_error"),
    call = call,
    .envir = .envir
  )
}

#' Signal a classed positively warning
#'
#' A thin wrapper around [cli::cli_warn()] that attaches the
#' `positively_warning` condition class to every warning the package raises. A
#' caller may layer an additional subclass through `warning_class`.
#'
#' @param ... Passed to [cli::cli_warn()]. The first element is the message,
#'   which may use cli inline markup.
#' @param warning_class Optional character vector of extra condition subclasses
#'   layered ahead of `positively_warning`.
#' @param call The execution environment used to build the warning's call.
#' @param .envir The environment in which to evaluate cli glue expressions.
#'
#' @return Invisibly `NULL`; raises a warning as a side effect.
#' @keywords internal
#' @noRd
warn <- function(
  ...,
  warning_class = NULL,
  call = rlang::caller_env(),
  .envir = parent.frame()
) {
  cli::cli_warn(
    ...,
    class = c(warning_class, "positively_warning"),
    call = call,
    .envir = .envir
  )
}

#' Emit an informational message
#'
#' A wrapper around [cli::cli_alert_info()] that respects
#' `options(positively.quiet)`. When the option is `TRUE`, the message is
#' suppressed entirely.
#'
#' @param .message The message text, which may use cli inline markup.
#' @param .envir The environment in which to evaluate cli glue expressions.
#'
#' @return Invisibly `NULL`; emits a message as a side effect.
#' @keywords internal
#' @noRd
alert_info <- function(.message, .envir = parent.frame()) {
  if (!be_quiet()) {
    cli::cli_alert_info(text = .message, .envir = .envir)
  }
  invisible()
}
