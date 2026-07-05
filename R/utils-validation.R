# Validation helpers for the positively package. Each helper raises a classed
# `positively_error` on failure and returns its input invisibly on success, so
# that validation can be chained inside a pipeline without altering the value.

#' Validate that an input is a data frame
#'
#' @param .data The object to validate.
#' @param arg_name The argument name used in error messages.
#' @param call The calling environment, used to build the error's call.
#'
#' @return `.data`, invisibly, when it is a data frame.
#' @keywords internal
#' @noRd
validate_data_frame <- function(
  .data,
  arg_name = ".data",
  call = rlang::caller_env()
) {
  if (!is.data.frame(.data)) {
    data_class <- class(.data)[1]
    abort(
      "{.arg {arg_name}} must be a data frame, not a {.cls {data_class}}.",
      error_class = "positively_type_error",
      call = call
    )
  }
  invisible(.data)
}

#' Validate a resolved column selection
#'
#' Takes an already-resolved [tidyselect::eval_select()]-style selection, a
#' named integer vector of column positions, and aborts when it is empty.
#'
#' @param selection A named integer vector of resolved column positions.
#' @param arg_name The argument name used in error messages.
#' @param call The calling environment, used to build the error's call.
#'
#' @return `selection`, invisibly, when it selects at least one column.
#' @keywords internal
#' @noRd
validate_column_selection <- function(
  selection,
  arg_name = ".covariates",
  call = rlang::caller_env()
) {
  if (length(selection) == 0) {
    abort(
      "{.arg {arg_name}} must select at least one column.",
      error_class = "positively_empty_error",
      call = call
    )
  }
  invisible(selection)
}

#' Validate probabilities in the unit interval
#'
#' @param probs A numeric vector of probabilities.
#' @param arg_name The argument name used in error messages.
#' @param call The calling environment, used to build the error's call.
#'
#' @return `probs`, invisibly, when every element lies in `[0, 1]`.
#' @keywords internal
#' @noRd
validate_prob <- function(
  probs,
  arg_name = "probs",
  call = rlang::caller_env()
) {
  if (!is.numeric(probs)) {
    probs_class <- class(probs)[1]
    abort(
      "{.arg {arg_name}} must be numeric, not a {.cls {probs_class}}.",
      error_class = "positively_type_error",
      call = call
    )
  }
  if (any(probs < 0 | probs > 1, na.rm = TRUE)) {
    bad <- probs[!is.na(probs) & (probs < 0 | probs > 1)]
    abort(
      c(
        "{.arg {arg_name}} must be between {.val {0}} and {.val {1}}.",
        x = "Found {.val {bad}}."
      ),
      error_class = "positively_range_error",
      call = call
    )
  }
  invisible(probs)
}
