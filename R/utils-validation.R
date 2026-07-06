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
  if (length(probs) == 0) {
    abort(
      "{.arg {arg_name}} must contain at least one value.",
      error_class = "positively_empty_error",
      call = call
    )
  }
  if (anyNA(probs)) {
    abort(
      "{.arg {arg_name}} must not contain missing values.",
      error_class = "positively_missing_error",
      call = call
    )
  }
  if (any(probs < 0 | probs > 1)) {
    bad <- probs[probs < 0 | probs > 1]
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

#' Validate a single probability in the open unit interval
#'
#' @param x A candidate probability.
#' @param arg_name The argument name used in error messages.
#' @param call The calling environment, used to build the error's call.
#'
#' @return `x`, invisibly, when it is a single value in `(0, 1)`.
#' @keywords internal
#' @noRd
validate_probability <- function(
  x,
  arg_name,
  call = rlang::caller_env()
) {
  if (!is.numeric(x) || length(x) != 1 || is.na(x)) {
    abort(
      "{.arg {arg_name}} must be a single number.",
      error_class = "positively_type_error",
      call = call
    )
  }
  if (x <= 0 || x >= 1) {
    abort(
      c(
        "{.arg {arg_name}} must be between {.val {0}} and {.val {1}}, exclusive.",
        x = "Found {.val {x}}."
      ),
      error_class = "positively_range_error",
      call = call
    )
  }
  invisible(x)
}

#' Validate a single positive number
#'
#' @param x A candidate value.
#' @param arg_name The argument name used in error messages.
#' @param call The calling environment, used to build the error's call.
#'
#' @return `x`, invisibly, when it is a single finite positive number.
#' @keywords internal
#' @noRd
validate_positive_number <- function(
  x,
  arg_name,
  call = rlang::caller_env()
) {
  if (!is.numeric(x) || length(x) != 1 || is.na(x) || !is.finite(x)) {
    abort(
      "{.arg {arg_name}} must be a single finite number.",
      error_class = "positively_type_error",
      call = call
    )
  }
  if (x <= 0) {
    abort(
      c(
        "{.arg {arg_name}} must be positive.",
        x = "Found {.val {x}}."
      ),
      error_class = "positively_range_error",
      call = call
    )
  }
  invisible(x)
}

#' Validate a single positive whole number
#'
#' @param x A candidate value.
#' @param arg_name The argument name used in error messages.
#' @param call The calling environment, used to build the error's call.
#'
#' @return `x`, invisibly, when it is a single whole number of at least one.
#' @keywords internal
#' @noRd
validate_count <- function(
  x,
  arg_name,
  call = rlang::caller_env()
) {
  is_whole <- is.numeric(x) &&
    length(x) == 1 &&
    !is.na(x) &&
    is.finite(x) &&
    x == round(x)
  if (!is_whole || x < 1) {
    abort(
      c(
        "{.arg {arg_name}} must be a single whole number of at least {.val {1}}.",
        x = "Found {.val {x}}."
      ),
      error_class = "positively_range_error",
      call = call
    )
  }
  invisible(x)
}

#' Validate a history-window lag
#'
#' A lag is a single non-negative whole number of time points, or `Inf` for the
#' full history.
#'
#' @param x A candidate lag value.
#' @param arg_name The argument name used in error messages.
#' @param call The calling environment, used to build the error's call.
#'
#' @return `x`, invisibly, when it is a valid lag.
#' @keywords internal
#' @noRd
validate_lag <- function(
  x,
  arg_name = "lag",
  call = rlang::caller_env()
) {
  is_valid <- is.numeric(x) &&
    length(x) == 1 &&
    !is.na(x) &&
    x >= 0 &&
    (is.infinite(x) || x == round(x))
  if (!is_valid) {
    abort(
      c(
        "{.arg {arg_name}} must be a single non-negative whole number or {.code Inf}.",
        x = "Found {.val {x}}."
      ),
      error_class = "positively_range_error",
      call = call
    )
  }
  invisible(x)
}

#' Validate that selected columns are numeric and complete
#'
#' @param .data The data frame.
#' @param columns A character vector of column names to check.
#' @param arg_name The argument name used in error messages.
#' @param call The calling environment, used to build the error's call.
#'
#' @return `.data`, invisibly, when every named column is numeric and free of
#'   missing values.
#' @keywords internal
#' @noRd
validate_numeric_columns <- function(
  .data,
  columns,
  arg_name,
  call = rlang::caller_env()
) {
  non_numeric <- columns[
    !vapply(
      columns,
      function(column) is.numeric(.data[[column]]),
      logical(1)
    )
  ]
  if (length(non_numeric) > 0) {
    abort(
      c(
        "{.arg {arg_name}} must select numeric columns.",
        x = "{.val {non_numeric}} {?is/are} not numeric.",
        i = "Encode factor or character covariates as numeric indicators first."
      ),
      error_class = "positively_type_error",
      call = call
    )
  }
  missing_columns <- columns[vapply(
    columns,
    function(column) anyNA(.data[[column]]),
    logical(1)
  )]
  if (length(missing_columns) > 0) {
    abort(
      c(
        "{.arg {arg_name}} must not contain missing values.",
        x = "Missing values in {.val {missing_columns}}."
      ),
      error_class = "positively_missing_error",
      call = call
    )
  }
  invisible(.data)
}
