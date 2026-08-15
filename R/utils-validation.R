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

#' Match an argument against its menu, reporting a classed positively error
#'
#' Wraps an [rlang::arg_match()] call so that a value outside the menu carries
#' `positively_error` like every other failure the package raises, instead of
#' the bare `rlang_error` `arg_match()` throws. The matching itself is left to
#' `arg_match()`: its handling of an untouched formal default, of a supplied
#' value, and of a near miss is the contract each `check_*()` argument is
#' written against, and its message names every allowed value and, on a near
#' miss, the one that was meant. Only the condition class changes.
#'
#' `matched` is evaluated inside [rlang::try_fetch()] but stays a promise of the
#' calling function, so `arg_match()` still reads the menu from that function's
#' own formal and needs no second copy of the allowed values here.
#'
#' Forcing that promise also evaluates the caller's argument, so the handler is
#' narrowed to `rlang_error`, which is what `arg_match()` raises. An argument
#' that fails to evaluate at all raises a base R `simpleError`, and relabelling
#' that as a failure to match a menu would be untrue, so it passes through. An
#' argument expression that raises an `rlang_error` of its own is the residual
#' case this cannot separate; its message survives intact, so the report stays
#' accurate even where the class is not.
#'
#' The caught headline and bullets reach cli as values rather than as template
#' text, so an argument value containing braces stays a value instead of being
#' read as an expression, which would replace the argument error with an
#' internal cli failure and lose the package's condition classes with it.
#'
#' @param matched An [rlang::arg_match()] call on one of the caller's formals.
#' @param call The calling environment, used to build the error's call.
#'
#' @return The matched value.
#' @keywords internal
#' @noRd
resolve_arg_match <- function(matched, call = rlang::caller_env()) {
  rlang::try_fetch(
    matched,
    rlang_error = function(cnd) {
      bullets <- if (is.character(cnd$body)) cnd$body else character(0)
      templates <- sprintf("{bullets[[%d]]}", seq_along(bullets))
      abort(
        c("{cnd$message}", rlang::set_names(templates, names(bullets))),
        error_class = "positively_args_error",
        call = call
      )
    }
  )
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

#' Select exactly one column, reporting failures as a positively error
#'
#' Evaluates a single-column data-masked selection and translates a tidyselect
#' failure (an unknown column, say) into a classed `positively_selection_error`,
#' so that every user-facing failure carries the package's condition class.
#'
#' @param quo A quosure holding the data-masked selection.
#' @param .data The data frame to select from.
#' @param arg_name The argument name used in error messages.
#' @param call The calling environment, used to build the error's call.
#'
#' @return A named integer of length one, the selected column position.
#' @keywords internal
#' @noRd
eval_select_column <- function(
  quo,
  .data,
  arg_name,
  call = rlang::caller_env()
) {
  pos <- rlang::try_fetch(
    tidyselect::eval_select(quo, .data),
    error = function(cnd) {
      abort(
        "{.arg {arg_name}} must select a column that exists in {.arg .data}.",
        error_class = "positively_selection_error",
        call = call,
        parent = cnd
      )
    }
  )
  if (length(pos) != 1) {
    abort(
      "{.arg {arg_name}} must select exactly one column, not {length(pos)}.",
      error_class = "positively_selection_error",
      call = call
    )
  }
  pos
}

#' Select one or more columns, reporting failures as a positively error
#'
#' Evaluates a multi-column data-masked selection and translates a tidyselect
#' failure (an unknown column, say) into a classed `positively_selection_error`,
#' so that every user-facing failure carries the package's condition class. The
#' emptiness of the selection is left to `validate_column_selection()`, which
#' each caller applies with its own argument name.
#'
#' @param quo A quosure holding the data-masked selection.
#' @param .data The data frame to select from.
#' @param arg_name The argument name used in error messages.
#' @param call The calling environment, used to build the error's call.
#'
#' @return A named integer vector of the selected column positions.
#' @keywords internal
#' @noRd
eval_select_columns <- function(
  quo,
  .data,
  arg_name,
  call = rlang::caller_env()
) {
  rlang::try_fetch(
    tidyselect::eval_select(quo, .data),
    error = function(cnd) {
      abort(
        "{.arg {arg_name}} must select columns that exist in {.arg .data}.",
        error_class = "positively_selection_error",
        call = call,
        parent = cnd
      )
    }
  )
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

#' Validate a single logical flag
#'
#' @param x A candidate flag value.
#' @param arg_name The argument name used in error messages.
#' @param call The calling environment, used to build the error's call.
#'
#' @return `x`, invisibly, when it is a single non-missing logical.
#' @keywords internal
#' @noRd
validate_flag <- function(
  x,
  arg_name,
  call = rlang::caller_env()
) {
  if (!is.logical(x) || length(x) != 1 || is.na(x)) {
    abort(
      "{.arg {arg_name}} must be {.code TRUE} or {.code FALSE}.",
      error_class = "positively_type_error",
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
#' @param min The smallest whole number accepted. Defaults to `1`.
#' @param call The calling environment, used to build the error's call.
#'
#' @return `x`, invisibly, when it is a single whole number of at least `min`.
#' @keywords internal
#' @noRd
validate_count <- function(
  x,
  arg_name,
  min = 1,
  call = rlang::caller_env()
) {
  is_whole <- is.numeric(x) &&
    length(x) == 1 &&
    !is.na(x) &&
    is.finite(x) &&
    x == round(x)
  if (!is_whole || x < min) {
    abort(
      c(
        "{.arg {arg_name}} must be a single whole number of at least {.val {min}}.",
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

#' Validate that selected columns hold a type a model formula can encode
#'
#' Accepts numeric, logical, factor, and character columns, the types a model
#' formula expands into terms on its own. Missingness is a separate question, so
#' a caller that rejects `NA` pairs this with `validate_complete_columns()`.
#'
#' @param .data The data frame.
#' @param columns A character vector of column names to check.
#' @param arg_name The argument name used in error messages.
#' @param call The calling environment, used to build the error's call.
#'
#' @return `.data`, invisibly, when every named column holds an accepted type.
#' @keywords internal
#' @noRd
validate_encodable_columns <- function(
  .data,
  columns,
  arg_name,
  call = rlang::caller_env()
) {
  unsupported <- columns[
    !vapply(
      columns,
      function(column) {
        values <- .data[[column]]
        is.numeric(values) ||
          is.logical(values) ||
          is.factor(values) ||
          is.character(values)
      },
      logical(1)
    )
  ]
  if (length(unsupported) > 0) {
    abort(
      c(
        "{.arg {arg_name}} must select numeric, logical, factor, or character columns.",
        x = "{.val {unsupported}} {?is/are} another type.",
        i = "Convert date or other ordered columns to numeric first."
      ),
      error_class = "positively_type_error",
      call = call
    )
  }
  invisible(.data)
}

#' Validate that selected columns contain no missing values
#'
#' @param .data The data frame.
#' @param columns A character vector of column names to check.
#' @param arg_name The argument name used in error messages.
#' @param call The calling environment, used to build the error's call.
#'
#' @return `.data`, invisibly, when every named column is complete.
#' @keywords internal
#' @noRd
validate_complete_columns <- function(
  .data,
  columns,
  arg_name,
  call = rlang::caller_env()
) {
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

#' Validate that selected numeric columns contain no non-finite values
#'
#' Non-finite exposure or covariate values poison every kernel they enter: the
#' default standard-deviation half-distance becomes `NaN`, and even under a
#' supplied half-distance the self-difference `Inf - Inf` evaluates to `NaN`.
#' Non-numeric columns are skipped.
#'
#' @param .data The data frame.
#' @param columns A character vector of column names to check.
#' @param arg_name The argument name used in error messages.
#' @param call The calling environment, used to build the error's call.
#'
#' @return `.data`, invisibly, when every numeric column named is finite.
#' @keywords internal
#' @noRd
validate_finite_columns <- function(
  .data,
  columns,
  arg_name,
  call = rlang::caller_env()
) {
  non_finite_columns <- columns[vapply(
    columns,
    function(column) {
      is.numeric(.data[[column]]) && !all(is.finite(.data[[column]]))
    },
    logical(1)
  )]
  if (length(non_finite_columns) > 0) {
    abort(
      c(
        "{.arg {arg_name}} must not contain non-finite values.",
        x = "Non-finite values in {.val {non_finite_columns}}."
      ),
      error_class = "positively_range_error",
      call = call
    )
  }
  invisible(.data)
}

#' Validate a single truncation bound
#'
#' A truncation bound is a length-two numeric vector `c(lower, upper)`, with both
#' values in `[0, 1]` and `lower` strictly below `upper`. The propensity score is
#' bounded into this interval inside the estimator.
#'
#' @param truncation A candidate truncation bound.
#' @param arg_name The argument name used in error messages.
#' @param call The calling environment, used to build the error's call.
#'
#' @return `truncation`, invisibly, when it is a valid bound.
#' @keywords internal
#' @noRd
validate_truncation <- function(
  truncation,
  arg_name = "truncation",
  call = rlang::caller_env()
) {
  if (!is.numeric(truncation) || length(truncation) != 2) {
    abort(
      "{.arg {arg_name}} must be a length-two numeric vector {.code c(lower, upper)}.",
      error_class = "positively_type_error",
      call = call
    )
  }
  if (anyNA(truncation)) {
    abort(
      "{.arg {arg_name}} must not contain missing values.",
      error_class = "positively_missing_error",
      call = call
    )
  }
  if (any(truncation < 0 | truncation > 1)) {
    abort(
      c(
        "{.arg {arg_name}} must lie between {.val {0}} and {.val {1}}.",
        x = "Found {.val {truncation}}."
      ),
      error_class = "positively_range_error",
      call = call
    )
  }
  if (truncation[[1]] >= truncation[[2]]) {
    abort(
      c(
        "{.arg {arg_name}} must have its lower bound below its upper bound.",
        x = "Found lower {.val {truncation[[1]]}} and upper {.val {truncation[[2]]}}."
      ),
      error_class = "positively_range_error",
      call = call
    )
  }
  invisible(truncation)
}

#' Validate a truncation sweep grid
#'
#' A truncation grid is a numeric vector of lower bounds, each in `[0, 1/k)` for
#' an exposure with `k` levels, so that bounding every level from below leaves
#' room on the simplex. Two levels give the familiar `[0, 0.5)`, where the paired
#' upper bound `1 - lower` always stays above `lower`. The sweep shares one set
#' of bootstrap draws across every grid point.
#'
#' @param grid A candidate grid of lower bounds.
#' @param k The number of exposure levels.
#' @param arg_name The argument name used in error messages.
#' @param call The calling environment, used to build the error's call.
#'
#' @return `grid`, invisibly, when every lower bound is valid.
#' @keywords internal
#' @noRd
validate_truncation_grid <- function(
  grid,
  k = 2L,
  arg_name = "truncation_grid",
  call = rlang::caller_env()
) {
  if (!is.numeric(grid)) {
    grid_class <- class(grid)[1]
    abort(
      "{.arg {arg_name}} must be numeric, not a {.cls {grid_class}}.",
      error_class = "positively_type_error",
      call = call
    )
  }
  if (length(grid) == 0) {
    abort(
      "{.arg {arg_name}} must contain at least one lower bound.",
      error_class = "positively_empty_error",
      call = call
    )
  }
  if (anyNA(grid)) {
    abort(
      "{.arg {arg_name}} must not contain missing values.",
      error_class = "positively_missing_error",
      call = call
    )
  }
  ceiling <- 1 / k
  if (any(grid < 0 | grid >= ceiling)) {
    bad <- grid[grid < 0 | grid >= ceiling]
    abort(
      c(
        "{.arg {arg_name}} lower bounds must lie in {.code [0, {signif(ceiling, 3)})}.",
        i = "The ceiling is {.code 1/k} for an exposure with {k} level{?s}.",
        x = "Found {.val {bad}}."
      ),
      error_class = "positively_range_error",
      call = call
    )
  }
  invisible(grid)
}
