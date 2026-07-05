# Exposure-type detection. positively mirrors the propensity idiom rather than
# importing it: the resolution logic, condition text, and informational message
# match so that users experience one r-causal ecosystem.

#' Does a vector have exactly two distinct values?
#'
#' @param x A vector.
#'
#' @return A single logical value.
#' @keywords internal
#' @noRd
has_two_levels <- function(x) {
  length(unique(x)) == 2
}

#' Is a numeric vector categorical by the unique-value heuristic?
#'
#' A numeric exposure is treated as categorical when the proportion of unique
#' values to non-missing observations falls below 20 percent.
#'
#' @param .exposure A numeric exposure vector.
#'
#' @return A single logical value.
#' @keywords internal
#' @noRd
is_categorical <- function(.exposure) {
  n_non_na <- sum(!is.na(.exposure))
  if (n_non_na == 0) {
    return(FALSE)
  }

  ratio <- length(unique(.exposure)) / n_non_na
  if (is.nan(ratio)) {
    return(FALSE)
  }

  ratio < 0.2
}

#' Detect the type of an exposure vector
#'
#' Classifies an exposure as `"binary"`, `"categorical"`, or `"continuous"`. A
#' vector with exactly two distinct values is binary; a factor or character
#' vector with more than two values is categorical; a numeric vector is
#' categorical when its share of unique values is small, and continuous
#' otherwise. The detected type is announced through `alert_info()` unless
#' `options(positively.quiet)` is `TRUE`.
#'
#' @param .exposure The exposure vector.
#'
#' @return A single string: `"binary"`, `"categorical"`, or `"continuous"`.
#' @keywords internal
#' @noRd
detect_exposure_type <- function(.exposure) {
  exposure_type <- if (has_two_levels(.exposure)) {
    "binary"
  } else if (is.factor(.exposure) || is.character(.exposure)) {
    if (length(unique(.exposure)) > 2) {
      "categorical"
    } else {
      "binary"
    }
  } else if (is_categorical(.exposure)) {
    "categorical"
  } else {
    "continuous"
  }

  alert_info("Treating {.arg .exposure} as {exposure_type}")

  exposure_type
}

#' Resolve an exposure type, detecting it when requested
#'
#' Matches `exposure_type` against the permitted values with
#' [rlang::arg_match()]. When it is `"auto"`, the type is inferred from the data
#' through `detect_exposure_type()`; otherwise the supplied type is honoured
#' without detection.
#'
#' @param exposure_type One of `"auto"`, `"binary"`, `"categorical"`, or
#'   `"continuous"`.
#' @param .exposure The exposure vector, used only when `exposure_type` is
#'   `"auto"`.
#' @param valid_types The permitted values for `exposure_type`.
#'
#' @return A single string: `"binary"`, `"categorical"`, or `"continuous"`.
#' @keywords internal
#' @noRd
match_exposure_type <- function(
  exposure_type = c("auto", "binary", "categorical", "continuous"),
  .exposure,
  valid_types = c("auto", "binary", "categorical", "continuous")
) {
  .exposure_type <- rlang::arg_match(exposure_type, valid_types)
  if (.exposure_type == "auto") {
    detect_exposure_type(.exposure)
  } else {
    .exposure_type
  }
}
