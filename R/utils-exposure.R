# Exposure-type detection. positively mirrors the propensity idiom rather than
# importing it: has_two_levels(), is_categorical(), and detect_exposure_type()
# classify a column exactly as propensity does, down to the informational
# message, so that users experience one r-causal ecosystem.
#
# The resolution layer built on top of detection deliberately diverges.
# Detection supplies a default here, never an override: an explicit
# exposure_type wins outright, and resolve_exposure_type() then validates the
# resolved type against the column on both paths. An explicit type therefore
# fails only when that type's math cannot run on the data as given, never
# because the unique-value heuristic in is_categorical() read the column
# differently than the user did. Under "auto" the heuristic is the only source
# of a type, so a detected type outside the diagnostic's supported set is still
# an abort; that gate is deliberate, and removing it would silently run a
# diagnostic on a column no one declared it could handle.

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
#' @param announce Whether to announce the detected type through `alert_info()`.
#'   Defaults to `TRUE`; callers that detect several exposures at once pass
#'   `FALSE` to stay silent.
#'
#' @return A single string: `"binary"`, `"categorical"`, or `"continuous"`.
#' @keywords internal
#' @noRd
detect_exposure_type <- function(.exposure, announce = TRUE) {
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

  if (announce) {
    alert_info("Treating {.arg .exposure} as {exposure_type}")
  }

  exposure_type
}

#' Validate that an exposure vector can carry its resolved type
#'
#' Checks the structural requirements a resolved exposure type places on the
#' column itself. A continuous exposure must be numeric; `is.numeric()` rather
#' than coercibility, because `as.double(factor(c("a", "b", "c")))` silently
#' returns the level codes `1 2 3` and would admit a factor as a dose. A binary
#' exposure must have exactly two distinct non-missing values. A categorical
#' exposure imposes nothing: the diagnostics that need more than one level check
#' for it themselves.
#'
#' @param .exposure The exposure vector.
#' @param exposure_type The resolved exposure type.
#' @param fn The name of the calling diagnostic, used in error messages.
#' @param arg The argument name used in error messages.
#' @param call The calling environment, used to build the error's call.
#'
#' @return `.exposure`, invisibly, when it can carry `exposure_type`.
#' @keywords internal
#' @noRd
validate_exposure_structure <- function(
  .exposure,
  exposure_type,
  fn,
  arg = ".exposure",
  call = rlang::caller_env()
) {
  if (exposure_type == "continuous" && !is.numeric(.exposure)) {
    abort(
      c(
        "{.fn {fn}} needs a numeric {.arg {arg}} for a continuous exposure.",
        x = "{.arg {arg}} is {.cls {class(.exposure)}}."
      ),
      error_class = "positively_exposure_type_error",
      call = call
    )
  }

  # Missing values are counted out so that a binary column with NAs reports the
  # missingness through the diagnostic's own check rather than as a spurious
  # third level here.
  if (exposure_type == "binary") {
    n_levels <- length(unique(.exposure[!is.na(.exposure)]))
    if (n_levels != 2) {
      abort(
        c(
          "{.fn {fn}} needs exactly two distinct values in {.arg {arg}} for a
           binary exposure.",
          x = "{.arg {arg}} has {n_levels} distinct value{?s}."
        ),
        error_class = "positively_exposure_type_error",
        call = call
      )
    }
  }

  invisible(.exposure)
}

#' Resolve an exposure type, detecting it when requested
#'
#' Matches `exposure_type` against `"auto"` and the types the calling diagnostic
#' supports with [rlang::arg_match()]. When it is `"auto"`, the type is inferred
#' from the data through `detect_exposure_type()` and must land in `supported`;
#' otherwise the supplied type is honored without detection. Either way the
#' resolved type is checked against the column through
#' `validate_exposure_structure()`, so a detected type and a declared type fail
#' for the same reasons.
#'
#' @param exposure_type `"auto"` or one of `supported`.
#' @param .exposure The exposure vector.
#' @param supported The exposure types the calling diagnostic can compute.
#' @param fn The name of the calling diagnostic, used in error messages.
#' @param arg The argument name used in error messages.
#' @param call The calling environment, used to build the error's call.
#'
#' @return A single string: `"binary"`, `"categorical"`, or `"continuous"`.
#' @keywords internal
#' @noRd
resolve_exposure_type <- function(
  exposure_type,
  .exposure,
  supported,
  fn,
  arg = ".exposure",
  call = rlang::caller_env()
) {
  exposure_type <- rlang::arg_match(
    exposure_type,
    values = c("auto", supported),
    error_arg = "exposure_type",
    error_call = call
  )

  if (exposure_type == "auto") {
    exposure_type <- detect_exposure_type(.exposure)
    if (!exposure_type %in% supported) {
      type_code <- paste0("exposure_type = \"", supported, "\"")
      abort(
        c(
          "{.fn {fn}} supports {.or {supported}} exposures only.",
          i = "{.arg {arg}} was detected as {.val {exposure_type}}.",
          i = "If {.arg {arg}} is {.or {supported}}, set
               {.or {.code {type_code}}}."
        ),
        error_class = "positively_exposure_type_error",
        call = call
      )
    }
  }

  validate_exposure_structure(
    .exposure,
    exposure_type,
    fn = fn,
    arg = arg,
    call = call
  )

  exposure_type
}
