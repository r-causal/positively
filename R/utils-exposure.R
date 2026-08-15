# Exposure-type detection is causalgenerics'. Reading the type from there is
# what makes a column classify the same way in every r-causal package. The
# announcement stays here, because it goes through alert_info() and so
# options(positively.quiet) silences it along with every other informational
# message the package emits; causalgenerics is asked to stay silent, so one
# detection reports itself once.
#
# The resolution layer built on top of detection deliberately diverges.
# Detection supplies a default here, never an override: an explicit
# exposure_type wins outright, and resolve_exposure_type() then validates the
# resolved type against the column on both paths. An explicit type therefore
# fails only when that type's math cannot run on the data as given, never
# because the unique-value heuristic in causalgenerics::is_categorical() read
# the column differently than the user did. Under "auto" the heuristic is the
# only source of a type, so a detected type outside the diagnostic's supported
# set is still an abort; that gate is deliberate, and removing it would
# silently run a diagnostic on a column no one declared it could handle.
#
# Whether the reading is announced diverges too. The announcement reports which
# branch a call took, so it is worth making only for a diagnostic that supports
# more than one exposure type. A diagnostic supporting one type reports nothing
# by resolving successfully, and an unsupported reading is named in the abort
# either way.

#' Detect the type of an exposure vector
#'
#' Classifies an exposure as `"binary"`, `"categorical"`, or `"continuous"`
#' through [causalgenerics::detect_exposure_type()]. A vector with exactly two
#' observed values is binary; a factor or character vector with more than two is
#' categorical; a numeric vector is categorical when its share of unique values
#' is small, and continuous otherwise. Missing values are not levels on either
#' side of that heuristic, so a 0/1 column with an `NA` is binary.
#'
#' The announcement is made here rather than by causalgenerics, which reads an
#' option of its own, so that the message a user silences with
#' `options(positively.quiet)` covers this one too. causalgenerics is asked to
#' stay silent, leaving exactly one announcement per detection.
#'
#' @param .exposure The exposure vector.
#' @param arg The argument name the announcement points at. Defaults to
#'   `".exposure"`; a diagnostic whose exposure argument is named otherwise
#'   passes its own, so the announcement never names an argument the calling
#'   function does not have.
#' @param announce Whether to announce the detected type through `alert_info()`.
#'   Defaults to `TRUE`; a caller passes `FALSE` when the reading is not worth
#'   reporting, such as one detecting several exposures at once.
#'
#' @return A single string: `"binary"`, `"categorical"`, or `"continuous"`.
#' @keywords internal
#' @noRd
detect_exposure_type <- function(
  .exposure,
  arg = ".exposure",
  announce = TRUE
) {
  exposure_type <- causalgenerics::detect_exposure_type(
    .exposure,
    arg = arg,
    announce = FALSE
  )

  if (announce) {
    alert_info("Treating {.arg {arg}} as {exposure_type}")
  }

  exposure_type
}

#' Count the distinct values an exposure vector takes
#'
#' Missing values are counted out so that a binary column with NAs reports the
#' missingness through the diagnostic's own check rather than as a spurious
#' third level here.
#'
#' @param .exposure The exposure vector.
#'
#' @return A single integer.
#' @keywords internal
#' @noRd
n_exposure_levels <- function(.exposure) {
  length(unique(.exposure[!is.na(.exposure)]))
}

#' Can an exposure vector carry an exposure type?
#'
#' The structural requirements an exposure type places on the column itself,
#' as a single predicate. A continuous exposure must be numeric; `is.numeric()`
#' rather than coercibility, because `as.double(factor(c("a", "b", "c")))`
#' silently returns the level codes `1 2 3` and would admit a factor as a dose.
#' A binary exposure must have exactly two distinct non-missing values. A
#' categorical exposure imposes nothing: the diagnostics that need more than one
#' level check for it themselves.
#'
#' Both the gate that rejects a type the column cannot carry and the advice that
#' names a type worth switching to read the rule from here, so every type the
#' advice offers is one the gate accepts.
#'
#' @param .exposure The exposure vector.
#' @param exposure_type The exposure type to test.
#'
#' @return A single logical value.
#' @keywords internal
#' @noRd
exposure_carries_type <- function(.exposure, exposure_type) {
  switch(
    exposure_type,
    continuous = is.numeric(.exposure),
    binary = n_exposure_levels(.exposure) == 2,
    TRUE
  )
}

#' The exposure types each diagnostic supports
#'
#' The one place that records what each diagnostic can compute, keyed by suffix
#' so that entry `x` belongs to `check_x()`. Every diagnostic that takes an
#' `exposure_type` argument has an entry; `check_density_ratios()` has none,
#' because it takes no exposure type.
#'
#' Two statements have to agree for every diagnostic: the `exposure_type`
#' formal, which is the menu the user is offered and what [rlang::arg_match()]
#' matches against, and the `supported` set handed to `resolve_exposure_type()`,
#' which is what actually gates the call. Reading the second from here leaves
#' the formal as the only one of the two that can drift. The formal has to stay
#' a literal, because roxygen2 derives the usage block by deparsing it and a
#' call to an unexported function would be printed there instead.
#'
#' @return A named list of character vectors of exposure types.
#' @keywords internal
#' @noRd
diagnostic_supported_types <- function() {
  list(
    edp = c("binary", "categorical", "continuous"),
    port = c("binary", "categorical", "continuous"),
    hat_values = "continuous",
    hdr = "continuous",
    extrapolation = "binary",
    eta_bias = c("binary", "categorical"),
    port_seq = c("binary", "categorical", "continuous"),
    hdr_seq = "continuous",
    positivity = c("binary", "categorical", "continuous")
  )
}

#' Validate that an exposure vector can carry its resolved type
#'
#' Applies `exposure_carries_type()` and explains a failure. The rule itself
#' lives in the predicate; only the wording of the abort is chosen here.
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
  if (exposure_carries_type(.exposure, exposure_type)) {
    return(invisible(.exposure))
  }

  # Only the continuous and binary rules in exposure_carries_type() can fail, so
  # the type alone selects which of the two explains the failure.
  if (exposure_type == "continuous") {
    abort(
      c(
        "{.fn {fn}} needs a numeric {.arg {arg}} for a continuous exposure.",
        x = "{.arg {arg}} is {.cls {class(.exposure)}}."
      ),
      error_class = "positively_exposure_type_error",
      call = call
    )
  }

  n_levels <- n_exposure_levels(.exposure)
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

#' Decide an exposure type without raising
#'
#' The deciding half of `resolve_exposure_type()`. It asks the same questions in
#' the same order, and returns what it found rather than aborting on it. A
#' caller resolving one column hands the answer straight to the reporting half;
#' a caller resolving several collects every answer first and reports once, so a
#' user with three offending exposure columns is told about all three instead of
#' correcting one and rerunning to meet the next.
#'
#' The menu gate is not applied here. `rlang::arg_match()` is a statement about
#' the argument rather than about any one column, so a caller resolving several
#' columns matches once before its loop and a caller resolving one matches in
#' `resolve_exposure_type()`.
#'
#' @param exposure_type `"auto"` or one of `supported`, already matched against
#'   the menu.
#' @param .exposure The exposure vector.
#' @param supported The exposure types the calling diagnostic can compute.
#' @param arg The argument name the announcement points at.
#' @param announce Whether a detected type is announced through `alert_info()`.
#'   Defaults to `TRUE`; the caller decides, and one classifying several
#'   exposures at once passes `FALSE` to stay silent.
#'
#' @return A list with `type`, the resolved exposure type; `detected`, the type
#'   detection read from the column, or `NA_character_` when a type was
#'   declared; and `problem`, one of `"none"`, `"unsupported"` for a detected
#'   type outside `supported`, and `"structure"` for a resolved type the column
#'   cannot carry.
#' @keywords internal
#' @noRd
classify_exposure_type <- function(
  exposure_type,
  .exposure,
  supported,
  arg = ".exposure",
  announce = TRUE
) {
  detected <- NA_character_
  type <- exposure_type

  if (exposure_type == "auto") {
    detected <- detect_exposure_type(.exposure, arg = arg, announce = announce)
    type <- detected
    if (!detected %in% supported) {
      return(list(type = type, detected = detected, problem = "unsupported"))
    }
  }

  problem <- if (exposure_carries_type(.exposure, type)) "none" else "structure"
  list(type = type, detected = detected, problem = problem)
}

#' Resolve an exposure type, detecting it when requested
#'
#' Matches `exposure_type` against `"auto"` and the types the calling diagnostic
#' supports with [rlang::arg_match()], then reports whatever
#' `classify_exposure_type()` decided. When it is `"auto"`, the type is inferred
#' from the data and must land in `supported`; otherwise the supplied type is
#' honored without detection. Either way the resolved type is checked against
#' the column, so a detected type and a declared type fail for the same reasons.
#'
#' @param exposure_type `"auto"` or one of `supported`.
#' @param .exposure The exposure vector.
#' @param supported The exposure types the calling diagnostic can compute.
#' @param fn The name of the calling diagnostic, used in error messages.
#' @param arg The argument name used in error messages.
#' @param announce Whether a detected type may be announced through
#'   `alert_info()`. Defaults to `TRUE`; callers that resolve several exposures
#'   at once pass `FALSE` to stay silent. A single-element `supported` silences
#'   the announcement regardless, since the reading names no branch a caller
#'   supporting one type could have taken differently.
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
  announce = TRUE,
  call = rlang::caller_env()
) {
  exposure_type <- resolve_arg_match(
    rlang::arg_match(
      exposure_type,
      values = c("auto", supported),
      error_arg = "exposure_type",
      error_call = call
    ),
    call = call
  )

  classified <- classify_exposure_type(
    exposure_type,
    .exposure,
    supported = supported,
    arg = arg,
    announce = announce && length(supported) > 1
  )

  if (classified$problem == "unsupported") {
    detected <- classified$detected
    type_code <- paste0("exposure_type = \"", supported, "\"")
    abort(
      c(
        "{.fn {fn}} supports {.or {supported}} exposures only.",
        i = "{.arg {arg}} was detected as {.val {detected}}.",
        i = "If {.arg {arg}} is {.or {supported}}, set
             {.or {.code {type_code}}}."
      ),
      error_class = "positively_exposure_type_error",
      call = call
    )
  }

  # The classifier has already applied exposure_carries_type(), so the gate is
  # entered only to raise; the wording of that failure belongs to it.
  if (classified$problem == "structure") {
    validate_exposure_structure(
      .exposure,
      classified$type,
      fn = fn,
      arg = arg,
      call = call
    )
  }

  classified$type
}
