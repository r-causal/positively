# S7 methods defined here are wired up at load time by S7::methods_register().

#' The HDR non-overlap diagnostic result class
#'
#' `hdr_result` is the S7 class returned by [check_hdr()] and [check_hdr_seq()].
#' It extends [positivity_diagnostic] with the HDR probability mass and a label
#' for the conditional-density estimator. It is created internally and is not
#' constructed directly by users.
#'
#' @usage NULL
#' @keywords internal
#' @noRd
hdr_result <- new_class(
  "hdr_result",
  parent = positivity_diagnostic,
  properties = list(
    mass = class_double,
    density_estimator = class_character
  ),
  validator = function(self) {
    if (length(self@mass) != 1) {
      "@mass must be a single value"
    } else if (length(self@density_estimator) != 1) {
      "@density_estimator must be a single label"
    }
  }
)

#' Diagnose positivity for continuous exposures with HDR non-overlap
#'
#' `check_hdr()` computes the highest-density-region (HDR) non-overlap ratio of
#' Bao and Schomaker (2025) for continuous exposures. For a common target
#' exposure value, it reports the fraction of covariate profiles whose
#' highest-density region excludes that value, that is, the share of the
#' population for which the target dose is not supported.
#'
#' @details
#' For a continuous exposure \eqn{A} with covariates \eqn{L}, fix a probability
#' mass `mass` (the paper's \eqn{\alpha}). At covariate profile \eqn{l} the HDR
#' is the smallest set of exposure values capturing that mass of the conditional
#' density,
#' \deqn{A_\alpha(l) = \{ a : f(a \mid l) \ge f_\alpha(l) \},
#' \qquad P(A \in A_\alpha(l) \mid L = l) = \mathrm{mass}.}
#' The non-overlap ratio at a target `a` is the fraction of profiles whose HDR
#' excludes it,
#' \deqn{\hat{\tau}(a) = \frac{1}{n} \sum_{j=1}^{n}
#' \mathbf{1}\{ a \notin A_\alpha(l_j) \},}
#' a value in \eqn{[0, 1]}: zero when `a` is supported everywhere, one when it is
#' supported nowhere.
#'
#' The conditional density is supplied by `density_estimator`. The default
#' [hdr_density_normal()] fits `lm(exposure ~ covariates)`, treats the density
#' as Gaussian with the residual standard deviation \eqn{\hat{\sigma}}, and uses
#' the closed-form cutoff \eqn{\mathrm{dnorm}(z) / \hat{\sigma}} with
#' \eqn{z = \Phi^{-1}((1 + \mathrm{mass}) / 2)}. Membership then reduces to the
#' interval test \eqn{|a - \hat{\mu}(l)| \le z\,\hat{\sigma}}, so
#' \eqn{\hat{\tau}(a)} is the fraction of fitted means more than
#' \eqn{z\,\hat{\sigma}} from `a`.
#'
#' The non-overlap ratio is a population-level, common-target quantity. It is not
#' small merely because the density model fits well: when the exposure is
#' strongly covariate-driven, even the best central target has a nonzero floor,
#' because setting the whole population to one common dose is infeasible. The default normal estimator detects mean-shift support gaps,
#' where a stratum's supported dose moves away from a target. It does not detect
#' multimodal gaps: a hole between two modes of the true conditional density is
#' filled by the fitted normal and reported as supported. Supply a flexible
#' estimator through [new_hdr_density()] when multimodal structure is expected.
#'
#' When the exposure is close to a deterministic function of the covariates, the
#' fitted conditional density is narrow enough that no profile's region reaches
#' another profile's dose, and the ratio is one everywhere inside the observed
#' exposure range. That reading is literally correct, since a dose fixed by the
#' covariates is supported at no other profile, but it is set by the width of the
#' fitted density rather than by any comparison between profiles, so
#' `check_hdr()` warns when it happens. The warning is decided on a dense grid
#' spanning the observed range, not on `values`, so it reports a property of the
#' fit and does not depend on which doses were asked about: a ratio of one at a
#' few probes aimed into a genuine support gap is a finding rather than an
#' artifact, and stays unremarked.
#'
#' @param .data A data frame.
#' @param .exposure The continuous exposure column, selected with data-masking.
#'   Under the default `exposure_type = "auto"`, `check_hdr()` aborts when the
#'   exposure is detected as binary or categorical; declaring
#'   `exposure_type = "continuous"` accepts any numeric column that takes more
#'   than one value.
#' @param .covariates The covariate columns, selected with tidyselect.
#' @param mass The HDR probability mass, a single value strictly between 0 and 1.
#'   Defaults to `0.95`.
#' @param values A numeric vector of target exposure values at which to evaluate
#'   the non-overlap ratio. `NULL` (the default) uses a 100-point grid spanning
#'   the observed exposure range.
#' @param density_estimator An `hdr_density` conditional-density estimator, built
#'   with [hdr_density_normal()] (the default) or [new_hdr_density()].
#' @param exposure_type One of `"auto"` (detect from the data, the default) or
#'   `"continuous"`. A supplied type is authoritative and detection is not
#'   consulted, so `"continuous"` is rejected only when the exposure column is
#'   not numeric. Declaring it on a numeric column with few distinct values, a
#'   dose recorded at a handful of milligram levels for instance, is exactly the
#'   supported use: the unique-value heuristic reads such a column as
#'   categorical, so under `"auto"` the call aborts. Because numeric is the whole
#'   of what the type asks of the column, a two-valued numeric column is then
#'   accepted as well, which is the price of an authoritative declaration. A
#'   constant column is refused whichever way the type was decided, since the
#'   ratio compares covariate profiles at a common target dose and one observed
#'   dose leaves nothing to compare.
#'
#' @return An `hdr_result` object, an S7 subclass of [positivity_diagnostic].
#'   Its `@results` tibble has one row per target value with columns `value` (the
#'   target `a`) and `nonoverlap` (the ratio \eqn{\hat{\tau}(a)}). It also
#'   carries the properties `@mass` and `@density_estimator`.
#'
#'   [generics::glance()] returns a one-row tibble with `n` (the sample size),
#'   `mass`, `density_estimator`, `n_values` (the number of target values), and
#'   `nonoverlap_min` and `nonoverlap_max`, the range of the ratio across those
#'   targets.
#'
#' @references
#' Bao Y, Schomaker M (2025). Feasible Dose-Response Curves for Continuous
#' Treatments Under Positivity Violations.
#'
#' @examples
#' set.seed(1)
#' n <- 300
#' l <- rnorm(n)
#' dose <- rnorm(n, mean = l)
#' df <- data.frame(dose = dose, l = l)
#'
#' result <- check_hdr(df, dose, l, values = c(-2, 0, 2))
#' result
#'
#' @export
check_hdr <- function(
  .data,
  .exposure,
  .covariates,
  mass = 0.95,
  values = NULL,
  density_estimator = hdr_density_normal(),
  exposure_type = c("auto", "continuous")
) {
  validate_data_frame(.data)
  validate_probability(mass, arg_name = "mass")
  validate_hdr_estimator(density_estimator)

  exposure_name <- select_single_exposure(rlang::enquo(.exposure), .data)
  exposure_vec <- .data[[exposure_name]]

  exposure_type <- resolve_exposure_type(
    exposure_type,
    exposure_vec,
    supported = diagnostic_supported_types()[["hdr"]],
    fn = "check_hdr"
  )

  covariate_pos <- eval_select_columns(
    rlang::enquo(.covariates),
    .data,
    ".covariates"
  )
  validate_column_selection(covariate_pos, ".covariates")
  covariate_names <- names(covariate_pos)

  # Conditioning the exposure on itself makes the overlap check degenerate. The
  # comparison is by position so a renamed exposure inside the selection is also
  # caught.
  exposure_pos <- match(exposure_name, names(.data))
  if (length(intersect(covariate_pos, exposure_pos)) > 0) {
    abort(
      c(
        "{.arg .covariates} must not include the exposure column.",
        x = "{.val {exposure_name}} is also selected by {.arg .exposure}."
      ),
      error_class = "positively_selection_error"
    )
  }

  validate_numeric_columns(.data, exposure_name, ".exposure")
  validate_numeric_columns(.data, covariate_names, ".covariates")

  n <- nrow(.data)
  if (n < 2) {
    abort(
      "{.arg .data} must have at least two observations, not {n}.",
      error_class = "positively_size_error"
    )
  }

  validate_exposure_variation(
    .data,
    exposure_name,
    ".exposure",
    fn = "check_hdr"
  )

  targets <- resolve_targets(values, exposure_vec)

  analysis <- .data[c(exposure_name, covariate_names)]
  exposure <- as.double(exposure_vec)
  fitted <- hdr_fit(
    density_estimator,
    exposure_name = exposure_name,
    covariate_names = covariate_names,
    analysis = analysis,
    exposure = exposure,
    mass = mass
  )
  nonoverlap <- hdr_ratios(density_estimator, fitted, analysis, targets)
  if (
    nonoverlap_is_flat_at_one(density_estimator, fitted, analysis, exposure)
  ) {
    warn_flat_nonoverlap()
  }

  results <- tibble::tibble(value = targets, nonoverlap = nonoverlap)

  hdr_result(
    results = results,
    exposure = exposure_name,
    exposure_type = exposure_type,
    n = as.integer(n),
    params = list(mass = mass, values = values),
    call = rlang::current_call(),
    mass = mass,
    density_estimator = density_estimator@label
  )
}

# ---- Shared internals -----------------------------------------------------

#' Validate that a density estimator is an `hdr_density`
#'
#' @param density_estimator The object supplied to `density_estimator`.
#' @param call The calling environment, used to build the error's call.
#'
#' @return `density_estimator`, invisibly, when it is an `hdr_density`.
#' @keywords internal
#' @noRd
validate_hdr_estimator <- function(
  density_estimator,
  call = rlang::caller_env()
) {
  if (!S7::S7_inherits(density_estimator, hdr_density)) {
    supplied <- class(density_estimator)[1]
    abort(
      c(
        "{.arg density_estimator} must be an {.cls hdr_density} object.",
        i = "Build one with {.fn hdr_density_normal} or {.fn new_hdr_density}.",
        x = "Received a {.cls {supplied}}."
      ),
      error_class = "positively_type_error",
      call = call
    )
  }
  invisible(density_estimator)
}

#' Validate that every named exposure column varies
#'
#' A constant exposure is a precondition failure of the HDR diagnostic rather
#' than a property the exposure type is missing. `"continuous"` says the column
#' carries a dose the density model can be fitted to, and a constant numeric
#' column carries one; what it cannot do is support the comparison the ratio
#' makes. Folding the rule into `exposure_carries_type()` would therefore refuse
#' a declared `"continuous"` in [check_edp()] and [check_port()], which
#' categorize the exposure and ask nothing of its spread, and would pre-empt the
#' rank check in [check_hat_values()] with a vaguer message.
#'
#' Aborting rather than warning is what the returned number requires. With one
#' distinct value `resolve_targets()` returns a single target instead of the
#' documented 100-point grid, the fitted residual standard deviation is zero or
#' within rounding of it, and the ratio comes back as zero, which reads as
#' perfectly supported: not an imprecise answer but the opposite of the truth.
#' The sibling diagnostics already hold this line, [check_hat_values()] through
#' its rank check and [check_extrapolation()] and [check_eta_bias()] through
#' their exposure-type checks.
#'
#' @param .data The data frame.
#' @param columns The exposure column names.
#' @param arg_name The argument name used in error messages.
#' @param fn The name of the calling diagnostic, used in error messages.
#' @param call The calling environment, used to build the error's call.
#'
#' @return `.data`, invisibly, when every named column takes more than one
#'   value.
#' @keywords internal
#' @noRd
validate_exposure_variation <- function(
  .data,
  columns,
  arg_name,
  fn,
  call = rlang::caller_env()
) {
  constant <- columns[vapply(
    columns,
    function(column) n_exposure_levels(.data[[column]]) < 2,
    logical(1)
  )]
  if (length(constant) == 0) {
    return(invisible(.data))
  }

  abort(
    c(
      "{.fn {fn}} needs {.arg {arg_name}} to take more than one value.",
      x = "{.arg {constant}} {?is/are} constant.",
      i = "The non-overlap ratio compares covariate profiles at a common target
           dose, and the target grid spans the observed exposure range. One
           observed dose leaves one target and a fitted density with no width,
           so the ratio would report that dose as supported at every profile."
    ),
    error_class = "positively_variance_error",
    call = call
  )
}

#' A dense grid of probe values strictly inside the observed exposure range
#'
#' The grid the flat-at-one test is read on. It is built here rather than taken
#' from `values` because how many targets a user asked for, and where, is not
#' something the test can be allowed to depend on: three probes aimed into a
#' genuine support gap all read one, and a criterion reading them would call a
#' real finding an artifact, while the same fit judged on a dense grid does not.
#' Generating the grid makes the reading a property of the fit.
#'
#' The endpoints are dropped. A ratio of one beyond the observed range is what a
#' deliberate probe at an unobserved dose is asking for, and each range endpoint
#' is a dose some profile actually received, so the ratio there is held below one
#' by that profile alone. Keeping them would leave the test unreachable.
#'
#' @param exposure The observed exposure vector, which varies by the time this
#'   is reached.
#' @param n_probe The number of interior probe values.
#'
#' @return A numeric vector of `n_probe` values strictly inside the range.
#' @keywords internal
#' @noRd
interior_probe_grid <- function(exposure, n_probe = 100L) {
  limits <- range(exposure)
  grid <- seq(limits[1], limits[2], length.out = n_probe + 2L)
  grid[-c(1L, length(grid))]
}

#' Is the non-overlap curve flat at one inside the observed exposure range?
#'
#' Reads the degeneracy off the estimator's output rather than its internals, so
#' the test holds for any estimator built through [new_hdr_density()] and not
#' only for the fitted residual standard deviation of the default one. The
#' fitted state is reused, so the extra ratios cost no second fit. Exact equality
#' is safe because each ratio is the `mean()` of a logical vector.
#'
#' A missing ratio is not a flat reading. A saturated or over-parameterized fit
#' has no residual standard deviation to speak of, which makes every ratio `NA`,
#' and a comparison over those is `NA` rather than `FALSE`. Returning it would
#' put `NA` where the caller branches and raise an unclassed base R error.
#'
#' @param density_estimator An `hdr_density` estimator.
#' @param fitted The fitted state and per-row cutoff from `hdr_fit()`.
#' @param analysis A data frame holding the exposure and covariate columns.
#' @param exposure The observed exposure vector as a double.
#'
#' @return A single logical value, never `NA`.
#' @keywords internal
#' @noRd
nonoverlap_is_flat_at_one <- function(
  density_estimator,
  fitted,
  analysis,
  exposure
) {
  probes <- interior_probe_grid(exposure)
  ratios <- hdr_ratios(density_estimator, fitted, analysis, probes)
  length(ratios) > 0 && !anyNA(ratios) && all(ratios == 1)
}

#' Report a non-overlap curve flat at one inside the observed range
#'
#' The exposure varies, so the reading is not wrong and the diagnostic returns
#' it. What the warning adds is where the number came from: the fitted density
#' is too narrow for any profile's region to reach another profile's dose, so
#' the ratio is set by that width rather than by a comparison between profiles.
#'
#' @param times The offending time points, or `NULL` for the point diagnostic.
#' @param call The calling environment, used to build the warning's call.
#'
#' @return Invisibly `NULL`; raises a warning as a side effect.
#' @keywords internal
#' @noRd
warn_flat_nonoverlap <- function(times = NULL, call = rlang::caller_env()) {
  headline <- if (is.null(times)) {
    "The non-overlap ratio is 1 throughout the observed exposure range."
  } else {
    "The non-overlap ratio is 1 throughout the observed exposure range of
     {cli::qty(length(times))}time point{?s} {times}."
  }

  warn(
    c(
      headline,
      i = "The exposure is close to a deterministic function of the columns it
           is conditioned on, so no covariate profile's highest-density region
           reaches another profile's dose.",
      i = "The reading is correct, but it is set by the width of the fitted
           conditional density rather than by any comparison between profiles."
    ),
    warning_class = "positively_degenerate_density_warning",
    call = call
  )
}

#' Resolve the exposure type of every named exposure column
#'
#' Classifies each exposure column with `classify_exposure_type()`, so the
#' sequential variant holds the same continuous-only line as [check_hdr()] and
#' honors a declared type the same way, then reports every offender in one
#' abort. Collecting before reporting is what keeps a user with three offending
#' waves from correcting one column, rerunning, and being told about the next.
#'
#' The menu gate runs once, before any column is read. A type outside the
#' supported menu is a fact about the argument rather than about any one column,
#' so matching it first keeps it an argument error instead of one offender among
#' several.
#'
#' @param .data The data frame.
#' @param exposure_names The exposure column names.
#' @param exposure_type The `exposure_type` argument.
#' @param call The calling environment, used to build the error's call.
#'
#' @return `.data`, invisibly, when every exposure resolves to continuous.
#' @keywords internal
#' @noRd
resolve_seq_exposure_types <- function(
  .data,
  exposure_names,
  exposure_type,
  call = rlang::caller_env()
) {
  supported <- diagnostic_supported_types()[["hdr_seq"]]
  exposure_type <- resolve_arg_match(
    rlang::arg_match(
      exposure_type,
      values = c("auto", supported),
      error_arg = "exposure_type",
      error_call = call
    ),
    call = call
  )

  classified <- lapply(exposure_names, function(name) {
    classify_exposure_type(
      exposure_type,
      .data[[name]],
      supported = supported,
      arg = ".exposures",
      announce = FALSE
    )
  })
  problems <- vapply(classified, function(x) x$problem, character(1))
  unsupported <- problems == "unsupported"
  cannot_carry <- problems == "structure"

  # The two failures are reported in the order the single-column path raises
  # them, detection before structure, so the message a user meets first is the
  # same one either way.
  if (any(unsupported)) {
    abort_seq_detected_types(
      exposure_names[unsupported],
      vapply(classified[unsupported], function(x) x$detected, character(1)),
      supported = supported,
      call = call
    )
  }
  if (any(cannot_carry)) {
    abort_seq_exposure_structure(
      .data,
      exposure_names[cannot_carry],
      call = call
    )
  }

  invisible(.data)
}

#' Report every exposure column detection read outside the supported set
#'
#' One info bullet per distinct detected type, so columns read the same way are
#' named together and columns read differently are named apart, followed by a
#' single hint naming every offender. Detection is a guess, so the bullets stay
#' informational and the hint names the escape hatch; that keeps the
#' one-offender message the shape [check_hdr()] raises, differing only in the
#' function and the argument it names.
#'
#' How many bullets there are depends on the data, so the templates are built by
#' index. Interpolating a column name into a template instead would send user
#' text through cli twice, once here and once in `abort()`: a name containing
#' braces would be read as an expression, replacing the diagnostic error with an
#' internal cli failure and losing the package's own condition classes with it.
#'
#' @param offenders The exposure column names whose detected type is
#'   unsupported.
#' @param detected The detected type of each offender, aligned with `offenders`.
#' @param supported The exposure types `check_hdr_seq()` can compute.
#' @param call The calling environment, used to build the error's call.
#'
#' @return Never returns; always raises an error.
#' @keywords internal
#' @noRd
abort_seq_detected_types <- function(
  offenders,
  detected,
  supported,
  call = rlang::caller_env()
) {
  types <- unique(detected)
  named <- lapply(types, function(type) offenders[detected == type])
  type_code <- paste0("exposure_type = \"", supported, "\"")

  readings <- sprintf(
    "{.arg {named[[%d]]}} {?was/were} detected as {.val {types[[%d]]}}.",
    seq_along(named),
    seq_along(named)
  )
  hint <- "If {.arg {offenders}} {?is/are} {.or {supported}}, set
           {.or {.code {type_code}}}."

  bullets <- c(readings, hint)
  abort(
    c(
      "{.fn check_hdr_seq} supports {.or {supported}} exposures only.",
      rlang::set_names(bullets, rep("i", length(bullets)))
    ),
    error_class = "positively_exposure_type_error",
    call = call
  )
}

#' Report every exposure column that cannot carry the resolved type
#'
#' One danger bullet per distinct column class, matching the shape
#' `validate_exposure_structure()` raises for a single column. No escape hatch
#' is offered: reaching here under a declared type means the user has already
#' used it, and reaching here under detection means the column is not numeric
#' whatever it was read as. `check_hdr_seq()` supports continuous exposures
#' alone, so the numeric requirement is the only one that can fail and the
#' wording states it outright.
#'
#' As in `abort_seq_detected_types()`, the templates are built by index so that
#' a column name reaches cli once, as a value, rather than twice as text.
#'
#' @param .data The data frame.
#' @param offenders The exposure column names that cannot carry the type.
#' @param call The calling environment, used to build the error's call.
#'
#' @return Never returns; always raises an error.
#' @keywords internal
#' @noRd
abort_seq_exposure_structure <- function(
  .data,
  offenders,
  call = rlang::caller_env()
) {
  classes <- lapply(.data[offenders], class)
  keys <- vapply(classes, paste, character(1), collapse = "/")
  distinct <- unique(keys)
  named <- lapply(distinct, function(key) offenders[keys == key])
  shapes <- lapply(distinct, function(key) classes[[match(key, keys)]])

  bullets <- sprintf(
    "{.arg {named[[%d]]}} {?is/are} {.cls {shapes[[%d]]}}.",
    seq_along(named),
    seq_along(named)
  )

  abort(
    c(
      "{.fn check_hdr_seq} needs a numeric {.arg .exposures} for a continuous
       exposure.",
      rlang::set_names(bullets, rep("x", length(bullets)))
    ),
    error_class = "positively_exposure_type_error",
    call = call
  )
}

#' Resolve a data-masked selection to a single exposure column name
#'
#' @param exposure_quo A quosure capturing the `.exposure` selection.
#' @param .data The data frame.
#' @param call The calling environment, used to build the error's call.
#'
#' @return The name of the single selected exposure column.
#' @keywords internal
#' @noRd
select_single_exposure <- function(
  exposure_quo,
  .data,
  call = rlang::caller_env()
) {
  names(eval_select_column(exposure_quo, .data, ".exposure", call = call))
}

#' Resolve the target grid of exposure values
#'
#' Returns a validated user grid, or a default 100-point grid spanning the
#' observed exposure range when `values` is `NULL`.
#'
#' @param values The user-supplied `values`, or `NULL`.
#' @param exposure The observed exposure vector.
#' @param call The calling environment, used to build the error's call.
#'
#' @return A numeric vector of target values.
#' @keywords internal
#' @noRd
resolve_targets <- function(values, exposure, call = rlang::caller_env()) {
  if (is.null(values)) {
    return(unique(seq(min(exposure), max(exposure), length.out = 100)))
  }
  if (!is.numeric(values)) {
    values_class <- class(values)[1]
    abort(
      "{.arg values} must be numeric, not a {.cls {values_class}}.",
      error_class = "positively_type_error",
      call = call
    )
  }
  if (length(values) == 0) {
    abort(
      "{.arg values} must contain at least one value.",
      error_class = "positively_empty_error",
      call = call
    )
  }
  if (anyNA(values)) {
    abort(
      "{.arg values} must not contain missing values.",
      error_class = "positively_missing_error",
      call = call
    )
  }
  if (!all(is.finite(values))) {
    abort(
      "{.arg values} must be finite.",
      error_class = "positively_range_error",
      call = call
    )
  }
  as.double(values)
}

#' Fit the conditional-density estimator and resolve its HDR cutoff
#'
#' The fit and the per-row density cutoff are the whole of what the ratios are
#' read from, so they are produced once and handed to every reading taken off
#' them: the ratios at the user's target values, and the ratios on the interior
#' probe grid the flat-at-one test uses.
#'
#' @param density_estimator An `hdr_density` estimator.
#' @param exposure_name The exposure column name.
#' @param covariate_names The covariate column names.
#' @param analysis A data frame holding the exposure and covariate columns.
#' @param exposure The observed exposure vector as a double.
#' @param mass The HDR probability mass.
#' @param call The calling environment, used to build the error's call.
#'
#' @return A list with `state`, the estimator's fitted state, and `cutoff`, the
#'   per-row HDR density cutoff.
#' @keywords internal
#' @noRd
hdr_fit <- function(
  density_estimator,
  exposure_name,
  covariate_names,
  analysis,
  exposure,
  mass,
  call = rlang::caller_env()
) {
  formula <- stats::reformulate(covariate_names, response = exposure_name)
  state <- density_estimator@fit(formula, analysis)
  cutoff <- hdr_thresholds(
    density_estimator,
    state = state,
    newdata = analysis,
    mass = mass,
    exposure = exposure,
    call = call
  )
  list(state = state, cutoff = cutoff)
}

#' The HDR non-overlap ratio at each target value
#'
#' For every target value, the fraction of rows whose density at that value falls
#' below the cutoff, that is, whose HDR excludes it.
#'
#' @param density_estimator An `hdr_density` estimator.
#' @param fitted The fitted state and per-row cutoff from `hdr_fit()`.
#' @param analysis A data frame holding the exposure and covariate columns.
#' @param targets The target exposure values.
#'
#' @return A numeric vector of non-overlap ratios, one per target value.
#' @keywords internal
#' @noRd
hdr_ratios <- function(density_estimator, fitted, analysis, targets) {
  vapply(
    targets,
    function(a) {
      density <- density_estimator@density(fitted$state, a, analysis)
      mean(density < fitted$cutoff)
    },
    numeric(1)
  )
}

#' Diagnose sequential positivity for continuous exposures with HDR non-overlap
#'
#' `check_hdr_seq()` applies the HDR non-overlap ratio of Bao and Schomaker
#' (2025) to a time-varying continuous exposure, one time point at a time. At
#' each time point it fits the conditional density of that exposure on the
#' history entering its conditioning set, following the sequential-support logic
#' of sPoRT (Chatton et al.), and reports the non-overlap ratio at the target
#' values.
#'
#' @details
#' Data are supplied in wide form, one row per subject with one exposure column
#' per time point. For time point \eqn{t} the diagnostic fits the conditional
#' density of the exposure at \eqn{t} on its conditioning set, then computes the
#' point-in-time non-overlap ratio \eqn{\hat{\tau}_t(a)} exactly as
#' [check_hdr()] does. The conditioning set is the baseline covariates, the
#' time-varying covariates for each time point within the `lag` history window,
#' and the exposures at earlier time points within that window. With the default
#' `lag = Inf` the full history enters.
#'
#' Because each time point is treated separately, a support gap that moves a
#' stratum's supported dose away from a target is isolated to the time point
#' where it occurs. The default normal estimator sees mean-shift gaps only; a
#' multimodal gap at a single time point is invisible to it, as for
#' [check_hdr()].
#'
#' A time point whose exposure is close to a deterministic function of its
#' conditioning set reads a ratio of one throughout that time point's observed
#' exposure range, for the reason [check_hdr()] describes. Each time point is
#' judged on a dense grid spanning its own range, since the shared target grid
#' spans the pooled range and can leave a time point with a single target of its
#' own. One warning names every such time point.
#'
#' @param .data A data frame in wide form, one row per subject.
#' @param .exposures An ordered tidyselect of exposure columns, one per time
#'   point.
#' @param .covariates A list of tidyselect expressions, one per time point, of
#'   the time-varying covariates. A length-one list is recycled across time
#'   points. Write it as a literal `list()` call, for example
#'   `list(l0, l1, l2)`, rather than a pre-built list held in a variable.
#' @param .baseline A tidyselect of baseline covariates always included in every
#'   conditioning set. Defaults to `NULL`.
#' @param mass The HDR probability mass, a single value strictly between 0 and 1.
#'   Defaults to `0.95`.
#' @param values A numeric vector of target exposure values. `NULL` (the
#'   default) uses a 100-point grid spanning the pooled exposure range.
#' @param density_estimator An `hdr_density` conditional-density estimator, built
#'   with [hdr_density_normal()] (the default) or [new_hdr_density()].
#' @param lag The history window: the number of earlier time points whose
#'   covariates and exposures enter each conditioning set. Defaults to `Inf`, the
#'   full history.
#' @param exposure_type One of `"auto"` (detect from the data, the default) or
#'   `"continuous"`, applied to every exposure column in turn. A supplied type is
#'   authoritative and detection is not consulted, so `"continuous"` is rejected
#'   only when an exposure column is not numeric. Either way the error names
#'   every offending column at once, rather than the first one it meets.
#'   Declaring it on numeric columns with few distinct values, doses recorded at
#'   a handful of milligram levels for instance, is exactly the supported use:
#'   the unique-value heuristic reads such a column as categorical, so under
#'   `"auto"` the call aborts. Because numeric is the whole of what the type asks
#'   of a column, a two-valued numeric column is then accepted as well, which is
#'   the price of an authoritative declaration. A constant column is refused
#'   whichever way the type was decided, and the error names every constant
#'   column at once.
#'
#' @return An `hdr_result` object, an S7 subclass of [positivity_diagnostic].
#'   Its `@results` tibble has one row per time point and target value with
#'   columns `time`, `value`, and `nonoverlap`.
#'
#'   [generics::glance()] returns the columns [check_hdr()] reports plus
#'   `n_times`, the number of time points. The non-overlap range spans every
#'   wave rather than any one of them.
#'
#' @references
#' Bao Y, Schomaker M (2025). Feasible Dose-Response Curves for Continuous
#' Treatments Under Positivity Violations.
#'
#' Chatton A, Schomaker M, Luque-Fernandez MA, Platt RW, Schnitzer ME (2025).
#' Is checking for sequential positivity violations getting you down? Try
#' sPoRT!
#'
#' @examples
#' set.seed(1)
#' n <- 200
#' df <- data.frame(
#'   l0 = rnorm(n),
#'   a1 = rnorm(n),
#'   l1 = rnorm(n),
#'   a2 = rnorm(n),
#'   l2 = rnorm(n),
#'   a3 = rnorm(n)
#' )
#'
#' result <- check_hdr_seq(
#'   df,
#'   c(a1, a2, a3),
#'   list(l0, l1, l2),
#'   values = c(-2, 0, 2)
#' )
#' result
#'
#' @export
check_hdr_seq <- function(
  .data,
  .exposures,
  .covariates,
  .baseline = NULL,
  mass = 0.95,
  values = NULL,
  density_estimator = hdr_density_normal(),
  lag = Inf,
  exposure_type = c("auto", "continuous")
) {
  validate_data_frame(.data)
  validate_probability(mass, arg_name = "mass")
  validate_hdr_estimator(density_estimator)
  validate_lag(lag)

  exposure_pos <- eval_select_columns(
    rlang::enquo(.exposures),
    .data,
    ".exposures"
  )
  validate_column_selection(exposure_pos, ".exposures")
  exposure_names <- names(exposure_pos)
  n_times <- length(exposure_names)

  resolve_seq_exposure_types(.data, exposure_names, exposure_type)

  covariate_sets <- parse_covariate_list(
    rlang::enquo(.covariates),
    .data,
    n_times = n_times
  )
  baseline_names <- eval_optional_selection(rlang::enquo(.baseline), .data)

  validate_numeric_columns(.data, exposure_names, ".exposures")
  validate_numeric_columns(
    .data,
    unique(unlist(covariate_sets)),
    ".covariates"
  )
  if (length(baseline_names) > 0) {
    validate_numeric_columns(.data, baseline_names, ".baseline")
  }

  n <- nrow(.data)
  if (n < 2) {
    abort(
      "{.arg .data} must have at least two observations, not {n}.",
      error_class = "positively_size_error"
    )
  }

  validate_exposure_variation(
    .data,
    exposure_names,
    ".exposures",
    fn = "check_hdr_seq"
  )

  pooled_exposure <- as.double(unlist(.data[exposure_names], use.names = FALSE))
  targets <- resolve_targets(values, pooled_exposure)

  # Assemble the per-time conditioning sets outside the mapping loop so that a
  # selection error is reported against check_hdr_seq() rather than wrapped in
  # purrr's indexed context.
  conditioning_sets <- lapply(seq_len(n_times), function(t) {
    conditioning_set(
      time = t,
      covariate_sets = covariate_sets,
      exposure_names = exposure_names,
      baseline_names = baseline_names,
      lag = lag
    )
  })
  validate_conditioning_sets(
    conditioning_sets,
    exposure_names,
    call = rlang::current_env()
  )

  # Each wave is fitted on its own conditioning set, so its flat-at-one reading
  # is taken on a probe grid spanning that wave's own observed range. The shared
  # target grid spans the pooled range instead, which can leave a wave with one
  # target of its own, or none, too few points to tell a degenerate fit from a
  # genuine support gap.
  seq_call <- rlang::current_env()
  per_time <- purrr::map(seq_len(n_times), function(t) {
    conditioning <- conditioning_sets[[t]]
    exposure_name <- exposure_names[[t]]
    analysis <- .data[c(exposure_name, conditioning)]
    exposure <- as.double(.data[[exposure_name]])
    fitted <- hdr_fit(
      density_estimator,
      exposure_name = exposure_name,
      covariate_names = conditioning,
      analysis = analysis,
      exposure = exposure,
      mass = mass,
      call = seq_call
    )
    list(
      results = tibble::tibble(
        time = rep(t, length(targets)),
        value = targets,
        nonoverlap = hdr_ratios(density_estimator, fitted, analysis, targets)
      ),
      flat = nonoverlap_is_flat_at_one(
        density_estimator,
        fitted,
        analysis,
        exposure
      )
    )
  })
  results <- vctrs::vec_rbind(!!!lapply(per_time, `[[`, "results"))

  flat <- vapply(per_time, function(wave) wave$flat, logical(1))
  if (any(flat)) {
    warn_flat_nonoverlap(times = which(flat))
  }

  # resolve_seq_exposure_types() supports continuous exposures only, so every
  # column has already resolved to "continuous" and the type is stated here
  # rather than threaded back. Widening that supported set requires returning
  # the resolved types and carrying them through, or this property becomes a
  # silent lie.
  hdr_result(
    results = results,
    exposure = exposure_names,
    exposure_type = "continuous",
    n = as.integer(n),
    params = list(mass = mass, values = values, lag = lag),
    call = rlang::current_call(),
    mass = mass,
    density_estimator = density_estimator@label
  )
}

#' Parse a list of tidyselect expressions into per-time covariate names
#'
#' Accepts the `list(...)` form of `.covariates`, one selection per time point,
#' and returns a list of character vectors of resolved column names. A single
#' selection (length-one list or a bare selection) is recycled across all time
#' points.
#'
#' @param covariates_quo A quosure capturing `.covariates`.
#' @param .data The data frame.
#' @param n_times The number of time points.
#' @param call The calling environment, used to build the error's call.
#'
#' @return A list of length `n_times` of character vectors of column names.
#' @keywords internal
#' @noRd
parse_covariate_list <- function(
  covariates_quo,
  .data,
  n_times,
  call = rlang::caller_env()
) {
  expr <- rlang::quo_get_expr(covariates_quo)
  env <- rlang::quo_get_env(covariates_quo)
  is_list_call <- rlang::is_call(expr, "list")
  elements <- if (is_list_call) rlang::call_args(expr) else list(expr)

  sets <- lapply(elements, function(element) {
    element_quo <- rlang::new_quosure(element, env)
    selection <- eval_select_columns(
      element_quo,
      .data,
      ".covariates",
      call = call
    )
    names(selection)
  })

  if (length(sets) == 1 && n_times > 1) {
    sets <- rep(sets, n_times)
  }
  if (length(sets) != n_times) {
    abort(
      c(
        "{.arg .covariates} must supply one selection per time point.",
        i = "Found {length(sets)} selection{?s} for {n_times} exposure{?s}."
      ),
      error_class = "positively_selection_error",
      call = call
    )
  }
  sets
}

#' Evaluate an optional data-masked selection to column names
#'
#' @param selection_quo A quosure capturing the selection.
#' @param .data The data frame.
#' @param arg_name The argument name used in error messages.
#' @param call The calling environment, used to build the error's call.
#'
#' @return A character vector of column names, empty when the selection is
#'   `NULL`.
#' @keywords internal
#' @noRd
eval_optional_selection <- function(
  selection_quo,
  .data,
  arg_name = ".baseline",
  call = rlang::caller_env()
) {
  if (rlang::quo_is_null(selection_quo)) {
    return(character(0))
  }
  names(eval_select_columns(selection_quo, .data, arg_name, call = call))
}

#' The conditioning set for one time point
#'
#' Assembles the baseline covariates, the time-varying covariates for each time
#' point within the `lag` window, and the exposures at earlier time points
#' within that window.
#'
#' @param time The current time point index.
#' @param covariate_sets The per-time covariate name list.
#' @param exposure_names The exposure column names, time-ordered.
#' @param baseline_names The baseline covariate names.
#' @param lag The history window.
#'
#' @return A character vector of unique conditioning-set column names.
#' @keywords internal
#' @noRd
conditioning_set <- function(
  time,
  covariate_sets,
  exposure_names,
  baseline_names,
  lag
) {
  earliest <- if (is.infinite(lag)) 1L else max(1L, time - lag)
  covariate_names <- unlist(
    covariate_sets[earliest:time],
    use.names = FALSE
  )
  prior_exposures <- if (time > earliest) {
    exposure_names[earliest:(time - 1)]
  } else {
    character(0)
  }
  unique(c(baseline_names, covariate_names, prior_exposures))
}

#' Validate the assembled per-time conditioning sets
#'
#' Rejects a conditioning set that contains the exposure at its own time point
#' (which would model the response on itself) and a conditioning set left empty
#' after baseline covariates and lagged history are folded in. Prior exposures
#' at earlier times remain legitimate conditioning columns.
#'
#' @param conditioning_sets The per-time conditioning-set name list.
#' @param exposure_names The exposure column names, time-ordered.
#' @param call The calling environment, used to build the error's call.
#'
#' @return `conditioning_sets`, invisibly, when every set is valid.
#' @keywords internal
#' @noRd
validate_conditioning_sets <- function(
  conditioning_sets,
  exposure_names,
  call = rlang::caller_env()
) {
  for (t in seq_along(conditioning_sets)) {
    own_exposure <- exposure_names[[t]]
    if (own_exposure %in% conditioning_sets[[t]]) {
      abort(
        c(
          "{.arg .covariates} and {.arg .baseline} must not include an exposure at its own time point.",
          x = "{.val {own_exposure}} conditions its own model at time {t}."
        ),
        error_class = "positively_selection_error",
        call = call
      )
    }
    if (length(conditioning_sets[[t]]) == 0) {
      abort(
        c(
          "Every time point must have at least one conditioning column.",
          x = "The conditioning set at time {t} is empty.",
          i = "Supply time-varying covariates in {.arg .covariates} or baseline covariates in {.arg .baseline}."
        ),
        error_class = "positively_empty_error",
        call = call
      )
    }
  }
  invisible(conditioning_sets)
}

# ---- Methods --------------------------------------------------------------

method(glance, hdr_result) <- function(x, ...) {
  results <- x@results
  columns <- list(n = x@n)
  if ("time" %in% names(results)) {
    columns$n_times <- length(x@exposure)
  }
  columns$mass <- x@mass
  columns$density_estimator <- x@density_estimator
  columns$n_values <- length(unique(results$value))
  # A sequential run reports the same targets at every wave, so the range spans
  # the whole history rather than any one wave.
  columns$nonoverlap_min <- min(results$nonoverlap)
  columns$nonoverlap_max <- max(results$nonoverlap)
  tibble::as_tibble(columns)
}

method(diagnostic_label, hdr_result) <- function(x) {
  "HDR non-overlap"
}

method(diagnostic_headline, hdr_result) <- function(x) {
  nonoverlap <- x@results$nonoverlap
  n_target <- length(unique(x@results$value))
  mass <- signif(x@mass * 100, 3)
  span <- round(range(nonoverlap), 3)
  if ("time" %in% names(x@results)) {
    n_times <- length(unique(x@results$time))
    cli::format_inline(
      "{mass}% HDR, {x@density_estimator} density; non-overlap {span[[1]]} to {span[[2]]} over {n_target} target{?s} and {n_times} time point{?s}"
    )
  } else {
    cli::format_inline(
      "{mass}% HDR, {x@density_estimator} density; non-overlap {span[[1]]} to {span[[2]]} over {n_target} target{?s}"
    )
  }
}

method(print, hdr_result) <- function(x, ...) {
  nonoverlap <- x@results$nonoverlap
  n_target <- length(unique(x@results$value))
  cat_cli({
    cli::cli_h1("{diagnostic_label(x)}")
    cli::cli_text("Exposure: {.val {x@exposure}} ({x@exposure_type})")
    cli::cli_text("Observations: {x@n}")
    cli::cli_text("HDR mass: {x@mass}")
    cli::cli_text("Density estimator: {x@density_estimator}")
    if ("time" %in% names(x@results)) {
      cli::cli_text("Time points: {length(unique(x@results$time))}")
    }
    cli::cli_text(
      "Non-overlap over {n_target} target{?s}: {round(min(nonoverlap), 3)} to {round(max(nonoverlap), 3)}"
    )
  })
  invisible(x)
}
