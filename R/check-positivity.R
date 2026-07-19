# check_positivity() is the top-level entry point. It resolves the exposure type
# once, announces it once, composes the applicable check_*() diagnostics with
# their defaults (overridable through `args`), and returns a positivity_check
# container. It renders no verdicts. The container's tidy(), glance(), and
# print() methods live here alongside it; the container class itself is defined
# in aaa-classes.R.

# The diagnostics check_positivity() can dispatch, in design-doc order.
# check_eta_bias() and check_density_ratios() are deliberately absent: they need
# an outcome or user-supplied ratios, so they are never composed here.
composed_diagnostics <- function() {
  c("edp", "port", "hat_values", "hdr", "extrapolation")
}

# Diagnostics check_positivity() knows about but will not compose, with the
# reason surfaced to users who request them.
excluded_diagnostics <- function() {
  c("eta_bias", "density_ratios")
}

# The exposure types each composable diagnostic accepts.
diagnostic_exposure_types <- function() {
  list(
    edp = c("binary", "categorical", "continuous"),
    port = c("binary", "categorical", "continuous"),
    hat_values = "continuous",
    hdr = "continuous",
    extrapolation = "binary"
  )
}

# The default diagnostic set for an exposure type, in design-doc order.
default_diagnostics <- function(exposure_type) {
  switch(
    exposure_type,
    binary = c("edp", "port", "extrapolation"),
    categorical = c("edp", "port"),
    continuous = c("edp", "port", "hat_values", "hdr")
  )
}

#' Diagnose positivity across the applicable methods
#'
#' `check_positivity()` is the general entry point to positively. It detects the
#' exposure type, runs the diagnostics that apply to that type with their
#' defaults, and returns a [positivity_check] container whose print method shows
#' each diagnostic's summary in its own section. Like every positively function,
#' it is diagnostic only: it reports numbers, not verdicts.
#'
#' @details
#' `check_positivity()` resolves the exposure type a single time and announces it
#' once through an informational message (suppressed by
#' `options(positively.quiet = TRUE)`). [check_edp()] and [check_port()] receive
#' the resolved type directly; [check_hat_values()], [check_hdr()], and
#' [check_extrapolation()] re-detect it internally, but their repeated
#' announcement is muffled so that the type is reported only once. Every other
#' informational message a diagnostic emits still reaches you. When `exposure_type`
#' is supplied explicitly, it is checked against the detected type up front and a
#' genuine disagreement aborts before any diagnostic runs. The default diagnostic
#' set depends on the exposure type:
#'
#' - **binary**: [check_edp()], [check_port()], [check_extrapolation()]
#' - **categorical**: [check_edp()], [check_port()]
#' - **continuous**: [check_edp()], [check_port()] (with the exposure
#'   categorized), [check_hat_values()], [check_hdr()]
#'
#' Pass `diagnostics` to run an explicit subset instead; the requested order is
#' honored. Every name must be a diagnostic that `check_positivity()` composes
#' and that applies to the resolved exposure type, otherwise the call aborts with
#' the valid options listed.
#'
#' [check_eta_bias()] and [check_density_ratios()] are never composed here.
#' The first needs an outcome and the second needs user-supplied density ratios,
#' so both are called directly rather than through `check_positivity()`.
#'
#' Per-method options are threaded through `args`, a named list of option lists,
#' for example `args = list(port = list(alpha = 0.1))`. Every name in `args` must
#' be a diagnostic that is actually being run.
#'
#' @param .data A data frame.
#' @param .exposure The exposure column, selected with data-masking.
#' @param .covariates The covariate columns, selected with tidyselect. Required,
#'   with no default, so that outcome columns are not swept in by accident. The
#'   exposure column must not be selected.
#' @param diagnostics The diagnostics to run. `NULL` (the default) selects the
#'   applicable set for the exposure type; otherwise a subset of `"edp"`,
#'   `"port"`, `"hat_values"`, `"hdr"`, and `"extrapolation"`, run in the order
#'   given.
#' @param exposure_type One of `"auto"` (detect from the data, the default),
#'   `"binary"`, `"categorical"`, or `"continuous"`. An explicit value is passed
#'   to [check_edp()] and [check_port()], which honor it directly. The three
#'   diagnostics that re-detect the type ([check_extrapolation()],
#'   [check_hat_values()], and [check_hdr()]) are instead checked up front
#'   against the detected type; if a requested one of them cannot run on the
#'   detected type, the call aborts before any diagnostic runs.
#' @param args A named list of per-diagnostic option lists, for example
#'   `list(port = list(alpha = 0.1))`. Each name must be a diagnostic being run.
#'
#' @return A [positivity_check] object bundling one [positivity_diagnostic] child
#'   per diagnostic, aligned with the `@diagnostics` names.
#'
#' @seealso [check_edp()], [check_port()], [check_hat_values()], [check_hdr()],
#'   and [check_extrapolation()] for the individual diagnostics, and
#'   [check_eta_bias()] and [check_density_ratios()] for the two that
#'   `check_positivity()` does not call.
#'
#' @examples
#' set.seed(1)
#' n <- 150
#' x1 <- rnorm(n)
#' x2 <- rnorm(n)
#' ps <- 0.2 + 0.6 * plogis(0.5 * x1 - 0.5 * x2)
#' df <- data.frame(exposure = rbinom(n, 1, ps), x1 = x1, x2 = x2)
#'
#' # The binary default set: edp, port, and extrapolation.
#' check_positivity(df, exposure, c(x1, x2))
#'
#' # Run a single diagnostic with a tuned option.
#' check_positivity(
#'   df,
#'   exposure,
#'   c(x1, x2),
#'   diagnostics = "port",
#'   args = list(port = list(alpha = 0.1))
#' )
#'
#' @export
check_positivity <- function(
  .data,
  .exposure,
  .covariates,
  diagnostics = NULL,
  exposure_type = c("auto", "binary", "categorical", "continuous"),
  args = list()
) {
  validate_data_frame(.data)

  exposure_pos <- tidyselect::eval_select(rlang::enquo(.exposure), .data)
  if (length(exposure_pos) != 1) {
    abort(
      "{.arg .exposure} must select exactly one column, not {length(exposure_pos)}.",
      error_class = "positively_selection_error"
    )
  }
  exposure_name <- names(exposure_pos)
  exposure_vec <- .data[[exposure_pos]]

  if (anyNA(exposure_vec)) {
    abort(
      "{.arg .exposure} must not contain missing values.",
      error_class = "positively_missing_error"
    )
  }

  covariate_pos <- tidyselect::eval_select(rlang::enquo(.covariates), .data)
  validate_column_selection(covariate_pos, ".covariates")
  covariate_names <- names(covariate_pos)

  if (exposure_name %in% covariate_names) {
    abort(
      c(
        "{.arg .covariates} must not include the exposure column {.val {exposure_name}}.",
        i = "Exclude it from the selection, for example {.code c(everything(), -{exposure_name})}."
      ),
      error_class = "positively_selection_error"
    )
  }

  # Resolve the exposure type once, before children run. When the type is auto it
  # is detected and announced here; when it is explicit the data are still
  # detected silently so the feasibility gate below can catch a re-detecting
  # child that the forced type would send off a cliff.
  type <- resolve_composed_exposure_type(exposure_type, exposure_vec)
  resolved_type <- type$resolved

  diagnostics <- diagnostics %||% default_diagnostics(resolved_type)
  validate_composed_diagnostics(diagnostics, resolved_type)
  if (type$is_explicit) {
    validate_composed_feasibility(diagnostics, resolved_type, type$detected)
  }
  validate_composed_args(args, diagnostics)

  # Run children one at a time. The repeated exposure-detection message from the
  # three children that re-detect internally is muffled so the type is announced
  # only once, while every other informational alert (hull skip, chunked Gower)
  # still reaches the user. A child failure is rethrown as a classed error that
  # names the diagnostic and chains the child's condition.
  checks <- vector("list", length(diagnostics))
  for (i in seq_along(diagnostics)) {
    name <- diagnostics[[i]]
    child_fn <- paste0("check_", name)
    checks[[i]] <- tryCatch(
      muffle_detection_message(run_composed_diagnostic(
        name,
        .data,
        exposure_name,
        covariate_names,
        resolved_type,
        args[[name]] %||% list()
      )),
      error = function(cnd) {
        abort(
          "{.fn {child_fn}} failed while composing diagnostics.",
          error_class = "positively_composition_error",
          parent = cnd,
          .envir = environment()
        )
      }
    )
  }

  positivity_check(checks = checks, diagnostics = diagnostics)
}

#' Resolve the composed exposure type
#'
#' Detects the exposure type from the data. When `exposure_type` is `"auto"`, the
#' detected type is announced and returned as the resolved type. When it is
#' explicit, the explicit value is the resolved type but the detected type is
#' still returned so that the feasibility gate can check the re-detecting
#' children.
#'
#' @param exposure_type The `exposure_type` argument.
#' @param exposure_vec The exposure vector.
#'
#' @return A list with `resolved` (the type children see), `detected` (the type
#'   the data imply), and `is_explicit` (whether the caller forced the type).
#' @keywords internal
#' @noRd
resolve_composed_exposure_type <- function(exposure_type, exposure_vec) {
  explicit <- rlang::arg_match(
    exposure_type,
    c("auto", "binary", "categorical", "continuous")
  )
  detected <- detect_exposure_type(exposure_vec, announce = FALSE)

  if (explicit == "auto") {
    alert_info("Treating {.arg .exposure} as {detected}")
    list(resolved = detected, detected = detected, is_explicit = FALSE)
  } else {
    list(resolved = explicit, detected = detected, is_explicit = TRUE)
  }
}

# The re-detecting children and the detected exposure type each one requires.
# check_edp() and check_port() are absent: they receive `exposure_type` and
# honor it, so a forced type is never sent off a cliff inside them.
redetecting_requirements <- function() {
  list(
    extrapolation = "binary",
    hat_values = "continuous",
    hdr = "continuous"
  )
}

#' Gate a forced exposure type against the re-detecting children
#'
#' The three diagnostics that re-detect the exposure type internally
#' ([check_extrapolation()], [check_hat_values()], [check_hdr()]) each require a
#' particular detected type. When the caller forces `exposure_type`, any of those
#' children among the requested set whose requirement the data do not meet would
#' abort mid-run, so this catches them up front and names them.
#'
#' @param diagnostics The requested diagnostic names.
#' @param resolved_type The forced exposure type children see.
#' @param detected The type the data imply.
#' @param call The calling environment, used to build the error's call.
#'
#' @return `diagnostics`, invisibly.
#' @keywords internal
#' @noRd
validate_composed_feasibility <- function(
  diagnostics,
  resolved_type,
  detected,
  call = rlang::caller_env()
) {
  requirements <- redetecting_requirements()
  requested <- intersect(diagnostics, names(requirements))
  infeasible <- requested[vapply(
    requested,
    function(name) requirements[[name]] != detected,
    logical(1)
  )]

  if (length(infeasible) > 0) {
    abort(
      c(
        "{.arg exposure_type} was set to {.val {resolved_type}}, but {.arg .exposure} is detected as {.val {detected}}.",
        x = "{.val {infeasible}} cannot run on a detected {.val {detected}} exposure.",
        i = "Drop {cli::qty(infeasible)}{?it/them} from {.arg diagnostics}, or call {?it/them} directly."
      ),
      error_class = "positively_exposure_type_error",
      call = call
    )
  }

  invisible(diagnostics)
}

#' Muffle only the repeated exposure-detection message
#'
#' Evaluates `expr` while suppressing the `detect_exposure_type()` announcement
#' ("Treating ...") that children emit when they re-detect the exposure type. All
#' other messages pass through untouched.
#'
#' @param expr An expression to evaluate.
#'
#' @return The value of `expr`.
#' @keywords internal
#' @noRd
muffle_detection_message <- function(expr) {
  withCallingHandlers(
    expr,
    message = function(cnd) {
      # Match the announcement as cli actually emits it: an info-symbol prefix
      # followed by "Treating `.exposure` as <type>". Anchoring on the backticked
      # column reference keeps an unrelated child alert that merely contains the
      # word "Treating" from being muffled.
      if (
        grepl("Treating `.exposure` as ", conditionMessage(cnd), fixed = TRUE)
      ) {
        rlang::cnd_muffle(cnd)
      }
    }
  )
}

#' Dispatch a single composed diagnostic
#'
#' Builds and evaluates a symbolic `check_*()` call for one diagnostic. The
#' selections are forwarded as resolved column names and the call is evaluated in
#' an environment binding `.data`, so the child's `@call` property prints as a
#' readable `check_*(.data = .data, .exposure = "...", ...)` rather than an
#' anonymous call with the data frame inlined. Only `check_edp()` and
#' `check_port()` accept `exposure_type`; the other three re-detect it, quietly,
#' under `muffle_detection_message()`.
#'
#' @param name The diagnostic name.
#' @param .data The data frame.
#' @param exposure_name The resolved exposure column name.
#' @param covariate_names The resolved covariate column names.
#' @param exposure_type The resolved exposure type.
#' @param method_args A list of per-method options.
#'
#' @return A [positivity_diagnostic] subclass object.
#' @keywords internal
#' @noRd
run_composed_diagnostic <- function(
  name,
  .data,
  exposure_name,
  covariate_names,
  exposure_type,
  method_args
) {
  call_args <- c(
    list(
      .data = rlang::sym(".data"),
      .exposure = exposure_name,
      .covariates = covariate_names
    ),
    if (name %in% c("edp", "port")) list(exposure_type = exposure_type),
    method_args
  )

  child_call <- rlang::call2(paste0("check_", name), !!!call_args)
  eval_env <- rlang::env(rlang::ns_env("positively"), .data = .data)
  rlang::eval_bare(child_call, eval_env)
}

#' Validate a requested set of diagnostics
#'
#' Rejects non-character input, duplicated entries, unrecognised names,
#' diagnostics that `check_positivity()` never composes, and diagnostics that do
#' not apply to the resolved exposure type.
#'
#' @param diagnostics The requested diagnostic names.
#' @param exposure_type The resolved exposure type.
#' @param call The calling environment, used to build the error's call.
#'
#' @return `diagnostics`, invisibly.
#' @keywords internal
#' @noRd
validate_composed_diagnostics <- function(
  diagnostics,
  exposure_type,
  call = rlang::caller_env()
) {
  if (!is.character(diagnostics) || length(diagnostics) == 0) {
    abort(
      "{.arg diagnostics} must be a character vector of diagnostic names.",
      error_class = "positively_diagnostic_error",
      call = call
    )
  }

  if (anyDuplicated(diagnostics) > 0) {
    duplicated_names <- unique(diagnostics[duplicated(diagnostics)])
    abort(
      c(
        "{.arg diagnostics} must not repeat a diagnostic.",
        x = "Duplicated: {.val {duplicated_names}}."
      ),
      error_class = "positively_diagnostic_error",
      call = call
    )
  }

  known <- c(composed_diagnostics(), excluded_diagnostics())
  unknown <- setdiff(diagnostics, known)
  if (length(unknown) > 0) {
    abort(
      c(
        "{.arg diagnostics} names {?an/} unrecognised diagnostic{?s}: {.val {unknown}}.",
        i = "Valid diagnostics are {.val {composed_diagnostics()}}."
      ),
      error_class = "positively_diagnostic_error",
      call = call
    )
  }

  excluded <- intersect(diagnostics, excluded_diagnostics())
  if (length(excluded) > 0) {
    abort(
      c(
        "{.fn check_positivity} does not run {.val {excluded}}.",
        i = "{.fn check_eta_bias} needs an outcome and {.fn check_density_ratios} needs user-supplied ratios.",
        i = "Call {cli::qty(excluded)}{?it/them} directly instead."
      ),
      error_class = "positively_diagnostic_error",
      call = call
    )
  }

  valid_types <- diagnostic_exposure_types()
  applies <- vapply(
    diagnostics,
    function(name) exposure_type %in% valid_types[[name]],
    logical(1)
  )
  if (!all(applies)) {
    invalid <- diagnostics[!applies]
    valid_here <- names(Filter(
      function(types) exposure_type %in% types,
      valid_types
    ))
    abort(
      c(
        "{.val {invalid}} {?does/do} not apply to a {.val {exposure_type}} exposure.",
        i = "Valid diagnostics for a {.val {exposure_type}} exposure are {.val {valid_here}}."
      ),
      error_class = "positively_diagnostic_error",
      call = call
    )
  }

  invisible(diagnostics)
}

#' Validate the per-method `args` list against the diagnostics being run
#'
#' @param args The `args` list supplied to `check_positivity()`.
#' @param diagnostics The diagnostics being run.
#' @param call The calling environment, used to build the error's call.
#'
#' @return `args`, invisibly.
#' @keywords internal
#' @noRd
validate_composed_args <- function(
  args,
  diagnostics,
  call = rlang::caller_env()
) {
  if (length(args) == 0) {
    return(invisible(args))
  }

  arg_names <- names(args)
  if (!is.list(args) || is.null(arg_names) || any(arg_names == "")) {
    abort(
      "{.arg args} must be a named list of per-diagnostic option lists.",
      error_class = "positively_args_error",
      call = call
    )
  }

  is_option_list <- vapply(args, is.list, logical(1))
  if (!all(is_option_list)) {
    bad <- arg_names[!is_option_list]
    abort(
      c(
        "Each element of {.arg args} must be a list of options.",
        x = "{.val {bad}} {?is/are} not {?a list/lists}."
      ),
      error_class = "positively_args_error",
      call = call
    )
  }

  extra <- setdiff(arg_names, diagnostics)
  if (length(extra) > 0) {
    abort(
      c(
        "{.arg args} names {?a diagnostic/diagnostics} that {?is/are} not being run: {.val {extra}}.",
        i = "The diagnostics being run are {.val {diagnostics}}."
      ),
      error_class = "positively_args_error",
      call = call
    )
  }

  invisible(args)
}

# ---- Container methods -----------------------------------------------------

#' @rdname positivity_check
#' @usage NULL
#' @order 2
method(tidy, positivity_check) <- function(x, diagnostic = NULL, ...) {
  if (is.null(diagnostic)) {
    return(combined_container_summary(x))
  }

  idx <- match(diagnostic, x@diagnostics)
  if (length(diagnostic) != 1 || is.na(idx)) {
    abort(
      c(
        "{.arg diagnostic} must name a diagnostic in this container.",
        i = "Available diagnostics are {.val {x@diagnostics}}."
      ),
      error_class = "positively_diagnostic_error"
    )
  }
  x@checks[[idx]]@results
}

#' The long combined summary across a container's children
#'
#' Stacks each child's [generics::glance()] row into `diagnostic`, `statistic`,
#' and `value` rows, coercing values to character so that statistics of mixed
#' type share one column.
#'
#' @param x A [positivity_check] object.
#'
#' @return A tibble with columns `diagnostic`, `statistic`, and `value`.
#' @keywords internal
#' @noRd
combined_container_summary <- function(x) {
  if (length(x@checks) == 0) {
    return(tibble::tibble(
      diagnostic = character(0),
      statistic = character(0),
      value = character(0)
    ))
  }
  rows <- purrr::map2(x@diagnostics, x@checks, function(name, check) {
    glanced <- generics::glance(check)
    tibble::tibble(
      diagnostic = name,
      statistic = names(glanced),
      value = vapply(
        glanced,
        function(column) as.character(column[[1]]),
        character(1)
      )
    )
  })
  purrr::list_rbind(rows)
}

#' @rdname positivity_check
#' @usage NULL
#' @order 3
method(glance, positivity_check) <- function(x, ...) {
  if (length(x@checks) == 0) {
    return(tibble::tibble(diagnostic = character(0)))
  }
  rows <- purrr::map2(x@diagnostics, x@checks, function(name, check) {
    tibble::tibble(diagnostic = name, generics::glance(check))
  })
  purrr::list_rbind(rows)
}
