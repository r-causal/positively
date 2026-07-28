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
#' `options(positively.quiet = TRUE)`). Every diagnostic it composes is handed
#' that resolved type, so none of them detects a type of its own and the type is
#' reported once. Every other informational message a diagnostic emits still
#' reaches you. The resolved type, whether you supplied it or it was detected, is
#' checked against the exposure column itself before any diagnostic runs: a
#' continuous type needs a numeric column, and a binary type needs exactly two
#' distinct values. The default diagnostic set depends on the exposure type:
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
#'   `"binary"`, `"categorical"`, or `"continuous"`. The resolved type is
#'   forwarded to every diagnostic that runs, so a supplied type is honored
#'   throughout rather than re-derived from the data by each diagnostic in turn.
#'   It must be a type the exposure column can carry, and it must apply to every
#'   entry in `diagnostics`; both are checked before any diagnostic runs.
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

  exposure_pos <- eval_select_column(
    rlang::enquo(.exposure),
    .data,
    ".exposure"
  )
  exposure_name <- names(exposure_pos)
  exposure_vec <- .data[[exposure_pos]]

  if (anyNA(exposure_vec)) {
    abort(
      "{.arg .exposure} must not contain missing values.",
      error_class = "positively_missing_error"
    )
  }

  covariate_pos <- eval_select_columns(
    rlang::enquo(.covariates),
    .data,
    ".covariates"
  )
  validate_column_selection(covariate_pos, ".covariates")
  covariate_names <- names(covariate_pos)

  if (length(intersect(covariate_pos, exposure_pos)) > 0) {
    abort(
      c(
        "{.arg .covariates} must not include the exposure column {.val {exposure_name}}.",
        i = "Exclude it from the selection, for example {.code c(everything(), -{exposure_name})}."
      ),
      error_class = "positively_selection_error"
    )
  }

  # Resolve the exposure type once, before children run, then check that the
  # exposure column can carry it. Both the default diagnostic set and the
  # applicability check take the resolved type as their premise, so the premise
  # is validated first: advice derived from a type the column cannot carry would
  # send the user to another type it cannot carry either.
  type <- resolve_composed_exposure_type(exposure_type, exposure_vec)
  resolved_type <- type$resolved
  validate_exposure_structure(
    exposure_vec,
    resolved_type,
    fn = "check_positivity"
  )

  diagnostics <- diagnostics %||% default_diagnostics(resolved_type)
  validate_composed_diagnostics(
    diagnostics,
    resolved_type,
    exposure_vec,
    type$is_explicit
  )
  validate_composed_args(args, diagnostics)

  # Run children one at a time. Each is handed the resolved exposure type and so
  # detects none of its own, which is why the type is announced once, while every
  # other informational alert (hull skip, chunked Gower) still reaches the user.
  # A child failure is rethrown as a classed error that names the diagnostic and
  # chains the child's condition.
  checks <- vector("list", length(diagnostics))
  for (i in seq_along(diagnostics)) {
    name <- diagnostics[[i]]
    child_fn <- paste0("check_", name)
    checks[[i]] <- tryCatch(
      run_composed_diagnostic(
        name,
        .data,
        exposure_name,
        covariate_names,
        resolved_type,
        args[[name]] %||% list()
      ),
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
#' Under `"auto"` the type is detected from the data and announced. A supplied
#' type is returned as given: detection is not consulted at all, so the
#' unique-value heuristic can never contradict a declaration. Whether the type
#' was supplied is reported back, because the advice that follows a rejected
#' diagnostic differs between a guess and a premise.
#'
#' @param exposure_type The `exposure_type` argument.
#' @param exposure_vec The exposure vector.
#' @param call The calling environment, used to build the error's call.
#'
#' @return A list with `resolved` (the type every diagnostic sees) and
#'   `is_explicit` (whether the caller supplied it).
#' @keywords internal
#' @noRd
resolve_composed_exposure_type <- function(
  exposure_type,
  exposure_vec,
  call = rlang::caller_env()
) {
  exposure_type <- rlang::arg_match(
    exposure_type,
    c("auto", "binary", "categorical", "continuous"),
    error_call = call
  )

  if (exposure_type == "auto") {
    detected <- detect_exposure_type(exposure_vec, announce = TRUE)
    list(resolved = detected, is_explicit = FALSE)
  } else {
    list(resolved = exposure_type, is_explicit = TRUE)
  }
}

#' Dispatch a single composed diagnostic
#'
#' Builds and evaluates a symbolic `check_*()` call for one diagnostic. The
#' selections are forwarded as resolved column names and the call is evaluated in
#' an environment binding `.data`, so the child's `@call` property prints as a
#' readable `check_*(.data = .data, .exposure = "...", ...)` rather than an
#' anonymous call with the data frame inlined. The resolved exposure type is
#' forwarded to every diagnostic, so none of them detects a type of its own.
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
  # The exposure type is forwarded to every diagnostic without exception, so
  # check_positivity() may only ever forward a type the diagnostic accepts. What
  # decides that is the `supported` set each diagnostic hands to
  # resolve_exposure_type(), since that is what gates the arg_match(); a
  # diagnostic whose `supported` set is narrower than its entry in
  # diagnostic_exposure_types() would turn an argument error into a composition
  # failure. validate_composed_diagnostics() has already rejected every
  # diagnostic the resolved type does not apply to.
  call_args <- c(
    list(
      .data = rlang::sym(".data"),
      .exposure = exposure_name,
      .covariates = covariate_names,
      exposure_type = exposure_type
    ),
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
#' @param exposure_vec The exposure vector, which decides which exposure types
#'   the advice may offer.
#' @param is_explicit Whether the caller supplied the exposure type. A detected
#'   type is a guess, so the rejection also names the types that would run the
#'   rejected diagnostics; a supplied type is a premise, and the advice that
#'   follows is to fix the requested set instead.
#' @param call The calling environment, used to build the error's call.
#'
#' @return `diagnostics`, invisibly.
#' @keywords internal
#' @noRd
validate_composed_diagnostics <- function(
  diagnostics,
  exposure_type,
  exposure_vec,
  is_explicit,
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
        i = "Valid diagnostics for a {.val {exposure_type}} exposure are {.val {valid_here}}.",
        if (!is_explicit) exposure_type_advice(invalid, exposure_vec)
      ),
      error_class = "positively_diagnostic_error",
      call = call
    )
  }

  invisible(diagnostics)
}

#' Name the exposure types a rejected diagnostic could be run under
#'
#' Builds the advice that follows an applicability rejection under a detected
#' exposure type. One bullet is emitted per type the rejected diagnostics need,
#' naming the diagnostics that type would run, so diagnostics needing the same
#' type are named together and diagnostics needing different types are named
#' apart. Only types the exposure column can structurally carry are offered:
#' `exposure_carries_type()` is the same rule the entry point's own gate applies,
#' so a caller who acts on a bullet is never sent into that gate. When no needed
#' type survives the filter there is nothing worth suggesting and no bullet is
#' emitted. The bullets are formatted here rather than handed to `abort()` as cli
#' templates because how many of them there are depends on the request.
#'
#' The rule is written for any number of surviving types, but today's diagnostic
#' set can yield at most one. Only `check_extrapolation()` needs a binary type,
#' and it is rejected only when the exposure is not binary, which is exactly when
#' a binary type is not on offer.
#'
#' @param invalid The rejected diagnostic names.
#' @param exposure_vec The exposure vector, which decides which needed types are
#'   available to suggest.
#'
#' @return A character vector of info bullets, one per feasible needed exposure
#'   type, empty when none is feasible.
#' @keywords internal
#' @noRd
exposure_type_advice <- function(invalid, exposure_vec) {
  valid_types <- diagnostic_exposure_types()[invalid]
  needed <- unique(unlist(valid_types, use.names = FALSE))
  feasible <- needed[vapply(
    needed,
    function(type) exposure_carries_type(exposure_vec, type),
    logical(1)
  )]

  bullets <- vapply(
    feasible,
    function(type) {
      runs <- invalid[vapply(
        valid_types,
        function(types) type %in% types,
        logical(1)
      )]
      type_code <- paste0("exposure_type = \"", type, "\"")
      cli::format_inline("Set {.code {type_code}} to run {.val {runs}}.")
    },
    character(1)
  )

  rlang::set_names(bullets, rep("i", length(bullets)))
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

  if (anyDuplicated(arg_names) > 0) {
    duplicated_names <- unique(arg_names[duplicated(arg_names)])
    abort(
      c(
        "{.arg args} must not repeat a diagnostic name.",
        x = "Duplicated: {.val {duplicated_names}}."
      ),
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

  # The resolved exposure type reaches every diagnostic as a named argument, so
  # a second one inside an option list would arrive as a duplicate formal and
  # surface as an unclassed base R error blaming the diagnostic.
  sets_type <- vapply(
    args,
    function(options) "exposure_type" %in% names(options),
    logical(1)
  )
  if (any(sets_type)) {
    collided <- arg_names[sets_type]
    abort(
      c(
        "{.arg args} must not set {.arg exposure_type}.",
        x = "{.val {collided}} {?sets/set} it.",
        i = "{.fn check_positivity} forwards its own {.arg exposure_type} to every
             diagnostic, so declare the type there instead."
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
