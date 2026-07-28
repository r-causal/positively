# check_positivity() is the top-level entry point: it detects the exposure type, runs
# the default diagnostics for that type (overridable through `diagnostics`),
# threads per-method options from `args`, and returns a positivity_check
# container. The default sets are fixed by exposure type: binary runs edp, port,
# and extrapolation; categorical runs edp and port; continuous runs edp, port,
# hat_values, and hdr. check_eta_bias() and check_density_ratios() are never run
# because they need an outcome or user-supplied ratios. These specs assert the
# plumbing (selection, passthrough, container shape, errors), so they lean on the
# cheapest diagnostics and shrink the expensive null resampling through `args`.

# The child produced by a named diagnostic, located through the aligned
# @diagnostics vector rather than list names.
child_named <- function(res, name) {
  res@checks[[match(name, res@diagnostics)]]
}

# Three-level categorical exposure with two numeric covariates.
sim_pos_categorical <- function(n = 150, seed = 1) {
  withr::local_seed(seed)
  x1 <- stats::rnorm(n)
  x2 <- stats::rnorm(n)
  level <- sample(c("low", "mid", "high"), n, replace = TRUE)
  tibble::tibble(
    exposure = factor(level, levels = c("low", "mid", "high")),
    x1 = x1,
    x2 = x2
  )
}

# ---- Default diagnostic sets by exposure type -----------------------------

test_that("a binary exposure runs edp, port, and extrapolation by default", {
  local_quiet()
  data <- dgp_good_positivity(n = 200)
  res <- check_positivity(data, exposure, c(x1, x2))

  expect_true(S7::S7_inherits(res, positivity_check))
  expect_identical(res@diagnostics, c("edp", "port", "extrapolation"))
  expect_length(res@checks, 3L)
  expect_true(all(vapply(
    res@checks,
    function(check) S7::S7_inherits(check, positivity_diagnostic),
    logical(1)
  )))
})

test_that("a binary default run excludes eta_bias and density_ratios", {
  local_quiet()
  data <- dgp_good_positivity(n = 200)
  res <- check_positivity(data, exposure, c(x1, x2))
  expect_false(any(c("eta_bias", "density_ratios") %in% res@diagnostics))
})

test_that("each binary child carries the resolved exposure name and type", {
  local_quiet()
  data <- dgp_good_positivity(n = 200)
  res <- check_positivity(data, exposure, c(x1, x2))
  expect_true(all(
    vapply(
      res@checks,
      function(check) check@exposure_type,
      character(1)
    ) ==
      "binary"
  ))
  expect_true(all(vapply(
    res@checks,
    function(check) identical(check@exposure, "exposure"),
    logical(1)
  )))
})

test_that("a categorical exposure runs edp and port by default", {
  local_quiet()
  data <- sim_pos_categorical(150)
  res <- check_positivity(data, exposure, c(x1, x2))

  expect_true(S7::S7_inherits(res, positivity_check))
  expect_identical(res@diagnostics, c("edp", "port"))
  expect_length(res@checks, 2L)
  expect_true(all(
    vapply(
      res@checks,
      function(check) check@exposure_type,
      character(1)
    ) ==
      "categorical"
  ))
})

test_that("a continuous exposure runs edp, port, hat_values, and hdr by default", {
  local_quiet()
  data <- dgp_continuous_support_gap(n = 150)
  # null_reps is shrunk so the hat-values null resampling stays cheap.
  res <- check_positivity(
    data,
    exposure,
    x1,
    args = list(hat_values = list(null_reps = 25))
  )

  expect_true(S7::S7_inherits(res, positivity_check))
  expect_identical(res@diagnostics, c("edp", "port", "hat_values", "hdr"))
  expect_length(res@checks, 4L)
  expect_true(all(
    vapply(
      res@checks,
      function(check) check@exposure_type,
      character(1)
    ) ==
      "continuous"
  ))
})

# ---- Explicit method selection --------------------------------------------

test_that("an explicit diagnostics vector is honoured in order", {
  local_quiet()
  data <- dgp_good_positivity(n = 200)
  res <- check_positivity(
    data,
    exposure,
    c(x1, x2),
    diagnostics = c("port", "extrapolation")
  )
  expect_identical(res@diagnostics, c("port", "extrapolation"))
  expect_length(res@checks, 2L)
})

test_that("a single requested diagnostic yields a one-child container", {
  local_quiet()
  data <- dgp_good_positivity(n = 200)
  res <- check_positivity(data, exposure, c(x1, x2), diagnostics = "port")
  expect_identical(res@diagnostics, "port")
  expect_length(res@checks, 1L)
})

test_that("an explicit exposure_type overrides auto detection", {
  local_quiet()
  data <- dgp_good_positivity(n = 200)
  res <- check_positivity(
    data,
    exposure,
    c(x1, x2),
    diagnostics = "port",
    exposure_type = "categorical"
  )
  expect_identical(child_named(res, "port")@exposure_type, "categorical")
})

# A coarse numeric exposure: few unique values, so detection calls it
# categorical, yet the column is numeric and so can carry a declared continuous
# type without complaint.
sim_pos_coarse_numeric <- function(n = 150, seed = 1) {
  withr::local_seed(seed)
  tibble::tibble(
    exposure = sample(1:8, n, replace = TRUE),
    x1 = stats::rnorm(n),
    x2 = stats::rnorm(n)
  )
}

test_that("a declared type the exposure cannot carry aborts before any child runs", {
  local_quiet()
  data <- sim_pos_categorical(150)
  # The exposure is a three-level factor, so a declared binary type is
  # structurally impossible: binary needs exactly two distinct values. The abort
  # must come from the entry point's own structural check, before any diagnostic
  # is dispatched, so that no child burns a full computation on a request that
  # cannot succeed. Recording the dispatches pins that directly; the condition
  # class pins the source, since anything a child raised would arrive wrapped as
  # a composition error.
  dispatched <- character()
  testthat::local_mocked_bindings(
    check_port = function(...) {
      dispatched <<- c(dispatched, "port")
      NULL
    },
    check_extrapolation = function(...) {
      dispatched <<- c(dispatched, "extrapolation")
      NULL
    }
  )

  cnd <- rlang::catch_cnd(
    check_positivity(
      data,
      exposure,
      c(x1, x2),
      diagnostics = c("port", "extrapolation"),
      exposure_type = "binary"
    ),
    classes = "error"
  )
  expect_s3_class(cnd, "positively_exposure_type_error")
  expect_false(inherits(cnd, "positively_composition_error"))
  expect_identical(dispatched, character())

  # port accepts every exposure type, so nothing about the requested set signals
  # the problem. Only the column's own structure does, and the gate must catch it
  # just the same.
  cnd_port <- rlang::catch_cnd(
    check_positivity(
      data,
      exposure,
      c(x1, x2),
      diagnostics = "port",
      exposure_type = "binary"
    ),
    classes = "error"
  )
  expect_s3_class(cnd_port, "positively_exposure_type_error")
  expect_false(inherits(cnd_port, "positively_composition_error"))
  expect_identical(dispatched, character())
})

test_that("the structural gate runs before the applicability check", {
  local_quiet()
  data <- sim_pos_categorical(150)
  # Two guards would fire on this call: the three-level factor cannot carry the
  # declared binary type, and hat_values does not apply to a binary exposure
  # either. Only the structural one gives advice worth acting on. Reported the
  # other way round, the user is told which type would run hat_values, and the
  # column cannot carry that type either, so following the advice aborts again.
  # The condition class pins which guard spoke first.
  dispatched <- character()
  testthat::local_mocked_bindings(
    check_hat_values = function(...) {
      dispatched <<- c(dispatched, "hat_values")
      NULL
    }
  )

  cnd <- rlang::catch_cnd(
    check_positivity(
      data,
      exposure,
      c(x1, x2),
      diagnostics = "hat_values",
      exposure_type = "binary"
    ),
    classes = "error"
  )
  expect_s3_class(cnd, "positively_exposure_type_error")
  expect_false(inherits(cnd, "positively_diagnostic_error"))
  expect_false(inherits(cnd, "positively_composition_error"))
  expect_identical(dispatched, character())
})

test_that("a declared type the exposure can carry is honoured", {
  local_quiet()
  data <- sim_pos_coarse_numeric(150)
  # Detection calls the coarse numeric exposure categorical, but the column is
  # numeric and so carries a declared continuous type, which check_port()
  # computes on.
  res <- check_positivity(
    data,
    exposure,
    c(x1, x2),
    diagnostics = "port",
    exposure_type = "continuous"
  )
  expect_identical(child_named(res, "port")@exposure_type, "continuous")
})

test_that("a declared continuous type runs the continuous set on a coarse dose", {
  local_quiet()
  data <- dgp_coarse_dose(n = 150)
  # Eight distinct milligram levels over 150 rows put the unique-value ratio
  # below the is_categorical() cutoff, so detection reads a genuinely continuous
  # dose as categorical. Declaring the type must be enough on its own: the full
  # continuous default set runs, hat_values and hdr included, and every child
  # computes on the declared type rather than on its own reading of the column.
  expect_identical(
    detect_exposure_type(data$exposure, announce = FALSE),
    "categorical"
  )

  # null_reps is shrunk so the hat-values null resampling stays cheap.
  res <- check_positivity(
    data,
    exposure,
    x1,
    exposure_type = "continuous",
    args = list(hat_values = list(null_reps = 25))
  )

  expect_identical(res@diagnostics, c("edp", "port", "hat_values", "hdr"))
  expect_identical(
    vapply(res@checks, function(check) check@exposure_type, character(1)),
    rep("continuous", 4L)
  )
})

# ---- Forwarding the resolved exposure type --------------------------------

test_that("every composed diagnostic accepts each type it can be forwarded", {
  # check_positivity() forwards the resolved exposure type to every child, so a
  # type listed for a diagnostic here that the child's own exposure_type
  # argument does not offer would reach that child as an arg_match() failure and
  # surface as a composition error rather than as the argument error it is.
  # Deriving both sides from the package keeps the invariant true of diagnostics
  # added later, not only of the five composed today.
  offered <- diagnostic_exposure_types()
  expect_setequal(names(offered), composed_diagnostics())

  unsupported <- vapply(
    composed_diagnostics(),
    function(name) {
      child <- rlang::env_get(
        rlang::ns_env("positively"),
        paste0("check_", name)
      )
      accepted <- eval(formals(child)$exposure_type)
      paste(setdiff(offered[[name]], accepted), collapse = ", ")
    },
    character(1)
  )
  expect_identical(
    unsupported,
    rlang::set_names(
      rep("", length(composed_diagnostics())),
      composed_diagnostics()
    )
  )
})

# ---- Child alerts and failures --------------------------------------------

test_that("a child informational alert survives the entry point", {
  withr::local_options(
    positively.quiet = FALSE,
    positively.gower_chunk_threshold = 20
  )
  data <- dgp_good_positivity(n = 60)
  # The chunked-Gower notice comes from within check_extrapolation(). The child
  # is handed the resolved exposure type and so announces no type of its own,
  # but every other informational alert it raises must still reach the user.
  expect_message(
    expect_message(
      check_positivity(
        data,
        exposure,
        c(x1, x2),
        diagnostics = "extrapolation"
      ),
      "row chunks"
    ),
    "Treating `.exposure` as binary"
  )
})

test_that("the detection message is announced exactly once at default verbosity", {
  withr::local_options(positively.quiet = FALSE)
  data <- dgp_good_positivity(n = 60)
  # The entry point resolves the exposure type once and forwards it to every
  # child, so extrapolation receives a concrete type and never detects one for
  # itself. The announcement therefore has exactly one source.
  messages <- character()
  withCallingHandlers(
    check_positivity(data, exposure, c(x1, x2), diagnostics = "extrapolation"),
    message = function(cnd) {
      messages <<- c(messages, conditionMessage(cnd))
      invokeRestart("muffleMessage")
    }
  )
  detection <- grepl("Treating `.exposure` as", messages, fixed = TRUE)
  expect_equal(sum(detection), 1L)
})

test_that("a failing child aborts with a classed error naming the diagnostic", {
  local_quiet()
  data <- dgp_good_positivity(n = 200)
  # alpha = 1.5 is rejected inside check_port(); check_positivity() must rethrow a
  # classed positively error rather than a purrr_error_indexed.
  err <- expect_error(
    check_positivity(
      data,
      exposure,
      c(x1, x2),
      diagnostics = "port",
      args = list(port = list(alpha = 1.5))
    ),
    class = "positively_error"
  )
  expect_s3_class(err, "positively_composition_error")
  expect_false(inherits(err, "purrr_error_indexed"))
})

test_that("a misspelled port control option surfaces as a classed error", {
  local_quiet()
  data <- dgp_good_positivity(n = 200)
  # A typo in an rpart.control option threaded through args must not be silently
  # dropped; check_port() aborts and check_positivity() rethrows it wrapped.
  err <- expect_error(
    check_positivity(
      data,
      exposure,
      c(x1, x2),
      diagnostics = "port",
      args = list(port = list(alpa = 0.4))
    ),
    class = "positively_error"
  )
  expect_s3_class(err, "positively_composition_error")
  expect_s3_class(err$parent, "positively_args_error")
})

# ---- Duplicate and malformed arguments ------------------------------------

test_that("duplicate diagnostics are rejected", {
  local_quiet()
  data <- dgp_good_positivity(n = 200)
  expect_error(
    check_positivity(
      data,
      exposure,
      c(x1, x2),
      diagnostics = c("port", "port")
    ),
    class = "positively_error"
  )
})

test_that("duplicate args names are rejected", {
  local_quiet()
  data <- dgp_good_positivity(n = 200)
  # A named list with a repeated name silently keeps only the first entry, so the
  # second set of options would vanish without warning.
  expect_error(
    # jarl-ignore duplicated_arguments: the repeated name is the condition under test
    check_positivity(
      data,
      exposure,
      c(x1, x2),
      diagnostics = "port",
      args = list(port = list(alpha = 0.2), port = list(alpha = 0.4))
    ),
    class = "positively_args_error"
  )
  expect_snapshot(
    # jarl-ignore duplicated_arguments: the repeated name is the condition under test
    check_positivity(
      data,
      exposure,
      c(x1, x2),
      diagnostics = "port",
      args = list(port = list(alpha = 0.2), port = list(alpha = 0.4))
    ),
    error = TRUE
  )
})

test_that("an args entry that is not a list is rejected", {
  local_quiet()
  data <- dgp_good_positivity(n = 200)
  expect_error(
    check_positivity(
      data,
      exposure,
      c(x1, x2),
      diagnostics = "port",
      args = list(port = "oops")
    ),
    class = "positively_error"
  )
})

test_that("args naming exposure_type for a diagnostic is rejected", {
  local_quiet()
  data <- dgp_good_positivity(n = 200)
  # check_positivity() forwards its own exposure_type to every child, so an
  # exposure_type inside a per-diagnostic option list collides with it. Left to
  # R that surfaces as an unclassed "matched by multiple actual arguments"
  # error, wrapped as a composition failure that blames the child rather than
  # the argument the caller wrote.
  cnd <- rlang::catch_cnd(
    check_positivity(
      data,
      exposure,
      c(x1, x2),
      diagnostics = "port",
      args = list(port = list(exposure_type = "categorical"))
    ),
    classes = "error"
  )
  expect_s3_class(cnd, "positively_args_error")
  expect_false(inherits(cnd, "positively_composition_error"))
})

test_that("tidy() is reachable without qualifying the generic", {
  local_quiet()
  data <- dgp_good_positivity(n = 200)
  res <- check_positivity(data, exposure, c(x1, x2), diagnostics = "port")
  # tidy() and glance() are re-exported, so the bare generics resolve.
  expect_s3_class(tidy(res), "tbl_df")
  expect_s3_class(glance(res), "tbl_df")
})

# ---- Per-method argument passthrough --------------------------------------

test_that("args threads a per-method option into the child call", {
  local_quiet()
  data <- dgp_good_positivity(n = 200)
  res <- check_positivity(
    data,
    exposure,
    c(x1, x2),
    diagnostics = "port",
    args = list(port = list(alpha = 0.2))
  )
  expect_identical(child_named(res, "port")@alpha, 0.2)
})

test_that("args reaches a diagnostic that only some exposure types run", {
  local_quiet()
  data <- dgp_continuous_support_gap(n = 150)
  res <- check_positivity(
    data,
    exposure,
    x1,
    diagnostics = "hat_values",
    args = list(hat_values = list(null_reps = 25))
  )
  # null_dist has one entry per null replicate, so its length proves the
  # passthrough of null_reps.
  expect_length(child_named(res, "hat_values")@null_dist, 25L)
})

# ---- Container tidy(), glance(), and print() ------------------------------

test_that("tidy() on the container returns a long combined summary", {
  local_quiet()
  data <- dgp_good_positivity(n = 200)
  res <- check_positivity(
    data,
    exposure,
    c(x1, x2),
    diagnostics = c("port", "extrapolation")
  )
  combined <- generics::tidy(res)
  expect_s3_class(combined, "tbl_df")
  expect_setequal(names(combined), c("diagnostic", "statistic", "value"))
  expect_setequal(unique(combined$diagnostic), c("port", "extrapolation"))
})

test_that("tidy() with a diagnostic name returns that child's results", {
  local_quiet()
  data <- dgp_good_positivity(n = 200)
  res <- check_positivity(
    data,
    exposure,
    c(x1, x2),
    diagnostics = c("port", "extrapolation")
  )
  expect_identical(
    generics::tidy(res, diagnostic = "port"),
    child_named(res, "port")@results
  )
})

test_that("tidy() rejects an unknown diagnostic name", {
  local_quiet()
  data <- dgp_good_positivity(n = 200)
  res <- check_positivity(data, exposure, c(x1, x2), diagnostics = "port")
  expect_error(
    generics::tidy(res, diagnostic = "nonexistent"),
    class = "positively_error"
  )
})

test_that("glance() on the container has one row per diagnostic", {
  local_quiet()
  data <- dgp_good_positivity(n = 200)
  res <- check_positivity(
    data,
    exposure,
    c(x1, x2),
    diagnostics = c("port", "extrapolation")
  )
  glanced <- generics::glance(res)
  expect_s3_class(glanced, "tbl_df")
  expect_identical(nrow(glanced), 2L)
  expect_true("diagnostic" %in% names(glanced))
  expect_setequal(glanced$diagnostic, c("port", "extrapolation"))
})

test_that("printing the container names each diagnostic section", {
  local_quiet()
  data <- dgp_good_positivity(n = 200)
  res <- check_positivity(
    data,
    exposure,
    c(x1, x2),
    diagnostics = c("port", "extrapolation")
  )
  expect_output(print(res), "Positivity check")
  expect_output(print(res), "port")
  expect_output(print(res), "extrapolation")
})

# ---- Extracting a child diagnostic ----------------------------------------

test_that("[[ extracts the identical child diagnostic object by name", {
  local_quiet()
  data <- dgp_good_positivity(n = 200)
  res <- check_positivity(
    data,
    exposure,
    c(x1, x2),
    diagnostics = c("port", "extrapolation")
  )
  expect_identical(res[["port"]], child_named(res, "port"))
  expect_identical(res[["extrapolation"]], child_named(res, "extrapolation"))
})

test_that("autoplot works on an extracted child", {
  local_quiet()
  data <- dgp_good_positivity(n = 200)
  res <- check_positivity(data, exposure, c(x1, x2), diagnostics = "port")
  expect_s3_class(ggplot2::autoplot(res[["port"]]), "ggplot")
})

test_that("[[ aborts on an unknown diagnostic name", {
  local_quiet()
  data <- dgp_good_positivity(n = 200)
  res <- check_positivity(data, exposure, c(x1, x2), diagnostics = "port")
  expect_error(res[["nonexistent"]], class = "positively_error")
})

test_that("names() returns the composed diagnostic names", {
  local_quiet()
  data <- dgp_good_positivity(n = 200)
  res <- check_positivity(
    data,
    exposure,
    c(x1, x2),
    diagnostics = c("port", "extrapolation")
  )
  expect_identical(names(res), c("port", "extrapolation"))
})

test_that("names() on an empty container is character(0)", {
  empty <- positivity_check(checks = list(), diagnostics = character(0))
  expect_identical(names(empty), character(0))
})

test_that("[[ extracts a child by numeric position", {
  local_quiet()
  data <- dgp_good_positivity(n = 200)
  res <- check_positivity(
    data,
    exposure,
    c(x1, x2),
    diagnostics = c("port", "extrapolation")
  )
  expect_identical(res[[1]], child_named(res, "port"))
  expect_identical(res[[2]], child_named(res, "extrapolation"))
})

test_that("[[ aborts on an out-of-bounds numeric index", {
  local_quiet()
  data <- dgp_good_positivity(n = 200)
  res <- check_positivity(data, exposure, c(x1, x2), diagnostics = "port")
  expect_error(res[[2]], class = "positively_bounds_error")
})

test_that("[[ aborts on a non-scalar index", {
  local_quiet()
  data <- dgp_good_positivity(n = 200)
  res <- check_positivity(
    data,
    exposure,
    c(x1, x2),
    diagnostics = c("port", "extrapolation")
  )
  expect_error(res[[c("port", "extrapolation")]], class = "positively_error")
})

# ---- Classed errors -------------------------------------------------------

test_that("an unknown diagnostic name aborts", {
  local_quiet()
  data <- dgp_good_positivity(n = 200)
  expect_error(
    check_positivity(data, exposure, c(x1, x2), diagnostics = "bogus"),
    class = "positively_error"
  )
})

test_that("a diagnostic that needs an outcome cannot be requested", {
  local_quiet()
  data <- dgp_good_positivity(n = 200)
  # eta_bias needs an outcome and density_ratios needs user-supplied ratios, so
  # neither is a valid composable diagnostic.
  expect_error(
    check_positivity(data, exposure, c(x1, x2), diagnostics = "eta_bias"),
    class = "positively_error"
  )
  expect_error(
    check_positivity(data, exposure, c(x1, x2), diagnostics = "density_ratios"),
    class = "positively_error"
  )
})

test_that("covariates overlapping the exposure are rejected", {
  local_quiet()
  data <- dgp_good_positivity(n = 200)
  # everything() sweeps the exposure column into the covariate selection.
  expect_error(
    check_positivity(data, exposure, tidyselect::everything()),
    class = "positively_selection_error"
  )
  # A covariate selection disjoint from the exposure still succeeds.
  res <- check_positivity(data, exposure, c(x1, x2))
  expect_true(S7::S7_inherits(res, positivity_check))
})

test_that("a renamed exposure inside the covariate selection is rejected", {
  local_quiet()
  data <- dgp_good_positivity(n = 200)
  # Renaming the exposure column inside the selection hides it from a name-based
  # overlap check, so the guard must compare resolved positions like
  # check_eta_bias() does. Otherwise the renamed column slips through and the
  # failure surfaces later, deep inside a composed child, as a confusing missing
  # column. The top-level condition, not merely something in its parent chain,
  # must be the overlap selection error raised before any child runs.
  cnd <- rlang::catch_cnd(
    check_positivity(data, exposure, c(foo = exposure, x1)),
    classes = "error"
  )
  expect_s3_class(cnd, "positively_selection_error")
  expect_false(inherits(cnd, "positively_composition_error"))
})

test_that("hat_values is invalid for a binary exposure", {
  local_quiet()
  data <- dgp_good_positivity(n = 200)
  expect_error(
    check_positivity(data, exposure, c(x1, x2), diagnostics = "hat_values"),
    class = "positively_error"
  )
})

test_that("extrapolation is invalid for a continuous exposure", {
  local_quiet()
  data <- dgp_continuous_support_gap(n = 150)
  expect_error(
    check_positivity(data, exposure, x1, diagnostics = "extrapolation"),
    class = "positively_error"
  )
})

# The lines of a condition's message that mention every one of `terms`, or the
# whole message when no line mentions them all, so that a failed expectation
# prints what the message actually said. Matching a line at a time keeps a
# diagnostic and the type it needs read as one statement rather than as two
# words that happen to appear somewhere in the message.
lines_mentioning <- function(cnd, terms) {
  message <- conditionMessage(cnd)
  lines <- strsplit(message, "\n", fixed = TRUE)[[1]]
  mentions <- vapply(
    lines,
    function(line) {
      all(vapply(terms, grepl, logical(1), x = line, fixed = TRUE))
    },
    logical(1)
  )
  if (any(mentions)) lines[mentions] else message
}

test_that("inapplicable diagnostics that share a type name it once", {
  local_quiet()
  # cli wraps to the console width and the advice is read one line at a time, so
  # the width is widened to keep a long sentence intact.
  withr::local_options(cli.width = 500)
  data <- dgp_good_positivity(n = 200)
  # Under auto detection the exposure type is a heuristic guess rather than
  # something the caller asserted, and the abort is the only place a user learns
  # the guess can be overridden. hat_values and hdr both need a continuous
  # exposure and the 0/1 exposure column is numeric, so that type is one the
  # column can carry and a single sentence names it and both diagnostics.
  cnd <- rlang::catch_cnd(
    check_positivity(
      data,
      exposure,
      c(x1, x2),
      diagnostics = c("hat_values", "hdr")
    ),
    classes = "error"
  )
  expect_s3_class(cnd, "positively_diagnostic_error")
  expect_match(conditionMessage(cnd), "exposure_type", fixed = TRUE)
  expect_match(
    lines_mentioning(cnd, c("hat_values", "hdr")),
    "continuous",
    fixed = TRUE,
    all = FALSE
  )
})

test_that("no type is suggested when the column can carry none of them", {
  local_quiet()
  withr::local_options(cli.width = 500)
  data <- sim_pos_categorical(150)
  # The exposure is a three-level factor. hat_values needs a continuous type,
  # which a factor cannot carry, and extrapolation needs a binary one, which a
  # column with three distinct values cannot carry. Either suggestion would abort
  # on the structural gate the moment the caller acted on it, so neither is
  # offered and the list of valid diagnostics is the whole of the advice.
  cnd <- rlang::catch_cnd(
    check_positivity(
      data,
      exposure,
      c(x1, x2),
      diagnostics = c("hat_values", "extrapolation")
    ),
    classes = "error"
  )
  expect_s3_class(cnd, "positively_diagnostic_error")
  expect_no_match(conditionMessage(cnd), "exposure_type", fixed = TRUE)
})

test_that("a suggested exposure type runs the diagnostics the advice names", {
  local_quiet()
  withr::local_options(cli.width = 500)
  data <- dgp_coarse_dose(n = 150)
  # Eight distinct milligram levels over 150 rows put the unique-value ratio
  # below the is_categorical() cutoff, so detection reads a genuinely continuous
  # dose as categorical and rejects both continuous diagnostics. The column is
  # numeric, so the continuous type the advice names is one it can carry, and
  # acting on the advice is what the advice is for: the same call with that type
  # declared runs both diagnostics through.
  cnd <- rlang::catch_cnd(
    check_positivity(data, exposure, x1, diagnostics = c("hat_values", "hdr")),
    classes = "error"
  )
  expect_s3_class(cnd, "positively_diagnostic_error")
  expect_match(
    lines_mentioning(cnd, c("hat_values", "hdr")),
    "exposure_type = \"continuous\"",
    fixed = TRUE,
    all = FALSE
  )

  # null_reps is shrunk so the hat-values null resampling stays cheap.
  res <- check_positivity(
    data,
    exposure,
    x1,
    diagnostics = c("hat_values", "hdr"),
    exposure_type = "continuous",
    args = list(hat_values = list(null_reps = 25))
  )
  expect_identical(res@diagnostics, c("hat_values", "hdr"))
})

test_that("a declared type is answered with the diagnostic set, not a type", {
  local_quiet()
  withr::local_options(cli.width = 500)
  data <- dgp_coarse_dose(n = 150)
  # The two calls resolve to the same categorical type and reject the same
  # diagnostic; only the source of the type differs. A detected type is a guess,
  # so naming the type that would run hat_values tells the caller the guess is
  # overridable. A declared type is a premise, and the advice that follows from a
  # premise is to fix the requested set, which the list of valid diagnostics
  # already gives. The coarse dose is numeric, so the continuous type is on offer
  # in the first call and its absence from the second is the difference under
  # test.
  detected <- rlang::catch_cnd(
    check_positivity(data, exposure, x1, diagnostics = "hat_values"),
    classes = "error"
  )
  declared <- rlang::catch_cnd(
    check_positivity(
      data,
      exposure,
      x1,
      diagnostics = "hat_values",
      exposure_type = "categorical"
    ),
    classes = "error"
  )
  expect_s3_class(declared, "positively_diagnostic_error")
  expect_match(conditionMessage(detected), "continuous", fixed = TRUE)
  expect_no_match(conditionMessage(declared), "continuous", fixed = TRUE)
})

test_that("args naming a diagnostic that is not run aborts", {
  local_quiet()
  data <- dgp_good_positivity(n = 200)
  expect_error(
    check_positivity(
      data,
      exposure,
      c(x1, x2),
      diagnostics = "port",
      args = list(bogus = list(alpha = 0.2))
    ),
    class = "positively_error"
  )
})

# ---- Quiet-mode behavior --------------------------------------------------

test_that("auto detection announces the resolved exposure type", {
  data <- dgp_good_positivity(n = 200)
  withr::local_options(positively.quiet = FALSE)
  expect_message(
    check_positivity(data, exposure, c(x1, x2), diagnostics = "port"),
    "Treating"
  )
})

test_that("quiet mode suppresses the detection message", {
  local_quiet()
  data <- dgp_good_positivity(n = 200)
  expect_no_message(
    check_positivity(data, exposure, c(x1, x2), diagnostics = "port")
  )
})

# ---- Snapshots: container print and classed errors ------------------------

test_that("the container print method is stable", {
  local_quiet()
  data <- dgp_good_positivity(n = 200)
  res <- check_positivity(
    data,
    exposure,
    c(x1, x2),
    diagnostics = c("port", "extrapolation")
  )
  expect_snapshot(print(res))
})

test_that("the classed errors are stable", {
  local_quiet()
  data <- dgp_good_positivity(n = 200)
  continuous <- dgp_continuous_support_gap(n = 150)

  # An unrecognised diagnostic name.
  expect_snapshot(
    check_positivity(data, exposure, c(x1, x2), diagnostics = "bogus"),
    error = TRUE
  )
  # Diagnostics that need an outcome or user-supplied ratios.
  expect_snapshot(
    check_positivity(data, exposure, c(x1, x2), diagnostics = "eta_bias"),
    error = TRUE
  )
  expect_snapshot(
    check_positivity(data, exposure, c(x1, x2), diagnostics = "density_ratios"),
    error = TRUE
  )
  # A diagnostic that does not apply to the exposure type. The 0/1 exposure is
  # numeric, so the continuous type hat_values needs is on offer; the continuous
  # exposure has far more than two distinct values, so the binary type
  # extrapolation needs is not.
  expect_snapshot(
    check_positivity(data, exposure, c(x1, x2), diagnostics = "hat_values"),
    error = TRUE
  )
  expect_snapshot(
    check_positivity(continuous, exposure, x1, diagnostics = "extrapolation"),
    error = TRUE
  )
  # An args name that is not being run.
  expect_snapshot(
    check_positivity(
      data,
      exposure,
      c(x1, x2),
      diagnostics = "port",
      args = list(bogus = list(alpha = 0.2))
    ),
    error = TRUE
  )
  # A duplicated diagnostic name.
  expect_snapshot(
    check_positivity(
      data,
      exposure,
      c(x1, x2),
      diagnostics = c("port", "port")
    ),
    error = TRUE
  )
  # An args entry that is not a list.
  expect_snapshot(
    check_positivity(
      data,
      exposure,
      c(x1, x2),
      diagnostics = "port",
      args = list(port = "oops")
    ),
    error = TRUE
  )
  # A declared exposure type the exposure column cannot carry.
  expect_snapshot(
    check_positivity(
      sim_pos_categorical(150),
      exposure,
      c(x1, x2),
      diagnostics = c("port", "extrapolation"),
      exposure_type = "binary"
    ),
    error = TRUE
  )
  # A failing child names the diagnostic and chains the child's condition.
  expect_snapshot(
    check_positivity(
      data,
      exposure,
      c(x1, x2),
      diagnostics = "port",
      args = list(port = list(alpha = 1.5))
    ),
    error = TRUE
  )
  # A single-column exposure selection is required.
  expect_snapshot(
    check_positivity(data, c(exposure, x1), x2),
    error = TRUE
  )
  # A covariate selection that includes the exposure column.
  expect_snapshot(
    check_positivity(data, exposure, tidyselect::everything()),
    error = TRUE
  )
  # tidy() rejects an unknown diagnostic name.
  res <- check_positivity(data, exposure, c(x1, x2), diagnostics = "port")
  expect_snapshot(
    generics::tidy(res, diagnostic = "nonexistent"),
    error = TRUE
  )
  # [[ rejects an unknown diagnostic name.
  expect_snapshot(
    res[["nonexistent"]],
    error = TRUE
  )
  # [[ rejects an out-of-bounds numeric index.
  expect_snapshot(
    res[[5]],
    error = TRUE
  )
  # [[ rejects a non-scalar index.
  expect_snapshot(
    res[[c("port", "extrapolation")]],
    error = TRUE
  )
  # Inapplicable diagnostics that all need the same type name it once.
  expect_snapshot(
    check_positivity(
      data,
      exposure,
      c(x1, x2),
      diagnostics = c("hat_values", "hdr")
    ),
    error = TRUE
  )
  # Inapplicable diagnostics whose needed types a three-level factor cannot
  # carry, so no type is suggested.
  expect_snapshot(
    check_positivity(
      sim_pos_categorical(150),
      exposure,
      c(x1, x2),
      diagnostics = c("hat_values", "extrapolation")
    ),
    error = TRUE
  )
  # An exposure_type inside a per-diagnostic option list.
  expect_snapshot(
    check_positivity(
      data,
      exposure,
      c(x1, x2),
      diagnostics = "port",
      args = list(port = list(exposure_type = "categorical"))
    ),
    error = TRUE
  )
})

# ---- Missing exposure ------------------------------------------------------

test_that("a missing exposure value aborts before any diagnostic runs", {
  local_quiet()
  data <- dgp_good_positivity(n = 200, seed = 1)
  data$exposure[1] <- NA
  # Under auto detection the missing value would otherwise surface only when a
  # child fails, wrapped as a composition error.
  expect_error(
    check_positivity(data, exposure, c(x1, x2)),
    class = "positively_missing_error"
  )
  # An explicit binary type must reach the same missing-value guard rather than
  # misdetecting the column as categorical off the extra NA level.
  expect_error(
    check_positivity(data, exposure, c(x1, x2), exposure_type = "binary"),
    class = "positively_missing_error"
  )
})

test_that("the missing-exposure message is stable", {
  local_quiet()
  data <- dgp_good_positivity(n = 200, seed = 1)
  data$exposure[1] <- NA
  expect_snapshot(check_positivity(data, exposure, c(x1, x2)), error = TRUE)
})
