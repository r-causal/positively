# check_positivity() is the top-level entry point: it detects the exposure type, runs
# the default diagnostics for that type (overridable through `diagnostics`),
# threads per-method options from `args`, and returns a positivity_check
# container. The default sets are fixed by exposure type: binary runs edp, port,
# and extrapolation; categorical runs edp and port; continuous runs edp, port,
# hat_values, and hdr. check_eta_bias() and check_density_ratios() are never run
# because they need an outcome or user-supplied ratios. These specs assert the
# plumbing (selection, passthrough, container shape, errors), so they lean on the
# cheapest diagnostics and shrink the expensive null resampling through `args`.

# The child produced by a named diagnostic, read straight off the named @checks
# list rather than through `[[` or `$`, which several of these blocks test.
child_named <- function(res, name) {
  res@checks[[name]]
}

# ---- Default diagnostic sets by exposure type -----------------------------

test_that("a binary exposure runs edp, port, and extrapolation by default", {
  local_quiet()
  data <- dgp_good_positivity(n = 200)
  res <- check_positivity(data, exposure, c(x1, x2))

  expect_true(S7::S7_inherits(res, positivity_check))
  expect_identical(names(res), c("edp", "port", "extrapolation"))
  expect_identical(names(res@checks), c("edp", "port", "extrapolation"))
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
  expect_false(any(c("eta_bias", "density_ratios") %in% names(res)))
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
  data <- dgp_categorical(150)
  res <- check_positivity(data, exposure, c(x1, x2))

  expect_true(S7::S7_inherits(res, positivity_check))
  expect_identical(names(res), c("edp", "port"))
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
  expect_identical(names(res), c("edp", "port", "hat_values", "hdr"))
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
  expect_identical(names(res), c("port", "extrapolation"))
  expect_length(res@checks, 2L)
})

test_that("a single requested diagnostic yields a one-child container", {
  local_quiet()
  data <- dgp_good_positivity(n = 200)
  res <- check_positivity(data, exposure, c(x1, x2), diagnostics = "port")
  expect_identical(names(res), "port")
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

test_that("a declared type the exposure cannot carry aborts before any child runs", {
  local_quiet()
  data <- dgp_categorical(150)
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
  data <- dgp_categorical(150)
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
  data <- dgp_coarse_dose(n = 150)
  # Detection calls the coarse numeric exposure categorical, but the column is
  # numeric and so carries a declared continuous type. The declaration has to
  # reach the arithmetic and not only the child's metadata, so the shape of the
  # result is what says which type was computed on: check_port() bins a
  # continuous exposure into quantile ranges and reports one row per subgroup and
  # range, where a categorical exposure reports a row per distinct dose. Fewer
  # reported levels than distinct doses is therefore only reachable through the
  # continuous path.
  res <- check_positivity(
    data,
    exposure,
    x1,
    diagnostics = "port",
    exposure_type = "continuous"
  )
  port <- child_named(res, "port")
  expect_identical(port@exposure_type, "continuous")
  expect_lt(
    length(unique(port@results$exposure_level)),
    length(unique(data$exposure))
  )
  # The covariate tracks the dose, so extreme doses are implausible at some
  # covariate values and the run has structure to find. Both outcomes must
  # appear: an empty result, or one where every subgroup reads the same way,
  # would satisfy the type assertion above just as well while proving nothing was
  # computed.
  expect_true(any(port@results$low_support))
  expect_true(!all(port@results$low_support))
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

  expect_identical(names(res), c("edp", "port", "hat_values", "hdr"))
  expect_identical(
    unname(vapply(
      res@checks,
      function(check) check@exposure_type,
      character(1)
    )),
    rep("continuous", 4L)
  )
})

# ---- Forwarding the resolved exposure type --------------------------------

test_that("diagnostic_exposure_types() projects the registry in composition order", {
  # Subsetting a named list by a character vector is what keeps the projection
  # faithful, so names, order, and values are each asserted rather than assumed.
  # validate_composed_diagnostics() reads names, and exposure_type_advice()
  # reads values, so a projection that lost either would misreport rather than
  # fail loudly.
  offered <- diagnostic_exposure_types()
  registry <- diagnostic_supported_types()

  expect_identical(names(offered), composed_diagnostics())
  expect_identical(
    unname(offered),
    unname(registry[composed_diagnostics()])
  )
})

# ---- Container metadata ---------------------------------------------------

test_that("the container carries the metadata check_positivity() resolved", {
  local_quiet()
  data <- dgp_good_positivity(n = 200)
  res <- check_positivity(
    data,
    exposure,
    c(x1, x2),
    diagnostics = c("port", "extrapolation")
  )

  expect_identical(res@exposure, "exposure")
  # No type was declared, so the recorded type is the one detection resolved
  # rather than the "auto" the caller passed.
  expect_identical(res@exposure_type, "binary")
  expect_identical(res@covariates, c("x1", "x2"))
  expect_identical(res@n, 200L)
  expect_true(rlang::is_call(res@call, "check_positivity"))
})

test_that("a declared exposure type is what the container records", {
  local_quiet()
  data <- dgp_good_positivity(n = 200)
  res <- check_positivity(
    data,
    exposure,
    c(x1, x2),
    diagnostics = "port",
    exposure_type = "categorical"
  )
  expect_identical(res@exposure_type, "categorical")
})

test_that("the container records a lone covariate as a length-one vector", {
  local_quiet()
  data <- dgp_continuous_support_gap(n = 150)
  res <- check_positivity(data, exposure, x1, diagnostics = "port")

  expect_identical(res@covariates, "x1")
  expect_identical(res@exposure_type, "continuous")
  expect_identical(res@n, 150L)
})

test_that("the container metadata agrees with every child's", {
  local_quiet()
  data <- dgp_good_positivity(n = 200)
  res <- check_positivity(
    data,
    exposure,
    c(x1, x2),
    diagnostics = c("port", "extrapolation")
  )

  agrees <- function(reader) {
    vapply(
      res@checks,
      function(check) identical(reader(check), reader(res)),
      logical(1)
    )
  }
  expect_true(all(agrees(function(x) x@exposure)))
  expect_true(all(agrees(function(x) x@exposure_type)))
  expect_true(all(agrees(function(x) x@n)))
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
  expect_snapshot_abort(
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

test_that("tidy() and glance() are reachable without qualifying the generic", {
  local_quiet()
  data <- dgp_good_positivity(n = 200)
  res <- check_positivity(data, exposure, c(x1, x2), diagnostics = "port")
  # tidy() and glance() are re-exported, so the bare generics resolve on a
  # child. Neither is defined on the container.
  expect_s3_class(tidy(res$port), "tbl_df")
  expect_s3_class(glance(res$port), "tbl_df")
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

# ---- The container's removed broom methods --------------------------------

test_that("the container has no tidy() and no glance() method", {
  local_quiet()
  data <- dgp_good_positivity(n = 200)
  res <- check_positivity(data, exposure, c(x1, x2), diagnostics = "port")
  container_class <- class(res)[[1]]

  # Nothing is registered for the container, so both calls fail in dispatch
  # rather than in a method kept only to refuse. summary() is the container's
  # overview now.
  expect_null(utils::getS3method("tidy", container_class, optional = TRUE))
  expect_null(utils::getS3method("glance", container_class, optional = TRUE))
  expect_error(generics::tidy(res), "no applicable method")
  expect_error(generics::glance(res), "no applicable method")

  # The second shape, a wide results tibble for one named child, goes with it.
  expect_error(
    generics::tidy(res, diagnostic = "port"),
    "no applicable method"
  )

  # Only the container lost them. The children keep both.
  expect_s3_class(generics::tidy(res$port), "tbl_df")
  expect_s3_class(generics::glance(res$port), "tbl_df")
})

test_that("tidy() on an extracted child returns its results unchanged", {
  local_quiet()
  data <- dgp_good_positivity(n = 200)
  res <- check_positivity(
    data,
    exposure,
    c(x1, x2),
    diagnostics = c("port", "extrapolation")
  )
  expect_identical(generics::tidy(res$port), child_named(res, "port")@results)
  expect_identical(
    generics::tidy(res$extrapolation),
    child_named(res, "extrapolation")@results
  )
})

# ---- Container summary() ---------------------------------------------------

test_that("summary() on the container reports the statistics of each child", {
  local_quiet()
  data <- dgp_practical_violation(n = 200)
  res <- check_positivity(
    data,
    exposure,
    c(x1, x2),
    diagnostics = c("port", "extrapolation", "edp")
  )
  summarized <- summary(res)

  expect_s3_class(summarized, "tbl_df")
  expect_identical(
    names(summarized),
    c("diagnostic", "statistic", "value", "threshold")
  )
  expect_type(summarized$diagnostic, "character")
  expect_type(summarized$statistic, "character")
  expect_type(summarized$value, "double")
  expect_type(summarized$threshold, "double")

  # Children appear in the order they were requested.
  expect_identical(
    unique(summarized$diagnostic),
    c("port", "extrapolation", "edp")
  )

  # One row per statistic the diagnostic computed, in the order glance() states
  # them. The sample size and the resolved method parameters are not statistics:
  # the container states n once and the parameters reach the threshold column.
  expect_identical(
    summarized$statistic[summarized$diagnostic == "port"],
    c("n_subgroups", "n_low_support")
  )
  # The geometric variability is absent because prop_supported reports it as
  # its threshold, and a column that reaches the threshold column would
  # otherwise be stated twice.
  expect_identical(
    summarized$statistic[summarized$diagnostic == "extrapolation"],
    c("mean_frac_nearby", "prop_supported", "prop_in_hull")
  )
  expect_identical(
    summarized$statistic[summarized$diagnostic == "edp"],
    c("n_values", "edp_min", "edp_max")
  )
  expect_false(any(
    c("n", "n_rows", "exposure", "exposure_type", "alpha", "beta", "gamma") %in%
      summarized$statistic
  ))
})

test_that("the container's summary values come from the children's results", {
  local_quiet()
  data <- dgp_practical_violation(n = 200)
  res <- check_positivity(
    data,
    exposure,
    c(x1, x2),
    diagnostics = c("port", "extrapolation", "edp")
  )
  summarized <- summary(res)
  value_of <- function(diagnostic, statistic) {
    row <- summarized$diagnostic == diagnostic &
      summarized$statistic == statistic
    summarized$value[row]
  }

  # Every value is pinned against the child's own results tibble or scalar
  # properties rather than against glance(), so the assertion does not travel
  # with whatever glance() happens to report.
  port_results <- child_named(res, "port")@results
  expect_identical(
    value_of("port", "n_subgroups"),
    as.numeric(nrow(port_results))
  )
  expect_identical(
    value_of("port", "n_low_support"),
    as.numeric(sum(port_results$low_support))
  )

  extrapolation <- child_named(res, "extrapolation")

  # The geometric variability reaches the reader in the threshold column of the
  # row it cuts rather than as a value of its own, so it is pinned there.
  gv_row <- summarized$diagnostic == "extrapolation" &
    summarized$statistic == "prop_supported"
  expect_false(
    "geometric_variability" %in%
      summarized$statistic[summarized$diagnostic == "extrapolation"]
  )
  expect_equal(
    summarized$threshold[gv_row],
    extrapolation@geometric_variability
  )

  expect_equal(
    value_of("extrapolation", "mean_frac_nearby"),
    mean(extrapolation@results$frac_nearby)
  )
  expect_equal(
    value_of("extrapolation", "prop_supported"),
    mean(!extrapolation@results$low_support)
  )
  expect_equal(
    value_of("extrapolation", "prop_in_hull"),
    mean(extrapolation@results$in_hull)
  )

  edp <- child_named(res, "edp")
  expect_identical(
    value_of("edp", "n_values"),
    as.numeric(length(edp@params$values))
  )
  expect_equal(value_of("edp", "edp_min"), min(edp@results$edp))
  expect_equal(value_of("edp", "edp_max"), max(edp@results$edp))
})

test_that("the container's summary pairs a statistic with the cut behind it", {
  local_quiet()
  data <- dgp_practical_violation(n = 200)
  res <- check_positivity(
    data,
    exposure,
    c(x1, x2),
    diagnostics = c("port", "extrapolation", "edp")
  )
  summarized <- summary(res)
  threshold_of <- function(diagnostic, statistic) {
    row <- summarized$diagnostic == diagnostic &
      summarized$statistic == statistic
    summarized$threshold[row]
  }

  # beta is the prevalence a subgroup is compared against, so it is the cut that
  # produced the count of low-support subgroups.
  expect_identical(
    threshold_of("port", "n_low_support"),
    child_named(res, "port")@beta
  )

  # alpha and gamma set which subgroups are reported at all rather than cutting
  # a reported value, and no single one of them produced the count, so the
  # number of subgroups has no threshold.
  expect_identical(threshold_of("port", "n_subgroups"), NA_real_)

  # prop_supported is the fraction of units whose nearest opposite unit lies
  # within one geometric variability, so the cut behind it is that distance.
  expect_identical(
    threshold_of("extrapolation", "prop_supported"),
    child_named(res, "extrapolation")@geometric_variability
  )

  # The radius governing frac_nearby is in different units from the fraction
  # itself, and hull membership is a containment test with no numeric cut, so
  # neither states one.
  expect_identical(threshold_of("extrapolation", "mean_frac_nearby"), NA_real_)
  expect_identical(threshold_of("extrapolation", "prop_in_hull"), NA_real_)

  # The edp statistics are raw magnitudes with nothing behind them.
  is_edp <- summarized$diagnostic == "edp"
  expect_true(all(is.na(summarized$threshold[is_edp])))
})

test_that("summary() keeps value numeric where the statistics are not", {
  local_quiet()
  data <- dgp_continuous_support_gap(n = 150)
  res <- check_positivity(
    data,
    exposure,
    x1,
    exposure_type = "continuous",
    diagnostics = c("port", "hat_values", "hdr", "edp"),
    args = list(hat_values = list(null_reps = 25))
  )
  summarized <- summary(res)

  # This run's glance() output includes two character statistics, variant and
  # density_estimator, and one logical, exceeds_null. A numeric value column
  # cannot carry the first two without the character coercion this overview
  # exists to remove, so none of the three appear here. All keep their types in
  # glance(), which is wide.
  expect_type(summarized$value, "double")
  expect_false(any(
    c("variant", "density_estimator", "exceeds_null", "mass") %in%
      summarized$statistic
  ))
  expect_identical(
    summarized$statistic[summarized$diagnostic == "hdr"],
    c("n_values", "nonoverlap_min", "nonoverlap_max")
  )
  expect_identical(
    summarized$statistic[summarized$diagnostic == "hat_values"],
    c("phi_hat", "n_high_leverage", "n_candidates", "p")
  )
})

test_that("phi_hat carries the null quantile it is compared against", {
  local_quiet()
  data <- dgp_continuous_support_gap(n = 150)
  res <- check_positivity(
    data,
    exposure,
    x1,
    exposure_type = "continuous",
    diagnostics = "hat_values",
    args = list(hat_values = list(null_reps = 25))
  )
  summarized <- summary(res)
  hat_values <- child_named(res, "hat_values")
  phi_row <- summarized$statistic == "phi_hat"

  expect_equal(summarized$value[phi_row], hat_values@phi_hat)
  expect_equal(summarized$threshold[phi_row], hat_values@null_quantile)

  # exceeds_null is dropped as a logical, but nothing it said is lost: it is
  # exactly this value against this threshold.
  exceeds <- summarized$value[phi_row] > summarized$threshold[phi_row]
  expect_identical(exceeds, hat_values@exceeds_null)

  # The null quantile is the cut beside phi_hat, so it is not also repeated as
  # a statistic of its own.
  expect_false("null_quantile" %in% summarized$statistic)

  # The count of high-leverage candidates carries the per-candidate leverage
  # cutoff, threshold * p / n, stated in the units of a hat value.
  count_row <- summarized$statistic == "n_high_leverage"
  cutoff <- hat_values@params$threshold * hat_values@p / hat_values@n
  expect_equal(summarized$threshold[count_row], cutoff)
  expect_identical(
    summarized$value[count_row],
    as.numeric(sum(hat_values@results$hat_value > cutoff))
  )
})

test_that("a child's summary is the container's rows for that child", {
  local_quiet()
  data <- dgp_practical_violation(n = 200)
  res <- check_positivity(
    data,
    exposure,
    c(x1, x2),
    diagnostics = c("port", "edp")
  )
  overview <- summary(res)

  # summary() is one operation at two levels, so the container's rows for a
  # child are that child's own summary with the diagnostic named.
  port_rows <- overview[overview$diagnostic == "port", , drop = FALSE]
  port_summary <- summary(res$port)
  expect_identical(
    names(port_summary),
    c("statistic", "value", "threshold")
  )
  expect_identical(port_summary$statistic, port_rows$statistic)
  expect_identical(port_summary$value, port_rows$value)
  expect_identical(port_summary$threshold, port_rows$threshold)

  edp_rows <- overview[overview$diagnostic == "edp", , drop = FALSE]
  edp_summary <- summary(res$edp)
  expect_identical(edp_summary$statistic, edp_rows$statistic)
  expect_identical(edp_summary$value, edp_rows$value)
  expect_identical(edp_summary$threshold, edp_rows$threshold)
})

# ---- Container print() -----------------------------------------------------

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

test_that("`$` extracts a child diagnostic instead of reaching for a property", {
  local_quiet()
  data <- dgp_good_positivity(n = 200)
  res <- check_positivity(
    data,
    exposure,
    c(x1, x2),
    diagnostics = c("port", "extrapolation")
  )
  expect_identical(res$port, child_named(res, "port"))
  expect_identical(res$extrapolation, res[["extrapolation"]])
  expect_error(res$nonexistent, class = "positively_diagnostic_error")
  # Property names are refused too, so `$` never silently returns metadata in
  # place of a child.
  expect_error(res$checks, class = "positively_diagnostic_error")
  expect_error(res$exposure, class = "positively_diagnostic_error")
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
  empty <- positivity_check(
    checks = list(),
    exposure = "exposure",
    exposure_type = "binary",
    covariates = c("x1", "x2"),
    n = 200L,
    call = quote(check_positivity())
  )
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
  data <- dgp_categorical(150)
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
  expect_identical(names(res), c("hat_values", "hdr"))
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
  skip_if_not_installed("lpSolve")
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
  expect_snapshot_abort(check_positivity(
    data,
    exposure,
    c(x1, x2),
    diagnostics = "bogus"
  ))
  # Diagnostics that need an outcome or user-supplied ratios.
  expect_snapshot_abort(check_positivity(
    data,
    exposure,
    c(x1, x2),
    diagnostics = "eta_bias"
  ))
  expect_snapshot_abort(check_positivity(
    data,
    exposure,
    c(x1, x2),
    diagnostics = "density_ratios"
  ))
  # A diagnostic that does not apply to the exposure type. The 0/1 exposure is
  # numeric, so the continuous type hat_values needs is on offer; the continuous
  # exposure has far more than two distinct values, so the binary type
  # extrapolation needs is not.
  expect_snapshot_abort(check_positivity(
    data,
    exposure,
    c(x1, x2),
    diagnostics = "hat_values"
  ))
  expect_snapshot_abort(check_positivity(
    continuous,
    exposure,
    x1,
    diagnostics = "extrapolation"
  ))
  # An args name that is not being run.
  expect_snapshot_abort(check_positivity(
    data,
    exposure,
    c(x1, x2),
    diagnostics = "port",
    args = list(bogus = list(alpha = 0.2))
  ))
  # A duplicated diagnostic name.
  expect_snapshot_abort(check_positivity(
    data,
    exposure,
    c(x1, x2),
    diagnostics = c("port", "port")
  ))
  # An args entry that is not a list.
  expect_snapshot_abort(check_positivity(
    data,
    exposure,
    c(x1, x2),
    diagnostics = "port",
    args = list(port = "oops")
  ))
  # A declared exposure type the exposure column cannot carry.
  expect_snapshot_abort(check_positivity(
    dgp_categorical(150),
    exposure,
    c(x1, x2),
    diagnostics = c("port", "extrapolation"),
    exposure_type = "binary"
  ))
  # A failing child names the diagnostic and chains the child's condition.
  expect_snapshot_abort(check_positivity(
    data,
    exposure,
    c(x1, x2),
    diagnostics = "port",
    args = list(port = list(alpha = 1.5))
  ))
  # A single-column exposure selection is required.
  expect_snapshot_abort(check_positivity(data, c(exposure, x1), x2))
  # A covariate selection that includes the exposure column.
  expect_snapshot_abort(check_positivity(
    data,
    exposure,
    tidyselect::everything()
  ))
  # [[ rejects an unknown diagnostic name.
  res <- check_positivity(data, exposure, c(x1, x2), diagnostics = "port")
  expect_snapshot_abort(res[["nonexistent"]])
  # [[ rejects an out-of-bounds numeric index.
  expect_snapshot_abort(res[[5]])
  # [[ rejects a non-scalar index.
  expect_snapshot_abort(res[[c("port", "extrapolation")]])
  # Inapplicable diagnostics that all need the same type name it once.
  expect_snapshot_abort(check_positivity(
    data,
    exposure,
    c(x1, x2),
    diagnostics = c("hat_values", "hdr")
  ))
  # Inapplicable diagnostics whose needed types a three-level factor cannot
  # carry, so no type is suggested.
  expect_snapshot_abort(check_positivity(
    dgp_categorical(150),
    exposure,
    c(x1, x2),
    diagnostics = c("hat_values", "extrapolation")
  ))
  # An exposure_type inside a per-diagnostic option list.
  expect_snapshot_abort(check_positivity(
    data,
    exposure,
    c(x1, x2),
    diagnostics = "port",
    args = list(port = list(exposure_type = "categorical"))
  ))
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
  expect_snapshot_abort(check_positivity(data, exposure, c(x1, x2)))
})

# ---- Container plotting ----------------------------------------------------

# A binary run holding three children, so the panel has more than one view to
# compose and a selected child is distinguishable from it.
plot_container <- function() {
  check_positivity(
    dgp_practical_violation(n = 200, seed = 2),
    exposure,
    c(x1, x2)
  )
}

test_that("autoplot() is part of the public surface", {
  expect_true("autoplot" %in% getNamespaceExports("positively"))
})

test_that("autoplot() panels every child", {
  local_quiet()
  # The block draws and also skips on a package, so the figure is announced
  # before the skip. An unannounced file looks unused to testthat's cleanup pass
  # and is deleted.
  announce_doppelganger("positivity check panel")
  skip_if_not_installed("patchwork")
  expect_doppelganger("positivity check panel", autoplot(plot_container()))
})

test_that("the panel is a ggplot", {
  local_quiet()
  skip_if_not_installed("patchwork")
  panel <- autoplot(plot_container())

  # Composing children must not cost the reader the ability to keep building on
  # the result.
  expect_true(ggplot2::is_ggplot(panel))
})

test_that("autoplot() selects a single child by name", {
  local_quiet()
  container <- plot_container()
  selected <- autoplot(container, "port")

  # The name is the one `names()` lists and `$` extracts by, so a reader can go
  # from a section of the report straight to its plot. Compared on what the
  # plots draw rather than on the objects, because an aesthetic carries the
  # environment it was captured in and the two paths capture different ones.
  expect_equal(
    ggplot2::ggplot_build(selected)$data,
    ggplot2::ggplot_build(ggplot2::autoplot(container$port))$data
  )
  expect_false(inherits(selected, "patchwork"))
})

test_that("autoplot() rejects a name the container does not hold", {
  local_quiet()
  container <- plot_container()
  expect_snapshot_abort(
    autoplot(container, "hdr"),
    class = "positively_diagnostic_error"
  )
})

test_that("plot() draws what autoplot() returns and gives back the container", {
  local_quiet()
  local_null_device()
  skip_if_not_installed("patchwork")
  container <- plot_container()

  expect_identical(plot(container), container)
  expect_identical(plot(container, "port"), container)
})

test_that("a container holding one child still plots", {
  local_quiet()
  local_null_device()
  skip_if_not_installed("patchwork")
  container <- check_positivity(
    dgp_practical_violation(n = 200, seed = 2),
    exposure,
    c(x1, x2),
    diagnostics = "port"
  )
  panel <- autoplot(container)

  expect_true(ggplot2::is_ggplot(panel))
  expect_no_error(print(panel))
})

test_that("autoplot() says what it needs when patchwork is absent", {
  local_quiet()
  # The container is built before the mock, because resolving the exposure type
  # and the hull test both ask whether a package is installed.
  container <- plot_container()

  # Mocked rather than skipped, so the reader without patchwork sees a message
  # this suite has read, rather than one that fires only on machines it never
  # runs on.
  local_mocked_bindings(is_installed = function(...) FALSE, .package = "rlang")
  expect_snapshot_abort(
    autoplot(container),
    class = "positively_missing_package_error"
  )

  # Naming a diagnostic composes nothing, so it is the half of the behavior that
  # still works, which is what the refusal above points the reader at.
  expect_true(ggplot2::is_ggplot(autoplot(container, "port")))
})

test_that("an empty container says there is nothing to plot", {
  # The class permits a container with no children, and both entry points reach
  # the same refusal rather than failing on an absent first child.
  container <- positivity_check(
    checks = list(),
    exposure = "exposure",
    exposure_type = "binary",
    covariates = c("x1", "x2"),
    n = 200L,
    call = quote(check_positivity())
  )

  expect_snapshot_abort(autoplot(container), class = "positively_empty_error")
  expect_snapshot_abort(plot(container), class = "positively_empty_error")
})

test_that("a named diagnostic takes its own arguments", {
  local_quiet()
  container <- plot_container()

  # The argument reaches the child rather than being dropped, so the view it
  # selects differs from the child's default.
  hull <- autoplot(container, "extrapolation", type = "hull")
  expect_true(ggplot2::is_ggplot(hull))
  expect_false(identical(
    ggplot2::ggplot_build(hull)$data,
    ggplot2::ggplot_build(autoplot(container, "extrapolation"))$data
  ))
})

test_that("the panel passes its arguments to every child", {
  local_quiet()
  skip_if_not_installed("patchwork")
  container <- plot_container()

  # An argument only extrapolation accepts reaches EDP as well, which refuses
  # it. This is why a per-diagnostic argument belongs with a named diagnostic.
  expect_error(
    autoplot(container, type = "hull"),
    class = "positively_args_error"
  )
})
