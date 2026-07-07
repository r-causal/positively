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
# categorical, yet a continuous override is a legitimate request check_port()
# honours directly.
sim_pos_coarse_numeric <- function(n = 150, seed = 1) {
  withr::local_seed(seed)
  tibble::tibble(
    exposure = sample(1:8, n, replace = TRUE),
    x1 = stats::rnorm(n),
    x2 = stats::rnorm(n)
  )
}

test_that("a forced type infeasible for a re-detecting child aborts up front", {
  local_quiet()
  data <- sim_pos_categorical(150)
  # The exposure is a three-level factor (detected categorical). Forcing binary
  # lets extrapolation pass the applicability check, but extrapolation re-detects
  # categorical, so the request must abort before any child runs.
  expect_error(
    check_positivity(
      data,
      exposure,
      c(x1, x2),
      diagnostics = c("port", "extrapolation"),
      exposure_type = "binary"
    ),
    class = "positively_error"
  )
})

test_that("a forced type a non-re-detecting child accepts is allowed", {
  local_quiet()
  data <- sim_pos_coarse_numeric(150)
  # Detection calls the coarse numeric exposure categorical, but check_port()
  # accepts a continuous override directly, so the request must succeed.
  res <- check_positivity(
    data,
    exposure,
    c(x1, x2),
    diagnostics = "port",
    exposure_type = "continuous"
  )
  expect_identical(child_named(res, "port")@exposure_type, "continuous")
})

# ---- Child alerts and failures --------------------------------------------

test_that("a child informational alert survives the entry point", {
  withr::local_options(
    positively.quiet = FALSE,
    positively.gower_chunk_threshold = 20
  )
  data <- dgp_good_positivity(n = 60)
  # The chunked-Gower notice comes from within check_extrapolation(); it must
  # still reach the user even though the repeated detection message is muffled.
  expect_message(
    check_positivity(data, exposure, c(x1, x2), diagnostics = "extrapolation"),
    "row chunks"
  )
})

test_that("the detection message is announced exactly once at default verbosity", {
  withr::local_options(positively.quiet = FALSE)
  data <- dgp_good_positivity(n = 60)
  # extrapolation re-detects the exposure type internally, so without muffling
  # the announcement would appear twice: once from the entry point and once from
  # the child.
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
  # A diagnostic that does not apply to the exposure type.
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
  # A forced exposure_type that sends a requested re-detecting child off a cliff.
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
})
