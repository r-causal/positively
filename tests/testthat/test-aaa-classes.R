# The abstract positivity_diagnostic parent cannot be instantiated, so a
# throwaway concrete subclass lets these specs exercise the inherited
# properties, validator, and shared tidy()/glance()/print wiring.
make_test_diagnostic <- function(
  results = tibble::tibble(.id = 1:3, value = c(0.1, 0.2, 0.3)),
  exposure = "a",
  exposure_type = "binary",
  n = 3L,
  params = list(alpha = 0.05),
  call = quote(check_test())
) {
  test_diagnostic <- S7::new_class(
    "test_diagnostic",
    parent = positivity_diagnostic
  )
  test_diagnostic(
    results = results,
    exposure = exposure,
    exposure_type = exposure_type,
    n = n,
    params = params,
    call = call
  )
}

# positivity_check carries its own metadata alongside a named list of children,
# so a hand-built container needs all six values. These defaults stand in for
# the ones a block is not exercising; every block overrides what it asserts on.
make_test_container <- function(
  checks = list(test = make_test_diagnostic()),
  exposure = "a",
  exposure_type = "binary",
  covariates = c("x1", "x2"),
  n = 3L,
  call = quote(check_positivity())
) {
  positivity_check(
    checks = checks,
    exposure = exposure,
    exposure_type = exposure_type,
    covariates = covariates,
    n = n,
    call = call
  )
}

test_that("positivity_diagnostic is an S7 class", {
  expect_true(S7::S7_inherits(make_test_diagnostic(), positivity_diagnostic))
})

test_that("positivity_diagnostic is abstract and cannot be instantiated", {
  expect_error(
    positivity_diagnostic(
      results = tibble::tibble(),
      exposure = "a",
      exposure_type = "binary",
      n = 0L,
      params = list(),
      call = quote(f())
    ),
    regexp = "abstract"
  )
})

test_that("positivity_diagnostic declares the fixed property set", {
  expect_setequal(
    names(positivity_diagnostic@properties),
    c("results", "exposure", "exposure_type", "n", "params", "call")
  )
})

test_that("the parent validator accepts a valid exposure_type", {
  obj <- make_test_diagnostic(exposure_type = "continuous")
  expect_identical(obj@exposure_type, "continuous")
})

test_that("the parent validator rejects an unknown exposure_type", {
  expect_error(
    make_test_diagnostic(exposure_type = "ordinal"),
    regexp = "exposure_type"
  )
})

test_that("tidy() returns the results tibble of a diagnostic", {
  results <- tibble::tibble(.id = 1:2, value = c(0.4, 0.6))
  obj <- make_test_diagnostic(results = results)
  expect_identical(generics::tidy(obj), results)
})

test_that("glance() returns a one-row tibble for a diagnostic", {
  glanced <- generics::glance(make_test_diagnostic())
  expect_s3_class(glanced, "tbl_df")
  expect_identical(nrow(glanced), 1L)
})

test_that("the parent validator rejects a non-scalar exposure_type", {
  expect_error(
    make_test_diagnostic(exposure_type = c("binary", "continuous")),
    regexp = "exposure_type"
  )
})

test_that("the parent validator rejects an empty exposure_type", {
  expect_error(
    make_test_diagnostic(exposure_type = character(0)),
    regexp = "exposure_type"
  )
})

test_that("the parent validator rejects a non-scalar n", {
  expect_error(
    make_test_diagnostic(n = c(1L, 2L)),
    regexp = "@n"
  )
})

test_that("the parent validator rejects a missing n", {
  expect_error(
    make_test_diagnostic(n = NA_integer_),
    regexp = "@n"
  )
})

test_that("the parent validator rejects an empty exposure", {
  expect_error(
    make_test_diagnostic(exposure = character(0)),
    regexp = "@exposure"
  )
})

# ---- The container's property set and metadata ----------------------------

test_that("positivity_check declares its property set", {
  expect_setequal(
    names(positivity_check@properties),
    c("checks", "exposure", "exposure_type", "covariates", "n", "call")
  )
})

test_that("positivity_check stores validated diagnostic children under their names", {
  obj <- make_test_diagnostic()
  container <- make_test_container(checks = list(test = obj))
  expect_length(container@checks, 1)
  expect_identical(names(container@checks), "test")
  expect_identical(container@checks[["test"]], obj)
})

test_that("positivity_check carries the metadata it was constructed with", {
  originating_call <- quote(check_positivity(df, dose, c(x1, x2, region)))
  container <- make_test_container(
    exposure = "dose",
    exposure_type = "continuous",
    covariates = c("x1", "x2", "region"),
    n = 500L,
    call = originating_call
  )
  expect_identical(container@exposure, "dose")
  expect_identical(container@exposure_type, "continuous")
  expect_identical(container@covariates, c("x1", "x2", "region"))
  expect_identical(container@n, 500L)
  expect_identical(container@call, originating_call)
})

# ---- The container validator ----------------------------------------------

test_that("positivity_check rejects non-diagnostic children", {
  expect_error(
    make_test_container(checks = list(bad = 1L)),
    regexp = "positivity_diagnostic"
  )
})

test_that("positivity_check rejects an unnamed checks list", {
  obj <- make_test_diagnostic()
  expect_error(
    make_test_container(checks = list(obj, obj)),
    regexp = "named"
  )
})

test_that("positivity_check rejects an empty or missing diagnostic name", {
  obj <- make_test_diagnostic()

  blank <- list(obj, obj)
  names(blank) <- c("edp", "")
  expect_error(make_test_container(checks = blank), regexp = "named")

  missing_name <- list(obj, obj)
  names(missing_name) <- c("edp", NA_character_)
  expect_error(make_test_container(checks = missing_name), regexp = "named")
})

test_that("positivity_check rejects duplicate diagnostic names", {
  obj <- make_test_diagnostic()
  repeated <- list(obj, obj)
  names(repeated) <- c("port", "port")
  expect_error(make_test_container(checks = repeated), regexp = "repeat")
})

test_that("the container validator rejects an unknown exposure_type", {
  expect_error(
    make_test_container(exposure_type = "ordinal"),
    regexp = "@exposure_type"
  )
})

test_that("the container validator rejects a non-scalar exposure_type", {
  expect_error(
    make_test_container(exposure_type = c("binary", "continuous")),
    regexp = "@exposure_type"
  )
})

test_that("the container validator rejects an empty exposure", {
  expect_error(
    make_test_container(exposure = character(0)),
    regexp = "@exposure"
  )
})

test_that("the container validator rejects empty covariates", {
  expect_error(
    make_test_container(covariates = character(0)),
    regexp = "@covariates"
  )
})

test_that("the container validator rejects a non-scalar n", {
  expect_error(
    make_test_container(n = c(1L, 2L)),
    regexp = "@n"
  )
})

test_that("the container validator rejects a missing n", {
  expect_error(
    make_test_container(n = NA_integer_),
    regexp = "@n"
  )
})

# ---- Naming and extracting children ---------------------------------------

test_that("names() returns the names of the checks list", {
  res <- make_test_container(
    checks = list(edp = make_test_diagnostic(), port = make_test_diagnostic())
  )
  expect_identical(names(res), c("edp", "port"))
  expect_identical(names(res), names(res@checks))
})

test_that("names() on a container with no children is character(0)", {
  expect_identical(names(make_test_container(checks = list())), character(0))
})

test_that("`[[` extracts a diagnostic by name and by whole-number position", {
  edp <- make_test_diagnostic()
  port <- make_test_diagnostic()
  res <- make_test_container(checks = list(edp = edp, port = port))
  expect_identical(res[["port"]], port)
  expect_identical(res[[2]], port)
  expect_identical(res[[2]], res[["port"]])
})

test_that("`[[` rejects invalid non-character indices", {
  res <- make_test_container(
    checks = list(edp = make_test_diagnostic(), port = make_test_diagnostic())
  )
  expect_error(res[[NA]], class = "positively_error")
  expect_error(res[[TRUE]], class = "positively_error")
  expect_error(res[[1.9]], class = "positively_error")
})

test_that("`[[` rejects a non-scalar index", {
  res <- make_test_container(
    checks = list(edp = make_test_diagnostic(), port = make_test_diagnostic())
  )
  expect_error(res[[c("edp", "port")]], class = "positively_diagnostic_error")
})

test_that("`[[` rejects an unknown diagnostic name", {
  res <- make_test_container(
    checks = list(edp = make_test_diagnostic(), port = make_test_diagnostic())
  )
  expect_error(res[["nonexistent"]], class = "positively_diagnostic_error")
})

test_that("`[[` rejects an out-of-bounds position", {
  res <- make_test_container(
    checks = list(edp = make_test_diagnostic(), port = make_test_diagnostic())
  )
  expect_error(res[[3]], class = "positively_bounds_error")
})

test_that("`$` extracts a diagnostic by name", {
  edp <- make_test_diagnostic()
  port <- make_test_diagnostic()
  res <- make_test_container(checks = list(edp = edp, port = port))
  expect_identical(res$port, port)
  expect_identical(res$port, res[["port"]])
})

test_that("`$` rejects an unknown diagnostic name rather than returning NULL", {
  res <- make_test_container(
    checks = list(edp = make_test_diagnostic(), port = make_test_diagnostic())
  )
  expect_error(res$nonexistent, class = "positively_diagnostic_error")
})

test_that("`$` does not partially match a diagnostic name", {
  res <- make_test_container(
    checks = list(edp = make_test_diagnostic(), port = make_test_diagnostic())
  )
  expect_error(res$por, class = "positively_diagnostic_error")
})

# Completion after `check$` runs through .DollarNames.default, which matches the
# supplied pattern against `names(x)`. The container defines no .DollarNames
# method of its own, so what completion offers depends entirely on the
# container's `names()` method. Completion passes a `^`-anchored regular
# expression rather than a bare prefix, so the pattern is matched as a regexp.
test_that(".DollarNames offers the diagnostic names and honours a pattern", {
  res <- make_test_container(
    checks = list(edp = make_test_diagnostic(), port = make_test_diagnostic())
  )
  expect_identical(utils::.DollarNames(res, ""), c("edp", "port"))
  expect_identical(utils::.DollarNames(res, "^p"), "port")
  expect_identical(utils::.DollarNames(res, "^z"), character(0))
})

# ---- Printing and pinned messages -----------------------------------------

test_that("printing a positivity_check produces sectioned output", {
  container <- make_test_container(checks = list(test = make_test_diagnostic()))
  expect_output(print(container), "test")
})

test_that("the invalid-index message is stable", {
  res <- make_test_container(
    checks = list(edp = make_test_diagnostic(), port = make_test_diagnostic())
  )
  expect_snapshot_abort(res[[1.9]])
})

test_that("the unknown-diagnostic message is stable", {
  res <- make_test_container(
    checks = list(edp = make_test_diagnostic(), port = make_test_diagnostic())
  )
  expect_snapshot_abort(res[["nonexistent"]])
  expect_snapshot_abort(res$nonexistent)
})

test_that("printing a diagnostic is stable", {
  expect_snapshot(print(make_test_diagnostic()))
})

test_that("printing a container with two children is stable", {
  container <- make_test_container(
    checks = list(edp = make_test_diagnostic(), port = make_test_diagnostic())
  )
  expect_snapshot(print(container))
})
