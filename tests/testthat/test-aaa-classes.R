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

test_that("positivity_check stores validated diagnostic children", {
  obj <- make_test_diagnostic()
  container <- positivity_check(checks = list(obj), diagnostics = "test")
  expect_length(container@checks, 1)
  expect_identical(container@diagnostics, "test")
})

test_that("positivity_check rejects non-diagnostic children", {
  expect_error(
    positivity_check(checks = list(1L), diagnostics = "bad"),
    regexp = "positivity_diagnostic"
  )
})

test_that("printing a positivity_check produces sectioned output", {
  obj <- make_test_diagnostic()
  container <- positivity_check(checks = list(obj), diagnostics = "test")
  expect_output(print(container), "test")
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

test_that("positivity_check rejects a checks/diagnostics length mismatch", {
  obj <- make_test_diagnostic()
  expect_error(
    positivity_check(checks = list(obj, obj), diagnostics = "only_one"),
    regexp = "one name per element"
  )
})

test_that("printing a diagnostic is stable", {
  expect_snapshot(print(make_test_diagnostic()))
})

test_that("printing a container with two children is stable", {
  container <- positivity_check(
    checks = list(make_test_diagnostic(), make_test_diagnostic()),
    diagnostics = c("edp", "port")
  )
  expect_snapshot(print(container))
})
