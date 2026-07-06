# Exposure-type detection mirrors the propensity idiom: match_exposure_type()
# resolves "auto" through detect_exposure_type(), which announces the detected
# type unless options(positively.quiet) suppresses it.

test_that("detect_exposure_type() identifies binary exposures", {
  withr::local_options(positively.quiet = TRUE)
  expect_identical(detect_exposure_type(c(0, 1, 1, 0)), "binary")
  expect_identical(detect_exposure_type(factor(c("a", "b", "a"))), "binary")
})

test_that("detect_exposure_type() identifies categorical exposures", {
  withr::local_options(positively.quiet = TRUE)
  expect_identical(
    detect_exposure_type(factor(c("a", "b", "c"))),
    "categorical"
  )
  expect_identical(
    detect_exposure_type(c("low", "med", "high", "low")),
    "categorical"
  )
})

test_that("detect_exposure_type() identifies continuous exposures", {
  withr::local_options(positively.quiet = TRUE)
  expect_identical(
    detect_exposure_type(seq(0, 1, length.out = 50)),
    "continuous"
  )
})

test_that("detect_exposure_type() announces the detected type by default", {
  withr::local_options(positively.quiet = FALSE)
  expect_message(detect_exposure_type(c(0, 1, 0, 1)), "Treating")
})

test_that("detect_exposure_type() is silent when positively.quiet is TRUE", {
  withr::local_options(positively.quiet = TRUE)
  expect_silent(detect_exposure_type(c(0, 1, 0, 1)))
})

test_that("match_exposure_type() honours an explicit type without detecting", {
  withr::local_options(positively.quiet = TRUE)
  expect_identical(match_exposure_type("continuous", c(0, 1)), "continuous")
})

test_that("match_exposure_type() detects the type when 'auto'", {
  withr::local_options(positively.quiet = TRUE)
  expect_identical(match_exposure_type("auto", c(0, 1, 0)), "binary")
})

test_that("match_exposure_type() rejects an unknown type", {
  expect_error(
    match_exposure_type("nonsense", c(0, 1)),
    regexp = "must be one of"
  )
})

test_that("a single-level factor is detected as binary", {
  withr::local_options(positively.quiet = TRUE)
  # A factor with one observed level is not caught by the exactly-two-values
  # rule, so it falls through the factor branch to binary rather than
  # categorical.
  expect_identical(detect_exposure_type(factor(c("a", "a", "a"))), "binary")
})

test_that("a 0/1 exposure with NA is detected as categorical", {
  withr::local_options(positively.quiet = TRUE)
  # The NA becomes a third unique value, so at a realistic n the unique-value
  # ratio sends the vector down the categorical branch. This deliberately
  # mirrors propensity's detection behavior; do not change one without the
  # other.
  exposure <- c(rep(0, 250), rep(1, 249), NA)
  expect_identical(detect_exposure_type(exposure), "categorical")
})

test_that("an all-NA exposure is detected as continuous", {
  withr::local_options(positively.quiet = TRUE)
  # With no non-missing observations the categorical heuristic returns FALSE, so
  # detection falls through to continuous. Locked to mirror propensity.
  expect_identical(
    detect_exposure_type(c(NA_real_, NA_real_, NA_real_)),
    "continuous"
  )
})

test_that("detect_exposure_type() announcement is stable", {
  withr::local_options(positively.quiet = FALSE)
  expect_snapshot(exposure_type <- detect_exposure_type(c(0, 1, 0, 1)))
})
