# Exposure-type detection mirrors the propensity idiom: detect_exposure_type()
# announces the type it infers unless options(positively.quiet) suppresses it.
# resolve_exposure_type() wraps that detection in the package's policy:
# detection supplies a default, an explicit exposure_type is authoritative, and
# an explicit type fails only when that type's math cannot run on the column as
# given. Under "auto" the detected type must still land in the diagnostic's
# supported set, so an unsupported detection aborts. Structure is validated on
# both paths, so a detected type and a declared type fail for the same reasons.

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

test_that("resolve_exposure_type() detects and announces when given 'auto'", {
  withr::local_options(positively.quiet = FALSE)
  expect_message(
    resolved <- resolve_exposure_type(
      "auto",
      c(0, 1, 0, 1),
      supported = c("binary", "categorical", "continuous"),
      fn = "check_edp"
    ),
    "Treating"
  )
  expect_identical(resolved, "binary")
})

test_that("resolve_exposure_type() suppresses the announcement on request", {
  withr::local_options(positively.quiet = FALSE)
  # The announcement names `.exposure`, so a caller resolving several exposure
  # columns in one call would emit a misnamed message per column. Such callers
  # pass announce = FALSE; the default announces.
  dose <- seq(0, 1, length.out = 50)

  expect_message(
    announced <- resolve_exposure_type(
      "auto",
      dose,
      supported = "continuous",
      fn = "check_hdr_seq"
    ),
    "Treating"
  )
  expect_identical(announced, "continuous")

  silent <- expect_silent(
    resolve_exposure_type(
      "auto",
      dose,
      supported = "continuous",
      fn = "check_hdr_seq",
      announce = FALSE
    )
  )
  expect_identical(silent, "continuous")
})

test_that("resolve_exposure_type() honors an explicit type without announcing", {
  withr::local_options(positively.quiet = FALSE)
  resolved <- expect_silent(
    resolve_exposure_type(
      "continuous",
      c(0, 1, 0, 1),
      supported = c("binary", "categorical", "continuous"),
      fn = "check_edp"
    )
  )
  expect_identical(resolved, "continuous")
})

test_that("resolve_exposure_type() keeps a type the heuristic would not pick", {
  withr::local_options(positively.quiet = FALSE)
  resolved <- expect_silent(
    resolve_exposure_type(
      "categorical",
      seq(0, 1, length.out = 50),
      supported = c("binary", "categorical", "continuous"),
      fn = "check_edp"
    )
  )
  expect_identical(resolved, "categorical")
})

test_that("an explicit continuous type overrides the coarse-numeric heuristic", {
  withr::local_options(positively.quiet = FALSE)
  # A genuinely continuous dose recorded at eight distinct levels over 150
  # observations: the unique-value ratio is 0.053, so detection reads it as
  # categorical. Declaring the type is the only way to run continuous math on a
  # column measured this coarsely, and the declaration must win outright rather
  # than be checked against the heuristic that misreads it.
  dose <- rep(seq_len(8), length.out = 150)

  expect_message(
    detected <- resolve_exposure_type(
      "auto",
      dose,
      supported = c("binary", "categorical", "continuous"),
      fn = "check_edp"
    ),
    "Treating"
  )
  expect_identical(detected, "categorical")

  resolved <- expect_silent(
    resolve_exposure_type(
      "continuous",
      dose,
      supported = c("binary", "categorical", "continuous"),
      fn = "check_edp"
    )
  )
  expect_identical(resolved, "continuous")
})

test_that("resolve_exposure_type() rejects a type outside the supported set", {
  withr::local_options(positively.quiet = TRUE)
  err <- expect_error(
    resolve_exposure_type(
      "continuous",
      c(0, 1, 0, 1),
      supported = "binary",
      fn = "check_extrapolation"
    )
  )
  expect_match(conditionMessage(err), "exposure_type", fixed = TRUE)
  expect_match(conditionMessage(err), "must be one of", fixed = TRUE)
})

test_that("resolve_exposure_type() aborts when a detected type is unsupported", {
  withr::local_options(positively.quiet = TRUE)
  err <- expect_error(
    resolve_exposure_type(
      "auto",
      c(0, 1, 0, 1),
      supported = "continuous",
      fn = "check_hdr"
    ),
    class = "positively_exposure_type_error"
  )
  message <- conditionMessage(err)
  expect_match(message, "check_hdr", fixed = TRUE)
  expect_match(message, "exposures only", fixed = TRUE)
  expect_match(message, "detected as", fixed = TRUE)
  # The hint names the escape hatch, so a user whose data the heuristic reads
  # differently can declare the type instead of reshaping the column.
  expect_match(message, "exposure_type", fixed = TRUE)
})

test_that("an explicit continuous type aborts on a character exposure", {
  withr::local_options(positively.quiet = TRUE)
  expect_error(
    resolve_exposure_type(
      "continuous",
      c("a", "b", "c"),
      supported = c("binary", "categorical", "continuous"),
      fn = "check_edp"
    ),
    class = "positively_exposure_type_error"
  )
})

test_that("an explicit continuous type aborts on a factor exposure", {
  withr::local_options(positively.quiet = TRUE)
  # as.double(factor(c("a", "b", "c"))) silently returns 1 2 3, so a check that
  # coerces first would read a factor as a usable continuous exposure.
  expect_error(
    resolve_exposure_type(
      "continuous",
      factor(c("a", "b", "c")),
      supported = c("binary", "categorical", "continuous"),
      fn = "check_edp"
    ),
    class = "positively_exposure_type_error"
  )
})

test_that("an explicit binary type aborts on three distinct values", {
  withr::local_options(positively.quiet = TRUE)
  expect_error(
    resolve_exposure_type(
      "binary",
      c(0, 1, 2),
      supported = c("binary", "categorical", "continuous"),
      fn = "check_edp"
    ),
    class = "positively_exposure_type_error"
  )
})

test_that("an explicit binary type accepts a factor with an unused level", {
  withr::local_options(positively.quiet = FALSE)
  # binary_to_01() re-levels with factor() before modeling, so an unused level
  # never reaches the math; the two observed levels are all it needs.
  exposure <- factor(c("a", "b", "a"), levels = c("a", "b", "c"))
  resolved <- expect_silent(
    resolve_exposure_type(
      "binary",
      exposure,
      supported = c("binary", "categorical", "continuous"),
      fn = "check_edp"
    )
  )
  expect_identical(resolved, "binary")
})

test_that("the auto path validates structure as strictly as an explicit type", {
  withr::local_options(positively.quiet = TRUE)
  # A single-level factor is detected as binary, but binary math needs two
  # distinct values, so the detected type does not fit the column it came from.
  # Validating the auto path turns that into the same classed abort an explicit
  # "binary" raises, instead of leaving a diagnostic to compute against an empty
  # opposite group.
  expect_error(
    resolve_exposure_type(
      "auto",
      factor(rep("a", 60)),
      supported = c("binary", "categorical", "continuous"),
      fn = "check_extrapolation"
    ),
    class = "positively_exposure_type_error"
  )
})

test_that("the auto path rejects date-like and logical exposures", {
  withr::local_options(positively.quiet = TRUE)
  # Dates, date-times, differences, and logicals store numbers underneath but
  # are not is.numeric(), so detection reads them as continuous while the
  # numeric requirement turns them away. Requiring is.numeric() rather than
  # coercibility is what keeps a factor from passing as a dose, and these
  # columns are the cost of that rule.
  exotic <- list(
    as.Date("2020-01-01") + 0:49,
    as.POSIXct("2020-01-01", tz = "UTC") + seq(0, 4900, by = 100),
    as.difftime(seq_len(50), units = "days"),
    c(TRUE, FALSE, NA)
  )
  detected <- vapply(exotic, detect_exposure_type, character(1))
  expect_identical(detected, rep("continuous", length(exotic)))

  for (exposure in exotic) {
    expect_error(
      resolve_exposure_type(
        "auto",
        exposure,
        supported = c("binary", "categorical", "continuous"),
        fn = "check_edp"
      ),
      class = "positively_exposure_type_error"
    )
  }
})

test_that("validate_exposure_structure() requires a numeric continuous exposure", {
  expect_no_error(
    validate_exposure_structure(
      seq(0, 1, length.out = 10),
      "continuous",
      fn = "check_edp"
    )
  )
  expect_error(
    validate_exposure_structure(
      c("a", "b", "c"),
      "continuous",
      fn = "check_edp"
    ),
    class = "positively_exposure_type_error"
  )
  expect_error(
    validate_exposure_structure(
      factor(c("a", "b", "c")),
      "continuous",
      fn = "check_edp"
    ),
    class = "positively_exposure_type_error"
  )
})

test_that("validate_exposure_structure() requires two values for binary", {
  expect_no_error(
    validate_exposure_structure(c(0, 1, 1, 0), "binary", fn = "check_edp")
  )
  expect_no_error(
    validate_exposure_structure(
      factor(c("a", "b", "a"), levels = c("a", "b", "c")),
      "binary",
      fn = "check_edp"
    )
  )
  # Missing values are not a level: an NA-bearing binary column has two distinct
  # values here, and its missingness is the diagnostic's own check to report.
  expect_no_error(
    validate_exposure_structure(c(0, 1, NA, 1), "binary", fn = "check_edp")
  )
  expect_error(
    validate_exposure_structure(c(0, 1, 2), "binary", fn = "check_edp"),
    class = "positively_exposure_type_error"
  )
  expect_error(
    validate_exposure_structure(rep(1, 4), "binary", fn = "check_edp"),
    class = "positively_exposure_type_error"
  )
})

test_that("validate_exposure_structure() adds no requirement for categorical", {
  expect_no_error(
    validate_exposure_structure(
      c("a", "b", "c"),
      "categorical",
      fn = "check_edp"
    )
  )
  expect_no_error(
    validate_exposure_structure(
      factor(c("a", "b", "c")),
      "categorical",
      fn = "check_edp"
    )
  )
  # A numeric column the heuristic would read as continuous is still a legal
  # categorical exposure once the user declares it.
  expect_no_error(
    validate_exposure_structure(
      seq(0, 1, length.out = 50),
      "categorical",
      fn = "check_edp"
    )
  )
})

test_that("resolve_exposure_type() errors point at the calling function", {
  withr::local_options(positively.quiet = TRUE)
  check_fake <- function(.exposure, exposure_type = "auto") {
    resolve_exposure_type(
      exposure_type,
      .exposure,
      supported = "continuous",
      fn = "check_fake"
    )
  }
  expect_snapshot(error = TRUE, check_fake(c(0, 1, 0, 1)))
  expect_snapshot(
    error = TRUE,
    check_fake(c(0, 1, 0, 1), exposure_type = "binary")
  )
  expect_snapshot(
    error = TRUE,
    check_fake(factor(c("a", "b", "c")), exposure_type = "continuous")
  )
})
