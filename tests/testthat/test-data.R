# The lazy-loaded datasets ship known positivity violations so examples and
# tests can exercise the diagnostics against ground truth. pos_violations is a
# binary-exposure point dataset carrying one structural subgroup (region "b"
# with x2 above 1 is never exposed) and one practical near-violation (the x1
# tails push the fitted propensity below 0.01 and above 0.99 while both exposure
# levels remain observed). pos_violations_long is a wide longitudinal dataset,
# one row per subject over three time points, whose time-2 exposure a2 has a
# hard support gap on the interval (2, 4) that a1 and a3 populate.

local_quiet <- function(.env = parent.frame()) {
  withr::local_options(positively.quiet = TRUE, .local_envir = .env)
}

# The propensity model documented for the practical near-violation. Fitting it
# recovers scores that cross the planted extremes.
pos_violations_ps <- function(data) {
  model <- stats::glm(
    exposure ~ x1 + x2 + region,
    family = stats::binomial(),
    data = data
  )
  stats::predict(model, type = "response")
}

# The structural subgroup: region "b" and x2 above 1.
pos_violations_structural <- function(data) {
  data$region == "b" & data$x2 > 1
}

# ---- pos_violations: schema ------------------------------------------------

test_that("pos_violations is a data frame with the expected dimensions", {
  expect_s3_class(pos_violations, "data.frame")
  expect_identical(nrow(pos_violations), 1000L)
  expect_identical(ncol(pos_violations), 4L)
})

test_that("pos_violations carries the expected columns and types", {
  expect_identical(names(pos_violations), c("exposure", "x1", "x2", "region"))
  expect_type(pos_violations$exposure, "integer")
  expect_true(all(pos_violations$exposure %in% c(0L, 1L)))
  expect_type(pos_violations$x1, "double")
  expect_type(pos_violations$x2, "double")
  expect_s3_class(pos_violations$region, "factor")
  expect_identical(levels(pos_violations$region), c("a", "b"))
})

test_that("pos_violations has no missing values", {
  expect_false(anyNA(pos_violations))
})

# ---- pos_violations: structural violation ----------------------------------

test_that("the structural subgroup is present and non-empty", {
  structural <- pos_violations_structural(pos_violations)
  expect_gt(sum(structural), 0L)
})

test_that("exactly one exposure level occurs inside the structural subgroup", {
  structural <- pos_violations_structural(pos_violations)
  inside <- pos_violations$exposure[structural]

  expect_length(unique(inside), 1L)
  expect_true(all(inside == 0L))
})

test_that("both exposure levels occur outside the structural subgroup", {
  structural <- pos_violations_structural(pos_violations)
  outside <- pos_violations$exposure[!structural]

  expect_setequal(unique(outside), c(0L, 1L))
})

# ---- pos_violations: practical near-violation ------------------------------

test_that("the documented propensity model crosses the planted extremes", {
  fitted <- pos_violations_ps(pos_violations)

  expect_lt(min(fitted), 0.01)
  expect_gt(max(fitted), 0.99)
})

test_that("both exposure levels remain observed despite the near-violation", {
  # A practical violation differs from a structural one in that both levels are
  # observed even where the fitted score is near-deterministic.
  expect_setequal(unique(pos_violations$exposure), c(0L, 1L))
})

# ---- pos_violations: diagnostic smoke --------------------------------------

test_that("pos_violations runs through a binary diagnostic without error", {
  local_quiet()
  expect_no_error(
    check_edp(pos_violations, exposure, c(x1, x2, region), values = 1)
  )
})

# ---- pos_violations_long: schema -------------------------------------------

test_that("pos_violations_long is a data frame with the expected dimensions", {
  expect_s3_class(pos_violations_long, "data.frame")
  expect_identical(nrow(pos_violations_long), 500L)
  expect_identical(ncol(pos_violations_long), 7L)
})

test_that("pos_violations_long carries the expected columns and types", {
  expect_identical(
    names(pos_violations_long),
    c("id", "l0", "a1", "l1", "a2", "l2", "a3")
  )
  expect_type(pos_violations_long$id, "integer")
  expect_type(pos_violations_long$l0, "double")
  expect_type(pos_violations_long$a1, "double")
  expect_type(pos_violations_long$l1, "double")
  expect_type(pos_violations_long$a2, "double")
  expect_type(pos_violations_long$l2, "double")
  expect_type(pos_violations_long$a3, "double")
})

test_that("pos_violations_long identifies 500 distinct subjects", {
  expect_setequal(pos_violations_long$id, seq_len(500L))
})

test_that("pos_violations_long has no missing values", {
  expect_false(anyNA(pos_violations_long))
})

test_that("each time point's exposure is continuous rather than a few levels", {
  expect_gt(length(unique(pos_violations_long$a1)), 100L)
  expect_gt(length(unique(pos_violations_long$a2)), 100L)
  expect_gt(length(unique(pos_violations_long$a3)), 100L)
})

# ---- pos_violations_long: support gap at time 2 ----------------------------

test_that("the time-2 exposure leaves the interval (2, 4) unpopulated", {
  expect_false(any(pos_violations_long$a2 > 2 & pos_violations_long$a2 < 4))
})

test_that("the time-1 and time-3 exposures populate the interval (2, 4)", {
  expect_true(any(pos_violations_long$a1 > 2 & pos_violations_long$a1 < 4))
  expect_true(any(pos_violations_long$a3 > 2 & pos_violations_long$a3 < 4))
})

# ---- pos_violations_long: diagnostic smoke ---------------------------------

test_that("the longitudinal data runs through the continuous HDR diagnostic", {
  local_quiet()
  expect_no_error(check_hdr(pos_violations_long, a2, c(l0, l1)))
})

test_that("the longitudinal data runs through the sequential HDR diagnostic", {
  local_quiet()
  expect_no_error(
    check_hdr_seq(
      pos_violations_long,
      c(a1, a2, a3),
      list(l0, c(l0, l1), c(l0, l1, l2))
    )
  )
})

test_that("the longitudinal data runs through the sequential PoRT diagnostic", {
  local_quiet()
  expect_no_error(
    check_port_seq(
      pos_violations_long,
      c(a1, a2, a3),
      list(l0, l1, l2),
      .baseline = l0
    )
  )
})
