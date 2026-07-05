# The data-generating processes are test infrastructure, so these specs assert
# their documented properties and pass without any package implementation.

test_that("dgp_good_positivity() keeps the propensity in [0.2, 0.8]", {
  data <- dgp_good_positivity(n = 400, seed = 1)
  expect_s3_class(data, "tbl_df")
  expect_identical(nrow(data), 400L)
  expect_true(all(data$ps >= 0.2 & data$ps <= 0.8))
  expect_setequal(unique(data$exposure), c(0L, 1L))
})

test_that("dgp_practical_violation() produces near-deterministic tails", {
  data <- dgp_practical_violation(n = 500, seed = 2)
  expect_identical(nrow(data), 500L)
  expect_lt(min(data$ps), 0.05)
  expect_gt(max(data$ps), 0.95)
  expect_setequal(unique(data$exposure), c(0L, 1L))
})

test_that("dgp_structural_subgroup() plants an untreated subgroup", {
  data <- dgp_structural_subgroup(n = 1000, seed = 3)
  in_gap <- data$x1 > 2 & data$x2 == "b"
  expect_gt(sum(in_gap), 0)
  expect_true(all(data$exposure[in_gap] == 0L))
  # Treatment still occurs outside the planted rule.
  expect_gt(sum(data$exposure[!in_gap]), 0)
})

test_that("dgp_continuous_support_gap() leaves the interval (2, 4) empty", {
  data <- dgp_continuous_support_gap(n = 500, seed = 4)
  expect_identical(nrow(data), 500L)
  expect_identical(sum(data$exposure > 2 & data$exposure < 4), 0L)
  expect_true(any(data$exposure <= 2))
  expect_true(any(data$exposure >= 4))
  expect_gt(length(unique(data$exposure)), 2)
})

test_that("dgp_longitudinal() is wide with a time-2 support gap", {
  data <- dgp_longitudinal(n = 300, seed = 5)
  expect_identical(nrow(data), 300L)
  expect_true(all(c("a1", "a2", "a3") %in% names(data)))
  expect_true(all(c("l0", "l1", "l2") %in% names(data)))
  expect_identical(sum(data$a2 > 2 & data$a2 < 4), 0L)
})

test_that("each data-generating process is deterministic given its seed", {
  expect_identical(dgp_good_positivity(seed = 7), dgp_good_positivity(seed = 7))
  expect_identical(
    dgp_practical_violation(seed = 7),
    dgp_practical_violation(seed = 7)
  )
  expect_identical(
    dgp_structural_subgroup(seed = 7),
    dgp_structural_subgroup(seed = 7)
  )
  expect_identical(
    dgp_continuous_support_gap(seed = 7),
    dgp_continuous_support_gap(seed = 7)
  )
  expect_identical(dgp_longitudinal(seed = 7), dgp_longitudinal(seed = 7))
})
