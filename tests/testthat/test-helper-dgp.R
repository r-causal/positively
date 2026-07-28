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

test_that("dgp_coarse_dose() dispenses eight levels detection misreads", {
  data <- dgp_coarse_dose(n = 150, seed = 1)
  expect_identical(nrow(data), 150L)
  # The grid literal is pinned here and nowhere else, so coarse_dose_grid() has
  # exactly one place that fixes its values.
  expect_setequal(unique(data$exposure), seq(2.5, 20, by = 2.5))
  # The fixture exists to be misread, so the heuristic itself is the assertion.
  # A hard-coded ratio would keep this test green after a cutoff change that had
  # already stopped detection from misreading the dose.
  expect_true(is_categorical(data$exposure))
  # The covariate tracks the dose, so a diagnostic run on the declared type has
  # structure to find.
  expect_gt(stats::cor(data$exposure, data$x1), 0.5)
})

test_that("dgp_longitudinal() is wide with a time-2 support gap", {
  data <- dgp_longitudinal(n = 300, seed = 5)
  expect_identical(nrow(data), 300L)
  expect_true(all(c("a1", "a2", "a3") %in% names(data)))
  expect_true(all(c("l0", "l1", "l2") %in% names(data)))
  expect_identical(sum(data$a2 > 2 & data$a2 < 4), 0L)
})

test_that("dgp_longitudinal_binary() plants two never-initiate regions at t = 2", {
  data <- dgp_longitudinal_binary(n = 2000, seed = 6)
  expect_true(all(c("a1", "a2", "a3") %in% names(data)))
  expect_true(all(c("l0", "l1", "l2") %in% names(data)))
  # Monotone treatment: once treated, stays treated.
  expect_true(all(data$a2 >= data$a1))
  expect_true(all(data$a3 >= data$a2))
  # Among the followers at t = 2, treatment never initiates where the lagged
  # covariate l0 is above 1 or where the current covariate l1 is above 1.
  followers2 <- data$a1 == 0
  region_l0 <- followers2 & data$l0 > 1
  region_l1 <- followers2 & data$l1 > 1
  expect_gt(sum(region_l0), 0)
  expect_gt(sum(region_l1), 0)
  expect_true(all(data$a2[region_l0] == 0L))
  expect_true(all(data$a2[region_l1] == 0L))
  # Treatment still initiates among followers outside the planted regions.
  free <- followers2 & data$l0 <= 1 & data$l1 <= 1
  expect_gt(sum(data$a2[free]), 0)
  # The follower risk set shrinks across time.
  expect_gt(sum(data$a1 == 0), sum(data$a2 == 0))
})

test_that("dgp_longitudinal_binary_censoring() plants near-certain t-2 censoring", {
  data <- dgp_longitudinal_binary_censoring(n = 2000, seed = 7)
  expect_true(all(c("c1", "c2", "c3") %in% names(data)))
  # The censoring indicators are binary 0/1.
  for (nm in c("c1", "c2", "c3")) {
    expect_setequal(unique(data[[nm]]), c(0L, 1L))
  }
  # Among subjects uncensored through t = 1, censoring at t = 2 is near-certain
  # where l0 > 1 and light otherwise.
  risk2 <- data$c1 == 0
  region <- risk2 & data$l0 > 1
  expect_gt(sum(region), 0)
  expect_gt(mean(data$c2[region]), 0.9)
  expect_lt(mean(data$c2[risk2 & data$l0 <= 1]), 0.2)
  # The uncensored-so-far risk set shrinks across time.
  risk_sets <- c(
    nrow(data),
    sum(data$c1 == 0),
    sum(data$c1 == 0 & data$c2 == 0)
  )
  expect_true(all(diff(risk_sets) < 0))
  # A subject censored at an earlier wave carries 0 at every later wave, having
  # left the risk set.
  expect_true(all(data$c2[data$c1 == 1] == 0))
  expect_true(all(data$c3[data$c2 == 1] == 0))
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
  expect_identical(dgp_coarse_dose(seed = 7), dgp_coarse_dose(seed = 7))
  expect_identical(dgp_longitudinal(seed = 7), dgp_longitudinal(seed = 7))
  expect_identical(
    dgp_longitudinal_binary(seed = 7),
    dgp_longitudinal_binary(seed = 7)
  )
  expect_identical(
    dgp_longitudinal_binary_censoring(seed = 7),
    dgp_longitudinal_binary_censoring(seed = 7)
  )
})
