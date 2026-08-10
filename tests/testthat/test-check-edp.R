# EDP ("effective data points") sums a product kernel over the sample for each
# intervened-on observation o_int = (l_i, a*): covariates held at their own
# values, exposure set to the intervention value. Each dimension contributes a
# factor in [0, 1]:
#   continuous dimension, half-distance h : 0.5 ^ ((delta / h) ^ 2)
#                                           (1 at delta 0, exactly 0.5 at |delta| = h)
#   categorical dimension                 : 1 on a match, else categorical_similarity
# so 0 <= edp <= n. Binary and categorical exposure dimensions use the match /
# similarity kernel; only continuous exposures use the half-distance kernel.
# Numeric covariates are continuous dimensions, factor and character covariates
# are categorical dimensions. The fixed data-variant results columns are
# .id / value / edp; the estimator variant replaces edp with edp_outcome /
# edp_treatment / ideal_weight. EDP is meaningful relative across observations
# for a fixed covariate set, never against a universal threshold, so magnitude
# claims are anchored to hand computations or to stratum-versus-stratum ratios.

# ---- Scenario generators --------------------------------------------------

# Fixed-density Gaussian exposure with an independent Gaussian covariate. a* = 0
# is well supported; a* = 4 sits four standard deviations from all support.
sim_edp_gaussian <- function(n, seed = 1) {
  withr::local_seed(seed)
  tibble::tibble(exposure = stats::rnorm(n), x1 = stats::rnorm(n))
}

# The "divided" regime: the exposure clusters near 0.2 for l = 0 and near 0.8
# for l = 1, so a* = 0.2 is supported within the l = 0 stratum and unsupported
# within the l = 1 stratum.
sim_edp_divided <- function(n, seed = 1) {
  withr::local_seed(seed)
  l <- stats::rbinom(n, 1L, 0.3)
  a <- ifelse(
    l == 1L,
    stats::rnorm(n, 0.8, 0.1),
    stats::rnorm(n, 0.2, 0.1)
  )
  tibble::tibble(exposure = a, l = l)
}

# Binary exposure with a structural boundary on a continuous covariate: nobody
# with x1 > 1.5 is ever treated, so the deep region x1 > 2.5 has little support
# for a* = 1.
sim_edp_structural <- function(n, seed = 1) {
  withr::local_seed(seed)
  x1 <- stats::rnorm(n)
  ps <- stats::plogis(0.5 * x1)
  ps[x1 > 1.5] <- 0
  exposure <- stats::rbinom(n, 1L, ps)
  tibble::tibble(exposure = exposure, x1 = x1)
}

# Binary exposure with a categorical covariate subgroup (s == "s1") that is
# never treated. At a* = 1 this subgroup has zero support once non-matching
# categories contribute zero.
sim_edp_masking <- function(n, seed = 1) {
  withr::local_seed(seed)
  s <- factor(sample(c("s0", "s1"), n, replace = TRUE))
  ps <- ifelse(s == "s1", 0, 0.5)
  exposure <- stats::rbinom(n, 1L, ps)
  tibble::tibble(exposure = exposure, s = s)
}

# Three-level categorical exposure where level 2 is absent when z2 > 1.
sim_edp_categorical <- function(n, seed = 1) {
  withr::local_seed(seed)
  z2 <- stats::rnorm(n)
  base <- sample(0:2, n, replace = TRUE)
  base[z2 > 1 & base == 2L] <- 0L
  tibble::tibble(exposure = factor(base), z2 = z2)
}

# Continuous exposure that tracks a covariate: g(a* | l) shrinks as l moves away
# from a*, which drives ideal_weight up where support is thin.
sim_edp_estimator <- function(n, seed = 1) {
  withr::local_seed(seed)
  l <- stats::rnorm(n)
  a <- stats::rnorm(n, mean = l, sd = 0.5)
  tibble::tibble(exposure = a, l = l)
}

# A g-model-only covariate g with an isolated subgroup eight standard deviations
# out. Placed in the treatment set only, it starves edp_treatment there while
# edp_outcome (which excludes g) stays high.
sim_edp_asymmetric <- function(n, seed = 1) {
  withr::local_seed(seed)
  l <- stats::rnorm(n)
  a <- stats::rnorm(n, mean = l, sd = 0.5)
  g <- stats::rnorm(n)
  idx <- sample(n, 20L)
  g[idx] <- g[idx] + 8
  tibble::tibble(
    exposure = a,
    l = l,
    g = g,
    subgroup = seq_len(n) %in% idx
  )
}

# ---- Argument validation --------------------------------------------------

test_that("check_edp() rejects non-data-frame input", {
  expect_error(
    check_edp(1:10, exposure, x1),
    class = "positively_error"
  )
})

test_that("check_edp() rejects an empty covariate selection", {
  data <- sim_edp_gaussian(60)
  expect_error(
    check_edp(
      data,
      exposure,
      tidyselect::starts_with("zzz"),
      exposure_type = "continuous"
    ),
    class = "positively_error"
  )
})

test_that("check_edp() requires a single exposure column", {
  data <- sim_edp_gaussian(60)
  expect_error(
    check_edp(data, c(exposure, x1), x1, exposure_type = "continuous"),
    class = "positively_error"
  )
})

test_that("check_edp() rejects a categorical_similarity outside the unit interval", {
  data <- sim_edp_gaussian(60)
  expect_error(
    check_edp(
      data,
      exposure,
      x1,
      categorical_similarity = -0.1,
      exposure_type = "continuous"
    ),
    class = "positively_error"
  )
  expect_error(
    check_edp(
      data,
      exposure,
      x1,
      categorical_similarity = 1.5,
      exposure_type = "continuous"
    ),
    class = "positively_error"
  )
})

test_that("check_edp() rejects a non-numeric or non-scalar categorical_similarity", {
  data <- sim_edp_gaussian(60)
  expect_error(
    check_edp(
      data,
      exposure,
      x1,
      categorical_similarity = "x",
      exposure_type = "continuous"
    ),
    class = "positively_error"
  )
  expect_error(
    check_edp(
      data,
      exposure,
      x1,
      categorical_similarity = c(0, 1),
      exposure_type = "continuous"
    ),
    class = "positively_error"
  )
})

test_that("check_edp() rejects a negative, non-numeric, or non-scalar bw_exposure", {
  data <- sim_edp_gaussian(60)
  expect_error(
    check_edp(
      data,
      exposure,
      x1,
      bw_exposure = -1,
      exposure_type = "continuous"
    ),
    class = "positively_error"
  )
  expect_error(
    check_edp(
      data,
      exposure,
      x1,
      bw_exposure = "1",
      exposure_type = "continuous"
    ),
    class = "positively_error"
  )
  expect_error(
    check_edp(
      data,
      exposure,
      x1,
      bw_exposure = c(1, 2),
      exposure_type = "continuous"
    ),
    class = "positively_error"
  )
})

test_that("check_edp() rejects a negative bw_covariates", {
  data <- sim_edp_gaussian(60)
  expect_error(
    check_edp(
      data,
      exposure,
      x1,
      bw_covariates = -1,
      exposure_type = "continuous"
    ),
    class = "positively_error"
  )
})

test_that("check_edp() aborts on missing exposure or covariate values", {
  data <- sim_edp_gaussian(60)

  na_exposure <- data
  na_exposure$exposure[1] <- NA
  expect_error(
    check_edp(na_exposure, exposure, x1, exposure_type = "continuous"),
    class = "positively_error"
  )

  na_covariate <- data
  na_covariate$x1[1] <- NA
  expect_error(
    check_edp(na_covariate, exposure, x1, exposure_type = "continuous"),
    class = "positively_error"
  )
})

test_that("check_edp() aborts on fewer than two observations", {
  data <- tibble::tibble(exposure = 0.5, x1 = -0.2)
  expect_error(
    check_edp(data, exposure, x1, exposure_type = "continuous"),
    class = "positively_error"
  )
})

test_that("check_edp() rejects an unknown variant or kernel", {
  data <- sim_edp_gaussian(60)
  # An argument the user chose from a fixed menu fails like every other failure
  # the package raises, so it carries the package's condition class rather than
  # arriving as the bare rlang error `arg_match()` throws.
  expect_error(
    check_edp(
      data,
      exposure,
      x1,
      variant = "bogus",
      exposure_type = "continuous"
    ),
    class = "positively_args_error"
  )
  expect_error(
    check_edp(
      data,
      exposure,
      x1,
      kernel = "bogus",
      exposure_type = "continuous"
    ),
    class = "positively_args_error"
  )
})

test_that("a near miss on a menu argument keeps its suggestion", {
  data <- sim_edp_gaussian(60)
  # Classing the failure must not cost the user the message `arg_match()` wrote:
  # it names every allowed value and, on a near miss, the one that was meant.
  err <- expect_error(
    check_edp(
      data,
      exposure,
      x1,
      variant = "dat",
      exposure_type = "continuous"
    ),
    class = "positively_args_error"
  )
  expect_match(
    conditionMessage(err),
    'Did you mean "data"?',
    fixed = TRUE
  )

  # An argument value is user text, and it reaches the message as a value rather
  # than as template text, so a value cli would otherwise read as an expression
  # stays a value and the package's condition classes survive.
  braced <- expect_error(
    check_edp(
      data,
      exposure,
      x1,
      variant = "dat{a}",
      exposure_type = "continuous"
    ),
    class = "positively_args_error"
  )
  expect_match(conditionMessage(braced), 'not "dat{a}"', fixed = TRUE)

  # A failure raised while evaluating the argument is not a failure to match it,
  # so it must not be relabelled as one on its way out.
  failed <- expect_error(check_edp(data, exposure, x1, variant = no_such_thing))
  expect_false(inherits(failed, "positively_error"))
  expect_match(conditionMessage(failed), "no_such_thing", fixed = TRUE)
})

test_that("declaring continuous on a character exposure aborts before any coercion", {
  local_quiet()
  # The continuous path reaches as.double(exposure_vec) when it sizes the
  # default bandwidth. On a character column that returns all NAs with a
  # coercion warning and the run then dies on a missing condition, so the type
  # has to be checked against the column before any arithmetic runs.
  data <- sim_edp_gaussian(60)
  data$exposure <- rep(c("a", "b", "c"), length.out = nrow(data))

  expect_error(
    check_edp(data, exposure, x1, exposure_type = "continuous"),
    class = "positively_exposure_type_error"
  )
})

test_that("declaring binary on a three-level exposure aborts rather than mislabeling", {
  local_quiet()
  # Nothing downstream of the type resolution narrows a three-level exposure to
  # two, so a binary declaration used to compute a categorical edp and stamp
  # "binary" on the result, which then propagated to @exposure_type, glance(),
  # and print().
  data <- sim_edp_categorical(300)

  expect_error(
    check_edp(data, exposure, z2, exposure_type = "binary"),
    class = "positively_exposure_type_error"
  )
})

test_that("check_edp() announces the type it detected", {
  withr::local_options(positively.quiet = FALSE)
  # check_edp() computes a different quantity for each exposure type it
  # supports, so which one it read decides what the result means and has to be
  # said. A declared type is taken as given, so there is nothing to announce.
  data <- sim_edp_gaussian(60)

  expect_message(
    check_edp(data, exposure, x1, values = 0),
    "Treating `.exposure` as continuous",
    fixed = TRUE
  )
  expect_no_message(
    check_edp(data, exposure, x1, values = 0, exposure_type = "continuous")
  )
})

# ---- Result class and structure -------------------------------------------

test_that("check_edp() returns an edp_result diagnostic", {
  data <- sim_edp_gaussian(150)
  res <- check_edp(
    data,
    exposure,
    x1,
    values = c(0, 1),
    exposure_type = "continuous"
  )

  expect_true(S7::S7_inherits(res, positivity_diagnostic))
  expect_identical(S7::S7_class(res)@name, "edp_result")
  expect_identical(res@exposure_type, "continuous")
  expect_identical(res@n, 150L)
})

test_that("the data variant carries the fixed results columns and variant property", {
  data <- sim_edp_gaussian(150)
  res <- check_edp(
    data,
    exposure,
    x1,
    values = c(0, 1),
    exposure_type = "continuous"
  )

  expect_s3_class(res@results, "tbl_df")
  expect_setequal(names(res@results), c(".id", "value", "edp"))
  expect_type(res@results$edp, "double")
  expect_identical(nrow(res@results), 2L * 150L)
  expect_identical(res@variant, "data")
  expect_type(res@bandwidths, "list")
})

test_that("nrow scales as n times the number of intervention values", {
  data <- sim_edp_gaussian(120)
  res <- check_edp(
    data,
    exposure,
    x1,
    values = c(-1, 0, 1),
    exposure_type = "continuous"
  )
  expect_identical(nrow(res@results), 3L * 120L)
  expect_setequal(unique(res@results$value), c(-1, 0, 1))
})

test_that("the continuous default grid spans several intervention values", {
  data <- sim_edp_gaussian(100)
  res <- check_edp(data, exposure, x1, exposure_type = "continuous")

  expect_gt(length(unique(res@results$value)), 1)
  expect_identical(
    nrow(res@results),
    length(unique(res@results$value)) * 100L
  )
})

test_that("check_edp() accepts a binary exposure with all levels by default", {
  local_quiet()
  data <- dgp_good_positivity(n = 200, seed = 1)
  res <- check_edp(data, exposure, c(x1, x2))

  expect_identical(res@exposure_type, "binary")
  expect_identical(nrow(res@results), 2L * 200L)
})

test_that("check_edp() accepts a categorical exposure with all levels by default", {
  data <- sim_edp_categorical(300)
  res <- check_edp(data, exposure, z2, exposure_type = "categorical")

  expect_identical(res@exposure_type, "categorical")
  expect_identical(nrow(res@results), 3L * 300L)
})

test_that("the estimator variant swaps in the outcome and treatment columns", {
  data <- sim_edp_gaussian(150)
  res <- check_edp(
    data,
    exposure,
    x1,
    variant = "estimator",
    values = 0,
    exposure_type = "continuous"
  )

  expect_setequal(
    names(res@results),
    c(".id", "value", "edp_outcome", "edp_treatment", "ideal_weight")
  )
  expect_false("edp" %in% names(res@results))
  expect_identical(res@variant, "estimator")
})

test_that("tidy() returns the results tibble", {
  data <- sim_edp_gaussian(100)
  res <- check_edp(data, exposure, x1, values = 0, exposure_type = "continuous")

  expect_identical(generics::tidy(res), res@results)
})

# Three observations at exposures 0, 1, and 2 on a constant covariate. With both
# half-distances at 1 the covariate factor is 1 for every pair, so EDP is the
# exposure kernel alone and is identical for every row: at a* = 0 it is
# 1 + 0.5 + 0.5^4 = 1.5625 and at a* = 1 it is 0.5 + 1 + 0.5 = 2.
edp_hand_data <- function() {
  tibble::tibble(exposure = c(0, 1, 2), x1 = c(0, 0, 0))
}

test_that("glance() reports the data variant's EDP range", {
  res <- check_edp(
    edp_hand_data(),
    exposure,
    x1,
    bw_exposure = 1,
    bw_covariates = 1,
    values = c(0, 1),
    exposure_type = "continuous"
  )
  glanced <- generics::glance(res)

  expect_s3_class(glanced, "tbl_df")
  expect_identical(nrow(glanced), 1L)
  expect_setequal(
    names(glanced),
    c("n", "variant", "n_values", "edp_min", "edp_max")
  )

  expect_identical(glanced$n, 3L)
  expect_identical(glanced$variant, "data")
  expect_identical(glanced$n_values, 2L)
  expect_equal(glanced$edp_min, 1.5625)
  expect_equal(glanced$edp_max, 2)
})

test_that("glance() reports the estimator variant's three measures", {
  # edp_outcome conditions on the exposure and the covariate, so it repeats the
  # data-variant values. edp_treatment conditions on the constant covariate
  # alone and so counts all three observations at every point. ideal_weight is
  # their ratio: 3 / 1.5625 = 1.92 at a* = 0 and 3 / 2 = 1.5 at a* = 1.
  res <- check_edp(
    edp_hand_data(),
    exposure,
    x1,
    variant = "estimator",
    bw_exposure = 1,
    bw_covariates = 1,
    values = c(0, 1),
    exposure_type = "continuous"
  )
  glanced <- generics::glance(res)

  expect_identical(nrow(glanced), 1L)
  expect_setequal(
    names(glanced),
    c(
      "n",
      "variant",
      "n_values",
      "edp_outcome_min",
      "edp_outcome_max",
      "edp_treatment_min",
      "edp_treatment_max",
      "ideal_weight_min",
      "ideal_weight_max"
    )
  )

  # The estimator variant has no `edp` measure at all, so the data variant's
  # columns are absent rather than missing.
  expect_false("edp_min" %in% names(glanced))
  expect_false("edp_max" %in% names(glanced))

  expect_identical(glanced$variant, "estimator")
  expect_identical(glanced$n_values, 2L)
  expect_equal(glanced$edp_outcome_min, 1.5625)
  expect_equal(glanced$edp_outcome_max, 2)
  expect_equal(glanced$edp_treatment_min, 3)
  expect_equal(glanced$edp_treatment_max, 3)
  expect_equal(glanced$ideal_weight_min, 1.5)
  expect_equal(glanced$ideal_weight_max, 1.92)
})

# ---- Exact kernel identities ----------------------------------------------

test_that("the half-distance kernel gives k(0) = 1 and k(h) = 0.5", {
  # A constant covariate contributes a factor of 1 to every pair, isolating the
  # exposure kernel. obs a = c(0, 1); h = 1; a* = 0 gives edp = 1 + 0.5.
  data <- tibble::tibble(exposure = c(0, 1), x1 = c(0, 0))
  res <- check_edp(
    data,
    exposure,
    x1,
    bw_exposure = 1,
    bw_covariates = 1,
    values = 0,
    exposure_type = "continuous"
  )

  expect_equal(unique(res@results$edp), 1.5, tolerance = 1e-8)
})

test_that("edp matches the hand-computed single-dimension sum", {
  # obs a = c(0, 1, 2, 5); h = 1; a* = 0. Weights 1, 0.5, 0.0625, ~3e-8 sum to
  # 1.5625. The constant covariate contributes 1 throughout.
  data <- tibble::tibble(exposure = c(0, 1, 2, 5), x1 = c(0, 0, 0, 0))
  res <- check_edp(
    data,
    exposure,
    x1,
    bw_exposure = 1,
    bw_covariates = 1,
    values = 0,
    exposure_type = "continuous"
  )

  expect_equal(unique(res@results$edp), 1.5625, tolerance = 1e-6)
})

test_that("edp matches the hand-computed two-dimension product", {
  # covariate x1 with h = 1, exposure with h = 0.5. obs (x1, a) =
  # (0, 0), (0, 1), (1, 0); a* = 0.5. The two x1 = 0 rows give edp = 1.25, the
  # x1 = 1 row gives 1.0.
  data <- tibble::tibble(exposure = c(0, 1, 0), x1 = c(0, 0, 1))
  res <- check_edp(
    data,
    exposure,
    x1,
    bw_exposure = 0.5,
    bw_covariates = 1,
    values = 0.5,
    exposure_type = "continuous"
  )

  edp <- res@results$edp[order(res@results$.id)]
  expect_equal(edp[1], 1.25, tolerance = 1e-8)
  expect_equal(edp[2], 1.25, tolerance = 1e-8)
  expect_equal(edp[3], 1.0, tolerance = 1e-8)
})

test_that("edp is bounded in [0, n] for arbitrary input", {
  data <- sim_edp_gaussian(80, seed = 3)
  res <- check_edp(
    data,
    exposure,
    x1,
    values = c(0, 4),
    exposure_type = "continuous"
  )

  expect_true(all(res@results$edp >= 0))
  expect_true(all(res@results$edp <= res@n))
})

# ---- Default rule-of-thumb bandwidths -------------------------------------

test_that("default bandwidths follow the paper rule of thumb", {
  # The paper's rule of thumb sets the exposure half-distance to 0.5 * sd of the
  # exposure and each numeric covariate half-distance to 1 * sd of that
  # covariate.
  data <- sim_edp_gaussian(120, seed = 1)
  res <- check_edp(data, exposure, x1, values = 0, exposure_type = "continuous")

  expect_equal(res@bandwidths$exposure, 0.5 * stats::sd(data$exposure))
  expect_equal(
    res@bandwidths$covariates,
    c(x1 = stats::sd(data$x1))
  )
})

# ---- Bandwidth limits and monotonicity ------------------------------------

test_that("a very large bandwidth drives edp to n", {
  data <- sim_edp_gaussian(100, seed = 2)
  res <- check_edp(
    data,
    exposure,
    x1,
    bw_exposure = 1e6,
    bw_covariates = 1e6,
    values = 0,
    exposure_type = "continuous"
  )

  expect_equal(mean(res@results$edp), res@n, tolerance = 0.5)
})

test_that("a zero bandwidth counts exact matches", {
  # a* = 0 with a constant matching covariate counts observed exposures equal to
  # zero: one for the distinct grid, two once a duplicate is added.
  distinct <- tibble::tibble(exposure = c(0, 1, 2, 5), x1 = c(0, 0, 0, 0))
  res_distinct <- check_edp(
    distinct,
    exposure,
    x1,
    bw_exposure = 0,
    bw_covariates = 1,
    values = 0,
    exposure_type = "continuous"
  )
  expect_true(all(res_distinct@results$edp == 1))

  duplicated <- tibble::tibble(exposure = c(0, 0, 1, 2), x1 = c(0, 0, 0, 0))
  res_duplicated <- check_edp(
    duplicated,
    exposure,
    x1,
    bw_exposure = 0,
    bw_covariates = 1,
    values = 0,
    exposure_type = "continuous"
  )
  expect_true(all(res_duplicated@results$edp == 2))
})

test_that("check_edp() rejects non-finite exposure or covariate values", {
  inf_exposure <- data.frame(a = c(1, 2, Inf, 4), x1 = 0)
  expect_error(
    check_edp(inf_exposure, a, x1, values = 1, exposure_type = "continuous"),
    class = "positively_error"
  )

  inf_covariate <- data.frame(a = c(1, 2, 3, 4), x1 = c(0, -Inf, 0, 0))
  expect_error(
    check_edp(
      inf_covariate,
      a,
      x1,
      bw_exposure = 1,
      bw_covariates = 1,
      values = 1,
      exposure_type = "continuous"
    ),
    class = "positively_range_error"
  )
})

test_that("mean edp is nondecreasing as the exposure bandwidth grows", {
  data <- sim_edp_gaussian(200, seed = 1)
  bws <- c(0.1, 0.25, 0.5, 1, 2, 4)
  edp_means <- vapply(
    bws,
    function(h) {
      res <- check_edp(
        data,
        exposure,
        x1,
        bw_exposure = h,
        bw_covariates = 1,
        values = 0,
        exposure_type = "continuous"
      )
      mean(res@results$edp)
    },
    numeric(1)
  )

  expect_true(all(diff(edp_means) >= -1e-9))
})

test_that("a steeper kernel sharpens a far-off violation toward zero", {
  data <- sim_edp_gaussian(200, seed = 1)
  bws <- c(2, 1, 0.5, 0.25)
  edp_medians <- vapply(
    bws,
    function(h) {
      res <- check_edp(
        data,
        exposure,
        x1,
        bw_exposure = h,
        bw_covariates = 1,
        values = 4,
        exposure_type = "continuous"
      )
      stats::median(res@results$edp)
    },
    numeric(1)
  )

  expect_true(all(diff(edp_medians) < 0))
  expect_lt(edp_medians[length(edp_medians)], 0.5)
})

# ---- Behavioral separation -----------------------------------------------

test_that("good binary overlap keeps every subject away from zero support", {
  local_quiet()
  data <- dgp_good_positivity(n = 1000, seed = 1)
  res <- check_edp(data, exposure, c(x1, x2), values = 1)

  expect_gt(min(res@results$edp), 1)
  expect_gt(stats::quantile(res@results$edp, 0.05), 5)
})

test_that("a stratified continuous violation collapses the unsupported stratum", {
  data <- sim_edp_divided(350, seed = 1)
  res <- check_edp(
    data,
    exposure,
    l,
    values = 0.2,
    exposure_type = "continuous"
  )

  l_by_id <- data$l[res@results$.id]
  supported <- stats::median(res@results$edp[l_by_id == 0])
  unsupported <- stats::median(res@results$edp[l_by_id == 1])

  expect_lt(unsupported, 0.10 * supported)
})

test_that("a continuous structural boundary starves the deep region", {
  local_quiet()
  data <- sim_edp_structural(1000, seed = 1)
  res <- check_edp(data, exposure, x1, values = 1)

  x1_by_id <- data$x1[res@results$.id]
  deep <- res@results$edp[x1_by_id > 2.5]
  bulk <- res@results$edp[abs(x1_by_id) < 0.5]

  # A continuous boundary degrades support gradually, so the deep region sits
  # well below the bulk median rather than exactly at zero.
  expect_lt(max(deep), 0.3 * stats::median(bulk))
})

test_that("a never-treated categorical subgroup reads as zero at c = 0", {
  local_quiet()
  data <- sim_edp_masking(600, seed = 1)
  res <- check_edp(data, exposure, s, values = 1, categorical_similarity = 0)

  subgroup <- data$s[res@results$.id] == "s1"
  expect_equal(stats::median(res@results$edp[subgroup]), 0)
})

test_that("a positive similarity masks the categorical violation", {
  local_quiet()
  data <- sim_edp_masking(600, seed = 1)
  res <- check_edp(data, exposure, s, values = 1, categorical_similarity = 0.5)

  subgroup <- data$s[res@results$.id] == "s1"
  supported <- data$s[res@results$.id] == "s0"
  masked <- stats::median(res@results$edp[subgroup])

  expect_gt(masked, 100)
  expect_gt(masked, 0.5 * stats::median(res@results$edp[supported]))
})

# ---- n-scaling ------------------------------------------------------------

test_that("edp per n is stable and edp scales roughly linearly with n", {
  ns <- c(100L, 400L, 1600L)
  supported <- vapply(
    ns,
    function(nn) {
      res <- check_edp(
        sim_edp_gaussian(nn, seed = 1),
        exposure,
        x1,
        values = 0,
        exposure_type = "continuous"
      )
      stats::median(res@results$edp)
    },
    numeric(1)
  )

  edp_per_n <- supported / ns
  expect_lt(stats::sd(edp_per_n) / mean(edp_per_n), 0.15)
  expect_gt(supported[2] / supported[1], 3.0)
  expect_lt(supported[2] / supported[1], 5.0)
})

test_that("a violation point stays at negligible edp per n across n", {
  ns <- c(100L, 400L, 1600L)
  violated <- vapply(
    ns,
    function(nn) {
      res <- check_edp(
        sim_edp_gaussian(nn, seed = 1),
        exposure,
        x1,
        values = 4,
        exposure_type = "continuous"
      )
      stats::median(res@results$edp) / nn
    },
    numeric(1)
  )

  expect_true(all(violated < 0.01))
})

# ---- Estimator variant ----------------------------------------------------

test_that("equal covariate sets guarantee edp_outcome <= edp_treatment", {
  data <- sim_edp_gaussian(300, seed = 1)
  res <- check_edp(
    data,
    exposure,
    x1,
    variant = "estimator",
    values = 0,
    exposure_type = "continuous"
  )

  expect_true(all(
    res@results$edp_outcome <= res@results$edp_treatment + 1e-9
  ))
})

test_that("ideal_weight is inversely related to treatment support", {
  data <- sim_edp_estimator(500, seed = 1)
  res <- check_edp(
    data,
    exposure,
    l,
    variant = "estimator",
    values = 0,
    exposure_type = "continuous"
  )
  d <- res@results

  expect_lt(stats::cor(log(d$ideal_weight), d$edp_treatment), 0)

  cuts <- stats::quantile(d$edp_treatment, c(0.1, 0.9))
  low_support <- d$ideal_weight[d$edp_treatment <= cuts[1]]
  high_support <- d$ideal_weight[d$edp_treatment >= cuts[2]]
  expect_gt(stats::median(low_support), stats::median(high_support))
})

test_that("a g-model-only violation drives edp_treatment below edp_outcome", {
  data <- sim_edp_asymmetric(400, seed = 1)
  res <- check_edp(
    data,
    exposure,
    l,
    .outcome_covariates = l,
    .treatment_covariates = c(l, g),
    variant = "estimator",
    values = 0,
    exposure_type = "continuous"
  )

  subgroup <- data$subgroup[res@results$.id]
  expect_lt(
    stats::median(res@results$edp_treatment[subgroup]),
    0.2 * stats::median(res@results$edp_outcome[subgroup])
  )
})

test_that("ideal_weight is infinite where outcome support is zero", {
  local_quiet()
  # The s1 subgroup is never treated, so at the treated target no observation in
  # its stratum supports the outcome model: edp_outcome is exactly zero and
  # ideal_weight is one over zero.
  data <- sim_edp_masking(200, seed = 1)
  res <- check_edp(
    data,
    exposure,
    s,
    variant = "estimator",
    values = 1,
    categorical_similarity = 0
  )

  subgroup <- data$s[res@results$.id] == "s1"
  expect_true(all(res@results$edp_outcome[subgroup] == 0))
  expect_true(all(is.infinite(res@results$ideal_weight[subgroup])))
})

# ---- Autoplot contract ----------------------------------------------------

test_that("autoplot() returns a ggplot for each data-variant type", {
  data <- sim_edp_gaussian(150)
  res <- check_edp(
    data,
    exposure,
    x1,
    values = c(0, 1),
    exposure_type = "continuous"
  )

  expect_s3_class(ggplot2::autoplot(res, type = "histogram"), "ggplot")
  expect_s3_class(ggplot2::autoplot(res, type = "ecdf"), "ggplot")
})

test_that("autoplot() rejects an unknown type as a classed error", {
  data <- sim_edp_gaussian(150)
  res <- check_edp(
    data,
    exposure,
    x1,
    values = c(0, 1),
    exposure_type = "continuous"
  )

  # A view name is chosen from a fixed menu like any other argument, and the
  # call rendering a figure does not make its failure a lesser one.
  expect_error(
    ggplot2::autoplot(res, type = "bogus"),
    class = "positively_args_error"
  )
})

test_that("plot() draws the view and returns the result invisibly", {
  local_null_device()
  data <- sim_edp_gaussian(150)
  res <- check_edp(
    data,
    exposure,
    x1,
    values = c(0, 1),
    exposure_type = "continuous"
  )
  expect_identical(plot(res), res)
})

test_that("the estimator scatter view is a ggplot and aborts for the data variant", {
  data <- sim_edp_gaussian(150)
  estimator <- check_edp(
    data,
    exposure,
    x1,
    variant = "estimator",
    values = 0,
    exposure_type = "continuous"
  )
  expect_s3_class(ggplot2::autoplot(estimator, type = "scatter"), "ggplot")

  plain <- check_edp(
    data,
    exposure,
    x1,
    values = 0,
    exposure_type = "continuous"
  )
  expect_error(
    ggplot2::autoplot(plain, type = "scatter"),
    class = "positively_error"
  )
})

test_that("the scatter view-gate abort names the autoplot method", {
  local_quiet()
  data <- sim_edp_gaussian(150)
  plain <- check_edp(
    data,
    exposure,
    x1,
    values = 0,
    exposure_type = "continuous"
  )
  expect_snapshot_abort(ggplot2::autoplot(plain, type = "scatter"))
})

test_that("the scatter view separates infinite ideal weights into their own layer", {
  local_quiet()
  # Three strata under exact categorical matching: two treated at different rates
  # give distinct finite ideal weights, and a never-treated stratum has zero
  # outcome support, so its ideal weight is infinite. The finite rows drive the
  # continuous colour gradient; the infinite rows need a fixed distinguishable
  # aesthetic in a separate layer.
  data <- data.frame(
    exposure = c(
      rep(c(0L, 1L), c(4L, 4L)),
      rep(c(0L, 1L), c(4L, 2L)),
      rep(0L, 5L)
    ),
    s = factor(rep(c("a", "b", "c"), c(8L, 6L, 5L)))
  )
  res <- check_edp(
    data,
    exposure,
    s,
    variant = "estimator",
    values = 1,
    categorical_similarity = 0
  )
  weights <- res@results$ideal_weight
  n_finite <- sum(is.finite(weights))
  n_infinite <- sum(is.infinite(weights))
  # The fixture must exercise both the gradient and the fixed layer.
  expect_gt(length(unique(weights[is.finite(weights)])), 1)
  expect_gt(n_infinite, 0)

  built <- ggplot2::ggplot_build(ggplot2::autoplot(res, type = "scatter"))
  layer_rows <- vapply(built$data, nrow, integer(1))
  expect_length(built$data, 2)
  expect_setequal(layer_rows, c(n_finite, n_infinite))

  finite_layer <- built$data[[which(layer_rows == n_finite)]]
  infinite_layer <- built$data[[which(layer_rows == n_infinite)]]
  # The finite rows keep the continuous colour mapping.
  expect_gt(length(unique(finite_layer$colour)), 1)
  # The infinite rows carry a single fixed shape distinct from the finite layer.
  expect_length(unique(infinite_layer$shape), 1)
  expect_false(unique(infinite_layer$shape) %in% unique(finite_layer$shape))
})

test_that("the scatter view builds without a colour label when every weight is infinite", {
  local_quiet()
  # A single never-treated stratum has zero outcome support at the treated target,
  # so every ideal weight is infinite and no finite rows drive the colour mapping.
  # Naming the colour label with no colour aesthetic draws an unknown-label
  # message, so the label must be dropped in this case.
  data <- data.frame(
    exposure = rep(0L, 8L),
    s = factor(rep("a", 8L))
  )
  res <- check_edp(
    data,
    exposure,
    s,
    variant = "estimator",
    values = 1,
    categorical_similarity = 0
  )
  expect_true(all(is.infinite(res@results$ideal_weight)))

  plot <- ggplot2::autoplot(res, type = "scatter")
  expect_no_condition(ggplot2::ggplot_build(plot))

  built <- ggplot2::ggplot_build(plot)
  expect_length(built$data, 1)
  expect_length(unique(built$data[[1]]$shape), 1)
  expect_false(is.null(plot$labels$subtitle))
})

test_that("EDP autoplot views render as expected", {
  local_quiet()
  announce_doppelganger(
    "EDP histogram by intervention value",
    "EDP ECDF by intervention value",
    "EDP estimator scatter",
    "EDP estimator scatter with infinite weights"
  )
  data <- sim_edp_gaussian(150)
  res <- check_edp(
    data,
    exposure,
    x1,
    values = c(0, 1),
    exposure_type = "continuous"
  )
  expect_doppelganger(
    "EDP histogram by intervention value",
    ggplot2::autoplot(res, type = "histogram")
  )
  expect_doppelganger(
    "EDP ECDF by intervention value",
    ggplot2::autoplot(res, type = "ecdf")
  )

  estimator <- check_edp(
    data,
    exposure,
    x1,
    variant = "estimator",
    values = 0,
    exposure_type = "continuous"
  )
  expect_doppelganger(
    "EDP estimator scatter",
    ggplot2::autoplot(estimator, type = "scatter")
  )

  infinite <- data.frame(
    exposure = c(
      rep(c(0L, 1L), c(4L, 4L)),
      rep(c(0L, 1L), c(4L, 2L)),
      rep(0L, 5L)
    ),
    s = factor(rep(c("a", "b", "c"), c(8L, 6L, 5L)))
  )
  res_infinite <- check_edp(
    infinite,
    exposure,
    s,
    variant = "estimator",
    values = 1,
    categorical_similarity = 0
  )
  expect_doppelganger(
    "EDP estimator scatter with infinite weights",
    ggplot2::autoplot(res_infinite, type = "scatter")
  )
})

# ---- Value and estimator-selection validation -----------------------------

test_that("check_edp() validates a user-supplied intervention-value grid", {
  data <- sim_edp_gaussian(60)

  expect_error(
    check_edp(data, exposure, x1, values = "a", exposure_type = "continuous"),
    class = "positively_error"
  )
  expect_error(
    check_edp(
      data,
      exposure,
      x1,
      values = c(0, NA),
      exposure_type = "continuous"
    ),
    class = "positively_error"
  )
  expect_error(
    check_edp(
      data,
      exposure,
      x1,
      values = numeric(0),
      exposure_type = "continuous"
    ),
    class = "positively_error"
  )
})

test_that("check_edp() rejects empty estimator covariate selections", {
  data <- sim_edp_gaussian(60)

  expect_error(
    check_edp(
      data,
      exposure,
      x1,
      .outcome_covariates = tidyselect::starts_with("zzz"),
      variant = "estimator",
      values = 0,
      exposure_type = "continuous"
    ),
    class = "positively_error"
  )
  expect_error(
    check_edp(
      data,
      exposure,
      x1,
      .treatment_covariates = tidyselect::starts_with("zzz"),
      variant = "estimator",
      values = 0,
      exposure_type = "continuous"
    ),
    class = "positively_error"
  )
})

test_that("check_edp() warns when a supplied bandwidth reaches no continuous dimension", {
  local_quiet()

  binary <- dgp_good_positivity(n = 60, seed = 1)
  expect_warning(
    check_edp(binary, exposure, c(x1, x2), bw_exposure = 0.1),
    class = "positively_warning"
  )

  categorical <- sim_edp_masking(60)
  expect_warning(
    check_edp(categorical, exposure, s, bw_covariates = 0.5),
    class = "positively_warning"
  )
})

test_that("the data variant warns when estimator covariates are supplied", {
  local_quiet()
  data <- sim_edp_gaussian(60)
  data$x2 <- rev(data$x1)

  expect_warning(
    check_edp(
      data,
      exposure,
      x1,
      .outcome_covariates = c(x1, x2),
      values = 0,
      exposure_type = "continuous"
    ),
    class = "positively_unused_arg_warning"
  )
  expect_warning(
    check_edp(
      data,
      exposure,
      x1,
      .treatment_covariates = c(x1, x2),
      values = 0,
      exposure_type = "continuous"
    ),
    class = "positively_unused_arg_warning"
  )
})

test_that("the estimator variant does not warn on its estimator covariates", {
  local_quiet()
  data <- sim_edp_gaussian(60)
  data$x2 <- rev(data$x1)

  expect_no_warning(
    check_edp(
      data,
      exposure,
      x1,
      .outcome_covariates = c(x1, x2),
      .treatment_covariates = c(x1, x2),
      variant = "estimator",
      values = 0,
      exposure_type = "continuous"
    ),
    class = "positively_unused_arg_warning"
  )
})

# ---- Snapshots ------------------------------------------------------------

test_that("the argument validation messages are stable", {
  local_quiet()
  data <- sim_edp_gaussian(60)

  expect_snapshot_abort(check_edp(1:10, exposure, x1))
  expect_snapshot_abort(check_edp(
    data,
    exposure,
    tidyselect::starts_with("zzz"),
    exposure_type = "continuous"
  ))
  expect_snapshot_abort(check_edp(
    data,
    c(exposure, x1),
    x1,
    exposure_type = "continuous"
  ))
  expect_snapshot_abort(check_edp(
    data,
    exposure,
    x1,
    categorical_similarity = 1.5,
    exposure_type = "continuous"
  ))
  expect_snapshot_abort(check_edp(
    data,
    exposure,
    x1,
    categorical_similarity = "x",
    exposure_type = "continuous"
  ))
  expect_snapshot_abort(check_edp(
    data,
    exposure,
    x1,
    bw_exposure = -1,
    exposure_type = "continuous"
  ))
  expect_snapshot_abort(check_edp(
    data,
    exposure,
    x1,
    bw_exposure = c(1, 2),
    exposure_type = "continuous"
  ))
  expect_snapshot_abort(check_edp(
    data,
    exposure,
    x1,
    bw_covariates = -1,
    exposure_type = "continuous"
  ))
  expect_snapshot_abort(check_edp(
    data,
    exposure,
    x1,
    variant = "bogus",
    exposure_type = "continuous"
  ))
  expect_snapshot_abort(check_edp(
    data,
    exposure,
    x1,
    kernel = "bogus",
    exposure_type = "continuous"
  ))
})

test_that("the exposure-type error messages are stable", {
  local_quiet()
  data <- sim_edp_gaussian(60)
  character_exposure <- data
  character_exposure$exposure <- rep(c("a", "b", "c"), length.out = nrow(data))

  expect_snapshot_abort(check_edp(
    character_exposure,
    exposure,
    x1,
    exposure_type = "continuous"
  ))
  expect_snapshot_abort(check_edp(
    sim_edp_categorical(300),
    exposure,
    z2,
    exposure_type = "binary"
  ))
  expect_snapshot_abort(check_edp(
    character_exposure,
    exposure,
    x1,
    exposure_type = "bogus"
  ))
})

test_that("the data-integrity error messages are stable", {
  local_quiet()
  data <- sim_edp_gaussian(60)

  one_row <- tibble::tibble(exposure = 0.5, x1 = -0.2)
  expect_snapshot_abort(check_edp(
    one_row,
    exposure,
    x1,
    exposure_type = "continuous"
  ))

  na_exposure <- data
  na_exposure$exposure[1] <- NA
  expect_snapshot_abort(check_edp(
    na_exposure,
    exposure,
    x1,
    exposure_type = "continuous"
  ))

  na_covariate <- data
  na_covariate$x1[1] <- NA
  expect_snapshot_abort(check_edp(
    na_covariate,
    exposure,
    x1,
    exposure_type = "continuous"
  ))
})

test_that("the value and estimator-selection error messages are stable", {
  local_quiet()
  data <- sim_edp_gaussian(60)

  expect_snapshot_abort(check_edp(
    data,
    exposure,
    x1,
    values = "a",
    exposure_type = "continuous"
  ))
  expect_snapshot_abort(check_edp(
    data,
    exposure,
    x1,
    values = c(0, NA),
    exposure_type = "continuous"
  ))
  expect_snapshot_abort(check_edp(
    data,
    exposure,
    x1,
    values = numeric(0),
    exposure_type = "continuous"
  ))
  expect_snapshot_abort(check_edp(
    data,
    exposure,
    x1,
    .outcome_covariates = tidyselect::starts_with("zzz"),
    variant = "estimator",
    values = 0,
    exposure_type = "continuous"
  ))
  expect_snapshot_abort(check_edp(
    data,
    exposure,
    x1,
    .treatment_covariates = tidyselect::starts_with("zzz"),
    variant = "estimator",
    values = 0,
    exposure_type = "continuous"
  ))
})

test_that("the unused-bandwidth warnings are stable", {
  local_quiet()
  withr::local_options(warn = 0)

  binary <- dgp_good_positivity(n = 60, seed = 1)
  expect_snapshot(
    res <- check_edp(binary, exposure, c(x1, x2), bw_exposure = 0.1)
  )

  categorical <- sim_edp_masking(60)
  expect_snapshot(
    res <- check_edp(categorical, exposure, s, bw_covariates = 0.5)
  )
})

test_that("the data-variant unused estimator-covariate warning is stable", {
  local_quiet()
  withr::local_options(warn = 0)
  data <- sim_edp_gaussian(60)
  data$x2 <- rev(data$x1)
  expect_snapshot(
    res <- check_edp(
      data,
      exposure,
      x1,
      .outcome_covariates = c(x1, x2),
      values = 0,
      exposure_type = "continuous"
    )
  )
})

test_that("the print method is stable", {
  local_quiet()

  continuous <- check_edp(
    sim_edp_gaussian(150),
    exposure,
    x1,
    values = c(0, 1),
    exposure_type = "continuous"
  )
  expect_snapshot(print(continuous))

  binary <- check_edp(
    dgp_good_positivity(n = 200, seed = 1),
    exposure,
    c(x1, x2)
  )
  expect_snapshot(print(binary))

  estimator <- check_edp(
    sim_edp_gaussian(150),
    exposure,
    x1,
    variant = "estimator",
    values = 0,
    exposure_type = "continuous"
  )
  expect_snapshot(print(estimator))
})

# ---- Display methods -------------------------------------------------------

# Seven intervention values under the data variant, which is the variant whose
# single edp column the range is read from.
edp_display_result <- function() {
  check_edp(
    sim_edp_gaussian(150),
    exposure,
    x1,
    values = c(-3, -2, -1, 0, 1, 2, 3),
    exposure_type = "continuous"
  )
}

test_that("diagnostic_label() names EDP without naming its class", {
  local_quiet()
  res <- edp_display_result()
  label <- expect_readable_label(res)

  printed <- printed_text(res)
  expect_no_match(printed, S7::S7_class(res)@name, fixed = TRUE)
  expect_match(printed, label, fixed = TRUE)
})

test_that("diagnostic_headline() reads the variant and the value count", {
  local_quiet()
  res <- edp_display_result()
  headline <- expect_readable_headline(res)
  text <- rendered_text(headline)

  # The two variants measure different quantities over the same intervention
  # values, so a reading that omitted the variant would be ambiguous.
  expect_match(tolower(text), "data", fixed = TRUE)
  expect_match(text, "\\b7\\b")
})

test_that("the EDP label and headline are stable", {
  local_quiet()
  res <- edp_display_result()
  expect_snapshot({
    diagnostic_label(res)
    writeLines(diagnostic_headline(res))
  })
})
