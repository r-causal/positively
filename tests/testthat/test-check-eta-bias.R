# ETA.Bias is a parametric bootstrap: the exact numbers depend on R's glm fits
# and RNG, so every magnitude claim is written as a direction, ordering, or a
# multiple of the run's own Monte Carlo error rather than a hard-coded value.
# Each fixed-seed dataset yields one realization; thresholds sit well inside the
# across-sample spread reported in the simulation study so a single draw stays on
# the correct side of every inequality. n_boot is held at 200 (100 where only a
# sign or ordering matters) because that is where the study calibrated the
# tolerances below. The bootstrap treatment is drawn from the untruncated fitted
# propensity, so truncation changes only the estimator; the truncation sweep
# tests read that as increasing bias and shrinking bootstrap spread.

# ---- Local scenario generators --------------------------------------------
# The shared helpers in helper-dgp.R carry no outcome, and ETA.Bias needs a
# binary exposure, an outcome, and a tunable degree of positivity violation.
# These seeded generators supply all three. A linear outcome model with a
# constant treatment effect tau means the correctly specified G-computation
# estimate recovers tau, so `truth` lands near tau.

# Good overlap: the true propensity is bounded inside (0.2, 0.8), so no fitted
# score is extreme and every estimator's ETA.Bias sits near zero.
sim_eta_good <- function(n = 1000, tau = 0.5, seed = 1) {
  withr::local_seed(seed)
  x1 <- stats::rnorm(n)
  x2 <- stats::rnorm(n)
  ps <- 0.2 + 0.6 * stats::plogis(0.5 * x1 - 0.5 * x2)
  a <- stats::rbinom(n, 1L, ps)
  y <- tau * a + x1 + x2 + stats::rnorm(n)
  tibble::tibble(a = a, y = y, x1 = x1, x2 = x2)
}

# Freedman-Berk style practical violation: a steep propensity model pushes
# fitted scores toward 0 and 1. `steepness` dials the severity; the default is
# the study's moderate scenario, and the outcome model is linear with effect
# tau, so `truth` is near tau.
sim_eta_violation <- function(n = 1000, steepness = 2, tau = 1, seed = 1) {
  withr::local_seed(seed)
  x1 <- stats::rnorm(n)
  x2 <- stats::rnorm(n)
  ps <- stats::plogis(steepness * (x1 + x2))
  a <- stats::rbinom(n, 1L, ps)
  y <- tau * a + x1 + x2 + stats::rnorm(n)
  tibble::tibble(a = a, y = y, x1 = x1, x2 = x2)
}

# The same practical violation with a binary outcome, so all bias magnitudes
# shrink under the bounded scale while the estimator ordering is preserved.
sim_eta_binary_outcome <- function(n = 1000, seed = 1) {
  withr::local_seed(seed)
  x1 <- stats::rnorm(n)
  x2 <- stats::rnorm(n)
  ps <- stats::plogis(2 * (x1 + x2))
  a <- stats::rbinom(n, 1L, ps)
  py <- stats::plogis(-0.5 + 0.5 * a + x1 + x2)
  y <- stats::rbinom(n, 1L, py)
  tibble::tibble(a = a, y = y, x1 = x1, x2 = x2)
}

# pos_violations carries the factor covariate region, which defines that
# dataset's structural violation, but no outcome. A seeded linear outcome with a
# region effect gives check_eta_bias() a factor covariate to condition on while
# keeping truth near the treatment effect.
sim_eta_factor <- function(seed = 1) {
  withr::local_seed(seed)
  data <- pos_violations
  region_effect <- ifelse(data$region == "b", 0.5, -0.5)
  data$y <- data$exposure +
    data$x1 +
    data$x2 +
    region_effect +
    stats::rnorm(nrow(data))
  data
}

# A two-level covariate whose second level is held by only two of 200 rows.
# A bootstrap resample can miss both rows and leave region with a single level.
# The outcome carries a region effect so truth stays near the treatment effect.
# `character` returns region as a character vector rather than a factor.
sim_eta_rare_level <- function(seed = 2, n = 200, character = FALSE) {
  withr::local_seed(seed)
  x1 <- stats::rnorm(n)
  x2 <- stats::rnorm(n)
  region <- c(rep("a", n - 2), rep("b", 2))
  if (!character) {
    region <- factor(region, levels = c("a", "b"))
  }
  a <- stats::rbinom(n, 1L, stats::plogis(x1))
  y <- a + x1 + (region == "b") + stats::rnorm(n)
  tibble::tibble(a = a, y = y, x1 = x1, x2 = x2, region = region)
}

# A balanced two-level character covariate with no rare level, so every
# bootstrap resample holds both levels. The outcome carries a region effect.
sim_eta_character <- function(seed = 1, n = 400) {
  withr::local_seed(seed)
  x1 <- stats::rnorm(n)
  x2 <- stats::rnorm(n)
  region <- sample(c("a", "b"), n, replace = TRUE)
  a <- stats::rbinom(n, 1L, stats::plogis(0.5 * x1))
  y <- a + x1 + (region == "b") + stats::rnorm(n)
  tibble::tibble(a = a, y = y, x1 = x1, x2 = x2, region = region)
}

# ---- Local accessors ------------------------------------------------------
# Each check_eta_bias() call is a single estimator, so a single-run result is a
# one-row tibble. A fixed seed makes the bootstrap draw reproducible.

fit_eta <- function(data, estimator, ..., seed = 2024) {
  withr::local_seed(seed)
  check_eta_bias(
    data,
    a,
    y,
    c(x1, x2),
    estimator = estimator,
    ...
  )
}

bias_of <- function(res) res@results$bias[[1]]
mc_se_of <- function(res) res@results$mc_se[[1]]

# ---- Argument validation --------------------------------------------------

test_that("check_eta_bias() rejects non-data-frame input", {
  local_quiet()
  expect_error(
    check_eta_bias(1:10, a, y, c(x1, x2), n_boot = 10),
    class = "positively_error"
  )
})

test_that("check_eta_bias() aborts on a continuous exposure", {
  local_quiet()
  data <- sim_eta_good(n = 100, seed = 1)
  data$a <- stats::rnorm(nrow(data))
  expect_error(
    check_eta_bias(data, a, y, c(x1, x2), n_boot = 10),
    class = "positively_error"
  )
})

test_that("check_eta_bias() aborts on a categorical exposure", {
  local_quiet()
  data <- sim_eta_good(n = 120, seed = 1)
  data$a <- factor(sample(c("a", "b", "c"), nrow(data), replace = TRUE))
  expect_error(
    check_eta_bias(data, a, y, c(x1, x2), n_boot = 10),
    class = "positively_error"
  )
})

test_that("check_eta_bias() rejects an empty covariate selection", {
  local_quiet()
  data <- sim_eta_good(n = 100, seed = 1)
  expect_error(
    check_eta_bias(
      data,
      a,
      y,
      tidyselect::starts_with("zzz"),
      n_boot = 10
    ),
    class = "positively_error"
  )
})

test_that("check_eta_bias() rejects covariates overlapping the exposure or outcome", {
  local_quiet()
  # A small design keeps a low n_boot; the abort must fire ahead of any
  # bootstrap work, so the resampling size never matters here.
  data <- sim_eta_good(n = 100, seed = 1)
  expect_error(
    check_eta_bias(data, a, y, c(a, x1, x2), n_boot = 5),
    class = "positively_selection_error"
  )
  expect_error(
    check_eta_bias(data, a, y, c(y, x1, x2), n_boot = 5),
    class = "positively_selection_error"
  )
  expect_error(
    check_eta_bias(data, a, y, tidyselect::everything(), n_boot = 5),
    class = "positively_selection_error"
  )
})

test_that("check_eta_bias() aborts on a missing outcome column", {
  local_quiet()
  data <- sim_eta_good(n = 100, seed = 1)
  expect_error(
    check_eta_bias(data, a, not_a_column, c(x1, x2), n_boot = 10),
    class = "positively_error"
  )
})

test_that("check_eta_bias() aborts on missing exposure or outcome values", {
  local_quiet()
  data <- sim_eta_good(n = 100, seed = 1)

  na_exposure <- data
  na_exposure$a[1] <- NA
  expect_error(
    check_eta_bias(na_exposure, a, y, c(x1, x2), n_boot = 10),
    class = "positively_error"
  )

  na_outcome <- data
  na_outcome$y[1] <- NA
  expect_error(
    check_eta_bias(na_outcome, a, y, c(x1, x2), n_boot = 10),
    class = "positively_error"
  )
})

test_that("check_eta_bias() validates the estimator against a fixed set", {
  local_quiet()
  data <- sim_eta_good(n = 80, seed = 1)
  expect_error(
    check_eta_bias(data, a, y, c(x1, x2), estimator = "banana", n_boot = 10),
    class = "rlang_error"
  )
})

test_that("check_eta_bias() validates outcome_type against a fixed set", {
  local_quiet()
  data <- sim_eta_good(n = 80, seed = 1)
  expect_error(
    check_eta_bias(
      data,
      a,
      y,
      c(x1, x2),
      outcome_type = "count",
      n_boot = 10
    ),
    class = "rlang_error"
  )
})

test_that("check_eta_bias() validates error_dist against a fixed set", {
  local_quiet()
  data <- sim_eta_good(n = 80, seed = 1)
  expect_error(
    check_eta_bias(
      data,
      a,
      y,
      c(x1, x2),
      error_dist = "poisson",
      n_boot = 10
    ),
    class = "rlang_error"
  )
})

test_that("check_eta_bias() rejects a malformed truncation bound", {
  local_quiet()
  data <- sim_eta_good(n = 80, seed = 1)

  # Not length two.
  expect_error(
    check_eta_bias(data, a, y, c(x1, x2), truncation = 0.05, n_boot = 10),
    class = "positively_error"
  )
  expect_error(
    check_eta_bias(
      data,
      a,
      y,
      c(x1, x2),
      truncation = c(0.05, 0.5, 0.95),
      n_boot = 10
    ),
    class = "positively_error"
  )
  # Outside the unit interval.
  expect_error(
    check_eta_bias(
      data,
      a,
      y,
      c(x1, x2),
      truncation = c(-0.1, 0.95),
      n_boot = 10
    ),
    class = "positively_error"
  )
  # Lower not below upper.
  expect_error(
    check_eta_bias(
      data,
      a,
      y,
      c(x1, x2),
      truncation = c(0.6, 0.4),
      n_boot = 10
    ),
    class = "positively_error"
  )
})

test_that("check_eta_bias() rejects a malformed truncation_grid", {
  local_quiet()
  data <- sim_eta_good(n = 80, seed = 1)

  expect_error(
    check_eta_bias(
      data,
      a,
      y,
      c(x1, x2),
      truncation_grid = c("a", "b"),
      n_boot = 10
    ),
    class = "positively_error"
  )
  # Lower bounds must sit in (0, 0.5) so upper = 1 - lower stays above lower.
  expect_error(
    check_eta_bias(
      data,
      a,
      y,
      c(x1, x2),
      truncation_grid = c(0.1, 0.6),
      n_boot = 10
    ),
    class = "positively_error"
  )
  expect_error(
    check_eta_bias(
      data,
      a,
      y,
      c(x1, x2),
      truncation_grid = c(0.1, NA),
      n_boot = 10
    ),
    class = "positively_error"
  )
})

test_that("check_eta_bias() rejects a non-whole or non-positive n_boot", {
  local_quiet()
  data <- sim_eta_good(n = 80, seed = 1)
  expect_error(
    check_eta_bias(data, a, y, c(x1, x2), n_boot = 0),
    class = "positively_error"
  )
  expect_error(
    check_eta_bias(data, a, y, c(x1, x2), n_boot = 2.5),
    class = "positively_error"
  )
  expect_error(
    check_eta_bias(data, a, y, c(x1, x2), n_boot = -10),
    class = "positively_error"
  )
})

test_that("check_eta_bias() accepts explicit model formula overrides", {
  local_quiet()
  data <- sim_eta_violation(n = 400, seed = 1)
  res <- fit_eta(
    data,
    "gcomp",
    exposure_formula = a ~ x1 + x2,
    outcome_formula = y ~ a + x1 + x2,
    n_boot = 100
  )
  expect_identical(S7::S7_class(res)@name, "eta_bias_result")
})

test_that("check_eta_bias() accepts a two-level factor exposure", {
  local_quiet()
  data <- sim_eta_good(n = 300, seed = 1)
  factor_data <- data
  factor_data$a <- factor(
    ifelse(data$a == 1, "treated", "control"),
    levels = c("control", "treated")
  )

  numeric_result <- fit_eta(data, "gcomp", n_boot = 50)
  factor_result <- fit_eta(factor_data, "gcomp", n_boot = 50)

  expect_identical(factor_result@exposure_type, "binary")
  # "control" maps to 0 and "treated" to 1, so the factor and numeric fits
  # face the same 0/1 coding and agree.
  expect_equal(bias_of(factor_result), bias_of(numeric_result))
})

test_that("check_eta_bias() rejects an unsupported covariate type", {
  local_quiet()
  data <- sim_eta_good(n = 60, seed = 1)
  data$bad <- as.list(seq_len(nrow(data)))
  expect_error(
    check_eta_bias(data, a, y, c(x1, bad), n_boot = 10),
    class = "positively_error"
  )
  expect_snapshot(
    check_eta_bias(data, a, y, c(x1, bad), n_boot = 10),
    error = TRUE
  )
})

test_that("check_eta_bias() aborts on a non-0/1 binary outcome", {
  local_quiet()
  data <- sim_eta_good(n = 100, seed = 1)
  expect_error(
    check_eta_bias(
      data,
      a,
      y,
      c(x1, x2),
      outcome_type = "binary",
      n_boot = 10
    ),
    class = "positively_error"
  )
})

test_that("check_eta_bias() runs with a constant covariate", {
  local_quiet()
  data <- sim_eta_violation(n = 300, seed = 1)
  data$x3 <- 1

  matrix_path <- withr::with_seed(
    2024,
    check_eta_bias(data, a, y, c(x1, x2, x3), estimator = "aipw", n_boot = 50)
  )
  formula_path <- withr::with_seed(
    2024,
    check_eta_bias(
      data,
      a,
      y,
      c(x1, x2, x3),
      estimator = "aipw",
      exposure_formula = a ~ x1 + x2 + x3,
      outcome_formula = y ~ a + x1 + x2 + x3,
      n_boot = 50
    )
  )

  expect_identical(S7::S7_class(matrix_path)@name, "eta_bias_result")
  # An aliased constant covariate contributes nothing, so the default matrix
  # path and the explicit formula path agree. The tolerance absorbs the tiny
  # numerical differences between glm.fit and glm accumulated over the refits.
  expect_equal(bias_of(matrix_path), bias_of(formula_path), tolerance = 1e-2)
})

test_that("check_eta_bias() accepts a factor covariate and returns finite estimates", {
  local_quiet()
  data <- sim_eta_factor(seed = 1)
  for (estimator in c("ipw", "gcomp", "aipw")) {
    res <- withr::with_seed(
      2024,
      check_eta_bias(
        data,
        exposure,
        y,
        c(x1, x2, region),
        estimator = estimator,
        n_boot = 50
      )
    )
    expect_identical(S7::S7_class(res)@name, "eta_bias_result")
    expect_true(is.finite(res@results$bias[[1]]))
    expect_true(is.finite(res@results$mc_se[[1]]))
    expect_true(is.finite(res@truth))
  }
})

test_that("the truncation sweep runs with a factor covariate present", {
  local_quiet()
  data <- sim_eta_factor(seed = 1)
  grid <- c(0, 0.05, 0.1)
  res <- withr::with_seed(
    2024,
    check_eta_bias(
      data,
      exposure,
      y,
      c(x1, x2, region),
      estimator = "ipw",
      truncation_grid = grid,
      n_boot = 50
    )
  )
  expect_identical(nrow(res@results), length(grid))
  expect_setequal(res@results$truncation_lower, grid)
  expect_true(all(is.finite(res@results$bias)))
  expect_length(res@boot_estimates, length(grid))
})

test_that("a rare factor level survives the bootstrap for every estimator", {
  local_quiet()
  data <- sim_eta_rare_level(seed = 2)
  for (estimator in c("gcomp", "ipw", "aipw")) {
    res <- withr::with_seed(
      9,
      check_eta_bias(
        data,
        a,
        y,
        c(x1, x2, region),
        estimator = estimator,
        n_boot = 50
      )
    )
    expect_true(is.finite(res@results$bias[[1]]))
    expect_true(is.finite(res@results$mc_se[[1]]))
    expect_true(is.finite(res@truth))
  }
})

test_that("a rare character level survives the bootstrap for every estimator", {
  local_quiet()
  data <- sim_eta_rare_level(seed = 2, character = TRUE)
  for (estimator in c("gcomp", "ipw", "aipw")) {
    res <- withr::with_seed(
      9,
      check_eta_bias(
        data,
        a,
        y,
        c(x1, x2, region),
        estimator = estimator,
        n_boot = 50
      )
    )
    expect_true(is.finite(res@results$bias[[1]]))
    expect_true(is.finite(res@results$mc_se[[1]]))
    expect_true(is.finite(res@truth))
  }
})

test_that("a rare level survives the bootstrap under formula overrides", {
  local_quiet()
  data <- sim_eta_rare_level(seed = 2)
  for (estimator in c("gcomp", "ipw", "aipw")) {
    res <- withr::with_seed(
      9,
      check_eta_bias(
        data,
        a,
        y,
        c(x1, x2, region),
        estimator = estimator,
        exposure_formula = a ~ x1 + x2 + region,
        outcome_formula = y ~ a + x1 + x2 + region,
        n_boot = 50
      )
    )
    expect_true(is.finite(res@results$bias[[1]]))
    expect_true(is.finite(res@results$mc_se[[1]]))
    expect_true(is.finite(res@truth))
  }
})

test_that("a balanced character covariate returns finite estimates", {
  local_quiet()
  data <- sim_eta_character(seed = 1)
  for (estimator in c("gcomp", "ipw", "aipw")) {
    res <- withr::with_seed(
      2024,
      check_eta_bias(
        data,
        a,
        y,
        c(x1, x2, region),
        estimator = estimator,
        n_boot = 50
      )
    )
    expect_true(is.finite(res@results$bias[[1]]))
    expect_true(is.finite(res@results$mc_se[[1]]))
    expect_true(is.finite(res@truth))
  }
})

test_that("the numeric matrix path is unchanged under a fixed seed (regression)", {
  local_quiet()
  data <- sim_eta_violation(n = 400, seed = 1)
  res <- withr::with_seed(
    2024,
    check_eta_bias(data, a, y, c(x1, x2), estimator = "ipw", n_boot = 100)
  )
  # Baseline captured from the implementation before factor covariates were
  # supported. Numeric covariates still take the matrix path untouched, so it
  # reproduces these values exactly.
  expect_equal(bias_of(res), 0.17472955288621894)
  expect_equal(mc_se_of(res), 0.05203235693919487)
  expect_equal(res@truth, 1.0636370442788932)
  expect_equal(res@results$boot_mean[[1]], 1.2383665971651121)
})

# ---- Result class and properties ------------------------------------------

test_that("check_eta_bias() returns an eta_bias_result diagnostic", {
  local_quiet()
  data <- sim_eta_good(n = 300, seed = 1)
  res <- fit_eta(data, "ipw", n_boot = 50)

  expect_true(S7::S7_inherits(res, positivity_diagnostic))
  expect_identical(S7::S7_class(res)@name, "eta_bias_result")
  expect_identical(res@exposure_type, "binary")
  expect_identical(res@n, 300L)
})

test_that("eta_bias_result carries the estimator, truth, and boot properties", {
  local_quiet()
  data <- sim_eta_good(n = 300, seed = 1)
  res <- fit_eta(data, "aipw", n_boot = 50)

  expect_identical(res@estimator, "aipw")
  expect_type(res@truth, "double")
  expect_length(res@truth, 1)
  expect_type(res@boot_estimates, "list")
  # A single run holds one truncation level, so one bootstrap-estimate vector.
  expect_length(res@boot_estimates, 1)
  expect_length(res@boot_estimates[[1]], 50)
  expect_type(res@boot_estimates[[1]], "double")
})

test_that("eta_bias_result has the fixed results columns", {
  local_quiet()
  data <- sim_eta_good(n = 300, seed = 1)
  res <- fit_eta(data, "ipw", n_boot = 50)

  expect_s3_class(res@results, "tbl_df")
  expect_setequal(
    names(res@results),
    c("truncation_lower", "truncation_upper", "bias", "mc_se", "boot_mean")
  )
  expect_false("flagged" %in% names(res@results))
})

# ---- Results shape (single run versus truncation grid) --------------------

test_that("a single run yields one row spanning the full weight range", {
  local_quiet()
  data <- sim_eta_good(n = 300, seed = 1)
  res <- fit_eta(data, "ipw", n_boot = 50)

  expect_identical(nrow(res@results), 1L)
  # No truncation means the retained range is the whole unit interval, and the
  # bounds obey the sweep relation upper = 1 - lower.
  expect_equal(res@results$truncation_lower, 0)
  expect_equal(res@results$truncation_upper, 1)
})

test_that("an explicit truncation bound is reported verbatim in one row", {
  local_quiet()
  data <- sim_eta_violation(n = 400, seed = 1)
  res <- fit_eta(data, "ipw", truncation = c(0.05, 0.95), n_boot = 50)

  expect_identical(nrow(res@results), 1L)
  expect_equal(res@results$truncation_lower, 0.05)
  expect_equal(res@results$truncation_upper, 0.95)
})

test_that("a truncation grid yields one row per lower bound with upper = 1 - lower", {
  local_quiet()
  data <- sim_eta_violation(n = 500, seed = 1)
  grid <- c(0, 0.025, 0.05, 0.1)
  res <- fit_eta(data, "ipw", truncation_grid = grid, n_boot = 100)

  expect_identical(nrow(res@results), length(grid))
  expect_setequal(res@results$truncation_lower, grid)
  expect_equal(
    res@results$truncation_upper,
    1 - res@results$truncation_lower
  )
  # One bootstrap-estimate vector per swept level.
  expect_length(res@boot_estimates, length(grid))
  expect_true(all(lengths(res@boot_estimates) == 100))
})

# ---- Statistical behavior (simulation REPORT expectations) ----------------

test_that("good overlap leaves every estimator near zero (expectation 1)", {
  local_quiet()
  data <- sim_eta_good(n = 1000, seed = 1)
  for (estimator in c("ipw", "gcomp", "aipw")) {
    res <- fit_eta(data, estimator, n_boot = 200)
    expect_lt(abs(bias_of(res)), 0.05)
  }
})

test_that("under violation gcomp stays near zero and ipw is large (expectations 2-5)", {
  local_quiet()
  data <- sim_eta_violation(n = 1000, seed = 1)
  ipw <- fit_eta(data, "ipw", n_boot = 500)
  gcomp <- fit_eta(data, "gcomp", n_boot = 200)
  aipw <- fit_eta(data, "aipw", n_boot = 200)

  # Expectation 2: gcomp is the most robust invariant, since truth is its own
  # point estimate; its residual bias is Monte Carlo noise.
  expect_lt(abs(bias_of(gcomp)), 3 * mc_se_of(gcomp))
  expect_lt(abs(bias_of(gcomp)), 0.03)

  # Expectation 3: ipw carries substantial, clearly non-null bias. The
  # multiple is 3 rather than the simulation REPORT's 5, and this fit uses
  # n_boot = 500: for this DGP the expected bias-to-mc_se ratio at
  # n_boot = 200 is about 4, so a 5-multiple holds only for lucky seeds.
  # At n_boot = 500 the expected ratio is about 6.3 with an across-seed
  # spread near 1, leaving a wide margin for RNG and refactor changes.
  expect_gt(bias_of(ipw), 0.05)
  expect_gt(bias_of(ipw), 3 * mc_se_of(ipw))

  # Expectation 4: ipw dominates both consistent estimators by a wide margin.
  expect_gt(abs(bias_of(ipw)), 3 * abs(bias_of(aipw)))
  expect_gt(abs(bias_of(ipw)), 3 * abs(bias_of(gcomp)))

  # Expectation 5: the doubly robust estimator stays near zero, though looser
  # than gcomp because of its larger across-sample spread.
  expect_lt(abs(bias_of(aipw)), 0.05)
})

test_that("truncation tightens ipw bias upward and its spread downward (expectation 6)", {
  local_quiet()
  data <- sim_eta_violation(n = 1000, seed = 1)
  grid <- c(0, 0.025, 0.05, 0.1)
  res <- fit_eta(data, "ipw", truncation_grid = grid, n_boot = 200)

  ord <- order(res@results$truncation_lower)
  bias <- res@results$bias[ord]
  boot_sd <- vapply(res@boot_estimates[ord], stats::sd, numeric(1))

  expect_true(all(diff(abs(bias)) > 0))
  expect_true(all(diff(boot_sd) < 0))

  # mc_se is boot_sd / sqrt(n_boot); the two recoveries of boot_sd must agree.
  expect_equal(res@results$mc_se[ord] * sqrt(200), boot_sd)
})

test_that("gcomp bias is flat across the truncation grid (internal consistency)", {
  local_quiet()
  data <- sim_eta_violation(n = 1000, seed = 1)
  grid <- c(0, 0.025, 0.05, 0.1)
  res <- fit_eta(data, "gcomp", truncation_grid = grid, n_boot = 100)

  # G-computation uses no propensity, so truncating weights cannot move it.
  expect_equal(diff(res@results$bias), rep(0, length(grid) - 1))
})

test_that("ipw bias grows with violation severity (expectation 7)", {
  local_quiet()
  bias_null <- bias_of(fit_eta(
    sim_eta_good(n = 1000, seed = 1),
    "ipw",
    n_boot = 100
  ))
  bias_moderate <- bias_of(fit_eta(
    sim_eta_violation(n = 1000, steepness = 2, seed = 1),
    "ipw",
    n_boot = 100
  ))
  bias_severe <- bias_of(fit_eta(
    sim_eta_violation(n = 1000, steepness = 3.5, seed = 1),
    "ipw",
    n_boot = 100
  ))

  expect_true(bias_null < bias_moderate)
  expect_true(bias_moderate < bias_severe)
})

test_that("Monte Carlo error scales as one over root n_boot (expectation 8)", {
  local_quiet()
  data <- sim_eta_violation(n = 1000, seed = 1)
  mc_100 <- mc_se_of(fit_eta(data, "ipw", n_boot = 100))
  mc_400 <- mc_se_of(fit_eta(data, "ipw", n_boot = 400))

  expect_lt(mc_400, mc_100)
  ratio <- mc_400 / mc_100
  # Expected ratio is sqrt(100 / 400) = 0.5; allow the study's 25 percent band.
  expect_gt(ratio, 0.375)
  expect_lt(ratio, 0.625)
})

test_that("binary outcomes run and preserve the estimator ordering (expectation 9)", {
  local_quiet()
  data <- sim_eta_binary_outcome(n = 1000, seed = 1)
  ipw <- fit_eta(data, "ipw", outcome_type = "binary", n_boot = 200)
  gcomp <- fit_eta(data, "gcomp", outcome_type = "binary", n_boot = 200)
  aipw <- fit_eta(data, "aipw", outcome_type = "binary", n_boot = 200)

  expect_gt(abs(bias_of(ipw)), abs(bias_of(aipw)))
  expect_gt(abs(bias_of(aipw)), abs(bias_of(gcomp)))
  expect_lt(abs(bias_of(gcomp)), 0.02)
})

test_that("truth is a scalar equal to the gcomp estimate near tau (expectation 10)", {
  local_quiet()
  data <- sim_eta_violation(n = 1000, tau = 1, seed = 1)
  ipw <- fit_eta(data, "ipw", n_boot = 100)
  gcomp <- fit_eta(data, "gcomp", n_boot = 100)

  # truth depends on the outcome model only, so it is identical across
  # estimators and constant across truncation levels of one run.
  expect_equal(ipw@truth, gcomp@truth)

  swept <- fit_eta(
    data,
    "ipw",
    truncation_grid = c(0, 0.1),
    n_boot = 100
  )
  expect_length(swept@truth, 1)

  # A correctly specified linear outcome model with effect tau = 1 recovers it.
  expect_equal(gcomp@truth, 1, tolerance = 0.15)
})

test_that("normal and empirical error models agree under violation (expectation 6, error_dist)", {
  local_quiet()
  data <- sim_eta_violation(n = 1000, seed = 1)
  normal <- fit_eta(data, "ipw", error_dist = "normal", n_boot = 200)
  empirical <- fit_eta(data, "ipw", error_dist = "empirical", n_boot = 200)

  expect_gt(bias_of(normal), 0.05)
  expect_gt(bias_of(empirical), 0.05)
  expect_lt(abs(bias_of(normal) - bias_of(empirical)), 0.15)
})

# ---- tidy() and glance() --------------------------------------------------

test_that("tidy() and glance() follow the shared diagnostic contract", {
  local_quiet()
  data <- sim_eta_good(n = 300, seed = 1)
  res <- fit_eta(data, "ipw", n_boot = 50)

  expect_identical(generics::tidy(res), res@results)
  glanced <- generics::glance(res)
  expect_s3_class(glanced, "tbl_df")
  expect_identical(nrow(glanced), 1L)
})

# ---- Autoplot contract ----------------------------------------------------

test_that("autoplot() returns a ggplot for each type", {
  local_quiet()
  data <- sim_eta_violation(n = 500, seed = 1)
  single <- fit_eta(data, "ipw", n_boot = 100)
  swept <- fit_eta(
    data,
    "ipw",
    truncation_grid = c(0, 0.025, 0.05, 0.1),
    n_boot = 100
  )

  expect_s3_class(ggplot2::autoplot(single, type = "bootstrap"), "ggplot")
  expect_s3_class(ggplot2::autoplot(swept, type = "sweep"), "ggplot")
})

test_that("ETA bias autoplot views render as expected", {
  local_quiet()
  data <- sim_eta_violation(n = 500, seed = 1)
  single <- fit_eta(data, "ipw", n_boot = 100)
  swept <- fit_eta(
    data,
    "ipw",
    truncation_grid = c(0, 0.025, 0.05, 0.1),
    n_boot = 100
  )

  expect_doppelganger(
    "ETA bias bootstrap distribution",
    ggplot2::autoplot(single, type = "bootstrap")
  )
  expect_doppelganger(
    "ETA bias truncation sweep",
    ggplot2::autoplot(swept, type = "sweep")
  )
})

test_that("the sweep view aborts on a single-run result", {
  local_quiet()
  data <- sim_eta_violation(n = 500, seed = 1)
  single <- fit_eta(data, "ipw", n_boot = 50)
  expect_error(
    ggplot2::autoplot(single, type = "sweep"),
    class = "positively_sweep_absent_error"
  )
  expect_snapshot(
    ggplot2::autoplot(single, type = "sweep"),
    error = TRUE
  )
})

# ---- Print method ---------------------------------------------------------

test_that("the single-run print method is stable", {
  local_quiet()
  data <- sim_eta_good(n = 300, seed = 1)
  res <- fit_eta(data, "ipw", n_boot = 50)
  expect_snapshot(print(res))
})

test_that("the truncation-sweep print method is stable", {
  local_quiet()
  data <- sim_eta_violation(n = 400, seed = 1)
  res <- fit_eta(data, "ipw", truncation_grid = c(0, 0.05, 0.1), n_boot = 50)
  expect_snapshot(print(res))
})

# ---- Classed error messages -----------------------------------------------

test_that("the input and exposure validation messages are stable", {
  local_quiet()
  good <- sim_eta_good(n = 100, seed = 1)

  continuous <- good
  continuous$a <- stats::rnorm(nrow(continuous))

  categorical <- good
  categorical$a <- factor(sample(c("a", "b", "c"), nrow(categorical), TRUE))

  na_exposure <- good
  na_exposure$a[1] <- NA

  expect_snapshot(
    check_eta_bias(1:10, a, y, c(x1, x2), n_boot = 10),
    error = TRUE
  )
  expect_snapshot(
    check_eta_bias(continuous, a, y, c(x1, x2), n_boot = 10),
    error = TRUE
  )
  expect_snapshot(
    check_eta_bias(categorical, a, y, c(x1, x2), n_boot = 10),
    error = TRUE
  )
  expect_snapshot(
    check_eta_bias(good, a, y, tidyselect::starts_with("zzz"), n_boot = 10),
    error = TRUE
  )
  expect_snapshot(
    check_eta_bias(good, a, not_a_column, c(x1, x2), n_boot = 10),
    error = TRUE
  )
  expect_snapshot(
    check_eta_bias(na_exposure, a, y, c(x1, x2), n_boot = 10),
    error = TRUE
  )
  expect_snapshot(
    check_eta_bias(good, a, y, c(x1, x2), outcome_type = "binary", n_boot = 10),
    error = TRUE
  )
  expect_snapshot(
    check_eta_bias(good, a, y, c(a, x1, x2), n_boot = 10),
    error = TRUE
  )
})

test_that("the truncation and n_boot validation messages are stable", {
  local_quiet()
  data <- sim_eta_good(n = 80, seed = 1)

  expect_snapshot(
    check_eta_bias(data, a, y, c(x1, x2), truncation = 0.05, n_boot = 10),
    error = TRUE
  )
  expect_snapshot(
    check_eta_bias(
      data,
      a,
      y,
      c(x1, x2),
      truncation = c(-0.1, 0.95),
      n_boot = 10
    ),
    error = TRUE
  )
  expect_snapshot(
    check_eta_bias(
      data,
      a,
      y,
      c(x1, x2),
      truncation = c(0.6, 0.4),
      n_boot = 10
    ),
    error = TRUE
  )
  expect_snapshot(
    check_eta_bias(
      data,
      a,
      y,
      c(x1, x2),
      truncation_grid = c(0.1, 0.6),
      n_boot = 10
    ),
    error = TRUE
  )
  expect_snapshot(
    check_eta_bias(
      data,
      a,
      y,
      c(x1, x2),
      truncation_grid = c(0.1, NA),
      n_boot = 10
    ),
    error = TRUE
  )
  expect_snapshot(
    check_eta_bias(data, a, y, c(x1, x2), n_boot = 2.5),
    error = TRUE
  )
})
