# ETA.Bias is a parametric bootstrap: the exact numbers depend on R's glm fits
# and RNG, so every magnitude claim is written as a direction, ordering, or a
# multiple of the run's own Monte Carlo error rather than a hard-coded value.
# Each fixed-seed dataset yields one realization; thresholds sit well inside
# the across-sample spread so a single draw stays on the correct side of every
# inequality. The tolerances below are calibrated at n_boot = 200. Tests that
# read only a sign or an ordering drop to 100; a few raise it where a tighter
# Monte Carlo error is what the assertion rests on. The bootstrap treatment is
# drawn from the untruncated fitted propensity, so truncation changes only the
# estimator; the truncation sweep tests read that as increasing bias and
# shrinking bootstrap spread.

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
# fitted scores toward 0 and 1. `steepness` dials the severity; the default of
# 2 is moderate, well short of the 3.5 the severity test reaches for, and the
# outcome model is linear with effect tau, so `truth` is near tau.
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
    class = "positively_args_error"
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
    class = "positively_args_error"
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
    class = "positively_args_error"
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

test_that("check_eta_bias() requires at least two bootstrap draws", {
  local_quiet()
  data <- sim_eta_good(n = 80, seed = 1)
  # A single bootstrap draw leaves the Monte Carlo standard error undefined, so
  # n_boot must be at least two.
  expect_error(
    check_eta_bias(data, a, y, x1, n_boot = 1),
    class = "positively_range_error"
  )
})

test_that("the n_boot floor message is stable", {
  local_quiet()
  data <- sim_eta_good(n = 80, seed = 1)
  expect_snapshot_abort(check_eta_bias(data, a, y, x1, n_boot = 1))
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
  expect_snapshot_abort(check_eta_bias(data, a, y, c(x1, bad), n_boot = 10))
})

test_that("check_eta_bias() handles a non-syntactic covariate name", {
  local_quiet()
  data <- withr::with_seed(1, {
    d <- data.frame(
      a = stats::rbinom(200, 1L, 0.5),
      y = stats::rnorm(200),
      x = stats::rnorm(200)
    )
    names(d)[3] <- "my var"
    d
  })
  res <- withr::with_seed(
    2024,
    check_eta_bias(data, a, y, `my var`, n_boot = 10)
  )
  expect_identical(S7::S7_class(res)@name, "eta_bias_result")
  expect_true(all(is.finite(res@results$bias)))
  expect_true(all(is.finite(res@results$mc_se)))
  expect_true(is.finite(res@truth))
})

test_that("check_eta_bias() handles non-syntactic exposure and outcome names", {
  local_quiet()
  data <- withr::with_seed(1, {
    d <- data.frame(
      a = stats::rbinom(200, 1L, 0.5),
      y = stats::rnorm(200),
      x = stats::rnorm(200)
    )
    names(d) <- c("my a", "my y", "x")
    d
  })
  res <- withr::with_seed(
    2024,
    check_eta_bias(data, `my a`, `my y`, x, n_boot = 10)
  )
  expect_identical(S7::S7_class(res)@name, "eta_bias_result")
  expect_true(all(is.finite(res@results$bias)))
  expect_true(all(is.finite(res@results$mc_se)))
  expect_true(is.finite(res@truth))
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

test_that("the matrix and formula paths agree tightly on an aliased design", {
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

  # The constant covariate x3 is aliased against the intercept. The observed
  # G-computation estimate does not depend on how the aliased column is handled,
  # so truth already matches to machine precision. The bias does depend on it:
  # it inherits the residual standard deviation, which divides the residual sum
  # of squares by the residual degrees of freedom. The formula path uses
  # n - rank; the matrix path must match it rather than counting the aliased
  # column, or the two sigmas and therefore the two biases disagree.
  expect_equal(matrix_path@truth, formula_path@truth, tolerance = 1e-8)
  expect_equal(bias_of(matrix_path), bias_of(formula_path), tolerance = 1e-8)
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

test_that("the rare-level gcomp bootstrap pins a reference bias, mc_se, and truth", {
  local_quiet()
  # A numerical pin on the aliased-design refit path: the rare second level can
  # drop out of a resample, leaving a rank-deficient design whose aliased column
  # is zeroed. The reference values come from a live run under this exact seed
  # and must reproduce to machine precision.
  data <- sim_eta_rare_level(seed = 2)
  res <- withr::with_seed(
    9,
    check_eta_bias(
      data,
      a,
      y,
      c(x1, x2, region),
      estimator = "gcomp",
      n_boot = 50
    )
  )
  expect_equal(res@results$bias[[1]], 0.014950150587088373, tolerance = 1e-12)
  expect_equal(res@results$mc_se[[1]], 0.026229325116211968, tolerance = 1e-12)
  expect_equal(res@truth, 0.91311162635432241, tolerance = 1e-12)
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

# ---- Declared exposure type ------------------------------------------------

test_that("a declared binary type reproduces the auto path exactly", {
  local_quiet()
  data <- sim_eta_good(n = 300, seed = 1)
  auto <- fit_eta(data, "ipw", n_boot = 50)
  declared <- fit_eta(data, "ipw", exposure_type = "binary", n_boot = 50)

  # Both calls run under the same seed, and on a genuinely binary exposure the
  # declaration can only agree with detection, so the bootstrap draws and every
  # summary built from them must match, not merely the fact that both returned.
  expect_identical(declared@exposure_type, auto@exposure_type)
  expect_identical(declared@results, auto@results)
  expect_identical(declared@truth, auto@truth)
  expect_identical(declared@boot_estimates, auto@boot_estimates)
})

test_that("a declared binary type skips the detection alert", {
  data <- sim_eta_good(n = 300, seed = 1)
  # The auto path announces what it inferred. A declared type is taken as given,
  # so there is nothing to announce.
  expect_message(
    fit_eta(data, "gcomp", n_boot = 10),
    "Treating `.exposure` as binary"
  )
  expect_no_message(
    fit_eta(data, "gcomp", exposure_type = "binary", n_boot = 10)
  )
})

test_that("a declared binary type aborts on a three-level exposure", {
  local_quiet()
  data <- sim_eta_good(n = 120, seed = 1)
  data$a <- factor(rep(c("a", "b", "c"), length.out = nrow(data)))

  err <- expect_error(
    check_eta_bias(
      data,
      a,
      y,
      c(x1, x2),
      exposure_type = "binary",
      n_boot = 10
    ),
    class = "positively_exposure_type_error"
  )
  expect_match(conditionMessage(err), "two distinct values", fixed = TRUE)
})

test_that("a constant exposure aborts on both the auto and declared paths", {
  local_quiet()
  data <- sim_eta_good(n = 120, seed = 1)
  data$a <- factor(rep("a", nrow(data)))

  # A one-level factor is not two distinct values, yet detection reads it as
  # binary, so detection alone can never reject it. binary_to_01() then maps the
  # single level to 0, leaving a treatment mechanism with no treated arm and
  # estimates that are degenerate rather than informative. Only the structural
  # requirement of the resolved type rules this out, and it rules it out whether
  # the type was declared or inferred.
  expect_identical(detect_exposure_type(data$a), "binary")

  expect_error(
    check_eta_bias(data, a, y, c(x1, x2), n_boot = 10),
    class = "positively_exposure_type_error"
  )
  expect_error(
    check_eta_bias(
      data,
      a,
      y,
      c(x1, x2),
      exposure_type = "binary",
      n_boot = 10
    ),
    class = "positively_exposure_type_error"
  )
})

test_that("check_eta_bias() rejects a type outside its supported menu", {
  local_quiet()
  data <- sim_eta_good(n = 80, seed = 1)
  err <- expect_error(
    check_eta_bias(
      data,
      a,
      y,
      c(x1, x2),
      exposure_type = "continuous",
      n_boot = 10
    ),
    class = "positively_args_error"
  )
  expect_match(
    conditionMessage(err),
    '`exposure_type` must be one of "auto" or "binary", not "continuous".',
    fixed = TRUE
  )
})

# ---- Degenerate bootstrap draws -------------------------------------------

test_that("non-finite bootstrap draws are dropped with a classed warning", {
  local_quiet()
  # A severe practical violation at small n drives some inverse-probability
  # refits to a zero weight sum, so a handful of bootstrap draws are non-finite.
  data <- withr::with_seed(4, {
    n <- 25
    x1 <- stats::rnorm(n)
    a <- stats::rbinom(n, 1L, stats::plogis(2 + 2 * x1))
    y <- a + x1 + stats::rnorm(n)
    tibble::tibble(a = a, y = y, x1 = x1)
  })

  expect_warning(
    res <- withr::with_seed(
      1,
      check_eta_bias(data, a, y, x1, estimator = "ipw", n_boot = 50)
    ),
    class = "positively_degenerate_boot_warning"
  )

  # The summaries recover to finite values once the non-finite draws are dropped.
  expect_true(is.finite(res@results$bias[[1]]))
  expect_true(is.finite(res@results$mc_se[[1]]))
  expect_true(is.finite(res@results$boot_mean[[1]]))
  # Two of the fifty draws are non-finite, so forty-eight are retained.
  expect_length(res@boot_estimates[[1]], 48L)
})

test_that("the degenerate bootstrap warning wording is stable", {
  local_quiet()
  # The same severe practical violation the class-only test above uses. This
  # snapshot pins the exact warning wording so that the count of dropped draws
  # is reported in terms that stay meaningful across truncation levels, rather
  # than as the single per-level maximum.
  data <- withr::with_seed(4, {
    n <- 25
    x1 <- stats::rnorm(n)
    a <- stats::rbinom(n, 1L, stats::plogis(2 + 2 * x1))
    y <- a + x1 + stats::rnorm(n)
    tibble::tibble(a = a, y = y, x1 = x1)
  })
  withr::local_options(warn = 0)
  expect_snapshot(
    res <- withr::with_seed(
      1,
      check_eta_bias(data, a, y, x1, estimator = "ipw", n_boot = 50)
    )
  )
})

test_that("a healthy run keeps every bootstrap draw and warns for none", {
  local_quiet()
  data <- sim_eta_good(n = 300, seed = 1)
  res <- expect_no_condition(
    fit_eta(data, "ipw", n_boot = 50),
    class = "positively_degenerate_boot_warning"
  )
  expect_length(res@boot_estimates[[1]], 50L)
  expect_true(is.finite(bias_of(res)))
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
  expect_false("low_support" %in% names(res@results))
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

# ---- Statistical behavior -------------------------------------------------

test_that("good overlap leaves every estimator near zero", {
  local_quiet()
  data <- sim_eta_good(n = 1000, seed = 1)
  for (estimator in c("ipw", "gcomp", "aipw")) {
    res <- fit_eta(data, estimator, n_boot = 200)
    expect_lt(abs(bias_of(res)), 0.05)
  }
})

test_that("under violation gcomp stays near zero and ipw is large", {
  local_quiet()
  data <- sim_eta_violation(n = 1000, seed = 1)
  ipw <- fit_eta(data, "ipw", n_boot = 500)
  gcomp <- fit_eta(data, "gcomp", n_boot = 200)
  aipw <- fit_eta(data, "aipw", n_boot = 200)

  # gcomp is the most robust invariant, since truth is its own point estimate;
  # its residual bias is Monte Carlo noise.
  expect_lt(abs(bias_of(gcomp)), 3 * mc_se_of(gcomp))
  expect_lt(abs(bias_of(gcomp)), 0.03)

  # ipw carries substantial, clearly non-null bias. The multiple is 3 rather
  # than 5, and this fit uses n_boot = 500: for this DGP the expected
  # bias-to-mc_se ratio at n_boot = 200 is about 4, so a 5-multiple holds only
  # for lucky seeds. At n_boot = 500 the expected ratio is about 6.3 with an
  # across-seed spread near 1, leaving a wide margin for RNG and refactor
  # changes.
  expect_gt(bias_of(ipw), 0.05)
  expect_gt(bias_of(ipw), 3 * mc_se_of(ipw))

  # ipw dominates both consistent estimators by a wide margin.
  expect_gt(abs(bias_of(ipw)), 3 * abs(bias_of(aipw)))
  expect_gt(abs(bias_of(ipw)), 3 * abs(bias_of(gcomp)))

  # The doubly robust estimator stays near zero, though looser than gcomp
  # because of its larger across-sample spread.
  expect_lt(abs(bias_of(aipw)), 0.05)
})

test_that("truncation tightens ipw bias upward and its spread downward", {
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

test_that("ipw bias grows with violation severity", {
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

test_that("Monte Carlo error scales as one over root n_boot", {
  local_quiet()
  data <- sim_eta_violation(n = 1000, seed = 1)
  mc_100 <- mc_se_of(fit_eta(data, "ipw", n_boot = 100))
  mc_400 <- mc_se_of(fit_eta(data, "ipw", n_boot = 400))

  expect_lt(mc_400, mc_100)
  ratio <- mc_400 / mc_100
  # Expected ratio is sqrt(100 / 400) = 0.5, allowed a 25 percent band.
  expect_gt(ratio, 0.375)
  expect_lt(ratio, 0.625)
})

test_that("binary outcomes run and preserve the estimator ordering", {
  local_quiet()
  data <- sim_eta_binary_outcome(n = 1000, seed = 1)
  ipw <- fit_eta(data, "ipw", outcome_type = "binary", n_boot = 200)
  gcomp <- fit_eta(data, "gcomp", outcome_type = "binary", n_boot = 200)
  aipw <- fit_eta(data, "aipw", outcome_type = "binary", n_boot = 200)

  expect_gt(abs(bias_of(ipw)), abs(bias_of(aipw)))
  expect_gt(abs(bias_of(aipw)), abs(bias_of(gcomp)))
  expect_lt(abs(bias_of(gcomp)), 0.02)
})

test_that("truth is a scalar equal to the gcomp estimate near tau", {
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

test_that("normal and empirical error_dist models agree under violation", {
  local_quiet()
  data <- sim_eta_violation(n = 1000, seed = 1)
  normal <- fit_eta(data, "ipw", error_dist = "normal", n_boot = 200)
  empirical <- fit_eta(data, "ipw", error_dist = "empirical", n_boot = 200)

  expect_gt(bias_of(normal), 0.05)
  expect_gt(bias_of(empirical), 0.05)
  expect_lt(abs(bias_of(normal) - bias_of(empirical)), 0.15)
})

# ---- tidy() and glance() --------------------------------------------------

test_that("tidy() returns the results tibble", {
  local_quiet()
  data <- sim_eta_good(n = 300, seed = 1)
  res <- fit_eta(data, "ipw", n_boot = 50)

  expect_identical(generics::tidy(res), res@results)
})

test_that("glance() reports the bias and its Monte Carlo error at one level", {
  local_quiet()
  data <- sim_eta_good(n = 300, seed = 1)
  res <- fit_eta(data, "ipw", n_boot = 50)
  glanced <- generics::glance(res)

  expect_s3_class(glanced, "tbl_df")
  expect_identical(nrow(glanced), 1L)
  expect_setequal(
    names(glanced),
    c("n", "estimator", "n_boot", "truth", "bias", "mc_se")
  )

  expect_identical(glanced$n, 300L)
  expect_identical(glanced$estimator, "ipw")
  expect_identical(glanced$n_boot, 50L)
  expect_identical(glanced$truth, res@truth)

  # ETA.Bias is the mean bootstrap estimate less the truth and its Monte Carlo
  # error is the standard error of those same draws, so both follow from the
  # retained draws on @boot_estimates without going through @results.
  draws <- res@boot_estimates[[1]]
  expect_equal(glanced$bias, mean(draws) - res@truth)
  expect_equal(glanced$mc_se, stats::sd(draws) / sqrt(length(draws)))
})

test_that("glance() reports the bias range across a truncation sweep", {
  local_quiet()
  data <- sim_eta_good(n = 300, seed = 1)
  res <- fit_eta(data, "ipw", truncation_grid = c(0, 0.05, 0.1), n_boot = 50)
  glanced <- generics::glance(res)

  expect_identical(nrow(glanced), 1L)
  expect_setequal(
    names(glanced),
    c("n", "estimator", "n_boot", "truth", "n_levels", "bias_min", "bias_max")
  )

  # A sweep has one bias per truncation level, so there is no single bias to
  # report and no single Monte Carlo error beside it.
  expect_false("bias" %in% names(glanced))
  expect_false("mc_se" %in% names(glanced))

  expect_identical(glanced$n_levels, 3L)
  expect_identical(glanced$truth, res@truth)

  per_level <- vapply(res@boot_estimates, mean, numeric(1)) - res@truth
  expect_equal(glanced$bias_min, min(per_level))
  expect_equal(glanced$bias_max, max(per_level))
})

test_that("summary() reports the bias and its Monte Carlo error", {
  local_quiet()
  data <- sim_eta_good(n = 300, seed = 1)
  res <- fit_eta(data, "ipw", n_boot = 50)
  summarized <- summary(res)

  # The estimator name is a character statistic and the bootstrap draw count is
  # a setting, so neither takes a row. The truth is a quantity the run computed
  # rather than a setting the caller chose, so it does.
  expect_identical(summarized$statistic, c("truth", "bias", "mc_se"))
  expect_type(summarized$value, "double")
  expect_identical(summarized$value, c(res@truth, bias_of(res), mc_se_of(res)))

  # ETA.Bias is read against no cut. Its own reference point, the truth, is a
  # row of its own rather than a threshold beside the bias, because the bias is
  # already the signed distance between the two.
  expect_true(all(is.na(summarized$threshold)))
})

test_that("summary() follows the sweep's glance() shape", {
  local_quiet()
  data <- sim_eta_good(n = 300, seed = 1)
  res <- fit_eta(data, "ipw", truncation_grid = c(0, 0.05, 0.1), n_boot = 50)
  summarized <- summary(res)

  # A sweep reports a range in place of a single bias, so summary() states
  # whatever this run's glance() computed rather than a fixed statistic set.
  expect_identical(
    summarized$statistic,
    c("truth", "n_levels", "bias_min", "bias_max")
  )
  expect_identical(
    summarized$value,
    c(res@truth, 3, min(res@results$bias), max(res@results$bias))
  )
  expect_true(all(is.na(summarized$threshold)))
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

test_that("autoplot() rejects an unknown type as a classed error", {
  local_quiet()
  data <- sim_eta_violation(n = 500, seed = 1)
  res <- fit_eta(data, "ipw", n_boot = 10)

  # A view name is chosen from a fixed menu like any other argument, and the
  # call rendering a figure does not make its failure a lesser one.
  expect_error(
    ggplot2::autoplot(res, type = "bogus"),
    class = "positively_args_error"
  )
})

test_that("ETA bias autoplot views render as expected", {
  local_quiet()
  announce_doppelganger(
    "ETA bias bootstrap distribution",
    "ETA bias truncation sweep"
  )
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

test_that("plot() draws the view and returns the result invisibly", {
  local_quiet()
  local_null_device()
  data <- sim_eta_violation(n = 500, seed = 1)
  swept <- fit_eta(
    data,
    "ipw",
    truncation_grid = c(0, 0.025, 0.05, 0.1),
    n_boot = 100
  )
  # The type argument reaches the plot method through the dots.
  expect_identical(plot(swept, type = "sweep"), swept)
})

test_that("the sweep view aborts on a single-run result", {
  local_quiet()
  data <- sim_eta_violation(n = 500, seed = 1)
  single <- fit_eta(data, "ipw", n_boot = 50)
  expect_snapshot_abort(
    ggplot2::autoplot(single, type = "sweep"),
    class = "positively_sweep_absent_error"
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

  expect_snapshot_abort(check_eta_bias(1:10, a, y, c(x1, x2), n_boot = 10))
  expect_snapshot_abort(check_eta_bias(
    continuous,
    a,
    y,
    c(x1, x2),
    n_boot = 10
  ))
  expect_snapshot_abort(check_eta_bias(
    categorical,
    a,
    y,
    c(x1, x2),
    n_boot = 10
  ))
  expect_snapshot_abort(check_eta_bias(
    good,
    a,
    y,
    tidyselect::starts_with("zzz"),
    n_boot = 10
  ))
  expect_snapshot_abort(check_eta_bias(
    good,
    a,
    not_a_column,
    c(x1, x2),
    n_boot = 10
  ))
  expect_snapshot_abort(check_eta_bias(
    na_exposure,
    a,
    y,
    c(x1, x2),
    n_boot = 10
  ))
  expect_snapshot_abort(check_eta_bias(
    good,
    a,
    y,
    c(x1, x2),
    outcome_type = "binary",
    n_boot = 10
  ))
  expect_snapshot_abort(check_eta_bias(good, a, y, c(a, x1, x2), n_boot = 10))
})

test_that("the truncation and n_boot validation messages are stable", {
  local_quiet()
  data <- sim_eta_good(n = 80, seed = 1)

  expect_snapshot_abort(check_eta_bias(
    data,
    a,
    y,
    c(x1, x2),
    truncation = 0.05,
    n_boot = 10
  ))
  expect_snapshot_abort(check_eta_bias(
    data,
    a,
    y,
    c(x1, x2),
    truncation = c(-0.1, 0.95),
    n_boot = 10
  ))
  expect_snapshot_abort(check_eta_bias(
    data,
    a,
    y,
    c(x1, x2),
    truncation = c(0.6, 0.4),
    n_boot = 10
  ))
  expect_snapshot_abort(check_eta_bias(
    data,
    a,
    y,
    c(x1, x2),
    truncation_grid = c(0.1, 0.6),
    n_boot = 10
  ))
  expect_snapshot_abort(check_eta_bias(
    data,
    a,
    y,
    c(x1, x2),
    truncation_grid = c(0.1, NA),
    n_boot = 10
  ))
  expect_snapshot_abort(check_eta_bias(data, a, y, c(x1, x2), n_boot = 2.5))
})

# ---- Display methods -------------------------------------------------------

# A single IPW run rather than a truncation sweep, so the reading is one bias
# with its Monte Carlo error rather than a range over levels.
eta_bias_display_result <- function() {
  fit_eta(sim_eta_good(n = 300, seed = 1), "ipw", n_boot = 50)
}

test_that("diagnostic_label() names the ETA bias check, not its class", {
  local_quiet()
  res <- eta_bias_display_result()
  label <- expect_readable_label(res)

  printed <- printed_text(res)
  expect_no_match(printed, S7::S7_class(res)@name, fixed = TRUE)
  expect_match(printed, label, fixed = TRUE)
})

test_that("diagnostic_headline() names the estimator the bias belongs to", {
  local_quiet()
  res <- eta_bias_display_result()
  headline <- expect_readable_headline(res)
  text <- rendered_text(headline)

  # ETA bias is a property of an estimator under a truncation rule, not of the
  # data alone, so a reading that omitted the estimator would not identify what
  # was measured.
  expect_match(text, "ipw", fixed = TRUE)
})

test_that("the ETA bias label and headline are stable", {
  local_quiet()
  res <- eta_bias_display_result()
  expect_snapshot({
    diagnostic_label(res)
    writeLines(diagnostic_headline(res))
  })
})

test_that("a fitted ETA bias run has nothing to report", {
  local_quiet()
  # The bias is a magnitude rather than a reading against a cut the caller set,
  # so the class declares no findings and a fitted run reports none.
  data <- sim_eta_violation(n = 400, seed = 1)
  res <- fit_eta(data, "ipw", truncation_grid = c(0, 0.05, 0.1), n_boot = 50)

  expect_gt(max(abs(res@results$bias)), 0)
  expect_identical(nrow(sniff_violations(res)), 0L)
})
