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

# A three-level exposure drawn from a multinomial logit, with the planted
# violation on the reference level: level "a" carries the flat linear predictor,
# so wherever x1 or x2 is large the other two levels take nearly all the mass
# and "a" becomes rare. That leaves both contrasts against "a" exposed to the
# same sparsity. `steepness` dials the severity; the default of 2 is moderate,
# well short of the 4 the statistical block reaches for. The outcome is linear
# with a per-level effect, so the two truths land near tau_b and tau_c.
sim_eta_categorical <- function(
  n = 600,
  steepness = 2,
  tau_b = 1,
  tau_c = 2,
  seed = 1
) {
  withr::local_seed(seed)
  x1 <- stats::rnorm(n)
  x2 <- stats::rnorm(n)
  odds <- exp(cbind(0, steepness * x1, steepness * x2))
  probs <- odds / rowSums(odds)
  # An inverse-CDF draw, one uniform per row, as the bootstrap exposure draw
  # makes past two levels.
  drawn <- rowSums(matrix(stats::runif(n), n, 3) > t(apply(probs, 1, cumsum)))
  a <- factor(c("a", "b", "c")[pmin(drawn, 2L) + 1L], levels = c("a", "b", "c"))
  y <- tau_b * (a == "b") + tau_c * (a == "c") + x1 + x2 + stats::rnorm(n)
  tibble::tibble(a = a, y = y, x1 = x1, x2 = x2)
}

# A continuous exposure whose treatment mechanism is linear in the two
# covariates with residual standard deviation `conditional_sd`. That is the one
# dial: the marginal exposure variance is `2 + conditional_sd^2`, so the
# stabilized density ratio f(A) / f(A | W) has finite variance only while
# `conditional_sd` exceeds sqrt(2). A wide setting leaves the ratio near one
# everywhere; a narrow one makes the mechanism near-deterministic and hands a
# few observations nearly all the weight, which is the continuous reading of a
# practical positivity violation. The outcome is linear with effect tau, so the
# counterfactual surface is exactly linear in the exposure and the default
# working model's truth is tau.
sim_eta_continuous <- function(
  n = 600,
  conditional_sd = 1,
  tau = 1,
  seed = 1
) {
  withr::local_seed(seed)
  x1 <- stats::rnorm(n)
  x2 <- stats::rnorm(n)
  a <- x1 + x2 + stats::rnorm(n, sd = conditional_sd)
  y <- tau * a + x1 + x2 + stats::rnorm(n)
  tibble::tibble(a = a, y = y, x1 = x1, x2 = x2)
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
  expect_equal(unname(res@truth), 0.91311162635432241, tolerance = 1e-12)
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
  expect_equal(unname(res@truth), 1.0636370442788932)
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

test_that("check_eta_bias() announces the detected type and no declared one", {
  withr::local_options(positively.quiet = FALSE)
  data <- sim_eta_good(n = 300, seed = 1)

  # Two supported types make the reading worth reporting: the announcement names
  # which of the two branches detection sent the call down. A declared type is
  # the caller's own statement read back, so it is not announced.
  expect_message(fit_eta(data, "gcomp", n_boot = 10), "as binary")
  expect_no_message(
    fit_eta(data, "gcomp", exposure_type = "binary", n_boot = 10)
  )
  expect_no_message(
    fit_eta(data, "gcomp", exposure_type = "categorical", n_boot = 10)
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
  # binary, so detection alone can never reject it. The factor coding in
  # eta_spec() maps the single level to 0, leaving a treatment mechanism with no
  # treated arm and estimates that are degenerate rather than informative. Only
  # the structural requirement of the resolved type rules this out, and it rules
  # it out whether the type was declared or inferred. A categorical exposure
  # requires nothing structurally, so it is the level count that rejects it
  # there.
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
  expect_error(
    check_eta_bias(
      data,
      a,
      y,
      c(x1, x2),
      exposure_type = "categorical",
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
      exposure_type = "ordinal",
      n_boot = 10
    ),
    class = "positively_args_error"
  )
  expect_match(
    conditionMessage(err),
    '`exposure_type` must be one of "auto", "binary", "categorical", or "continuous", not "ordinal".',
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
  # The order is pinned, not the set: the dimension column leads whatever the
  # exposure type, so a binary run reads the same way a categorical one does.
  expect_identical(
    names(res@results),
    c(
      "term",
      "truncation_lower",
      "truncation_upper",
      "bias",
      "mc_se",
      "boot_mean"
    )
  )
  expect_false("low_support" %in% names(res@results))
})

test_that("the estimand term names the contrast and keys the truth", {
  local_quiet()
  data <- sim_eta_good(n = 300, seed = 1)
  res <- fit_eta(data, "ipw", n_boot = 50)

  # A binary exposure has one contrast, so the term names both levels and the
  # truth carries that one name.
  expect_identical(res@results$term, "1 - 0")
  expect_identical(names(res@truth), "1 - 0")
  expect_length(res@truth, 1L)
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

# ---- Probability truncation -----------------------------------------------

test_that("truncation past two levels clamps below and renormalizes", {
  gmat <- rbind(
    c(0.05, 0.15, 0.80),
    c(0.02, 0.03, 0.95)
  )
  # A missing upper bound is what a run past two levels supplies, so a finite
  # result is itself the proof that the upper bound goes unused.
  bounded <- truncate_probabilities(gmat, 0.1, NA_real_)

  first <- c(0.1, 0.15, 0.80)
  second <- c(0.1, 0.1, 0.95)
  expect_equal(
    bounded,
    rbind(first / sum(first), second / sum(second))
  )
  expect_type(bounded, "double")
  expect_equal(rowSums(bounded), c(1, 1))
})

test_that("a row already above the bound survives truncation unchanged", {
  gmat <- rbind(
    c(0.2, 0.3, 0.5),
    c(0.1, 0.4, 0.5)
  )
  # Every entry is at or above the bound, so the clamp is inert and the row
  # already sums to one, leaving the renormalization inert too.
  expect_equal(truncate_probabilities(gmat, 0.1, NA_real_), gmat)
})

test_that("the truncation grid ceiling is 1 over the level count", {
  # Two levels give the familiar `[0, 0.5)`; three give `[0, 1/3)`. The
  # acceptance set is exact at the ceiling, so a bound a hair below it is
  # admitted while the ceiling itself is not.
  expect_error(
    validate_truncation_grid(0.5, k = 2),
    class = "positively_range_error"
  )
  expect_identical(validate_truncation_grid(c(0, 0.4999), k = 2), c(0, 0.4999))

  expect_error(
    validate_truncation_grid(1 / 3, k = 3),
    class = "positively_range_error"
  )
  expect_identical(validate_truncation_grid(c(0, 0.333), k = 3), c(0, 0.333))
})

test_that("the continuous truncation grid ceiling is one", {
  # A continuous exposure has no simplex and no levels, so a grid point is the
  # quantile level at which the stabilized weight is capped and the ceiling is
  # one rather than 1/k. A level of 0.999 clears every discrete ceiling, so
  # admitting it holds the ceiling to the exposure type rather than to `k`. The
  # acceptance set is exact at the ceiling, so the ceiling itself is refused: a
  # level of one would cap every weight at the smallest one observed.
  expect_identical(
    validate_truncation_grid(0.999, exposure_type = "continuous"),
    0.999
  )
  expect_error(
    validate_truncation_grid(1, exposure_type = "continuous"),
    class = "positively_range_error"
  )
  expect_error(
    validate_truncation_grid(1.5, exposure_type = "continuous"),
    class = "positively_range_error"
  )
  # The floor is shared with the discrete branch, since a negative number names
  # no quantile.
  expect_error(
    validate_truncation_grid(-0.1, exposure_type = "continuous"),
    class = "positively_range_error"
  )
})

test_that("the continuous truncation grid messages are stable", {
  local_quiet()
  data <- sim_eta_continuous(n = 100, seed = 1)

  # Both messages call a grid point a quantile level, where the discrete
  # wording calls it a lower bound, so a message that drifts onto the other
  # branch's phrasing shows up here.
  expect_snapshot_abort(
    check_eta_bias(
      data,
      a,
      y,
      c(x1, x2),
      truncation_grid = c(0.1, 1.5),
      n_boot = 10
    ),
    class = "positively_range_error"
  )
  expect_snapshot_abort(
    check_eta_bias(
      data,
      a,
      y,
      c(x1, x2),
      truncation_grid = numeric(0),
      n_boot = 10
    ),
    class = "positively_empty_error"
  )
})

# ---- Oracles against propensity -------------------------------------------
# The bootstrap loop cannot afford a call into propensity per draw, so it
# spells the truncation rule and the weight formula out itself. These blocks
# hold both spellings to the package they mirror.

test_that("simplex truncation matches propensity::ps_trunc()", {
  skip_if_not_installed("propensity")
  gmat <- rbind(
    c(0.02, 0.18, 0.80),
    c(0.30, 0.35, 0.35),
    c(0.05, 0.05, 0.90)
  )
  colnames(gmat) <- c("a", "b", "c")
  exposure <- factor(colnames(gmat), levels = colnames(gmat))

  oracle <- propensity::ps_trunc(
    gmat,
    .exposure = exposure,
    method = "ps",
    lower = 0.1
  )
  # ps_trunc() records no upper bound for a matrix, which is the property the
  # categorical results carry through as a missing `truncation_upper`.
  expect_identical(attr(oracle, "ps_trunc_meta")$upper_bound, NA_real_)

  expected <- matrix(
    as.numeric(oracle),
    nrow = nrow(gmat),
    ncol = ncol(gmat),
    dimnames = dimnames(gmat)
  )
  expect_equal(truncate_probabilities(gmat, 0.1, NA_real_), expected)
})

test_that("the categorical estimator weights match propensity::wt_ate()", {
  skip_if_not_installed("propensity")
  withr::local_seed(7)
  n <- 40
  level_labels <- c("a", "b", "c")
  gmat <- matrix(stats::runif(n * 3), n, 3)
  gmat <- gmat / rowSums(gmat)
  colnames(gmat) <- level_labels
  exposure <- factor(
    sample(level_labels, n, replace = TRUE),
    levels = level_labels
  )
  y <- stats::rnorm(n)
  # G-computation plays no part in the Hajek estimate, so the outcome
  # predictions are free to be anything.
  qmat <- matrix(0, n, 3)

  weights <- as.numeric(propensity::wt_ate(
    gmat,
    exposure,
    exposure_type = "categorical"
  ))
  spec <- eta_spec(exposure, "categorical")

  # The per-level weight is 1{A = k} / g_k, which is the ATE weight on the rows
  # holding level k and zero elsewhere.
  for (level in seq_along(level_labels)) {
    internal <- (spec$a == spec$grid[[level]]) / gmat[, level]
    expect_equal(
      internal,
      ifelse(exposure == level_labels[[level]], weights, 0)
    )
  }

  # The contrast the estimator forms from those weights is the same difference
  # of weighted means built from wt_ate() directly.
  hajek <- vapply(
    level_labels,
    function(label) {
      keep <- exposure == label
      sum(weights[keep] * y[keep]) / sum(weights[keep])
    },
    numeric(1)
  )
  expect_equal(
    eta_estimate("ipw", spec$a, y, gmat, qmat, spec),
    c(
      "b - a" = hajek[["b"]] - hajek[["a"]],
      "c - a" = hajek[["c"]] - hajek[["a"]]
    )
  )
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
  expect_equal(unname(gcomp@truth), 1, tolerance = 0.15)
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
  expect_identical(glanced$truth, unname(res@truth))

  # ETA.Bias is the mean bootstrap estimate less the truth and its Monte Carlo
  # error is the standard error of those same draws, so both follow from the
  # retained draws on @boot_estimates without going through @results.
  draws <- res@boot_estimates[[1]]
  expect_equal(glanced$bias, mean(draws) - unname(res@truth))
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
  expect_identical(glanced$truth, unname(res@truth))

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
  expect_identical(
    summarized$value,
    c(unname(res@truth), bias_of(res), mc_se_of(res))
  )

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
    c(unname(res@truth), 3, min(res@results$bias), max(res@results$bias))
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

  # Every exposure type the diagnostic can detect is one it runs, so the
  # exposure-shaped abort left to pin is the one level that no type can make an
  # estimand out of.
  one_level <- good
  one_level$a <- factor(rep("a", nrow(one_level)))

  na_exposure <- good
  na_exposure$a[1] <- NA

  expect_snapshot_abort(check_eta_bias(1:10, a, y, c(x1, x2), n_boot = 10))
  expect_snapshot_abort(check_eta_bias(
    one_level,
    a,
    y,
    c(x1, x2),
    exposure_type = "categorical",
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

# ---- Categorical exposures -------------------------------------------------
# An exposure of more than two levels reports one contrast per non-reference
# level, so a run holds k - 1 terms at every truncation level. Two levels is the
# reduction case: the machinery reads the level count rather than the declared
# label, so declaring "categorical" on a two-level exposure returns the binary
# run unchanged, down to the bootstrap draws.

test_that("a declared categorical type runs a three-level exposure", {
  local_quiet()
  skip_if_not_installed("nnet")
  data <- sim_eta_categorical(n = 300, seed = 1)
  res <- fit_eta(data, "gcomp", exposure_type = "categorical", n_boot = 30)

  expect_identical(S7::S7_class(res)@name, "eta_bias_result")
  expect_identical(res@exposure_type, "categorical")
  expect_identical(res@n, 300L)
})

test_that("a three-level exposure is detected and announced as categorical", {
  withr::local_options(positively.quiet = FALSE)
  skip_if_not_installed("nnet")
  data <- sim_eta_categorical(n = 300, seed = 1)

  # Two supported types make the reading worth reporting, since the
  # announcement now names which of the two branches the call took.
  expect_message(
    res <- fit_eta(data, "gcomp", n_boot = 30),
    "categorical"
  )
  expect_identical(res@exposure_type, "categorical")
})

test_that("a three-level exposure without nnet aborts", {
  local_quiet()
  data <- sim_eta_categorical(n = 200, seed = 1)

  # Mocked rather than skipped, so the reader without nnet sees a refusal this
  # suite has read, rather than one that fires only on machines it never runs on.
  testthat::local_mocked_bindings(
    is_installed = function(...) FALSE,
    .package = "rlang"
  )
  expect_snapshot_abort(
    check_eta_bias(
      data,
      a,
      y,
      c(x1, x2),
      exposure_type = "categorical",
      n_boot = 10
    ),
    class = "positively_missing_package_error"
  )
})

test_that("a two-level exposure declared categorical reduces to the binary path", {
  local_quiet()
  data <- sim_eta_good(n = 300, seed = 1)
  binary <- fit_eta(data, "ipw", exposure_type = "binary", n_boot = 50)
  categorical <- fit_eta(
    data,
    "ipw",
    exposure_type = "categorical",
    n_boot = 50
  )

  # The declared label is recorded, and nothing else follows from it: both runs
  # face two levels, so both take the two-level treatment model and the two-level
  # exposure draw, leaving the bootstrap stream and every summary built from it
  # identical rather than merely close.
  expect_identical(categorical@exposure_type, "categorical")
  expect_identical(categorical@results, binary@results)
  expect_identical(categorical@truth, binary@truth)
  expect_identical(categorical@boot_estimates, binary@boot_estimates)
})

test_that("a two-level categorical run keeps the two-sided truncation bound", {
  local_quiet()
  data <- sim_eta_violation(n = 400, seed = 1)
  binary <- fit_eta(
    data,
    "ipw",
    exposure_type = "binary",
    truncation = c(0.05, 0.95),
    n_boot = 50
  )
  categorical <- fit_eta(
    data,
    "ipw",
    exposure_type = "categorical",
    truncation = c(0.05, 0.95),
    n_boot = 50
  )

  # The two-sided bound is a statement about a one-dimensional simplex, so it is
  # the level count that decides whether it means anything. At two levels it
  # does, whichever label the caller declared.
  expect_identical(categorical@results, binary@results)
  expect_equal(categorical@results$truncation_upper, 0.95)
})

test_that("a categorical run yields k - 1 terms per truncation level", {
  local_quiet()
  skip_if_not_installed("nnet")
  data <- sim_eta_categorical(n = 300, seed = 1)
  grid <- c(0, 0.05, 0.1)
  res <- fit_eta(
    data,
    "ipw",
    exposure_type = "categorical",
    truncation_grid = grid,
    n_boot = 30
  )

  expect_identical(nrow(res@results), 6L)
  # Terms are contrasts against the reference level, which defaults to the first.
  expect_identical(unique(res@results$term), c("b - a", "c - a"))
  expect_setequal(res@results$truncation_lower, grid)

  # One truth per term, keyed by term, and one bootstrap-estimate vector per
  # term and truncation level.
  expect_length(res@truth, 2L)
  expect_identical(names(res@truth), c("b - a", "c - a"))
  expect_length(res@boot_estimates, nrow(res@results))
  expect_true(all(lengths(res@boot_estimates) == 30))
})

test_that("term leads the categorical results columns", {
  local_quiet()
  skip_if_not_installed("nnet")
  data <- sim_eta_categorical(n = 300, seed = 1)
  res <- fit_eta(data, "ipw", exposure_type = "categorical", n_boot = 30)

  # The dimension column leads, so the order is pinned rather than the set.
  expect_identical(
    names(res@results),
    c(
      "term",
      "truncation_lower",
      "truncation_upper",
      "bias",
      "mc_se",
      "boot_mean"
    )
  )
})

test_that("every estimator returns finite categorical estimates", {
  local_quiet()
  skip_if_not_installed("nnet")
  data <- sim_eta_categorical(n = 300, seed = 1)
  for (estimator in c("gcomp", "ipw", "aipw")) {
    res <- fit_eta(data, estimator, exposure_type = "categorical", n_boot = 30)
    expect_identical(nrow(res@results), 2L)
    expect_true(all(is.finite(res@results$bias)))
    expect_true(all(is.finite(res@results$mc_se)))
    expect_true(all(is.finite(res@truth)))
  }
})

test_that("the formula path serves a categorical exposure", {
  local_quiet()
  skip_if_not_installed("nnet")
  data <- sim_eta_categorical(n = 300, seed = 1)
  matrix_path <- fit_eta(
    data,
    "ipw",
    exposure_type = "categorical",
    n_boot = 30
  )
  formula_path <- fit_eta(
    data,
    "ipw",
    exposure_type = "categorical",
    exposure_formula = a ~ x1 + x2,
    outcome_formula = y ~ a + x1 + x2,
    n_boot = 30
  )

  # Both paths hand the multinomial treatment model the same level indicators
  # and the outcome model the same saturated coding of the exposure, so they fit
  # the same two models on the same designs and agree exactly rather than
  # closely.
  expect_identical(formula_path@results, matrix_path@results)
  expect_identical(formula_path@truth, matrix_path@truth)

  # A factor covariate takes the formula path with no override to send it there.
  data$region <- factor(rep(c("north", "south"), length.out = nrow(data)))
  with_factor <- withr::with_seed(
    2024,
    check_eta_bias(
      data,
      a,
      y,
      c(x1, x2, region),
      estimator = "aipw",
      exposure_type = "categorical",
      n_boot = 30
    )
  )
  expect_identical(unique(with_factor@results$term), c("b - a", "c - a"))
  expect_true(all(is.finite(with_factor@results$bias)))
})

test_that("reference_level defaults to the first level and can be reset", {
  local_quiet()
  skip_if_not_installed("nnet")
  data <- sim_eta_categorical(n = 300, seed = 1)
  default <- fit_eta(data, "gcomp", exposure_type = "categorical", n_boot = 30)
  reset <- fit_eta(
    data,
    "gcomp",
    exposure_type = "categorical",
    reference_level = "b",
    n_boot = 30
  )

  expect_identical(unique(default@results$term), c("b - a", "c - a"))
  expect_identical(names(default@truth), c("b - a", "c - a"))
  expect_identical(unique(reset@results$term), c("a - b", "c - b"))
  expect_identical(names(reset@truth), c("a - b", "c - b"))

  # The saturated model is the same fit either way, so the two contrast sets
  # describe the same three counterfactual means and compose exactly:
  # (c - b) - (a - b) is c - a.
  expect_equal(
    unname(reset@truth[["c - b"]] - reset@truth[["a - b"]]),
    unname(default@truth[["c - a"]])
  )
})

test_that("reference_level rejects an unknown or non-scalar level", {
  local_quiet()
  data <- sim_eta_categorical(n = 200, seed = 1)

  expect_snapshot_abort(
    check_eta_bias(
      data,
      a,
      y,
      c(x1, x2),
      exposure_type = "categorical",
      reference_level = "zzz",
      n_boot = 10
    ),
    class = "positively_args_error"
  )
  expect_snapshot_abort(
    check_eta_bias(
      data,
      a,
      y,
      c(x1, x2),
      exposure_type = "categorical",
      reference_level = c("a", "b"),
      n_boot = 10
    ),
    class = "positively_args_error"
  )
})

test_that("a categorical run records no upper truncation bound", {
  local_quiet()
  skip_if_not_installed("nnet")
  data <- sim_eta_categorical(n = 300, seed = 1)

  single <- fit_eta(data, "ipw", exposure_type = "categorical", n_boot = 20)
  swept <- fit_eta(
    data,
    "ipw",
    exposure_type = "categorical",
    truncation_grid = c(0, 0.05, 0.1),
    n_boot = 20
  )

  # A lower pin on every level already bounds every weight denominator, so past
  # two levels there is no upper bound to apply and none to report.
  expect_identical(single@results$truncation_upper, rep(NA_real_, 2L))
  expect_identical(single@results$truncation_lower, rep(0, 2L))
  expect_identical(swept@results$truncation_upper, rep(NA_real_, 6L))
})

test_that("the two-sided truncation bound is rejected past two levels", {
  local_quiet()
  # The nnet guard resolves ahead of the truncation bounds, so this abort is
  # reachable only where the categorical path could have run.
  skip_if_not_installed("nnet")
  data <- sim_eta_categorical(n = 200, seed = 1)

  # `truncation` names a lower and an upper probability bound, which describes
  # the binary simplex and nothing else. Silently ignoring it would change the
  # number reported without saying so. The abort has to name the argument it
  # refuses, since the exposure type is not the argument at fault.
  err <- expect_error(
    check_eta_bias(
      data,
      a,
      y,
      c(x1, x2),
      exposure_type = "categorical",
      truncation = c(0.05, 0.95),
      n_boot = 10
    ),
    class = "positively_args_error"
  )
  expect_match(conditionMessage(err), "truncation", fixed = TRUE)

  expect_snapshot_abort(
    check_eta_bias(
      data,
      a,
      y,
      c(x1, x2),
      truncation = c(0.05, 0.95),
      n_boot = 10
    ),
    class = "positively_args_error"
  )
})

test_that("the categorical truncation grid stops at 1 over the level count", {
  local_quiet()
  skip_if_not_installed("nnet")
  data <- sim_eta_categorical(n = 300, seed = 1)

  err <- expect_error(
    check_eta_bias(
      data,
      a,
      y,
      c(x1, x2),
      exposure_type = "categorical",
      truncation_grid = c(0.1, 0.4),
      n_boot = 10
    ),
    class = "positively_range_error"
  )
  expect_match(conditionMessage(err), "0.333", fixed = TRUE)
  expect_match(conditionMessage(err), "1/k", fixed = TRUE)
  expect_match(conditionMessage(err), "3 levels", fixed = TRUE)

  # A bound below the ceiling leaves room on the simplex and is accepted.
  res <- fit_eta(
    data,
    "ipw",
    exposure_type = "categorical",
    truncation_grid = c(0, 0.3),
    n_boot = 20
  )
  expect_identical(nrow(res@results), 4L)
})

test_that("a categorical violation leaves gcomp flat and inflates ipw", {
  local_quiet()
  skip_if_not_installed("nnet")
  # The reference level is rare wherever either covariate is large, so both
  # contrasts against it lean on a sparse arm.
  data <- sim_eta_categorical(n = 600, steepness = 4, seed = 1)
  grid <- c(0, 0.05, 0.1, 0.2)
  gcomp <- fit_eta(
    data,
    "gcomp",
    exposure_type = "categorical",
    truncation_grid = grid,
    n_boot = 100
  )
  ipw <- fit_eta(
    data,
    "ipw",
    exposure_type = "categorical",
    truncation_grid = grid,
    n_boot = 100
  )

  for (term in c("b - a", "c - a")) {
    gcomp_term <- gcomp@results[gcomp@results$term == term, ]
    ipw_term <- ipw@results[ipw@results$term == term, ]
    gcomp_term <- gcomp_term[order(gcomp_term$truncation_lower), ]
    ipw_term <- ipw_term[order(ipw_term$truncation_lower), ]

    # G-computation reads no fitted probability, so truncation cannot move it
    # and its residual bias is Monte Carlo noise.
    expect_lt(max(abs(gcomp_term$bias)), 0.05)
    expect_equal(diff(gcomp_term$bias), rep(0, length(grid) - 1))

    # The untruncated weights already carry clear bias, and clamping every level
    # from below trades more of it for a tighter bootstrap spread.
    expect_gt(ipw_term$bias[[1]], 0.15)
    expect_gt(ipw_term$bias[[1]], 3 * ipw_term$mc_se[[1]])
    expect_true(all(diff(abs(ipw_term$bias)) > 0))
  }
})

test_that("a multi-term run prints one reading per term", {
  local_quiet()
  skip_if_not_installed("nnet")
  data <- sim_eta_categorical(n = 300, seed = 1)
  res <- fit_eta(data, "ipw", exposure_type = "categorical", n_boot = 50)
  printed <- printed_text(res)

  # One truncation level and two terms is two rows, so a reading keyed on the
  # row count alone would report a range over levels there is only one of. Each
  # term states its own bias and Monte Carlo error instead.
  expect_match(printed, "b - a", fixed = TRUE)
  expect_match(printed, "c - a", fixed = TRUE)
  expect_match(printed, "MC SE", fixed = TRUE)
  expect_no_match(printed, "Truncation levels", fixed = TRUE)

  expect_snapshot(print(res))
})

test_that("glance() reports the term count in place of a single truth", {
  local_quiet()
  skip_if_not_installed("nnet")
  data <- sim_eta_categorical(n = 300, seed = 1)
  res <- fit_eta(data, "ipw", exposure_type = "categorical", n_boot = 50)
  glanced <- generics::glance(res)

  expect_identical(nrow(glanced), 1L)
  # An estimand of more than one term has one truth per term, so no single truth
  # describes the run.
  expect_false("truth" %in% names(glanced))
  expect_setequal(
    names(glanced),
    c(
      "n",
      "estimator",
      "n_boot",
      "n_terms",
      "n_levels",
      "bias_min",
      "bias_max"
    )
  )

  expect_identical(glanced$n_terms, 2L)
  expect_identical(glanced$n_levels, 1L)
  expect_equal(glanced$bias_min, min(res@results$bias))
  expect_equal(glanced$bias_max, max(res@results$bias))
})

test_that("summary() follows the multi-term glance() shape", {
  local_quiet()
  skip_if_not_installed("nnet")
  data <- sim_eta_categorical(n = 300, seed = 1)
  res <- fit_eta(data, "ipw", exposure_type = "categorical", n_boot = 50)
  summarized <- summary(res)

  expect_identical(
    summarized$statistic,
    c("n_terms", "n_levels", "bias_min", "bias_max")
  )
  expect_identical(
    summarized$value,
    c(2, 1, min(res@results$bias), max(res@results$bias))
  )
  expect_true(all(is.na(summarized$threshold)))
})

# ---- Continuous exposures ---------------------------------------------------
# A continuous exposure has no levels to contrast, so the estimand is the
# coefficient set of a marginal structural model read over a grid of exposure
# deciles: one term by default, since the default working model is linear in the
# exposure. The treatment mechanism is a linear model and the estimator's weight
# is the stabilized density ratio f(A) / f(A | W), so truncation caps a weight
# rather than bounding a probability, and the arguments that name a probability
# bound or a reference level have nothing to act on.

test_that("check_eta_bias() lists continuous among its supported types", {
  expect_identical(
    diagnostic_supported_types()[["eta_bias"]],
    c("binary", "categorical", "continuous")
  )
})

test_that("a declared continuous type runs a numeric exposure", {
  local_quiet()
  data <- sim_eta_continuous(n = 300, seed = 1)
  res <- fit_eta(data, "gcomp", exposure_type = "continuous", n_boot = 30)

  expect_identical(S7::S7_class(res)@name, "eta_bias_result")
  expect_identical(res@exposure_type, "continuous")
  expect_identical(res@n, 300L)
})

test_that("a numeric exposure is detected and announced as continuous", {
  withr::local_options(positively.quiet = FALSE)
  data <- sim_eta_continuous(n = 300, seed = 1)

  # Three supported types leave three branches the call could have taken, so the
  # announcement is what names the one it did.
  expect_message(
    res <- fit_eta(data, "gcomp", n_boot = 30),
    "continuous"
  )
  expect_identical(res@exposure_type, "continuous")
})

test_that("the stabilized weights match propensity::wt_ate()", {
  skip_if_not_installed("propensity")
  data <- sim_eta_continuous(n = 200, seed = 1)
  a <- data$a
  design <- cbind(1, data$x1, data$x2)
  fitted <- stats::glm.fit(design, a, family = stats::gaussian())$fitted.values
  # The pooled maximum-likelihood residual standard deviation, which is what
  # propensity derives for itself when it is handed no sigma.
  sigma <- sqrt(mean((a - fitted)^2))

  internal <- eta_stabilized_weights(a, fitted, sigma)

  # propensity forms the reciprocal of the conditional density and then
  # multiplies by the marginal one, where the single quotient is the natural
  # spelling, so the two agree to floating-point tolerance rather than bit for
  # bit.
  expect_equal(
    internal,
    as.numeric(propensity::wt_ate(
      fitted,
      a,
      exposure_type = "continuous",
      stabilize = TRUE
    ))
  )

  # Both densities carry the maximum-likelihood variance rather than the n - 1
  # estimate, on the numerator's marginal fit as much as the denominator's
  # conditional one.
  centre <- mean(a)
  expect_equal(
    internal,
    stats::dnorm(a, centre, sqrt(mean((a - centre)^2))) /
      stats::dnorm(a, fitted, sigma)
  )
})

test_that("the default working model is linear in the exposure", {
  local_quiet()
  data <- sim_eta_continuous(n = 300, seed = 1)
  res <- fit_eta(data, "gcomp", exposure_type = "continuous", n_boot = 30)

  # One coefficient past the intercept, so one term and one truth.
  expect_identical(nrow(res@results), 1L)
  expect_identical(res@results$term, "a")
  expect_identical(names(res@truth), "a")

  # A term is a working-model coefficient name, so it follows the exposure
  # column rather than a fixed label.
  names(data)[names(data) == "a"] <- "dose"
  renamed <- withr::with_seed(
    2024,
    check_eta_bias(
      data,
      dose,
      y,
      c(x1, x2),
      estimator = "gcomp",
      exposure_type = "continuous",
      n_boot = 30
    )
  )
  expect_identical(renamed@results$term, "dose")
  expect_identical(names(renamed@truth), "dose")
})

test_that("truth is the working model's projection over the exposure deciles", {
  local_quiet()
  # A quadratic outcome model curves the counterfactual surface, so projecting
  # it onto a linear working model depends on where the surface is read. Under a
  # linear outcome model every grid returns the same slope and this would pin
  # nothing.
  data <- sim_eta_continuous(n = 400, conditional_sd = 1.5, seed = 3)
  data$y <- data$y + 0.8 * data$a^2
  res <- fit_eta(
    data,
    "gcomp",
    exposure_type = "continuous",
    outcome_formula = y ~ a + I(a^2) + x1 + x2,
    n_boot = 20
  )

  fit <- stats::glm(y ~ a + I(a^2) + x1 + x2, data = data)
  grid <- unique(unname(stats::quantile(
    data$a,
    probs = seq(0.1, 0.9, by = 0.1),
    type = 7
  )))
  surface <- vapply(
    grid,
    function(value) {
      at_value <- data
      at_value$a <- value
      mean(stats::predict(fit, newdata = at_value, type = "response"))
    },
    numeric(1)
  )
  design <- cbind(1, grid)
  # Equal weights over a quantile grid weight by the marginal density of the
  # exposure, which is the standard choice of working-model weighting.
  projection <- drop(solve(crossprod(design), crossprod(design, surface)))

  expect_equal(unname(res@truth), projection[[2]])
})

test_that("msm_formula replaces the working model and its terms", {
  local_quiet()
  data <- sim_eta_continuous(n = 400, conditional_sd = 1.5, seed = 3)
  data$y <- data$y + 0.8 * data$a^2
  res <- fit_eta(
    data,
    "gcomp",
    exposure_type = "continuous",
    outcome_formula = y ~ a + I(a^2) + x1 + x2,
    msm_formula = ~ a + I(a^2),
    n_boot = 20
  )

  expect_identical(nrow(res@results), 2L)
  expect_identical(res@results$term, c("a", "I(a^2)"))
  expect_identical(names(res@truth), c("a", "I(a^2)"))

  # The counterfactual surface is exactly quadratic in the exposure here, so a
  # quadratic working model reproduces it rather than approximating it and the
  # projection returns the outcome model's own exposure coefficients.
  fit <- stats::glm(y ~ a + I(a^2) + x1 + x2, data = data)
  expect_equal(
    unname(res@truth),
    unname(stats::coef(fit)[c("a", "I(a^2)")])
  )
})

test_that("a quadratic working model takes the multi-term glance shape", {
  local_quiet()
  data <- sim_eta_continuous(n = 300, seed = 1)
  res <- fit_eta(
    data,
    "ipw",
    exposure_type = "continuous",
    msm_formula = ~ a + I(a^2),
    n_boot = 30
  )
  glanced <- generics::glance(res)

  expect_identical(nrow(glanced), 1L)
  # An estimand of more than one term has one truth per term, so no single truth
  # describes the run.
  expect_false("truth" %in% names(glanced))
  expect_setequal(
    names(glanced),
    c(
      "n",
      "estimator",
      "n_boot",
      "n_terms",
      "n_levels",
      "bias_min",
      "bias_max"
    )
  )
  expect_identical(glanced$n_terms, 2L)
  expect_identical(glanced$n_levels, 1L)
})

test_that("a single-term continuous run keeps the scalar-truth glance shape", {
  local_quiet()
  data <- sim_eta_continuous(n = 300, seed = 1)
  res <- fit_eta(data, "ipw", exposure_type = "continuous", n_boot = 30)
  glanced <- generics::glance(res)

  expect_setequal(
    names(glanced),
    c("n", "estimator", "n_boot", "truth", "bias", "mc_se")
  )
  expect_equal(glanced$truth, unname(res@truth))
  expect_equal(glanced$bias, res@results$bias)
})

test_that("the continuous print method is stable", {
  local_quiet()
  data <- sim_eta_continuous(n = 300, seed = 1)
  res <- fit_eta(data, "ipw", exposure_type = "continuous", n_boot = 30)
  printed <- printed_text(res)

  expect_match(printed, "continuous", fixed = TRUE)
  # One term at one truncation level is one row, which reads as a single bias
  # beside its Monte Carlo error rather than a range over levels.
  expect_match(printed, "MC SE", fixed = TRUE)
  expect_no_match(printed, "Truncation levels", fixed = TRUE)

  expect_snapshot(print(res))
})

test_that("a continuous run caps the weight instead of bounding a probability", {
  local_quiet()
  data <- sim_eta_continuous(n = 400, seed = 1)
  grid <- c(0, 0.05, 0.1)
  single <- fit_eta(data, "ipw", exposure_type = "continuous", n_boot = 20)
  swept <- fit_eta(
    data,
    "ipw",
    exposure_type = "continuous",
    truncation_grid = grid,
    n_boot = 20
  )

  # There is no lower cap to place on a density ratio, so the lower column is
  # missing at every level and the upper column carries the cap that was
  # applied. An untruncated run is capped at infinity, which is the widest
  # admissible cap rather than an absent one.
  expect_identical(single@results$truncation_lower, NA_real_)
  expect_identical(single@results$truncation_upper, Inf)
  expect_identical(swept@results$truncation_lower, rep(NA_real_, 3L))

  # The cap is the (1 - g) quantile of the stabilized weights under the observed
  # fit, resolved once and held fixed across bootstrap draws. A grid point of
  # zero leaves the weights alone rather than capping them at their own maximum,
  # which would still clip the heavier weights a bootstrap refit produces.
  ps <- stats::lm(a ~ x1 + x2, data = data)
  weights <- eta_stabilized_weights(
    data$a,
    stats::fitted(ps),
    sqrt(mean(stats::residuals(ps)^2))
  )
  expect_equal(
    swept@results$truncation_upper,
    c(
      Inf,
      unname(stats::quantile(weights, 0.95, type = 7)),
      unname(stats::quantile(weights, 0.90, type = 7))
    )
  )

  # A larger grid point reads a lower quantile, so the caps fall across the
  # sweep, and the sweep reports the number of levels it swept even though the
  # lower bound is missing at every one of them.
  expect_true(all(diff(swept@results$truncation_upper) < 0))
  expect_identical(generics::glance(swept)$n_levels, 3L)
})

# The three arguments a continuous exposure leaves nothing for are refused
# rather than ignored, since each of them shapes the number that would be
# reported. Each guard is stated against a detected type, which is the call a
# user with a numeric exposure makes, and each abort has to name the argument it
# refuses rather than the exposure type, which is not the thing at fault.

test_that("the two-sided truncation bound is rejected for a continuous exposure", {
  local_quiet()
  data <- sim_eta_continuous(n = 200, seed = 1)

  # `truncation` names a lower and an upper probability bound, and a continuous
  # exposure has no fitted probability to bound.
  err <- expect_error(
    check_eta_bias(
      data,
      a,
      y,
      c(x1, x2),
      truncation = c(0.05, 0.95),
      n_boot = 10
    ),
    class = "positively_args_error"
  )
  expect_match(conditionMessage(err), "truncation", fixed = TRUE)
})

test_that("the augmented estimator is rejected for a continuous exposure", {
  local_quiet()
  data <- sim_eta_continuous(n = 200, seed = 1)

  # The augmentation term integrates the outcome residual over the exposure
  # grid, which is not computed on this path.
  err <- expect_error(
    check_eta_bias(data, a, y, c(x1, x2), estimator = "aipw", n_boot = 10),
    class = "positively_args_error"
  )
  expect_match(conditionMessage(err), "aipw", fixed = TRUE)

  # The refusal comes before the bootstrap, so a caller learns of it without
  # paying for a run that cannot finish.
  testthat::local_mocked_bindings(
    eta_bias_bootstrap = function(...) stop("The bootstrap must not run.")
  )
  expect_error(
    check_eta_bias(data, a, y, c(x1, x2), estimator = "aipw", n_boot = 10),
    class = "positively_args_error"
  )
})

test_that("reference_level is rejected for a continuous exposure", {
  local_quiet()
  data <- sim_eta_continuous(n = 200, seed = 1)

  # `reference_level` names the level every contrast is taken against, which is
  # a statement about a discrete exposure. A continuous estimand is a set of
  # working-model coefficients, so there is no level to hold the others against
  # and honouring the argument would report a number it did not shape.
  err <- expect_error(
    check_eta_bias(data, a, y, c(x1, x2), reference_level = 0, n_boot = 10),
    class = "positively_args_error"
  )
  expect_match(conditionMessage(err), "reference_level", fixed = TRUE)
})

test_that("msm_formula is rejected for a discrete exposure", {
  local_quiet()
  data <- sim_eta_good(n = 200, seed = 1)

  # The mirror of the `reference_level` guard. A discrete estimand contrasts the
  # levels against the reference and imposes no working model, so a working model
  # supplied there would shape nothing and is refused rather than ignored.
  err <- expect_error(
    check_eta_bias(data, a, y, c(x1, x2), msm_formula = ~a, n_boot = 10),
    class = "positively_args_error"
  )
  expect_match(conditionMessage(err), "msm_formula", fixed = TRUE)
})

test_that("the continuous guard messages are stable", {
  local_quiet()
  data <- sim_eta_continuous(n = 200, seed = 1)

  expect_snapshot_abort(
    check_eta_bias(data, a, y, c(x1, x2), estimator = "aipw", n_boot = 10),
    class = "positively_args_error"
  )
  expect_snapshot_abort(
    check_eta_bias(
      data,
      a,
      y,
      c(x1, x2),
      truncation = c(0.05, 0.95),
      n_boot = 10
    ),
    class = "positively_args_error"
  )
  expect_snapshot_abort(
    check_eta_bias(data, a, y, c(x1, x2), reference_level = 0, n_boot = 10),
    class = "positively_args_error"
  )
})

test_that("a wide treatment mechanism leaves continuous ETA.Bias near zero", {
  local_quiet()
  # A conditional standard deviation of 3 against a marginal one near 3.3 leaves
  # the density ratio close to one everywhere, which is the continuous reading
  # of good overlap.
  data <- sim_eta_continuous(n = 600, conditional_sd = 3, seed = 1)
  ipw <- fit_eta(data, "ipw", exposure_type = "continuous", n_boot = 100)
  gcomp <- fit_eta(data, "gcomp", exposure_type = "continuous", n_boot = 100)

  expect_lt(abs(bias_of(ipw)), 0.05)
  expect_lt(abs(bias_of(gcomp)), 0.02)
})

test_that("a sharp treatment mechanism inflates ipw and leaves gcomp flat", {
  local_quiet()
  # A conditional standard deviation of 0.6 puts the density ratio below the
  # sqrt(2) threshold where its variance stops being finite, so a handful of
  # draws carry most of the weight.
  data <- sim_eta_continuous(n = 600, conditional_sd = 0.6, seed = 1)
  grid <- c(0, 0.05)
  ipw <- fit_eta(
    data,
    "ipw",
    exposure_type = "continuous",
    truncation_grid = grid,
    n_boot = 100
  )
  gcomp <- fit_eta(
    data,
    "gcomp",
    exposure_type = "continuous",
    truncation_grid = grid,
    n_boot = 100
  )

  uncapped <- ipw@results$truncation_upper == Inf
  bias_uncapped <- ipw@results$bias[uncapped]
  bias_capped <- ipw@results$bias[!uncapped]

  expect_gt(bias_uncapped, 0.15)
  expect_gt(bias_uncapped, 3 * ipw@results$mc_se[uncapped])

  # G-computation reads no weight, so the cap cannot move it and its residual
  # bias is Monte Carlo noise.
  expect_lt(max(abs(gcomp@results$bias)), 0.05)
  expect_equal(diff(gcomp@results$bias), 0)
  expect_gt(abs(bias_uncapped), 10 * max(abs(gcomp@results$bias)))

  # Capping the weight trades more bias for a tighter bootstrap spread, the same
  # trade the binary truncation sweep reports.
  expect_gt(abs(bias_capped - bias_uncapped), 0.05)
})

# ---- Multi-term and continuous views ---------------------------------------
# Two things about a multi-term or continuous run reach the figures. The
# estimand gains a term dimension, so the views have to keep the terms apart
# rather than overlay them, and a continuous run reports no probability bound,
# so a truncation level is identified by the grid point that produced it rather
# than by a lower bound that is missing at every level. The blocks below read
# the built plot rather than the recorded image, so what they pin is the
# structure the figures are drawn from.

# The facet keys a built plot carries. ggplot2 spells its own layout bookkeeping
# columns in upper case, so dropping those leaves the variables the faceting was
# built on, whatever they are named.
facet_keys <- function(built) {
  layout <- built$layout$layout
  layout[!grepl("^[A-Z_]+$", names(layout))]
}

# Built layers are read by what draws them rather than by the order they were
# added in: points carry a shape, reference lines carry an intercept, and a
# histogram carries a count.
point_layer <- function(built) {
  Filter(function(layer) "shape" %in% names(layer), built$data)[[1]]
}

xintercept_layer <- function(built) {
  Filter(function(layer) "xintercept" %in% names(layer), built$data)[[1]]
}

histogram_layer <- function(built) {
  Filter(function(layer) "count" %in% names(layer), built$data)[[1]]
}

# The bootstrap draws each panel of a built plot holds.
draws_per_panel <- function(built) {
  counts <- histogram_layer(built)
  unname(vapply(split(counts$count, counts$PANEL), sum, numeric(1)))
}

# Every strip label a plot will draw, with the labeller already applied. The
# labeller is applied to the built layout by hand rather than through
# ggplot2::get_strip_labels(), which the declared ggplot2 floor of 3.5.0
# predates. A grid renders its columns and then its rows, a wrap its one
# dimension, which is the order and the content get_strip_labels() reports.
strip_labels <- function(plot) {
  built <- ggplot2::ggplot_build(plot)
  params <- built$layout$facet_params
  keys <- facet_keys(built)
  dimensions <- if (inherits(built$layout$facet, "FacetGrid")) {
    list(cols = names(params$cols), rows = names(params$rows))
  } else {
    list(cols = names(params$facets))
  }
  labeller <- match.fun(params$labeller)
  rendered <- lapply(names(dimensions), function(type) {
    vars <- intersect(names(keys), dimensions[[type]])
    if (length(vars) == 0) {
      return(NULL)
    }
    # A labeller built with .rows or .cols reads the margin it is handed off
    # this attribute, so a dimension is tagged before it is rendered.
    margin <- unique(keys[vars])
    attr(margin, "type") <- type
    labeller(margin)
  })
  as.character(unlist(rendered))
}

test_that("the categorical bootstrap view holds one facet per term and level", {
  local_quiet()
  skip_if_not_installed("nnet")
  data <- sim_eta_categorical(n = 200, seed = 1)
  res <- fit_eta(
    data,
    "ipw",
    exposure_type = "categorical",
    truncation_grid = c(0, 0.05, 0.1),
    n_boot = 20
  )

  built <- ggplot2::ggplot_build(ggplot2::autoplot(res, type = "bootstrap"))
  keys <- facet_keys(built)

  # Two terms across three levels is six cells, and each cell is a panel of its
  # own: pooling a level's terms into one panel puts the bootstrap
  # distributions of two different estimands in one histogram.
  expect_identical(nrow(keys), 6L)
  expect_identical(nrow(unique(keys)), 6L)
  expect_true("term" %in% names(keys))
  expect_setequal(as.character(keys$term), c("b - a", "c - a"))

  # Every panel holds one cell's draws, so no panel counts a draw twice.
  expect_identical(draws_per_panel(built), rep(20, 6L))
})

test_that("each categorical bootstrap facet marks its own term's truth", {
  local_quiet()
  skip_if_not_installed("nnet")
  data <- sim_eta_categorical(n = 200, seed = 1)
  res <- fit_eta(
    data,
    "ipw",
    exposure_type = "categorical",
    truncation_grid = c(0, 0.05, 0.1),
    n_boot = 20
  )

  built <- ggplot2::ggplot_build(ggplot2::autoplot(res, type = "bootstrap"))
  layout <- built$layout$layout
  truths <- xintercept_layer(built)

  # The truth is per term, so one line per panel. A truth drawn from the whole
  # named vector puts every term's line in every panel, which reads as a second
  # target the panel's estimator was never aimed at.
  expect_identical(nrow(truths), nrow(layout))
  expect_true("term" %in% names(layout))

  terms <- as.character(layout[["term"]])
  for (panel in seq_along(terms)) {
    drawn <- truths$xintercept[truths$PANEL == layout$PANEL[[panel]]]
    expect_equal(drawn, res@truth[[terms[[panel]]]])
  }
})

test_that("the categorical sweep view draws one line per term", {
  local_quiet()
  skip_if_not_installed("nnet")
  data <- sim_eta_categorical(n = 200, seed = 1)
  res <- fit_eta(
    data,
    "ipw",
    exposure_type = "categorical",
    truncation_grid = c(0, 0.05, 0.1),
    n_boot = 20
  )

  plot <- ggplot2::autoplot(res, type = "sweep")
  built <- ggplot2::ggplot_build(plot)
  points <- point_layer(built)

  # A categorical run sweeps a probability bound, so the axis is the one the
  # binary sweep reads.
  expect_identical(built$plot$labels$x, "Truncation lower bound")

  # Terms are told apart by color, with a key that names them, rather than
  # joined into one line that walks the grid once per term.
  guide <- ggplot2::get_guide_data(plot, "colour")
  expect_setequal(as.character(guide$.label), c("b - a", "c - a"))
  expect_length(unique(points$colour), 2L)
  expect_length(unique(points$group), 2L)

  # Each term's line carries that term's bias across the grid.
  for (term in unique(res@results$term)) {
    rows <- res@results[res@results$term == term, ]
    rows <- rows[order(rows$truncation_lower), ]
    drawn <- points[points$colour %in% guide$colour[guide$.label == term], ]
    drawn <- drawn[order(drawn$x), ]
    expect_equal(drawn$x, rows$truncation_lower)
    expect_equal(drawn$y, rows$bias)
  }
})

test_that("a one-term binary sweep keeps its axis and draws no key", {
  local_quiet()
  data <- sim_eta_good(n = 200, seed = 1)
  res <- fit_eta(
    data,
    "ipw",
    truncation_grid = c(0, 0.05, 0.1),
    n_boot = 20
  )

  plot <- ggplot2::autoplot(res, type = "sweep")
  built <- ggplot2::ggplot_build(plot)

  # A binary run sweeps a probability bound and holds one term, which needs no
  # key to tell terms apart, so the term dimension leaves this figure as it is.
  expect_identical(built$plot$labels$x, "Truncation lower bound")
  expect_null(ggplot2::get_guide_data(plot, "colour"))
  expect_equal(point_layer(built)$x, c(0, 0.05, 0.1))
})

test_that("the continuous sweep reads its axis from the truncation grid", {
  local_quiet()
  data <- sim_eta_continuous(n = 200, seed = 1)
  grid <- c(0, 0.05, 0.1)
  res <- fit_eta(
    data,
    "ipw",
    exposure_type = "continuous",
    truncation_grid = grid,
    n_boot = 20
  )

  plot <- ggplot2::autoplot(res, type = "sweep")
  built <- ggplot2::ggplot_build(plot)
  points <- point_layer(built)

  # A continuous run caps a weight rather than bounding a probability, so the
  # lower bound is missing at every level and the sweep is drawn against the
  # grid points the caller swept. Reading the missing column leaves the panel
  # empty.
  expect_identical(res@params$truncation_grid, grid)
  expect_false(anyNA(points$x))
  expect_equal(points$x, grid)
  expect_equal(points$y, res@results$bias)

  # The axis names what it carries, which is a quantile of the weights rather
  # than a probability bound.
  expect_no_match(built$plot$labels$x, "Truncation lower bound", fixed = TRUE)
  expect_match(built$plot$labels$x, "quantile", ignore.case = TRUE)
})

test_that("the continuous bootstrap view holds one facet per level", {
  local_quiet()
  data <- sim_eta_continuous(n = 200, seed = 1)
  res <- fit_eta(
    data,
    "ipw",
    exposure_type = "continuous",
    truncation_grid = c(0, 0.05, 0.1),
    n_boot = 20
  )

  plot <- ggplot2::autoplot(res, type = "bootstrap")
  built <- ggplot2::ggplot_build(plot)

  # Three levels are three panels. Keying the facet on a bound a continuous run
  # never reports collapses the three bootstrap distributions into one panel
  # labeled for a missing value, which hides that the levels differ at all.
  expect_identical(nrow(built$layout$layout), 3L)
  expect_false(anyNA(unlist(facet_keys(built))))
  expect_false(any(grepl("NA", strip_labels(plot), fixed = TRUE)))
  expect_identical(draws_per_panel(built), rep(20, 3L))
})

test_that("the sweep view refuses a single-level multi-term or continuous run", {
  local_quiet()
  skip_if_not_installed("nnet")
  categorical <- fit_eta(
    sim_eta_categorical(n = 200, seed = 1),
    "ipw",
    exposure_type = "categorical",
    n_boot = 20
  )
  continuous <- fit_eta(
    sim_eta_continuous(n = 200, seed = 1),
    "ipw",
    exposure_type = "continuous",
    n_boot = 20
  )

  # A multi-term run at one truncation level has one row per term and still no
  # sweep to draw, so the refusal reads the level count rather than the row
  # count, and it says what the binary refusal says.
  expect_error(
    ggplot2::autoplot(categorical, type = "sweep"),
    "more than one level",
    class = "positively_sweep_absent_error"
  )
  expect_error(
    ggplot2::autoplot(continuous, type = "sweep"),
    "more than one level",
    class = "positively_sweep_absent_error"
  )
})

test_that("plot() draws the multi-term and continuous views", {
  local_quiet()
  local_null_device()
  skip_if_not_installed("nnet")
  categorical <- fit_eta(
    sim_eta_categorical(n = 200, seed = 1),
    "ipw",
    exposure_type = "categorical",
    truncation_grid = c(0, 0.05, 0.1),
    n_boot = 20
  )
  continuous <- fit_eta(
    sim_eta_continuous(n = 200, seed = 1),
    "ipw",
    exposure_type = "continuous",
    truncation_grid = c(0, 0.05, 0.1),
    n_boot = 20
  )

  # Both entry points draw the same views, so the type argument reaches the
  # multi-term and continuous shapes through the dots as it does the binary one.
  expect_identical(plot(categorical, type = "bootstrap"), categorical)
  expect_identical(plot(categorical, type = "sweep"), categorical)
  expect_identical(plot(continuous, type = "bootstrap"), continuous)
  expect_identical(plot(continuous, type = "sweep"), continuous)
})

test_that("the categorical ETA bias views render as expected", {
  local_quiet()
  announce_doppelganger(
    "ETA bias categorical bootstrap grid",
    "ETA bias categorical truncation sweep"
  )
  skip_if_not_installed("nnet")
  res <- fit_eta(
    sim_eta_categorical(n = 200, seed = 1),
    "ipw",
    exposure_type = "categorical",
    truncation_grid = c(0, 0.05, 0.1),
    n_boot = 20
  )

  expect_doppelganger(
    "ETA bias categorical bootstrap grid",
    ggplot2::autoplot(res, type = "bootstrap")
  )
  expect_doppelganger(
    "ETA bias categorical truncation sweep",
    ggplot2::autoplot(res, type = "sweep")
  )
})

test_that("the continuous ETA bias sweep renders as expected", {
  local_quiet()
  res <- fit_eta(
    sim_eta_continuous(n = 200, seed = 1),
    "ipw",
    exposure_type = "continuous",
    truncation_grid = c(0, 0.05, 0.1),
    n_boot = 20
  )

  expect_doppelganger(
    "ETA bias continuous truncation sweep",
    ggplot2::autoplot(res, type = "sweep")
  )
})

# ---- Function-valued estimators ---------------------------------------------
# `estimator` also takes an estimator of the caller's own: a bare function, or a
# named length-one list whose name labels it. The function is handed one
# bootstrap dataset per draw, in the terms the caller supplied the data in, and
# returns one number per estimand term. Nothing else about the run changes: the
# truth is still G-computation under the observed fit, and the draws are
# generated before any estimate is formed, so a function that spells a built-in
# estimator out sees the same datasets that built-in sees and has to reproduce
# it.

# The three replicas below are those spellings. Each mirrors the arithmetic of
# one built-in path in terms of the data frame the function receives, so a
# disagreement with the built-in is a disagreement about the datasets, the
# truth, or the summaries rather than about the estimator.

# The stabilized-weight working-model slope, which is continuous `ipw`: a normal
# treatment mechanism with the pooled maximum-likelihood residual spread, the
# density ratio f(A) / f(A | W) under it, and the weighted least-squares fit of
# the outcome on the working model linear in the exposure.
replica_continuous_ipw <- function(.data) {
  ps <- stats::lm(a ~ x1 + x2, data = .data)
  fitted <- stats::fitted(ps)
  sigma <- sqrt(mean(stats::residuals(ps)^2))
  centre <- mean(.data$a)
  weights <- stats::dnorm(.data$a, centre, sqrt(mean((.data$a - centre)^2))) /
    stats::dnorm(.data$a, fitted, sigma)
  design <- cbind(1, .data$a)
  stats::lm.wfit(design, .data$y, weights)$coefficients[[2]]
}

# The per-level Hajek contrasts against the reference level, which is
# categorical `ipw`, over a multinomial-logit treatment mechanism. The return is
# unnamed, so it is read in term order.
replica_categorical_ipw <- function(.data) {
  fit <- nnet::multinom(a ~ x1 + x2, data = .data, trace = FALSE)
  probabilities <- stats::fitted(fit)
  mean_at <- function(label) {
    weights <- (.data$a == label) / probabilities[, label]
    sum(weights * .data$y) / sum(weights)
  }
  reference <- mean_at("a")
  c(mean_at("b") - reference, mean_at("c") - reference)
}

# The G-computation contrast at two levels: an outcome regression, counterfactual
# predictions with the exposure set to each level in turn, and the mean
# difference. Reading the levels off the delivered column is what makes this a
# statement about the exposure arriving as its own levels.
replica_binary_gcomp <- function(.data) {
  fit <- stats::lm(y ~ a + x1 + x2, data = .data)
  labels <- levels(.data$a)
  mean_at <- function(label) {
    frame <- .data
    frame$a <- factor(label, levels = labels)
    mean(stats::predict(fit, newdata = frame))
  }
  mean_at(labels[[2]]) - mean_at(labels[[1]])
}

test_that("a bare function is accepted and labeled custom", {
  local_quiet()
  data <- sim_eta_good(n = 200, seed = 1)
  # A constant estimator fixes every bootstrap estimate, so the summaries are
  # determined and each one pins the arithmetic it is formed by rather than a
  # realization of it.
  res <- fit_eta(data, function(.data) 0.25, n_boot = 10)

  expect_identical(res@estimator, "custom")
  expect_identical(res@params$estimator, "custom")
  expect_identical(generics::glance(res)$estimator, "custom")

  expect_equal(res@boot_estimates[[1]], rep(0.25, 10))
  expect_equal(res@results$boot_mean, 0.25)
  expect_equal(res@results$mc_se, 0)
  expect_equal(res@results$bias, 0.25 - unname(res@truth))

  # The truth is G-computation under the observed fit, which no estimator has a
  # say in, so it is the same number the built-in run reports.
  expect_identical(res@truth, fit_eta(data, "gcomp", n_boot = 10)@truth)
})

test_that("a named list labels the estimator by its list name", {
  local_quiet()
  data <- sim_eta_good(n = 200, seed = 1)
  res <- fit_eta(data, list(cbps = function(.data) 0.25), n_boot = 10)

  expect_identical(res@estimator, "cbps")
  expect_identical(res@params$estimator, "cbps")
  expect_identical(generics::glance(res)$estimator, "cbps")
  expect_match(printed_text(res), "Estimator: cbps", fixed = TRUE)

  # The label is what the result carries. Holding the function would keep the
  # caller's whole calling environment alive inside a result they saved for its
  # numbers.
  expect_false(any(vapply(res@params, is.function, logical(1))))
})

test_that("the estimator argument rejects every other shape", {
  local_quiet()
  data <- sim_eta_good(n = 80, seed = 1)
  estimate <- function(.data) 0.25
  # Every refusal below shares its condition class with the refusal the argument
  # gives today for naming a function at all, so the message is what separates
  # them. An argument that takes a function may not say it takes a character
  # vector, and each refusal has to describe the form it was looking for.
  refusal <- function(estimator) {
    err <- expect_error(
      check_eta_bias(
        data,
        a,
        y,
        c(x1, x2),
        estimator = estimator,
        n_boot = 10
      ),
      class = "positively_args_error"
    )
    message <- conditionMessage(err)
    expect_no_match(message, "character vector", fixed = TRUE)
    expect_match(message, "function", fixed = TRUE)
    message
  }

  # A list carries the label, so a list without one names nothing. The message
  # is bound before it is matched, because `expect_match()` evaluates its first
  # argument twice and a second refusal would double every expectation above.
  unnamed <- refusal(list(estimate))
  expect_match(unnamed, "name", fixed = TRUE)
  # One run diagnoses one estimator, so a list of two describes two runs.
  refusal(list(cbps = estimate, tmle = estimate))
  # A named list is a labeled function and not a labeled anything else.
  refusal(list(cbps = "ipw"))
  refusal(1)
})

test_that("the estimator function is called once per bootstrap draw", {
  local_quiet()
  data <- sim_eta_good(n = 200, seed = 1)
  calls <- 0L
  seen <- NULL
  capture <- function(.data) {
    calls <<- calls + 1L
    if (is.null(seen)) {
      seen <<- .data
    }
    0.25
  }
  fit_eta(data, capture, n_boot = 12)

  expect_identical(calls, 12L)
  expect_true(is.data.frame(seen))
  expect_identical(nrow(seen), 200L)
  # The covariates, the exposure, and the outcome, under the names the caller
  # supplied them under and with nothing else alongside them.
  expect_setequal(names(seen), c("a", "y", "x1", "x2"))

  # The covariate rows are a resample of the observed ones, and the exposure and
  # outcome are drawn from the fitted mechanisms rather than carried over.
  expect_true(all(seen$x1 %in% data$x1))
  expect_true(all(seen$x2 %in% data$x2))
  expect_false(identical(as.character(seen$a), as.character(data$a)))
  expect_false(identical(seen$y, data$y))

  # A two-level exposure arrives as its own levels, which is a factor of the two
  # observed values rather than the codes the treatment mechanism is fitted on.
  expect_s3_class(seen$a, "factor")
  expect_identical(levels(seen$a), c("0", "1"))
})

test_that("a categorical exposure reaches the function as its original levels", {
  local_quiet()
  skip_if_not_installed("nnet")
  data <- sim_eta_categorical(n = 200, seed = 1)
  seen <- NULL
  capture <- function(.data) {
    if (is.null(seen)) {
      seen <<- .data
    }
    c(0.25, 0.5)
  }
  fit_eta(data, capture, exposure_type = "categorical", n_boot = 5)

  # The labels, not the 0-based codes the multinomial mechanism is fitted on,
  # and every level the observed exposure holds even where a draw missed one.
  expect_s3_class(seen$a, "factor")
  expect_identical(levels(seen$a), c("a", "b", "c"))
  expect_true(all(as.character(seen$a) %in% c("a", "b", "c")))
})

test_that("a continuous exposure reaches the function as the drawn values", {
  local_quiet()
  data <- sim_eta_continuous(n = 200, seed = 1)
  seen <- NULL
  capture <- function(.data) {
    if (is.null(seen)) {
      seen <<- .data
    }
    0.25
  }
  fit_eta(data, capture, exposure_type = "continuous", n_boot = 5)

  expect_type(seen$a, "double")
  expect_length(seen$a, 200L)
  # A continuous draw lands wherever the normal mechanism puts it, so it is
  # neither the observed column nor a resample of its values.
  expect_false(identical(seen$a, data$a))
  expect_false(all(seen$a %in% data$a))
})

test_that("the estimator function's return value is checked on the first draw", {
  local_quiet()
  data <- sim_eta_good(n = 150, seed = 1)
  # A binary run has one estimand term, so one number is what the estimator owes
  # it. A return that cannot be read is a fault in the estimator rather than a
  # property of the data, so it is refused at the first draw rather than after
  # the whole bootstrap has run on it.
  refusal <- function(estimate) {
    calls <- 0L
    counted <- function(.data) {
      calls <<- calls + 1L
      estimate
    }
    expect_error(
      fit_eta(data, counted, n_boot = 10),
      class = "positively_estimator_error"
    )
    expect_identical(calls, 1L)
  }

  # Too many values, too few, names that are not the terms, and a return that is
  # not numeric at all.
  refusal(c(0.25, 0.5))
  refusal(numeric(0))
  refusal(c(zzz = 0.25))
  refusal("0.25")
})

test_that("a named return is read by name and an unnamed one in term order", {
  local_quiet()
  skip_if_not_installed("nnet")
  data <- sim_eta_categorical(n = 200, seed = 1)

  named <- fit_eta(
    data,
    function(.data) c("b - a" = 0.25, "c - a" = 0.5),
    exposure_type = "categorical",
    n_boot = 5
  )
  positional <- fit_eta(
    data,
    function(.data) c(0.25, 0.5),
    exposure_type = "categorical",
    n_boot = 5
  )

  expect_identical(named@results$term, c("b - a", "c - a"))
  expect_equal(named@results$boot_mean, c(0.25, 0.5))
  expect_equal(positional@results$boot_mean, c(0.25, 0.5))
  expect_equal(named@results$bias, positional@results$bias)
})

test_that("a named return is read by name whatever order it arrives in", {
  local_quiet()
  skip_if_not_installed("nnet")
  data <- sim_eta_categorical(n = 200, seed = 1)
  # Reading by name is what spares the caller the term order, so a return that
  # states the terms in the other order has to report the same two numbers
  # against the same two terms.
  permuted <- fit_eta(
    data,
    function(.data) c("c - a" = 0.5, "b - a" = 0.25),
    exposure_type = "categorical",
    n_boot = 5
  )

  expect_identical(permuted@results$term, c("b - a", "c - a"))
  expect_equal(permuted@results$boot_mean, c(0.25, 0.5))
})

# The three blocks below state that the reading is per draw. Each pins a named
# return against a control run returning the same numbers unnamed in term order
# under the same seed, so a reading that carried an order across draws would
# report the two terms' numbers swapped rather than merely differ in summary.

# An estimator that fails its first draw and returns `values` on every draw
# after it.
fails_first_draw <- function(values) {
  calls <- 0L
  function(.data) {
    calls <<- calls + 1L
    if (calls == 1L) {
      stop("no estimate on this draw")
    }
    values
  }
}

test_that("a failed first draw leaves later names read by name", {
  local_quiet()
  skip_if_not_installed("nnet")
  data <- sim_eta_categorical(n = 200, seed = 1)

  expect_warning(
    permuted <- fit_eta(
      data,
      fails_first_draw(c("c - a" = 0.5, "b - a" = 0.25)),
      exposure_type = "categorical",
      n_boot = 6
    ),
    class = "positively_degenerate_boot_warning"
  )
  expect_warning(
    control <- fit_eta(
      data,
      fails_first_draw(c(0.25, 0.5)),
      exposure_type = "categorical",
      n_boot = 6
    ),
    class = "positively_degenerate_boot_warning"
  )

  expect_identical(permuted@results$term, c("b - a", "c - a"))
  expect_equal(permuted@results$boot_mean, c(0.25, 0.5))
  expect_equal(permuted@results, control@results)
  expect_equal(permuted@boot_estimates, control@boot_estimates)
})

test_that("names that are not the terms are refused on whatever draw carries them", {
  local_quiet()
  skip_if_not_installed("nnet")
  data <- sim_eta_categorical(n = 200, seed = 1)
  # The names are the estimator's own claim about which term each value belongs
  # to, so a claim on a term the estimand does not hold is refused rather than
  # read positionally, whether it arrives on the first draw or a later one.
  calls <- 0L
  drifting <- function(.data) {
    calls <<- calls + 1L
    if (calls < 3L) {
      c("b - a" = 0.25, "c - a" = 0.5)
    } else {
      c(zzz = 0.25, "c - a" = 0.5)
    }
  }
  expect_error(
    fit_eta(data, drifting, exposure_type = "categorical", n_boot = 6),
    class = "positively_estimator_error"
  )

  # Refused where it appeared, so the draws behind it were never offered.
  expect_identical(calls, 3L)
})

test_that("named and unnamed returns may alternate across draws", {
  local_quiet()
  skip_if_not_installed("nnet")
  data <- sim_eta_categorical(n = 200, seed = 1)
  # Neither form settles anything for the draws behind it, so an estimator that
  # states the terms on some draws and leaves them to term order on the rest is
  # read correctly throughout.
  calls <- 0L
  alternating <- function(.data) {
    calls <<- calls + 1L
    if (calls %% 2L == 0L) {
      c("c - a" = 0.5, "b - a" = 0.25)
    } else {
      c(0.25, 0.5)
    }
  }
  mixed <- fit_eta(
    data,
    alternating,
    exposure_type = "categorical",
    n_boot = 6
  )
  control <- fit_eta(
    data,
    function(.data) c(0.25, 0.5),
    exposure_type = "categorical",
    n_boot = 6
  )

  expect_identical(calls, 6L)
  expect_equal(mixed@results, control@results)
  expect_equal(mixed@boot_estimates, control@boot_estimates)
})

test_that("truncation and truncation_grid are refused with a function estimator", {
  local_quiet()
  data <- sim_eta_good(n = 150, seed = 1)
  estimate <- function(.data) 0.25
  # Both arguments bound an internal nuisance the estimator never receives, so
  # honouring them would leave the caller a sweep whose every level reports the
  # same untruncated number, and a flat sweep is itself a finding.
  bounded <- expect_error(
    check_eta_bias(
      data,
      a,
      y,
      c(x1, x2),
      estimator = estimate,
      truncation = c(0.05, 0.95),
      n_boot = 10
    ),
    class = "positively_args_error"
  )
  expect_match(conditionMessage(bounded), "truncation", fixed = TRUE)
  expect_no_match(conditionMessage(bounded), "character vector", fixed = TRUE)

  swept <- expect_error(
    check_eta_bias(
      data,
      a,
      y,
      c(x1, x2),
      estimator = estimate,
      truncation_grid = c(0, 0.05),
      n_boot = 10
    ),
    class = "positively_args_error"
  )
  expect_match(conditionMessage(swept), "truncation_grid", fixed = TRUE)
  expect_no_match(conditionMessage(swept), "character vector", fixed = TRUE)
})

test_that("a function estimator reproduces the built-in continuous ipw run", {
  local_quiet()
  data <- sim_eta_continuous(n = 300, seed = 1)
  builtin <- fit_eta(data, "ipw", exposure_type = "continuous", n_boot = 20)
  # The augmented estimator is refused for a continuous exposure because the
  # built-in augmentation term is not defined over the grid. That is a statement
  # about the built-in and not about what a caller's own estimator may compute
  # here, so a function runs.
  custom <- fit_eta(
    data,
    list(wls = replica_continuous_ipw),
    exposure_type = "continuous",
    n_boot = 20
  )

  # Not identical: the internal path fits the treatment mechanism with
  # `glm.fit()` and the replica with `lm()`, which agree to the last few bits
  # rather than to all of them.
  expect_equal(custom@boot_estimates, builtin@boot_estimates)
  expect_equal(custom@results, builtin@results)
  expect_identical(custom@truth, builtin@truth)
  expect_identical(custom@estimator, "wls")
})

test_that("a function estimator reproduces the built-in categorical ipw run", {
  local_quiet()
  skip_if_not_installed("nnet")
  data <- sim_eta_categorical(n = 300, seed = 1)
  builtin <- fit_eta(data, "ipw", exposure_type = "categorical", n_boot = 20)
  custom <- fit_eta(
    data,
    list(hajek = replica_categorical_ipw),
    exposure_type = "categorical",
    n_boot = 20
  )

  expect_equal(custom@boot_estimates, builtin@boot_estimates)
  expect_equal(custom@results, builtin@results)
  expect_identical(custom@truth, builtin@truth)
})

test_that("a function estimator reproduces the built-in binary gcomp run", {
  local_quiet()
  data <- sim_eta_good(n = 300, seed = 1)
  builtin <- fit_eta(data, "gcomp", n_boot = 20)
  custom <- fit_eta(data, list(gformula = replica_binary_gcomp), n_boot = 20)

  expect_equal(custom@boot_estimates, builtin@boot_estimates)
  expect_equal(custom@results, builtin@results)
  expect_identical(custom@truth, builtin@truth)
})

test_that("a draw the estimator function fails on is dropped", {
  local_quiet()
  data <- sim_eta_good(n = 150, seed = 1)
  calls <- 0L
  flaky <- function(.data) {
    calls <<- calls + 1L
    if (calls %% 5L == 0L) {
      stop("no estimate on this draw")
    }
    0.25
  }
  # A draw the estimator failed on carries no more information than a draw whose
  # built-in estimate came back non-finite, so it is dropped the same way and
  # reported by the same warning.
  expect_warning(
    res <- fit_eta(data, flaky, n_boot = 20),
    class = "positively_degenerate_boot_warning"
  )

  # Every draw is still offered to the estimator; a failure on one says nothing
  # about the next.
  expect_identical(calls, 20L)
  expect_length(res@boot_estimates[[1]], 16L)
  expect_equal(res@results$boot_mean, 0.25)
  expect_identical(generics::glance(res)[["n_boot_used"]], 16L)
  expect_match(printed_text(res), "(4 dropped)", fixed = TRUE)
})

test_that("a run whose every draw fails aborts", {
  local_quiet()
  data <- sim_eta_good(n = 150, seed = 1)
  err <- expect_error(
    fit_eta(data, function(.data) stop("no estimate here"), n_boot = 5),
    class = "positively_error"
  )
  # There is no bias to report and no draw to report it from, so the run has to
  # say so rather than return a summary of nothing.
  expect_no_match(conditionMessage(err), "character vector", fixed = TRUE)
  expect_match(conditionMessage(err), "draw", fixed = TRUE)
})

# Dropped draws are reported the same way whatever dropped them, so the two
# blocks below state the reporting on a built-in run, where a severe practical
# violation is what leaves an estimate non-finite.

test_that("dropped draws are reported by print and glance", {
  local_quiet()
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
  printed <- printed_text(res)

  # The draw count asked for is still what the line reports, since the count
  # retained is a reading about the run rather than an argument to it.
  expect_match(printed, "Bootstrap draws: 50", fixed = TRUE)
  expect_match(printed, "(2 dropped)", fixed = TRUE)

  glanced <- generics::glance(res)
  expect_true("n_boot_used" %in% names(glanced))
  expect_identical(glanced[["n_boot_used"]], 48L)
  # The minimum across the cells, so that no cell is described by a count of
  # draws it did not have.
  expect_identical(glanced[["n_boot_used"]], min(lengths(res@boot_estimates)))

  expect_warning(
    swept <- withr::with_seed(
      1,
      check_eta_bias(
        data,
        a,
        y,
        x1,
        estimator = "ipw",
        truncation_grid = c(0, 0.05, 0.1),
        n_boot = 50
      )
    ),
    class = "positively_degenerate_boot_warning"
  )
  retained <- min(lengths(swept@boot_estimates))
  expect_identical(generics::glance(swept)[["n_boot_used"]], retained)
  expect_match(
    printed_text(swept),
    paste0("(", 50L - retained, " dropped)"),
    fixed = TRUE
  )
})

test_that("a run that dropped nothing reports no count", {
  local_quiet()
  data <- sim_eta_good(n = 300, seed = 1)
  res <- fit_eta(data, "ipw", n_boot = 50)

  # The reading is conditional on there being something to report, so the shapes
  # every other run takes are the shapes they already took.
  expect_length(res@boot_estimates[[1]], 50L)
  expect_no_match(printed_text(res), "dropped", fixed = TRUE)
  expect_false("n_boot_used" %in% names(generics::glance(res)))
})

test_that("the function-estimator argument messages are stable", {
  local_quiet()
  data <- sim_eta_good(n = 80, seed = 1)
  estimate <- function(.data) 0.25

  expect_snapshot_abort(
    check_eta_bias(
      data,
      a,
      y,
      c(x1, x2),
      estimator = list(estimate),
      n_boot = 10
    ),
    class = "positively_args_error"
  )
  expect_snapshot_abort(
    check_eta_bias(
      data,
      a,
      y,
      c(x1, x2),
      estimator = estimate,
      truncation = c(0.05, 0.95),
      n_boot = 10
    ),
    class = "positively_args_error"
  )
  expect_snapshot_abort(
    check_eta_bias(
      data,
      a,
      y,
      c(x1, x2),
      estimator = estimate,
      truncation_grid = c(0, 0.05),
      n_boot = 10
    ),
    class = "positively_args_error"
  )
})

test_that("the estimator return messages are stable", {
  local_quiet()
  data <- sim_eta_good(n = 80, seed = 1)

  expect_snapshot_abort(
    check_eta_bias(
      data,
      a,
      y,
      c(x1, x2),
      estimator = function(.data) c(0.25, 0.5),
      n_boot = 10
    ),
    class = "positively_estimator_error"
  )
  expect_snapshot_abort(
    check_eta_bias(
      data,
      a,
      y,
      c(x1, x2),
      estimator = function(.data) c(zzz = 0.25),
      n_boot = 10
    ),
    class = "positively_estimator_error"
  )
  expect_snapshot_abort(
    check_eta_bias(
      data,
      a,
      y,
      c(x1, x2),
      estimator = function(.data) stop("no estimate here"),
      n_boot = 5
    ),
    class = "positively_degenerate_boot_error"
  )
})

test_that("the custom estimator print method is stable", {
  local_quiet()
  data <- sim_eta_good(n = 120, seed = 1)
  res <- fit_eta(data, list(cbps = function(.data) 0.25), n_boot = 10)

  expect_match(printed_text(res), "cbps", fixed = TRUE)
  expect_snapshot(print(res))
})
