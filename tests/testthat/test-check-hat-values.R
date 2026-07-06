# Absolute phi_hat depends steeply on the candidate percentile grid, so every
# magnitude expectation is anchored to the default grid seq(0.05, 0.95, by =
# 0.05) and to the relative comparison against the null, never a hard-coded
# absolute value.

# ---- Local scenario generators --------------------------------------------
# The shared helpers in helper-dgp.R fix the exposure-covariate dependence. Hat
# values need to vary that dependence strength (beta) and the dose distribution,
# so these seeded generators live here. Using one seed across beta values holds
# the covariates fixed and varies only the dose, which isolates the effect of
# dependence strength.

sim_hat_linear <- function(n, beta, q = 1, seed = 1) {
  withr::local_seed(seed)
  x <- matrix(stats::rnorm(n * q), nrow = n, ncol = q)
  colnames(x) <- paste0("x", seq_len(q))
  out <- tibble::as_tibble(as.data.frame(x))
  # The dose depends on the first covariate with strength beta; any extra
  # covariates are noise, present only to exercise p = q + 2.
  out$dose <- stats::rnorm(n, mean = beta * x[, 1], sd = 1)
  out
}

sim_hat_skewed_null <- function(n, seed = 1) {
  withr::local_seed(seed)
  # A right-skewed dose drawn independently of the covariate. This is the
  # scenario in which a uniform null fabricates a violation; the shipped null
  # must stay quiet here.
  tibble::tibble(
    dose = stats::rlnorm(n),
    x1 = stats::rnorm(n)
  )
}

# Proportion of null replicates below the observed phi_hat, i.e. an estimate of
# Pr(phi_hat > phi_hat_0). Centered slightly below 0.5 under the null.
pr_gt_null <- function(x) {
  mean(x@null_dist < x@phi_hat)
}

local_quiet <- function(.env = parent.frame()) {
  withr::local_options(positively.quiet = TRUE, .local_envir = .env)
}

# ---- Signature and argument validation ------------------------------------

test_that("check_hat_values() aborts on a binary exposure", {
  local_quiet()
  # dgp_good_positivity() has a 0/1 exposure. Hat values are continuous-only.
  data <- dgp_good_positivity(n = 100, seed = 1)
  expect_error(
    check_hat_values(data, exposure, c(x1, x2)),
    class = "positively_error"
  )
})

test_that("check_hat_values() aborts on a categorical exposure", {
  local_quiet()
  data <- dgp_structural_subgroup(n = 200, seed = 1)
  data$exposure <- factor(sample(c("a", "b", "c"), nrow(data), replace = TRUE))
  expect_error(
    check_hat_values(data, exposure, x1),
    class = "positively_error"
  )
})

test_that("check_hat_values() rejects non-data-frame input", {
  local_quiet()
  expect_error(
    check_hat_values(1:10, dose, x1),
    class = "positively_error"
  )
})

test_that("check_hat_values() rejects an empty covariate selection", {
  local_quiet()
  data <- sim_hat_linear(80, beta = 1, seed = 1)
  expect_error(
    check_hat_values(data, dose, tidyselect::starts_with("zzz")),
    class = "positively_error"
  )
})

test_that("check_hat_values() rejects probs outside the unit interval", {
  local_quiet()
  data <- sim_hat_linear(80, beta = 1, seed = 1)
  expect_error(
    check_hat_values(data, dose, x1, probs = c(0.5, 1.5)),
    class = "positively_error"
  )
})

test_that("check_hat_values() rejects a conf_level outside the unit interval", {
  local_quiet()
  data <- sim_hat_linear(80, beta = 1, seed = 1)
  expect_error(
    check_hat_values(data, dose, x1, conf_level = 1.5),
    class = "positively_error"
  )
})

test_that("check_hat_values() requires a single exposure column", {
  local_quiet()
  data <- sim_hat_linear(80, beta = 1, seed = 1)
  expect_error(
    check_hat_values(data, c(dose, x1), x1),
    class = "positively_error"
  )
})

test_that("check_hat_values() rejects a non-scalar conf_level", {
  local_quiet()
  data <- sim_hat_linear(80, beta = 1, seed = 1)
  expect_error(
    check_hat_values(data, dose, x1, conf_level = c(0.5, 0.9)),
    class = "positively_error"
  )
})

test_that("check_hat_values() rejects a non-positive threshold", {
  local_quiet()
  data <- sim_hat_linear(80, beta = 1, seed = 1)
  expect_error(
    check_hat_values(data, dose, x1, threshold = 0),
    class = "positively_error"
  )
  expect_error(
    check_hat_values(data, dose, x1, threshold = "2"),
    class = "positively_error"
  )
})

test_that("check_hat_values() rejects a non-integer or zero null_reps", {
  local_quiet()
  data <- sim_hat_linear(80, beta = 1, seed = 1)
  expect_error(
    check_hat_values(data, dose, x1, null_reps = 0),
    class = "positively_error"
  )
  expect_error(
    check_hat_values(data, dose, x1, null_reps = 2.5),
    class = "positively_error"
  )
})

test_that("check_hat_values() rejects empty and NA probs", {
  local_quiet()
  data <- sim_hat_linear(80, beta = 1, seed = 1)
  expect_error(
    check_hat_values(data, dose, x1, probs = numeric(0)),
    class = "positively_error"
  )
  expect_error(
    check_hat_values(data, dose, x1, probs = c(0.5, NA)),
    class = "positively_error"
  )
})

test_that("check_hat_values() aborts on missing exposure or covariate values", {
  local_quiet()
  data <- sim_hat_linear(80, beta = 1, seed = 1)

  na_exposure <- data
  na_exposure$dose[1] <- NA
  expect_error(
    check_hat_values(na_exposure, dose, x1, null_reps = 2),
    class = "positively_error"
  )

  na_covariate <- data
  na_covariate$x1[1] <- NA
  expect_error(
    check_hat_values(na_covariate, dose, x1, null_reps = 2),
    class = "positively_error"
  )
})

test_that("check_hat_values() aborts on a non-numeric covariate", {
  local_quiet()
  data <- sim_hat_linear(80, beta = 1, seed = 1)
  data$group <- factor(sample(c("a", "b", "c"), nrow(data), replace = TRUE))
  expect_error(
    check_hat_values(data, dose, group, null_reps = 2),
    class = "positively_error"
  )
})

test_that("check_hat_values() aborts on a rank-deficient design", {
  local_quiet()
  data <- sim_hat_linear(80, beta = 1, seed = 1)
  data$constant <- 1
  expect_error(
    check_hat_values(data, dose, constant, null_reps = 2),
    class = "positively_error"
  )
})

# ---- Result class and properties ------------------------------------------

test_that("check_hat_values() returns a hat_values_result diagnostic", {
  local_quiet()
  data <- sim_hat_linear(150, beta = 1, seed = 1)
  res <- check_hat_values(data, dose, x1, null_reps = 20)

  expect_true(S7::S7_inherits(res, positivity_diagnostic))
  expect_identical(S7::S7_class(res)@name, "hat_values_result")
  expect_identical(res@exposure_type, "continuous")
  expect_identical(res@n, 150L)
})

test_that("hat_values_result carries the scalar null-comparison properties", {
  local_quiet()
  data <- sim_hat_linear(150, beta = 1, seed = 1)
  res <- check_hat_values(data, dose, x1, null_reps = 20)

  expect_type(res@phi_hat, "double")
  expect_length(res@phi_hat, 1)
  expect_type(res@null_dist, "double")
  expect_length(res@null_dist, 20)
  expect_type(res@null_quantile, "double")
  expect_length(res@null_quantile, 1)
  expect_type(res@exceeds_null, "logical")
  expect_length(res@exceeds_null, 1)
})

test_that("hat_values_result exposes the number of model parameters p", {
  local_quiet()
  data <- sim_hat_linear(150, beta = 1, q = 1, seed = 1)
  res <- check_hat_values(data, dose, x1, null_reps = 20)

  expect_type(res@p, "integer")
  expect_length(res@p, 1)
  # p = intercept + dose + 1 covariate, the ncol of the design matrix.
  expect_identical(res@p, 3L)
})

test_that("p equals the design matrix column count for several covariates", {
  local_quiet()
  data <- sim_hat_linear(160, beta = 1, q = 3, seed = 4)
  res <- check_hat_values(data, dose, c(x1, x2, x3), null_reps = 20)

  # p = intercept + dose + 3 covariates.
  expect_identical(res@p, 5L)
})

test_that("the leverage cutoff equals 2 * p / n through the p property", {
  local_quiet()
  data <- sim_hat_linear(120, beta = 1, q = 1, seed = 2)
  res <- check_hat_values(data, dose, x1, null_reps = 2)

  expect_identical(
    res@results$high_leverage,
    res@results$hat_value > 2 * res@p / res@n
  )
})

test_that("hat_values_result has the fixed results columns", {
  local_quiet()
  data <- sim_hat_linear(150, beta = 1, seed = 1)
  res <- check_hat_values(data, dose, x1, null_reps = 20)

  expect_s3_class(res@results, "tbl_df")
  expect_setequal(
    names(res@results),
    c(".id", "prob", "value", "hat_value", "high_leverage")
  )
  expect_type(res@results$high_leverage, "logical")
  expect_true(all(res@results$hat_value > 0))
})

# ---- Results shape and structure (exact, non-stochastic) ------------------

test_that("nrow(results) equals length(probs) times n (expectation 7)", {
  local_quiet()
  data <- sim_hat_linear(150, beta = 1, seed = 1)

  res <- check_hat_values(data, dose, x1, null_reps = 2)
  expect_identical(nrow(res@results), length(seq(0.05, 0.95, by = 0.05)) * 150L)

  res2 <- check_hat_values(
    data,
    dose,
    x1,
    probs = c(0.1, 0.5, 0.9),
    null_reps = 2
  )
  expect_identical(nrow(res2@results), 3L * 150L)
  expect_setequal(unique(res2@results$prob), c(0.1, 0.5, 0.9))
})

test_that("the candidate dose value is constant within each prob", {
  local_quiet()
  data <- sim_hat_linear(150, beta = 1, seed = 1)
  res <- check_hat_values(
    data,
    dose,
    x1,
    probs = c(0.1, 0.5, 0.9),
    null_reps = 2
  )

  by_prob <- split(res@results$value, res@results$prob)
  constant <- vapply(by_prob, function(v) length(unique(v)) == 1L, logical(1))
  expect_true(all(constant))
})

test_that("high_leverage equals hat_value > 2p/n with p = q + 2 (expectation 6)", {
  local_quiet()
  data <- sim_hat_linear(120, beta = 1, q = 1, seed = 2)
  res <- check_hat_values(data, dose, x1, null_reps = 2)

  # p = intercept + dose + 1 covariate.
  p <- 3L
  expect_identical(
    res@results$high_leverage,
    res@results$hat_value > 2 * p / res@n
  )
})

test_that("p is q + 2 with several covariates (expectations 6 and 10)", {
  local_quiet()
  data <- sim_hat_linear(160, beta = 1, q = 3, seed = 4)
  res <- check_hat_values(data, dose, c(x1, x2, x3), null_reps = 2)

  # p = intercept + dose + 3 covariates, exposed directly on the result.
  expect_identical(res@p, 5L)
  expect_identical(
    res@results$high_leverage,
    res@results$hat_value > 2 * res@p / res@n
  )
  expect_identical(nrow(res@results), 19L * 160L)
})

test_that("the threshold multiplier scales the leverage cutoff", {
  local_quiet()
  data <- sim_hat_linear(120, beta = 1, q = 1, seed = 2)
  res <- check_hat_values(data, dose, x1, threshold = 1, null_reps = 2)

  p <- 3L
  expect_identical(
    res@results$high_leverage,
    res@results$hat_value > 1 * p / res@n
  )
})

# ---- Statistical behavior --------------------------------------------------

test_that("the null stays quiet when the dose is independent of x", {
  local_quiet()
  # Expectations 1 and 4: exceeds_null is FALSE, median Pr(phi > phi0) < 0.55,
  # and phi_hat sits within 0.10 of the null median.
  fits <- lapply(1:3, function(s) {
    check_hat_values(
      sim_hat_linear(200, beta = 0, seed = s),
      dose,
      x1,
      null_reps = 100
    )
  })

  expect_false(any(vapply(fits, function(f) f@exceeds_null, logical(1))))
  expect_lt(stats::median(vapply(fits, pr_gt_null, numeric(1))), 0.55)
  gaps <- vapply(
    fits,
    function(f) abs(f@phi_hat - stats::median(f@null_dist)),
    numeric(1)
  )
  expect_true(all(gaps < 0.10))
})

test_that("a strong violation exceeds the null (expectation 2)", {
  local_quiet()
  data <- sim_hat_linear(200, beta = 2, seed = 3)
  res <- check_hat_values(data, dose, x1, null_reps = 100)

  expect_true(res@exceeds_null)
  expect_gte(pr_gt_null(res), 0.95)
  expect_gt(res@phi_hat, 3 * res@null_quantile)
})

test_that("median phi_hat increases with dependence strength (expectation 3)", {
  local_quiet()
  # null_reps is minimal because this test reads only phi_hat, not the null
  # distribution. Strict monotonicity plus a wide overall span is the robust
  # signal on the default grid.
  betas <- c(0, 0.5, 1, 2)
  med_phi <- vapply(
    betas,
    function(b) {
      stats::median(vapply(
        1:3,
        function(s) {
          res <- check_hat_values(
            sim_hat_linear(200, beta = b, seed = s),
            dose,
            x1,
            null_reps = 2
          )
          res@phi_hat
        },
        numeric(1)
      ))
    },
    numeric(1)
  )

  expect_true(all(diff(med_phi) > 0))
  expect_gt(med_phi[4] - med_phi[1], 0.3)
})

test_that("phi_hat is stable across n at a fixed violation (expectation 5)", {
  local_quiet()
  # A single realization per n uses a 0.10 band around a stable phi_hat.
  ns <- c(100, 250, 500)
  fits <- lapply(ns, function(nn) {
    check_hat_values(
      sim_hat_linear(nn, beta = 2, seed = 7),
      dose,
      x1,
      null_reps = 50
    )
  })

  phis <- vapply(fits, function(f) f@phi_hat, numeric(1))
  expect_lt(max(phis) - min(phis), 0.10)
  expect_true(all(vapply(fits, function(f) f@exceeds_null, logical(1))))
})

test_that("flags localize to where the extreme dose is unobserved (expectation 8)", {
  local_quiet()
  data <- sim_hat_linear(300, beta = 3, seed = 11)
  res <- check_hat_values(data, dose, x1, null_reps = 2)

  # At the top dose percentile, high-leverage candidates are those whose
  # covariate makes the extreme dose implausible: low x under positive
  # confounding.
  top <- res@results[res@results$prob == max(res@results$prob), ]
  cov <- data$x1[top$.id]

  flagged <- cov[top$high_leverage]
  unflagged <- cov[!top$high_leverage]
  expect_lt(mean(flagged), mean(unflagged) - stats::sd(cov))

  rank_corr <- stats::cor(top$hat_value, cov, method = "spearman")
  expect_lt(rank_corr, -0.5)
})

test_that("the shipped null does not fabricate a violation on a skewed dose (expectation 9)", {
  local_quiet()
  # A right-skewed dose drawn independently of the covariate must not read as a
  # violation under the default null.
  fits <- lapply(1:3, function(s) {
    check_hat_values(
      sim_hat_skewed_null(200, seed = s),
      dose,
      x1,
      null_reps = 100
    )
  })

  expect_false(any(vapply(fits, function(f) f@exceeds_null, logical(1))))
  expect_lte(stats::median(vapply(fits, pr_gt_null, numeric(1))), 0.6)
})

# ---- Magnitude expectations (loose, wide tolerances) ----------------------

test_that("a two-point {5th, 95th} grid gives a null baseline near 0.13 (expectation 11)", {
  local_quiet()
  # Anchored to the relative range [0.08, 0.20], never the paper's 0.29.
  med_phi <- stats::median(vapply(
    1:3,
    function(s) {
      res <- check_hat_values(
        sim_hat_linear(300, beta = 0, seed = s),
        dose,
        x1,
        probs = c(0.05, 0.95),
        null_reps = 2
      )
      res@phi_hat
    },
    numeric(1)
  ))

  expect_gt(med_phi, 0.08)
  expect_lt(med_phi, 0.20)
})

test_that("the default 19-point grid gives a low null baseline and a high violation (expectation 12)", {
  local_quiet()
  null_phi <- stats::median(vapply(
    1:3,
    function(s) {
      check_hat_values(
        sim_hat_linear(300, beta = 0, seed = s),
        dose,
        x1,
        null_reps = 2
      )@phi_hat
    },
    numeric(1)
  ))
  violation_phi <- stats::median(vapply(
    1:3,
    function(s) {
      check_hat_values(
        sim_hat_linear(300, beta = 2, seed = s),
        dose,
        x1,
        null_reps = 2
      )@phi_hat
    },
    numeric(1)
  ))

  expect_gt(null_phi, 0.02)
  expect_lt(null_phi, 0.10)
  # phi_hat is deterministic on this DGP; the violation median sits at about
  # 0.47 on the default grid.
  expect_gt(violation_phi, 0.4)
})

# ---- Null-resampling method -----------------------------------------------

test_that("null_method rejects an unknown scheme", {
  local_quiet()
  data <- sim_hat_linear(80, beta = 1, seed = 1)
  expect_error(
    check_hat_values(data, dose, x1, null_method = "uniform", null_reps = 2),
    class = "rlang_error"
  )
})

test_that("every null method stays quiet on a skewed independent dose", {
  local_quiet()
  methods <- c("permutation", "bootstrap", "gaussian")
  for (method in methods) {
    quiet <- vapply(
      1:3,
      function(s) {
        res <- check_hat_values(
          sim_hat_skewed_null(200, seed = s),
          dose,
          x1,
          null_method = method,
          null_reps = 50
        )
        res@exceeds_null
      },
      logical(1)
    )
    expect_false(any(quiet))
  }
})

test_that("every null method flags a strong violation", {
  local_quiet()
  methods <- c("permutation", "bootstrap", "gaussian")
  for (method in methods) {
    flagged <- vapply(
      1:3,
      function(s) {
        res <- check_hat_values(
          sim_hat_linear(200, beta = 2, seed = s),
          dose,
          x1,
          null_method = method,
          null_reps = 50
        )
        res@exceeds_null
      },
      logical(1)
    )
    expect_true(all(flagged))
  }
})

test_that("results have identical shape across null methods", {
  local_quiet()
  data <- sim_hat_linear(150, beta = 1, seed = 1)
  results <- lapply(
    c("permutation", "bootstrap", "gaussian"),
    function(method) {
      check_hat_values(
        data,
        dose,
        x1,
        null_method = method,
        null_reps = 2
      )@results
    }
  )

  shapes <- lapply(results, function(r) list(names = names(r), n = nrow(r)))
  expect_identical(shapes[[2]], shapes[[1]])
  expect_identical(shapes[[3]], shapes[[1]])
})

# ---- Autoplot contract ----------------------------------------------------

test_that("autoplot() returns a ggplot for each type", {
  local_quiet()
  data <- sim_hat_linear(150, beta = 1, seed = 1)
  res <- check_hat_values(data, dose, x1, null_reps = 50)

  expect_s3_class(ggplot2::autoplot(res, type = "null"), "ggplot")
  expect_s3_class(ggplot2::autoplot(res, type = "profile"), "ggplot")
})

test_that("hat-value autoplot views render as expected", {
  local_quiet()
  data <- sim_hat_linear(150, beta = 1, seed = 1)
  res <- check_hat_values(data, dose, x1, null_reps = 50)

  expect_doppelganger(
    "Hat values null distribution",
    ggplot2::autoplot(res, type = "null")
  )
  expect_doppelganger(
    "Hat values leverage profile",
    ggplot2::autoplot(res, type = "profile")
  )
})

# ---- Message and print snapshots ------------------------------------------

test_that("the continuous-only abort message is stable", {
  local_quiet()
  binary <- dgp_good_positivity(n = 100, seed = 1)
  categorical <- dgp_structural_subgroup(n = 200, seed = 1)
  categorical$exposure <- factor(
    sample(c("a", "b", "c"), nrow(categorical), replace = TRUE)
  )

  expect_snapshot(
    check_hat_values(binary, exposure, c(x1, x2)),
    error = TRUE
  )
  expect_snapshot(
    check_hat_values(categorical, exposure, x1),
    error = TRUE
  )
})

test_that("the argument validation messages are stable", {
  local_quiet()
  data <- sim_hat_linear(80, beta = 1, seed = 1)

  expect_snapshot(check_hat_values(1:10, dose, x1), error = TRUE)
  expect_snapshot(
    check_hat_values(data, dose, tidyselect::starts_with("zzz")),
    error = TRUE
  )
  expect_snapshot(
    check_hat_values(data, dose, x1, probs = c(0.5, 1.5)),
    error = TRUE
  )
  expect_snapshot(
    check_hat_values(data, dose, x1, conf_level = 1.5),
    error = TRUE
  )
  expect_snapshot(
    check_hat_values(data, dose, x1, null_method = "uniform"),
    error = TRUE
  )
  expect_snapshot(check_hat_values(data, c(dose, x1), x1), error = TRUE)
  expect_snapshot(
    check_hat_values(data, dose, x1, conf_level = c(0.5, 0.9)),
    error = TRUE
  )
  expect_snapshot(
    check_hat_values(data, dose, x1, threshold = 0),
    error = TRUE
  )
  expect_snapshot(
    check_hat_values(data, dose, x1, null_reps = 0),
    error = TRUE
  )
  expect_snapshot(
    check_hat_values(data, dose, x1, probs = numeric(0)),
    error = TRUE
  )
})

test_that("the data-integrity error messages are stable", {
  local_quiet()
  data <- sim_hat_linear(80, beta = 1, seed = 1)

  na_covariate <- data
  na_covariate$x1[1] <- NA
  expect_snapshot(
    check_hat_values(na_covariate, dose, x1, null_reps = 2),
    error = TRUE
  )

  factor_covariate <- data
  factor_covariate$group <- factor(
    sample(c("a", "b"), nrow(data), replace = TRUE)
  )
  expect_snapshot(
    check_hat_values(factor_covariate, dose, group, null_reps = 2),
    error = TRUE
  )

  constant_covariate <- data
  constant_covariate$constant <- 1
  expect_snapshot(
    check_hat_values(constant_covariate, dose, constant, null_reps = 2),
    error = TRUE
  )
})

test_that("the print method is stable", {
  local_quiet()
  data <- sim_hat_linear(150, beta = 1, seed = 1)
  # Seed the null draw so the printed null quantile is reproducible.
  withr::local_seed(2024)
  res <- check_hat_values(data, dose, x1, null_reps = 50)
  expect_snapshot(print(res))
})
