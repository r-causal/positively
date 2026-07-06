# Absolute Gower distances rescale with the pooled sample range, so magnitude
# expectations are anchored to relative separations and to hand-computable
# fixtures, never to a hard-coded third significant figure. Geometric
# variability follows the WhatIf definition: one half the mean of the full
# Gower distance matrix.

# ---- Local scenario generators --------------------------------------------
# The shared helpers in helper-dgp.R do not vary covariate dimension or the
# planted separation, both of which drive the extrapolation statistics, so
# these seeded generators live here. A balanced two-group Gaussian holds the
# design simple: each group is half the sample, and the treated group mean is
# shifted by `sep` on every axis.

sim_extrap_gaussian <- function(n, p = 4, sep = 0, seed = 1) {
  withr::local_seed(seed)
  n0 <- n %/% 2L
  n1 <- n - n0
  x0 <- matrix(stats::rnorm(n0 * p), nrow = n0, ncol = p)
  x1 <- matrix(stats::rnorm(n1 * p), nrow = n1, ncol = p) + sep
  x <- rbind(x0, x1)
  colnames(x) <- paste0("x", seq_len(p))
  out <- tibble::as_tibble(as.data.frame(x))
  out$exposure <- c(rep(0L, n0), rep(1L, n1))
  out[c("exposure", paste0("x", seq_len(p)))]
}

# Two overlapping numeric covariates plus a categorical covariate that nearly
# determines treatment in the categorical_gap scenario (category "c" is treated
# with probability 0.98), leaving the numeric axes fully overlapping.
sim_extrap_mixed <- function(
  n,
  scenario = c("overlap", "categorical_gap"),
  seed = 1
) {
  scenario <- match.arg(scenario)
  withr::local_seed(seed)
  x1 <- stats::rnorm(n)
  x2 <- stats::rnorm(n)
  g <- factor(sample(c("a", "b", "c"), n, replace = TRUE))
  prob <- if (scenario == "overlap") {
    rep(0.5, n)
  } else {
    ifelse(g == "c", 0.98, 0.4)
  }
  exposure <- stats::rbinom(n, 1L, prob)
  tibble::tibble(exposure = exposure, x1 = x1, x2 = x2, g = g)
}

# Fraction of units with an opposite-group neighbour inside the geometric
# variability radius, the recommended binary support indicator.
prop_supported <- function(res) {
  mean(res@results$gower_min <= res@geometric_variability)
}

prop_in_hull <- function(res) {
  mean(res@results$in_hull)
}

local_quiet <- function(.env = parent.frame()) {
  withr::local_options(positively.quiet = TRUE, .local_envir = .env)
}

# ---- Signature and argument validation ------------------------------------

test_that("check_extrapolation() aborts on a continuous exposure", {
  local_quiet()
  data <- dgp_continuous_support_gap(n = 100, seed = 1)
  expect_error(
    check_extrapolation(data, exposure, x1),
    class = "positively_error"
  )
})

test_that("check_extrapolation() aborts on a categorical exposure", {
  local_quiet()
  data <- dgp_good_positivity(n = 200, seed = 1)
  data$exposure <- factor(sample(c("a", "b", "c"), nrow(data), replace = TRUE))
  expect_error(
    check_extrapolation(data, exposure, c(x1, x2)),
    class = "positively_error"
  )
})

test_that("check_extrapolation() aborts on a single-level exposure", {
  local_quiet()
  data <- dgp_good_positivity(n = 100, seed = 1)
  data$exposure <- 1L
  expect_error(
    check_extrapolation(data, exposure, c(x1, x2)),
    class = "positively_error"
  )
})

test_that("check_extrapolation() rejects non-data-frame input", {
  local_quiet()
  expect_error(
    check_extrapolation(1:10, exposure, x1),
    class = "positively_error"
  )
})

test_that("check_extrapolation() rejects an empty covariate selection", {
  local_quiet()
  data <- dgp_good_positivity(n = 100, seed = 1)
  expect_error(
    check_extrapolation(data, exposure, tidyselect::starts_with("zzz")),
    class = "positively_error"
  )
})

test_that("check_extrapolation() rejects a multi-column exposure selection", {
  local_quiet()
  data <- dgp_good_positivity(n = 100, seed = 1)
  expect_error(
    check_extrapolation(data, c(exposure, x1), x2),
    class = "positively_error"
  )
})

test_that("check_extrapolation() rejects a non-positive or non-scalar nearby", {
  local_quiet()
  data <- dgp_good_positivity(n = 100, seed = 1)
  expect_error(
    check_extrapolation(data, exposure, c(x1, x2), nearby = 0),
    class = "positively_error"
  )
  expect_error(
    check_extrapolation(data, exposure, c(x1, x2), nearby = -1),
    class = "positively_error"
  )
  expect_error(
    check_extrapolation(data, exposure, c(x1, x2), nearby = "1"),
    class = "positively_error"
  )
  expect_error(
    check_extrapolation(data, exposure, c(x1, x2), nearby = c(1, 2)),
    class = "positively_error"
  )
})

test_that("check_extrapolation() rejects a non-logical hull argument", {
  local_quiet()
  data <- dgp_good_positivity(n = 100, seed = 1)
  expect_error(
    check_extrapolation(data, exposure, c(x1, x2), hull = "banana"),
    class = "positively_error"
  )
  expect_error(
    check_extrapolation(data, exposure, c(x1, x2), hull = c(TRUE, FALSE)),
    class = "positively_error"
  )
})

test_that("check_extrapolation() aborts on missing exposure or covariate values", {
  local_quiet()
  data <- dgp_good_positivity(n = 100, seed = 1)

  na_exposure <- data
  na_exposure$exposure[1] <- NA
  expect_error(
    check_extrapolation(na_exposure, exposure, c(x1, x2)),
    class = "positively_error"
  )

  na_covariate <- data
  na_covariate$x1[1] <- NA
  expect_error(
    check_extrapolation(na_covariate, exposure, c(x1, x2)),
    class = "positively_error"
  )
})

test_that("check_extrapolation() aborts on fewer than two observations", {
  local_quiet()
  data <- tibble::tibble(exposure = 1L, x1 = 0.5, x2 = -0.2)
  expect_error(
    check_extrapolation(data, exposure, c(x1, x2)),
    class = "positively_error"
  )
})

# ---- Result class and properties ------------------------------------------

test_that("check_extrapolation() returns an extrapolation_result diagnostic", {
  local_quiet()
  data <- sim_extrap_gaussian(200, p = 4, sep = 0, seed = 1)
  res <- check_extrapolation(data, exposure, tidyselect::starts_with("x"))

  expect_true(S7::S7_inherits(res, positivity_diagnostic))
  expect_identical(S7::S7_class(res)@name, "extrapolation_result")
  expect_identical(res@exposure_type, "binary")
  expect_identical(res@n, 200L)
})

test_that("extrapolation_result carries the geometric_variability and hull_run properties", {
  local_quiet()
  data <- sim_extrap_gaussian(200, p = 4, sep = 0, seed = 1)
  res <- check_extrapolation(data, exposure, tidyselect::starts_with("x"))

  expect_type(res@geometric_variability, "double")
  expect_length(res@geometric_variability, 1)
  expect_gt(res@geometric_variability, 0)
  expect_type(res@hull_run, "logical")
  expect_length(res@hull_run, 1)
})

test_that("extrapolation_result has the fixed results columns", {
  local_quiet()
  data <- sim_extrap_gaussian(200, p = 4, sep = 0, seed = 1)
  res <- check_extrapolation(
    data,
    exposure,
    tidyselect::starts_with("x"),
    hull = TRUE
  )

  expect_s3_class(res@results, "tbl_df")
  expect_setequal(
    names(res@results),
    c(".id", "exposure", "frac_nearby", "gower_min", "gower_mean", "in_hull")
  )
  expect_type(res@results$frac_nearby, "double")
  expect_type(res@results$gower_min, "double")
  expect_type(res@results$gower_mean, "double")
  expect_type(res@results$in_hull, "logical")
  expect_identical(nrow(res@results), 200L)
})

test_that("the per-group summary getter returns one row per exposure level", {
  local_quiet()
  data <- sim_extrap_gaussian(200, p = 4, sep = 0, seed = 1)
  res <- check_extrapolation(data, exposure, tidyselect::starts_with("x"))

  expect_s3_class(res@summary, "tbl_df")
  expect_identical(nrow(res@summary), 2L)
  expect_true("exposure" %in% names(res@summary))
})

test_that("tidy() and glance() follow the shared diagnostic contract", {
  local_quiet()
  data <- sim_extrap_gaussian(200, p = 4, sep = 0, seed = 1)
  res <- check_extrapolation(data, exposure, tidyselect::starts_with("x"))

  expect_identical(generics::tidy(res), res@results)
  glanced <- generics::glance(res)
  expect_s3_class(glanced, "tbl_df")
  expect_identical(nrow(glanced), 1L)
  expect_setequal(
    names(glanced),
    c(
      "exposure",
      "exposure_type",
      "n",
      "nearby",
      "geometric_variability",
      "prop_supported",
      "hull_run"
    )
  )
  expect_identical(glanced$nearby, 1)
})

# ---- Gower and geometric variability on hand-checkable fixtures -----------

test_that("gower_min and gower_mean match hand-computed mixed-type distances", {
  local_quiet()
  # Pooled ranges are 2 on x1 and 2 on x2. The treated unit (.id 1) sits at
  # Gower 1/6 from the matching-category control and 5/6 from the other.
  data <- tibble::tibble(
    exposure = c(1L, 0L, 0L),
    x1 = c(1, 0, 2),
    x2 = c(0, 0, 2),
    g = factor(c("a", "a", "b"))
  )
  res <- check_extrapolation(data, exposure, c(x1, x2, g), hull = FALSE)

  treated <- res@results[res@results$.id == 1, ]
  expect_equal(treated$gower_min, 1 / 6)
  expect_equal(treated$gower_mean, 0.5)
})

test_that("geometric_variability is one half the mean of the full Gower matrix", {
  local_quiet()
  # A single covariate with values 0, 1, 2, 3 has pooled range 3. The full
  # 4-by-4 Gower matrix sums to 20/3, so gv = (20/3) / (2 * 16) = 5/24.
  data <- tibble::tibble(exposure = c(0L, 1L, 0L, 1L), x = c(0, 1, 2, 3))
  res <- check_extrapolation(data, exposure, x, hull = FALSE)

  expect_equal(res@geometric_variability, 5 / 24)
})

# ---- Statistical separation -----------------------------------------------

test_that("Gower separates a planted support gap in both directions", {
  local_quiet()
  overlap <- check_extrapolation(
    sim_extrap_gaussian(600, p = 4, sep = 0, seed = 1),
    exposure,
    tidyselect::starts_with("x")
  )
  gap <- check_extrapolation(
    sim_extrap_gaussian(600, p = 4, sep = 3, seed = 1),
    exposure,
    tidyselect::starts_with("x")
  )

  expect_gt(
    mean(gap@results$gower_min),
    3 * mean(overlap@results$gower_min)
  )
  expect_lt(mean(gap@results$frac_nearby), 0.01)
  expect_gt(mean(overlap@results$frac_nearby), 0.03)
  expect_lt(mean(overlap@results$frac_nearby), 0.15)
})

test_that("frac_nearby falls and gower_min rises as the gap grows", {
  local_quiet()
  seps <- c(0, 1, 2, 3)
  stats <- lapply(seps, function(s) {
    res <- check_extrapolation(
      sim_extrap_gaussian(600, p = 4, sep = s, seed = 1),
      exposure,
      tidyselect::starts_with("x")
    )
    c(
      frac = mean(res@results$frac_nearby),
      gmin = mean(res@results$gower_min)
    )
  })

  frac <- vapply(stats, function(z) z[["frac"]], numeric(1))
  gmin <- vapply(stats, function(z) z[["gmin"]], numeric(1))
  expect_true(all(diff(frac) < 0))
  expect_true(all(diff(gmin) > 0))
})

test_that("frac_nearby is stable across n under overlap", {
  local_quiet()
  ns <- c(200, 800, 3200)
  fracs <- vapply(
    ns,
    function(nn) {
      res <- check_extrapolation(
        sim_extrap_gaussian(nn, p = 4, sep = 0, seed = 1),
        exposure,
        tidyselect::starts_with("x")
      )
      mean(res@results$frac_nearby)
    },
    numeric(1)
  )
  expect_true(all(fracs > 0.06 & fracs < 0.10))

  gap <- check_extrapolation(
    sim_extrap_gaussian(800, p = 4, sep = 3, seed = 1),
    exposure,
    tidyselect::starts_with("x")
  )
  expect_gt(fracs[2] - mean(gap@results$frac_nearby), 0.05)
})

test_that("geometric_variability sits in a loose band for Gaussian covariates", {
  local_quiet()
  for (p in c(2, 4, 6)) {
    res <- check_extrapolation(
      sim_extrap_gaussian(400, p = p, sep = 0, seed = 1),
      exposure,
      tidyselect::starts_with("x")
    )
    expect_gt(res@geometric_variability, 0.03)
    expect_lt(res@geometric_variability, 0.25)
  }
})

test_that("the support indicator separates overlap from a strong gap", {
  local_quiet()
  overlap <- check_extrapolation(
    sim_extrap_gaussian(400, p = 4, sep = 0, seed = 1),
    exposure,
    tidyselect::starts_with("x")
  )
  gap <- check_extrapolation(
    sim_extrap_gaussian(400, p = 4, sep = 3, seed = 1),
    exposure,
    tidyselect::starts_with("x")
  )

  expect_gt(prop_supported(overlap), 0.85)
  expect_lt(prop_supported(gap), 0.20)
})

test_that("a categorical-only gap is detected with numeric axes overlapping", {
  local_quiet()
  overlap <- check_extrapolation(
    sim_extrap_mixed(800, scenario = "overlap", seed = 1),
    exposure,
    c(x1, x2, g),
    hull = FALSE
  )
  gap <- check_extrapolation(
    sim_extrap_mixed(800, scenario = "categorical_gap", seed = 1),
    exposure,
    c(x1, x2, g),
    hull = FALSE
  )

  expect_lt(
    mean(gap@results$frac_nearby),
    mean(overlap@results$frac_nearby)
  )
  expect_gt(
    mean(gap@results$gower_min),
    mean(overlap@results$gower_min)
  )
})

# ---- Convex-hull membership -----------------------------------------------

test_that("hull membership is exact on a hand-built square reference", {
  local_quiet()
  # Controls (.id 1 to 4) are the corners of the square [0, 2] x [0, 2]. The
  # treated queries are an interior point, a coincident vertex, an edge
  # midpoint, and an exterior point.
  data <- tibble::tibble(
    exposure = c(0L, 0L, 0L, 0L, 1L, 1L, 1L, 1L),
    x1 = c(0, 2, 0, 2, 1, 0, 1, 5),
    x2 = c(0, 0, 2, 2, 1, 0, 0, 5)
  )
  res <- check_extrapolation(data, exposure, c(x1, x2), hull = TRUE)

  in_hull <- res@results$in_hull[match(5:8, res@results$.id)]
  expect_true(in_hull[1]) # interior
  expect_true(in_hull[2]) # coincident vertex
  expect_true(in_hull[3]) # edge midpoint
  expect_false(in_hull[4]) # exterior
})

test_that("the in-hull fraction collapses monotonically with dimension", {
  local_quiet()
  ps <- c(2, 3, 5, 15)
  props <- vapply(
    ps,
    function(p) {
      res <- suppressWarnings(check_extrapolation(
        sim_extrap_gaussian(400, p = p, sep = 0, seed = 1),
        exposure,
        tidyselect::starts_with("x"),
        hull = TRUE
      ))
      prop_in_hull(res)
    },
    numeric(1)
  )

  expect_gt(props[1], 0.80)
  expect_gt(props[2], 0.80)
  expect_lt(props[4], 0.05)
  expect_true(all(diff(props) < 0))
})

test_that("hull runs automatically for low p and reports membership", {
  local_quiet()
  data <- sim_extrap_gaussian(100, p = 4, sep = 0, seed = 1)
  res <- check_extrapolation(data, exposure, tidyselect::starts_with("x"))

  expect_true(res@hull_run)
  expect_false(all(is.na(res@results$in_hull)))
})

test_that("hull is skipped for high p and in_hull is NA", {
  local_quiet()
  data <- sim_extrap_gaussian(100, p = 13, sep = 0, seed = 1)
  res <- check_extrapolation(data, exposure, tidyselect::starts_with("x"))

  expect_false(res@hull_run)
  expect_true(all(is.na(res@results$in_hull)))
})

test_that("forcing hull off yields NA membership and hull_run FALSE", {
  local_quiet()
  data <- sim_extrap_gaussian(100, p = 4, sep = 0, seed = 1)
  res <- check_extrapolation(
    data,
    exposure,
    tidyselect::starts_with("x"),
    hull = FALSE
  )

  expect_false(res@hull_run)
  expect_true(all(is.na(res@results$in_hull)))
})

test_that("forcing hull on at high p warns but still runs", {
  local_quiet()
  data <- sim_extrap_gaussian(60, p = 13, sep = 0, seed = 1)
  expect_warning(
    check_extrapolation(
      data,
      exposure,
      tidyselect::starts_with("x"),
      hull = TRUE
    ),
    class = "positively_warning"
  )
  res <- suppressWarnings(check_extrapolation(
    data,
    exposure,
    tidyselect::starts_with("x"),
    hull = TRUE
  ))
  expect_true(res@hull_run)
})

test_that("forcing hull on with no numeric covariates warns and skips", {
  local_quiet()
  data <- sim_extrap_mixed(200, scenario = "overlap", seed = 1)
  expect_warning(
    check_extrapolation(data, exposure, g, hull = TRUE),
    class = "positively_warning"
  )
  res <- suppressWarnings(check_extrapolation(data, exposure, g, hull = TRUE))
  expect_false(res@hull_run)
  expect_true(all(is.na(res@results$in_hull)))
})

test_that("the auto hull is skipped when lpSolve is not installed", {
  local_quiet()
  data <- sim_extrap_gaussian(100, p = 4, sep = 0, seed = 1)
  testthat::local_mocked_bindings(
    is_installed = function(...) FALSE,
    .package = "rlang"
  )
  res <- check_extrapolation(data, exposure, tidyselect::starts_with("x"))
  expect_false(res@hull_run)
  expect_true(all(is.na(res@results$in_hull)))
})

test_that("forcing hull on without lpSolve aborts", {
  local_quiet()
  data <- sim_extrap_gaussian(100, p = 4, sep = 0, seed = 1)
  testthat::local_mocked_bindings(
    is_installed = function(...) FALSE,
    .package = "rlang"
  )
  expect_error(
    check_extrapolation(
      data,
      exposure,
      tidyselect::starts_with("x"),
      hull = TRUE
    ),
    class = "positively_error"
  )
})

test_that("all-constant covariates give gv zero and full nearby support", {
  local_quiet()
  data <- tibble::tibble(
    exposure = rep(c(0L, 1L), each = 5),
    x1 = 3,
    x2 = -1
  )
  res <- check_extrapolation(data, exposure, c(x1, x2), hull = FALSE)

  expect_equal(res@geometric_variability, 0)
  expect_true(all(res@results$gower_min == 0))
  expect_true(all(res@results$frac_nearby == 1))
})

test_that("the chunked large-n Gower path matches the unchunked result", {
  local_quiet()
  data <- sim_extrap_gaussian(120, p = 3, sep = 1, seed = 1)
  unchunked <- check_extrapolation(
    data,
    exposure,
    tidyselect::starts_with("x"),
    hull = FALSE
  )
  withr::local_options(positively.gower_chunk_threshold = 25)
  chunked <- check_extrapolation(
    data,
    exposure,
    tidyselect::starts_with("x"),
    hull = FALSE
  )

  expect_equal(chunked@results, unchunked@results)
  expect_equal(
    chunked@geometric_variability,
    unchunked@geometric_variability
  )
})

test_that("a large sample emits a chunking alert", {
  withr::local_options(positively.gower_chunk_threshold = 25)
  data <- sim_extrap_gaussian(60, p = 2, sep = 0, seed = 1)
  expect_message(
    check_extrapolation(
      data,
      exposure,
      tidyselect::starts_with("x"),
      hull = FALSE
    ),
    "row chunks"
  )
})

test_that("autoplot(type = \"hull\") aborts when the hull did not run", {
  local_quiet()
  data <- sim_extrap_gaussian(100, p = 4, sep = 0, seed = 1)
  res <- check_extrapolation(
    data,
    exposure,
    tidyselect::starts_with("x"),
    hull = FALSE
  )
  expect_error(
    ggplot2::autoplot(res, type = "hull"),
    class = "positively_error"
  )
})

# ---- Autoplot contract ----------------------------------------------------

test_that("autoplot() returns a ggplot for each type", {
  local_quiet()
  data <- sim_extrap_gaussian(200, p = 4, sep = 0, seed = 1)
  res <- check_extrapolation(
    data,
    exposure,
    tidyselect::starts_with("x"),
    hull = TRUE
  )

  expect_s3_class(ggplot2::autoplot(res, type = "distribution"), "ggplot")
  expect_s3_class(ggplot2::autoplot(res, type = "hull"), "ggplot")
})

# ---- Message and print snapshots ------------------------------------------

test_that("the binary-only abort messages are stable", {
  local_quiet()
  continuous <- dgp_continuous_support_gap(n = 100, seed = 1)
  categorical <- dgp_good_positivity(n = 200, seed = 1)
  categorical$exposure <- factor(
    sample(c("a", "b", "c"), nrow(categorical), replace = TRUE)
  )

  expect_snapshot(
    check_extrapolation(continuous, exposure, x1),
    error = TRUE
  )
  expect_snapshot(
    check_extrapolation(categorical, exposure, c(x1, x2)),
    error = TRUE
  )
})

test_that("the argument validation messages are stable", {
  local_quiet()
  data <- dgp_good_positivity(n = 100, seed = 1)

  expect_snapshot(check_extrapolation(1:10, exposure, x1), error = TRUE)
  expect_snapshot(
    check_extrapolation(data, exposure, tidyselect::starts_with("zzz")),
    error = TRUE
  )
  expect_snapshot(
    check_extrapolation(data, exposure, c(x1, x2), nearby = 0),
    error = TRUE
  )
  expect_snapshot(
    check_extrapolation(data, exposure, c(x1, x2), hull = "banana"),
    error = TRUE
  )
  expect_snapshot(
    check_extrapolation(data, c(exposure, x1), x2),
    error = TRUE
  )
})

test_that("the missing-value and sample-size abort messages are stable", {
  local_quiet()
  data <- dgp_good_positivity(n = 100, seed = 1)

  na_exposure <- data
  na_exposure$exposure[1] <- NA
  expect_snapshot(
    check_extrapolation(na_exposure, exposure, c(x1, x2)),
    error = TRUE
  )

  na_covariate <- data
  na_covariate$x1[1] <- NA
  expect_snapshot(
    check_extrapolation(na_covariate, exposure, c(x1, x2)),
    error = TRUE
  )

  too_small <- tibble::tibble(exposure = 1L, x1 = 0.5, x2 = -0.2)
  expect_snapshot(
    check_extrapolation(too_small, exposure, c(x1, x2)),
    error = TRUE
  )
})

test_that("the hull warning and error messages are stable", {
  local_quiet()
  mixed <- sim_extrap_mixed(200, scenario = "overlap", seed = 1)
  expect_snapshot(res <- check_extrapolation(mixed, exposure, g, hull = TRUE))

  gaussian <- sim_extrap_gaussian(100, p = 4, sep = 0, seed = 1)
  testthat::local_mocked_bindings(
    is_installed = function(...) FALSE,
    .package = "rlang"
  )
  expect_snapshot(
    check_extrapolation(
      gaussian,
      exposure,
      tidyselect::starts_with("x"),
      hull = TRUE
    ),
    error = TRUE
  )
})

test_that("the hull-view abort without a hull run is stable", {
  local_quiet()
  data <- sim_extrap_gaussian(100, p = 4, sep = 0, seed = 1)
  res <- check_extrapolation(
    data,
    exposure,
    tidyselect::starts_with("x"),
    hull = FALSE
  )
  expect_snapshot(ggplot2::autoplot(res, type = "hull"), error = TRUE)
})

test_that("the high-dimensional hull messages are stable", {
  data <- sim_extrap_gaussian(60, p = 13, sep = 0, seed = 1)

  expect_snapshot(
    check_extrapolation(data, exposure, tidyselect::starts_with("x"))
  )
  expect_snapshot(
    check_extrapolation(
      data,
      exposure,
      tidyselect::starts_with("x"),
      hull = TRUE
    )
  )
})

test_that("the print method is stable", {
  local_quiet()
  data <- sim_extrap_gaussian(200, p = 4, sep = 0, seed = 1)
  res <- check_extrapolation(
    data,
    exposure,
    tidyselect::starts_with("x"),
    hull = TRUE
  )
  expect_snapshot(print(res))
})
