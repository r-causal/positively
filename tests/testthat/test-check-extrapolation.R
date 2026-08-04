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

# Fraction of units with an opposite-group neighbor inside the geometric
# variability radius, the recommended binary support indicator.
prop_supported <- function(res) {
  mean(res@results$gower_min <= res@geometric_variability)
}

prop_in_hull <- function(res) {
  mean(res@results$in_hull)
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
  data <- dgp_categorical(n = 200, seed = 1)
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

test_that("check_extrapolation() rejects non-finite covariate values", {
  local_quiet()
  data <- withr::with_seed(1, {
    data.frame(
      exposure = rep(0:1, each = 5),
      x1 = c(stats::rnorm(9), Inf),
      x2 = stats::rnorm(10)
    )
  })
  expect_error(
    check_extrapolation(data, exposure, c(x1, x2), hull = FALSE),
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

# ---- Declared exposure type ------------------------------------------------

test_that("a declared binary type reproduces the auto path exactly", {
  local_quiet()
  data <- sim_extrap_gaussian(200, p = 4, sep = 0, seed = 1)
  auto <- check_extrapolation(data, exposure, tidyselect::starts_with("x"))
  declared <- check_extrapolation(
    data,
    exposure,
    tidyselect::starts_with("x"),
    exposure_type = "binary"
  )

  # On a genuinely binary exposure the declaration can only agree with
  # detection, so every reported quantity must match, not merely the fact that
  # both calls returned.
  expect_identical(declared@exposure_type, auto@exposure_type)
  expect_identical(declared@results, auto@results)
  expect_identical(
    declared@geometric_variability,
    auto@geometric_variability
  )
  expect_identical(declared@hull_run, auto@hull_run)
  expect_identical(declared@n, auto@n)
})

test_that("a declared binary type skips the detection alert", {
  data <- sim_extrap_gaussian(200, p = 4, sep = 0, seed = 1)
  # This test observes messages rather than muffling them, so the hull is turned
  # off explicitly: the automatic hull announces its own skip when lpSolve is not
  # installed, which is a supported configuration and unrelated to detection.
  # The auto path announces what it inferred. A declared type is taken as given,
  # so there is nothing to announce.
  expect_message(
    check_extrapolation(
      data,
      exposure,
      tidyselect::starts_with("x"),
      hull = FALSE
    ),
    "Treating `.exposure` as binary"
  )
  expect_no_message(
    check_extrapolation(
      data,
      exposure,
      tidyselect::starts_with("x"),
      hull = FALSE,
      exposure_type = "binary"
    )
  )
})

test_that("a declared binary type aborts on a three-level exposure", {
  local_quiet()
  data <- sim_extrap_gaussian(60, p = 2, sep = 0, seed = 1)
  data$exposure <- factor(rep(c("a", "b", "c"), length.out = nrow(data)))

  # The opposite-group comparison is a not-equal test against the current unit's
  # group, so three levels would quietly become one against all others: a
  # finite, well-formed, and meaningless answer. Declaring the type must not
  # buy passage to it.
  err <- expect_error(
    check_extrapolation(
      data,
      exposure,
      c(x1, x2),
      exposure_type = "binary"
    ),
    class = "positively_exposure_type_error"
  )
  expect_match(conditionMessage(err), "two distinct values", fixed = TRUE)
})

test_that("a constant exposure aborts on both the auto and declared paths", {
  local_quiet()
  data <- sim_extrap_gaussian(60, p = 2, sep = 0, seed = 1)
  data$exposure <- factor(rep("a", nrow(data)))

  # A one-level factor is not two distinct values, yet detection reads it as
  # binary, so detection alone can never reject it. The opposite group is then
  # empty and every per-unit statistic is taken over nothing. Only the
  # structural requirement of the resolved type rules this out, and it rules it
  # out whether the type was declared or inferred.
  expect_identical(detect_exposure_type(data$exposure), "binary")

  expect_error(
    check_extrapolation(data, exposure, c(x1, x2)),
    class = "positively_exposure_type_error"
  )
  expect_error(
    check_extrapolation(
      data,
      exposure,
      c(x1, x2),
      exposure_type = "binary"
    ),
    class = "positively_exposure_type_error"
  )
})

test_that("check_extrapolation() rejects a type outside its supported menu", {
  local_quiet()
  data <- sim_extrap_gaussian(60, p = 2, sep = 0, seed = 1)
  err <- expect_error(
    check_extrapolation(
      data,
      exposure,
      c(x1, x2),
      exposure_type = "continuous"
    ),
    class = "positively_args_error"
  )
  expect_match(
    conditionMessage(err),
    '`exposure_type` must be one of "auto" or "binary", not "continuous".',
    fixed = TRUE
  )
})

# ---- Unsupported covariate types ------------------------------------------

test_that("check_extrapolation() rejects a Date covariate", {
  local_quiet()
  data <- withr::with_seed(1, {
    tibble::tibble(
      exposure = rep(0:1, each = 50),
      x1 = stats::rnorm(100),
      d = as.Date("2020-01-01") + seq_len(100)
    )
  })
  expect_snapshot_abort(check_extrapolation(
    data,
    exposure,
    c(x1, d),
    hull = FALSE
  ))
})

test_that("check_extrapolation() rejects a list-column covariate", {
  local_quiet()
  data <- withr::with_seed(1, {
    tibble::tibble(
      exposure = rep(0:1, each = 50),
      x1 = stats::rnorm(100),
      lc = as.list(seq_len(100))
    )
  })
  expect_error(
    check_extrapolation(data, exposure, c(x1, lc), hull = FALSE),
    class = "positively_error"
  )
})

test_that("a logical covariate returns finite geometric variability and results", {
  local_quiet()
  data <- withr::with_seed(1, {
    tibble::tibble(
      exposure = rep(0:1, each = 50),
      x1 = stats::rnorm(100),
      flag = sample(c(TRUE, FALSE), 100, replace = TRUE)
    )
  })
  res <- check_extrapolation(data, exposure, c(x1, flag), hull = FALSE)

  expect_true(is.finite(res@geometric_variability))
  expect_true(all(is.finite(res@results$gower_min)))
  expect_true(all(is.finite(res@results$gower_mean)))
  expect_true(all(is.finite(res@results$frac_nearby)))
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
  skip_if_not_installed("lpSolve")
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
    c(
      ".id",
      "exposure",
      "frac_nearby",
      "gower_min",
      "gower_mean",
      "in_hull",
      "low_support"
    )
  )
  expect_type(res@results$frac_nearby, "double")
  expect_type(res@results$gower_min, "double")
  expect_type(res@results$gower_mean, "double")
  expect_type(res@results$in_hull, "logical")
  expect_type(res@results$low_support, "logical")
  expect_false(anyNA(res@results$low_support))
  expect_identical(nrow(res@results), 200L)
})

test_that("extrapolation carries a per-unit low_support column", {
  local_quiet()
  # The column has to be there whether or not the hull ran, since the hull is
  # optional and skipped outright above hull_max_p. The hull-run shape is
  # covered by the fixed-columns test above.
  data <- sim_extrap_gaussian(150, p = 4, sep = 0, seed = 1)
  res <- check_extrapolation(
    data,
    exposure,
    tidyselect::starts_with("x"),
    hull = FALSE
  )

  expect_true("low_support" %in% names(res@results))
  expect_type(res@results$low_support, "logical")
  expect_false(anyNA(res@results$low_support))
  expect_length(res@results$low_support, 150L)

  # The column means "the nearest opposite-group unit lies farther away than one
  # geometric variability", and every prop_supported in the package is its
  # complement. Pin that sign directly. The comparison is exact rather than
  # tolerant because the two sides are complements of a single inequality.
  supported <- res@results$gower_min <= res@geometric_variability
  expect_identical(res@results$low_support, !supported)
  expect_equal(mean(res@results$low_support), 1 - glance(res)$prop_supported)
})

test_that("summary() returns one row per exposure level", {
  local_quiet()
  data <- sim_extrap_gaussian(200, p = 4, sep = 0, seed = 1)
  res <- check_extrapolation(data, exposure, tidyselect::starts_with("x"))
  summarized <- summary(res)

  expect_s3_class(summarized, "tbl_df")
  expect_identical(nrow(summarized), 2L)
  expect_identical(
    names(summarized),
    c("exposure", "n", "mean_gower_min", "prop_supported", "prop_in_hull")
  )

  # The summary derives prop_supported from low_support instead of repeating
  # the Gower comparison, so pin the per-group values against that comparison.
  supported <- res@results$gower_min <= res@geometric_variability
  expect_equal(
    summarized$prop_supported,
    as.numeric(tapply(supported, res@results$exposure, mean))
  )
})

test_that("the computed summary property is gone", {
  local_quiet()
  data <- sim_extrap_gaussian(200, p = 4, sep = 0, seed = 1)
  res <- check_extrapolation(data, exposure, tidyselect::starts_with("x"))

  # The per-group aggregate is a summary() method now, so `@` access is not
  # taught for it and the property no longer exists to be reached for.
  expect_false("summary" %in% S7::prop_names(res))
  expect_error(res@summary, "Can't find property")
})

test_that("summary() drops mean_frac_nearby, whose group means are forced equal", {
  local_quiet()
  data <- sim_extrap_gaussian(200, p = 4, sep = 0, seed = 1)
  res <- check_extrapolation(data, exposure, tidyselect::starts_with("x"))
  results <- res@results
  summarized <- summary(res)

  expect_s3_class(summarized, "tbl_df")
  expect_false("mean_frac_nearby" %in% names(summarized))

  # frac_nearby counts the fraction of the opposite group within the radius, so
  # the group-0 mean is cross_pairs / (n0 * n1) and the group-1 mean is the same
  # quantity with its factors swapped. The equality is an identity rather than a
  # property of this data, which is why the column carries no group-level
  # information. gower_min, computed over the same pairs, does separate the
  # groups.
  frac_means <- tapply(results$frac_nearby, results$exposure, mean)
  expect_equal(frac_means[[1]], frac_means[[2]])
  gower_means <- tapply(results$gower_min, results$exposure, mean)
  gower_gap <- abs(gower_means[[1]] - gower_means[[2]])
  expect_gt(gower_gap, 1e-6)

  # The per-unit column stays, because unit by unit it is informative.
  expect_true("frac_nearby" %in% names(results))
})

test_that("summary() aggregates the hand-computed example by exposure group", {
  local_quiet()
  # The three-unit example used for glance() below. The geometric variability is
  # 2/9 and the Gower distances are d(1, 2) = 1/6, d(1, 3) = 5/6, d(2, 3) = 1.
  # The single treated unit sits 1/6 from its nearest control, so group 1 is
  # fully supported and its mean nearest distance is 1/6. Of the two controls,
  # unit 2 is 1/6 from the treated unit and unit 3 is 5/6 away, which exceeds one
  # geometric variability, so group 0 is half supported and its mean nearest
  # distance is (1/6 + 5/6) / 2.
  data <- tibble::tibble(
    exposure = c(1L, 0L, 0L),
    x1 = c(1, 0, 2),
    x2 = c(0, 0, 2),
    g = factor(c("a", "a", "b"))
  )
  res <- check_extrapolation(data, exposure, c(x1, x2, g), hull = FALSE)

  expect_equal(
    summary(res),
    tibble::tibble(
      exposure = c(0L, 1L),
      n = c(2L, 1L),
      mean_gower_min = c(0.5, 1 / 6),
      prop_supported = c(0.5, 1),
      # The hull test did not run, so the fraction inside it is unknown rather
      # than zero, in the per-group summary as in glance().
      prop_in_hull = c(NA_real_, NA_real_)
    )
  )
})

test_that("tidy() returns the results tibble", {
  local_quiet()
  data <- sim_extrap_gaussian(200, p = 4, sep = 0, seed = 1)
  res <- check_extrapolation(data, exposure, tidyselect::starts_with("x"))

  expect_identical(generics::tidy(res), res@results)
})

test_that("glance() matches the hand-computed extrapolation statistics", {
  local_quiet()
  # Pooled ranges are 2 on x1 and 2 on x2, so the three Gower distances are
  # d(1, 2) = (0.5 + 0 + 0) / 3, d(1, 3) = (0.5 + 1 + 1) / 3, and
  # d(2, 3) = (1 + 1 + 1) / 3. Half the mean of the full matrix, diagonal
  # included, is (2 * (1/6 + 5/6 + 1) / 9) / 2 = 2/9. Unit 1 has one of the two
  # controls inside that radius, unit 2 has the one treated unit inside it, and
  # unit 3 has none, so frac_nearby is 0.5, 1, 0 and two of three units are
  # supported.
  data <- tibble::tibble(
    exposure = c(1L, 0L, 0L),
    x1 = c(1, 0, 2),
    x2 = c(0, 0, 2),
    g = factor(c("a", "a", "b"))
  )
  res <- check_extrapolation(data, exposure, c(x1, x2, g), hull = FALSE)
  glanced <- generics::glance(res)

  expect_s3_class(glanced, "tbl_df")
  expect_identical(nrow(glanced), 1L)
  expect_setequal(
    names(glanced),
    c(
      "n",
      "nearby",
      "geometric_variability",
      "mean_frac_nearby",
      "prop_supported",
      "hull_run",
      "prop_in_hull"
    )
  )

  expect_identical(glanced$n, 3L)
  expect_identical(glanced$nearby, 1)
  expect_equal(glanced$geometric_variability, 2 / 9)
  expect_equal(glanced$mean_frac_nearby, 0.5)
  expect_equal(glanced$prop_supported, 2 / 3)

  # The hull test did not run, so its fraction is unknown rather than zero, and
  # the flag that says so keeps its type.
  expect_type(glanced$hull_run, "logical")
  expect_false(glanced$hull_run)
  expect_type(glanced$prop_in_hull, "double")
  expect_true(is.na(glanced$prop_in_hull))
})

test_that("glance() reports the hull fraction once the hull test runs", {
  local_quiet()
  skip_if_not_installed("lpSolve")
  data <- sim_extrap_gaussian(200, p = 4, sep = 0, seed = 1)
  res <- check_extrapolation(data, exposure, tidyselect::starts_with("x"))
  glanced <- generics::glance(res)

  expect_true(glanced$hull_run)
  expect_equal(glanced$prop_in_hull, prop_in_hull(res))

  # The two readings answer different questions and disagree here, which is the
  # whole reason both are reported.
  expect_equal(glanced$prop_supported, prop_supported(res))
  expect_gt(glanced$prop_supported, glanced$prop_in_hull)
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

test_that("frac_nearby grows with the nearby radius", {
  local_quiet()
  # A wider nearby radius counts more opposite-group neighbors as support, so the
  # mean frac_nearby rises with it. The Gower distances themselves do not depend
  # on the radius, so gower_min and gower_mean are identical across the sweep.
  # low_support is fixed at one geometric variability rather than the nearby
  # radius, so it is identical across the sweep as well.
  data <- sim_extrap_gaussian(100, p = 2, sep = 0, seed = 1)
  results <- lapply(c(0.5, 1, 2), function(nb) {
    check_extrapolation(
      data,
      exposure,
      tidyselect::starts_with("x"),
      nearby = nb,
      hull = FALSE
    )@results
  })

  fracs <- vapply(results, function(r) mean(r$frac_nearby), numeric(1))
  expect_true(all(diff(fracs) > 0))

  expect_identical(results[[2]]$gower_min, results[[1]]$gower_min)
  expect_identical(results[[3]]$gower_min, results[[1]]$gower_min)
  expect_identical(results[[2]]$gower_mean, results[[1]]$gower_mean)
  expect_identical(results[[3]]$gower_mean, results[[1]]$gower_mean)
  expect_identical(results[[2]]$low_support, results[[1]]$low_support)
  expect_identical(results[[3]]$low_support, results[[1]]$low_support)
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
  skip_if_not_installed("lpSolve")
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

test_that("hull membership is exact on the square reference under a large offset", {
  local_quiet()
  skip_if_not_installed("lpSolve")
  # The same square reference and treated queries, translated by 1e6 on both
  # axes. Membership is invariant to translation, so the four verdicts must not
  # change when the coordinates carry a large constant offset.
  data <- tibble::tibble(
    exposure = c(0L, 0L, 0L, 0L, 1L, 1L, 1L, 1L),
    x1 = c(0, 2, 0, 2, 1, 0, 1, 5) + 1e6,
    x2 = c(0, 0, 2, 2, 1, 0, 0, 5) + 1e6
  )
  res <- check_extrapolation(data, exposure, c(x1, x2), hull = TRUE)

  in_hull <- res@results$in_hull[match(5:8, res@results$.id)]
  expect_true(in_hull[1]) # interior
  expect_true(in_hull[2]) # coincident vertex
  expect_true(in_hull[3]) # edge midpoint
  expect_false(in_hull[4]) # exterior
})

test_that("hull membership is invariant to a large offset on a Gaussian design", {
  local_quiet()
  skip_if_not_installed("lpSolve")
  data <- sim_extrap_gaussian(100, p = 2, sep = 0, seed = 1)
  base <- check_extrapolation(
    data,
    exposure,
    tidyselect::starts_with("x"),
    hull = TRUE
  )

  shifted <- data
  shifted$x1 <- shifted$x1 + 1e6
  shifted$x2 <- shifted$x2 + 1e6
  offset <- check_extrapolation(
    shifted,
    exposure,
    tidyselect::starts_with("x"),
    hull = TRUE
  )

  expect_identical(offset@results$in_hull, base@results$in_hull)
})

test_that("the in-hull fraction collapses monotonically with dimension", {
  local_quiet()
  skip_if_not_installed("lpSolve")
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
  skip_if_not_installed("lpSolve")
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
  skip_if_not_installed("lpSolve")
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
  # The degenerate case where the strict inequality behind low_support carries
  # its weight: a distance of zero is within a radius of zero, so no row is
  # unsupported even though gower_min and the geometric variability are equal.
  expect_false(any(res@results$low_support))
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
    expect_message(
      check_extrapolation(
        data,
        exposure,
        tidyselect::starts_with("x"),
        hull = FALSE
      ),
      "row chunks"
    ),
    "Treating `.exposure` as binary"
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
  skip_if_not_installed("lpSolve")
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

test_that("autoplot() rejects an unknown type as a classed error", {
  local_quiet()
  data <- sim_extrap_gaussian(100, p = 4, sep = 0, seed = 1)
  res <- check_extrapolation(data, exposure, tidyselect::starts_with("x"))

  # A view name is chosen from a fixed menu like any other argument, and the
  # call rendering a figure does not make its failure a lesser one.
  expect_error(
    ggplot2::autoplot(res, type = "bogus"),
    class = "positively_args_error"
  )
})

test_that("plot() draws the view and returns the result invisibly", {
  local_quiet()
  local_null_device()
  data <- sim_extrap_gaussian(100, p = 4, sep = 0, seed = 1)
  res <- check_extrapolation(
    data,
    exposure,
    tidyselect::starts_with("x"),
    hull = FALSE
  )
  expect_identical(plot(res), res)
})

test_that("the extrapolation distribution view renders as expected", {
  local_quiet()
  # The distribution view reads only `frac_nearby` and `exposure`, so it renders
  # identically whether or not the hull ran, and this half needs no lpSolve.
  data <- sim_extrap_gaussian(200, p = 4, sep = 0, seed = 1)
  res <- check_extrapolation(
    data,
    exposure,
    tidyselect::starts_with("x"),
    hull = FALSE
  )

  expect_doppelganger(
    "Extrapolation frac nearby distribution",
    ggplot2::autoplot(res, type = "distribution")
  )
})

test_that("the extrapolation hull view renders as expected", {
  local_quiet()
  announce_doppelganger("Extrapolation convex hull membership")
  skip_if_not_installed("lpSolve")
  data <- sim_extrap_gaussian(200, p = 4, sep = 0, seed = 1)
  res <- check_extrapolation(
    data,
    exposure,
    tidyselect::starts_with("x"),
    hull = TRUE
  )

  expect_doppelganger(
    "Extrapolation convex hull membership",
    ggplot2::autoplot(res, type = "hull")
  )
})

# ---- Message and print snapshots ------------------------------------------

test_that("the binary-only abort messages are stable", {
  local_quiet()
  continuous <- dgp_continuous_support_gap(n = 100, seed = 1)
  categorical <- dgp_categorical(n = 200, seed = 1)

  expect_snapshot_abort(check_extrapolation(continuous, exposure, x1))
  expect_snapshot_abort(check_extrapolation(categorical, exposure, c(x1, x2)))
})

test_that("the argument validation messages are stable", {
  local_quiet()
  data <- dgp_good_positivity(n = 100, seed = 1)

  expect_snapshot_abort(check_extrapolation(1:10, exposure, x1))
  expect_snapshot_abort(check_extrapolation(
    data,
    exposure,
    tidyselect::starts_with("zzz")
  ))
  expect_snapshot_abort(check_extrapolation(
    data,
    exposure,
    c(x1, x2),
    nearby = 0
  ))
  expect_snapshot_abort(check_extrapolation(
    data,
    exposure,
    c(x1, x2),
    hull = "banana"
  ))
  expect_snapshot_abort(check_extrapolation(data, c(exposure, x1), x2))
})

test_that("the missing-value and sample-size abort messages are stable", {
  local_quiet()
  data <- dgp_good_positivity(n = 100, seed = 1)

  na_exposure <- data
  na_exposure$exposure[1] <- NA
  expect_snapshot_abort(check_extrapolation(na_exposure, exposure, c(x1, x2)))

  na_covariate <- data
  na_covariate$x1[1] <- NA
  expect_snapshot_abort(check_extrapolation(na_covariate, exposure, c(x1, x2)))

  too_small <- tibble::tibble(exposure = 1L, x1 = 0.5, x2 = -0.2)
  expect_snapshot_abort(check_extrapolation(too_small, exposure, c(x1, x2)))

  non_finite <- withr::with_seed(1, {
    data.frame(
      exposure = rep(0:1, each = 5),
      x1 = c(stats::rnorm(9), Inf),
      x2 = stats::rnorm(10)
    )
  })
  expect_snapshot_abort(check_extrapolation(
    non_finite,
    exposure,
    c(x1, x2),
    hull = FALSE
  ))
})

test_that("the hull warning and error messages are stable", {
  local_quiet()
  withr::local_options(warn = 0)
  mixed <- sim_extrap_mixed(200, scenario = "overlap", seed = 1)
  expect_snapshot(res <- check_extrapolation(mixed, exposure, g, hull = TRUE))

  gaussian <- sim_extrap_gaussian(100, p = 4, sep = 0, seed = 1)
  testthat::local_mocked_bindings(
    is_installed = function(...) FALSE,
    .package = "rlang"
  )
  expect_snapshot_abort(check_extrapolation(
    gaussian,
    exposure,
    tidyselect::starts_with("x"),
    hull = TRUE
  ))
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
  expect_snapshot_abort(ggplot2::autoplot(res, type = "hull"))
})

test_that("the hull-view abort points to the numeric-covariate precondition when hull was skipped", {
  local_quiet()
  data <- withr::with_seed(1, {
    tibble::tibble(
      exposure = rep(0:1, each = 40),
      g1 = factor(sample(c("a", "b"), 80, replace = TRUE)),
      g2 = factor(sample(c("x", "y"), 80, replace = TRUE))
    )
  })
  # With no numeric covariate the hull test cannot run even though hull = TRUE
  # was requested, so it is skipped with a warning and hull_run stays FALSE.
  expect_warning(
    check_extrapolation(data, exposure, c(g1, g2), hull = TRUE),
    class = "positively_hull_covariate_warning"
  )
  res <- suppressWarnings(
    check_extrapolation(data, exposure, c(g1, g2), hull = TRUE)
  )
  expect_false(res@hull_run)

  # Rerunning with hull = TRUE would change nothing; the abort must name the
  # numeric-covariate precondition rather than repeat the request.
  expect_error(
    ggplot2::autoplot(res, type = "hull"),
    regexp = "numeric covariate"
  )
  expect_snapshot_abort(ggplot2::autoplot(res, type = "hull"))
})

test_that("the hull-view abort names the numeric precondition when hull was left unset", {
  local_quiet()
  data <- withr::with_seed(1, {
    tibble::tibble(
      exposure = rep(0:1, each = 40),
      g1 = factor(sample(c("a", "b"), 80, replace = TRUE)),
      g2 = factor(sample(c("x", "y"), 80, replace = TRUE))
    )
  })
  # With hull left at its NULL default and no numeric covariate, the auto path
  # skips the hull test silently: hull was not requested, so no warning fires.
  res <- expect_no_warning(check_extrapolation(data, exposure, c(g1, g2)))
  expect_false(res@hull_run)
  expect_null(res@params$hull)

  # A rerun with hull = TRUE cannot help here either, since there is still no
  # numeric covariate, so the abort must name that precondition rather than
  # advise a rerun.
  expect_error(
    ggplot2::autoplot(res, type = "hull"),
    regexp = "numeric covariate"
  )
  expect_snapshot_abort(ggplot2::autoplot(res, type = "hull"))
})

test_that("the high-dimensional hull messages are stable", {
  withr::local_options(warn = 0)
  skip_if_not_installed("lpSolve")
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
  skip_if_not_installed("lpSolve")
  data <- sim_extrap_gaussian(200, p = 4, sep = 0, seed = 1)
  res <- check_extrapolation(
    data,
    exposure,
    tidyselect::starts_with("x"),
    hull = TRUE
  )
  expect_snapshot(print(res))
})
