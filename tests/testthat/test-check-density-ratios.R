# check_density_ratios() is bring-your-own-numbers: a numeric vector is one
# point-treatment time point, an n-by-T matrix is T time points whose cumulative
# products are summarized alongside the per-time ratios, and an lmtp fit is read
# for its density-ratio component. The tidy results tibble is long, with the
# fixed columns time / statistic / value. The statistic vocabulary these specs
# rely on:
#   per-time    : mean, quantile_50, quantile_90, quantile_95, quantile_99, max,
#                 prop_gt_10, prop_gt_50, ess, ess_fraction, prop_zero
#   cumulative  : the same names prefixed with cumulative_ (matrix input only)
# A point treatment occupies a single time point labelled 1. The maximum is the
# probs = 1 quantile and is reported under the name "max". ess_fraction is the
# Kish effective sample size divided by n.

gen_lognormal_ratios <- function(n, k, seed = 1) {
  withr::local_seed(seed)
  # The exact Gaussian shift-MTP density ratio: log r = k Z - k^2 / 2.
  exp(k * stats::rnorm(n) - k^2 / 2)
}

stat_value <- function(summ, statistic, time = 1) {
  summ$value[summ$statistic == statistic & summ$time == time]
}

fake_lmtp <- function(density_ratios) {
  # A structurally faithful stand-in for an lmtp fit: only the density-ratio
  # component is inspected, so the accessor must not require the lmtp package.
  structure(list(density_ratios = density_ratios), class = "lmtp")
}

# ---- Argument validation and dispatch -------------------------------------

test_that("check_density_ratios() rejects negative ratios", {
  expect_error(
    check_density_ratios(c(1, -0.5, 2)),
    class = "positively_error"
  )
})

test_that("check_density_ratios() rejects missing ratios", {
  expect_error(
    check_density_ratios(c(1, NA, 2)),
    class = "positively_error"
  )
})

test_that("check_density_ratios() rejects infinite ratios", {
  expect_error(
    check_density_ratios(c(1, Inf, 2)),
    class = "positively_error"
  )
  expect_error(
    check_density_ratios(matrix(c(1, 2, Inf, 3), nrow = 2)),
    class = "positively_error"
  )
  # -Inf is negative, so the negativity guard fires first and must keep firing.
  expect_error(
    check_density_ratios(c(1, -Inf)),
    class = "positively_range_error"
  )
})

test_that("check_density_ratios() rejects an empty numeric input", {
  expect_error(
    check_density_ratios(numeric(0)),
    class = "positively_error"
  )
})

test_that("check_density_ratios() rejects a non-numeric input", {
  expect_error(
    check_density_ratios(c("a", "b")),
    class = "positively_error"
  )
})

test_that("check_density_ratios() rejects a data frame input", {
  expect_error(
    check_density_ratios(tibble::tibble(a = 1:3)),
    class = "positively_error"
  )
})

test_that("check_density_ratios() rejects probs outside the unit interval", {
  expect_error(
    check_density_ratios(c(1, 1, 1), probs = c(0.5, 1.5)),
    class = "positively_error"
  )
})

test_that("check_density_ratios() rejects non-positive or non-numeric thresholds", {
  expect_error(
    check_density_ratios(c(1, 1, 1), thresholds = c(-1, 50)),
    class = "positively_error"
  )
  expect_error(
    check_density_ratios(c(1, 1, 1), thresholds = "10"),
    class = "positively_error"
  )
})

test_that("check_density_ratios() rejects a non-numeric matrix", {
  # A character matrix dispatches to the matrix method, so the non-numeric
  # guard fires inside validate_ratios rather than the default-method guard.
  m <- matrix(c("a", "b", "c", "d"), nrow = 2)
  expect_error(
    check_density_ratios(m),
    class = "positively_error"
  )
})

test_that("check_density_ratios() rejects empty or missing thresholds", {
  expect_error(
    check_density_ratios(c(1, 1, 1), thresholds = numeric(0)),
    class = "positively_error"
  )
  expect_error(
    check_density_ratios(c(1, 1, 1), thresholds = c(10, NA)),
    class = "positively_error"
  )
})

test_that("check_density_ratios() rejects a matrix with negative entries", {
  m <- matrix(c(1, 2, -1, 3), nrow = 2)
  expect_error(
    check_density_ratios(m),
    class = "positively_error"
  )
})

test_that("check_density_ratios() rejects a matrix with missing entries", {
  m <- matrix(c(1, 2, NA, 3), nrow = 2)
  expect_error(
    check_density_ratios(m),
    class = "positively_error"
  )
})

test_that("check_density_ratios() aborts on an lmtp fit without a ratio component", {
  malformed <- structure(list(other = 1), class = "lmtp")
  expect_error(
    check_density_ratios(malformed),
    class = "positively_error"
  )
})

# ---- Label collisions ------------------------------------------------------

test_that("thresholds that collapse to one statistic label are rejected", {
  # 1e-9 and 2e-9 are distinct values but round to the same eight-decimal label,
  # so both would name the statistic prop_gt_0.
  expect_snapshot_abort(
    check_density_ratios(c(0.5, 1, 2), thresholds = c(1e-9, 2e-9)),
    class = "positively_duplicate_error"
  )
})

test_that("probs that collapse to one quantile label are rejected", {
  # 0.5 and 0.5 + 1e-11 are distinct but both map to the quantile_50 label.
  expect_snapshot_abort(
    check_density_ratios(c(0.5, 1, 2), probs = c(0.5, 0.5 + 1e-11)),
    class = "positively_duplicate_error"
  )
})

# ---- Result class and properties ------------------------------------------

test_that("check_density_ratios() returns a density_ratios_result diagnostic", {
  res <- check_density_ratios(gen_lognormal_ratios(200, k = 1))

  expect_true(S7::S7_inherits(res, positivity_diagnostic))
  expect_identical(S7::S7_class(res)@name, "density_ratios_result")
  expect_identical(res@exposure_type, "continuous")
  expect_identical(res@n, 200L)
})

test_that("the raw per-time ratios are carried as a list property", {
  res <- check_density_ratios(gen_lognormal_ratios(50, k = 0.5))
  expect_type(res@ratios, "list")
  expect_length(res@ratios, 1)

  m <- matrix(gen_lognormal_ratios(60, k = 0.5), nrow = 20, ncol = 3)
  res_m <- check_density_ratios(m)
  expect_type(res_m@ratios, "list")
  expect_length(res_m@ratios, 3)
  expect_identical(res_m@n, 20L)
})

test_that("density_ratios_result has the fixed results columns", {
  res <- check_density_ratios(gen_lognormal_ratios(200, k = 1))

  expect_s3_class(res@results, "tbl_df")
  expect_setequal(names(res@results), c("time", "statistic", "value"))
  expect_type(res@results$statistic, "character")
  expect_type(res@results$value, "double")
})

test_that("a point treatment occupies a single time point", {
  res <- check_density_ratios(gen_lognormal_ratios(100, k = 1))
  expect_setequal(unique(res@results$time), 1)
})

test_that("a matrix reports one time point per column plus cumulative rows", {
  m <- matrix(gen_lognormal_ratios(300, k = 0.5), nrow = 100, ncol = 3)
  res <- check_density_ratios(m)

  expect_setequal(unique(res@results$time), 1:3)
  expect_true("cumulative_ess_fraction" %in% res@results$statistic)
  expect_true("ess_fraction" %in% res@results$statistic)
})

test_that("awkward probs give unique, precision-preserving statistic keys", {
  res <- check_density_ratios(
    gen_lognormal_ratios(200, k = 1),
    probs = c(0.123, 0.1234, 0.995, 1)
  )
  quantile_names <- res@results$statistic[
    grepl("^quantile_|^max$", res@results$statistic)
  ]
  expect_setequal(
    quantile_names,
    c("quantile_12.3", "quantile_12.34", "quantile_99.5", "max")
  )

  key <- paste(res@results$time, res@results$statistic)
  expect_identical(anyDuplicated(key), 0L)
})

# ---- tidy / glance --------------------------------------------------------

test_that("tidy() returns the results tibble", {
  res <- check_density_ratios(gen_lognormal_ratios(200, k = 1))

  expect_identical(generics::tidy(res), res@results)
})

test_that("glance() reports the headline weight statistics", {
  # Five ratios, two of them exactly zero. The Kish effective sample size is
  # 8^2 / 30, so its fraction is 64 / 150; the mean is 8 / 5 and two of the five
  # are zero. All four follow from the input alone.
  res <- check_density_ratios(c(0, 0, 1, 2, 5))
  glanced <- generics::glance(res)

  expect_s3_class(glanced, "tbl_df")
  expect_identical(nrow(glanced), 1L)
  expect_setequal(
    names(glanced),
    c("n", "n_times", "mean", "max", "prop_zero", "ess_fraction")
  )

  expect_identical(glanced$n, 5L)
  expect_identical(glanced$n_times, 1L)
  expect_equal(glanced$mean, 1.6)
  expect_equal(glanced$max, 5)
  expect_equal(glanced$prop_zero, 0.4)
  expect_equal(glanced$ess_fraction, 64 / 150)
})

test_that("glance() headline values are the final-time cumulative summaries", {
  # For a matrix the headline effective sample size and maximum are the
  # cumulative-product summaries at the last time point, and n_times counts the
  # columns.
  m <- matrix(gen_lognormal_ratios(300, k = 1, seed = 1), nrow = 100, ncol = 3)
  res <- check_density_ratios(m)
  glanced <- generics::glance(res)
  summ <- generics::tidy(res)

  expect_equal(
    glanced$ess_fraction,
    stat_value(summ, "cumulative_ess_fraction", time = 3)
  )
  expect_equal(glanced$max, stat_value(summ, "cumulative_max", time = 3))
  expect_equal(glanced$mean, stat_value(summ, "cumulative_mean", time = 3))
  expect_equal(
    glanced$prop_zero,
    stat_value(summ, "cumulative_prop_zero", time = 3)
  )
  expect_identical(glanced$n_times, 3L)

  # The headline reads the whole history, not the last wave on its own.
  cumulative <- Reduce(`*`, res@ratios)
  expect_equal(glanced$max, max(cumulative))
  expect_equal(glanced$mean, mean(cumulative))
})

test_that("glance() headline values match the point-treatment summaries", {
  # A point treatment is a single time point, so the headline maximum is the
  # maximum ratio and the headline ESS fraction is the per-time ESS fraction.
  ratios <- gen_lognormal_ratios(200, k = 1, seed = 2)
  res <- check_density_ratios(ratios)
  glanced <- generics::glance(res)
  summ <- generics::tidy(res)

  expect_equal(glanced$max, max(ratios))
  expect_equal(glanced$ess_fraction, stat_value(summ, "ess_fraction"))
})

test_that("summary() reports the headline weight statistics", {
  # The five ratios of the glance() block above: the mean is 8 / 5, two of the
  # five are zero, and the Kish effective sample size is 8^2 / 30, so its
  # fraction is 64 / 150.
  res <- check_density_ratios(c(0, 0, 1, 2, 5))
  summarized <- summary(res)

  # The sample size is metadata the container states once, but the number of
  # time points is a reading of the ratios supplied and stays.
  expect_type(summarized$value, "double")
  expect_equal(
    summarized,
    tibble::tibble(
      statistic = c("n_times", "mean", "max", "prop_zero", "ess_fraction"),
      value = c(1, 1.6, 5, 0.4, 64 / 150),
      # Weight concentration is reported as a magnitude. Nothing here was read
      # against a cut, so no row carries one.
      threshold = rep(NA_real_, 5)
    )
  )
})

test_that("summary() counts the time points of a multi-wave run", {
  m <- matrix(gen_lognormal_ratios(300, k = 1, seed = 1), nrow = 100, ncol = 3)
  summarized <- summary(check_density_ratios(m))

  expect_identical(
    summarized$statistic,
    c("n_times", "mean", "max", "prop_zero", "ess_fraction")
  )
  expect_identical(summarized$value[[1]], 3)
  expect_true(all(is.na(summarized$threshold)))
})

test_that("glance() stays one row when probs omits the maximum", {
  res <- check_density_ratios(
    gen_lognormal_ratios(200, k = 1),
    probs = c(0.5, 0.9)
  )
  glanced <- generics::glance(res)
  expect_s3_class(glanced, "tbl_df")
  expect_identical(nrow(glanced), 1L)
})

# ---- Exact identities ------------------------------------------------------

test_that("all-ones ratios give unit quantiles, zero exceedance, and full ESS", {
  res <- check_density_ratios(rep(1, 10))
  summ <- generics::tidy(res)

  expect_equal(stat_value(summ, "mean"), 1)
  expect_equal(stat_value(summ, "quantile_50"), 1)
  expect_equal(stat_value(summ, "quantile_95"), 1)
  expect_equal(stat_value(summ, "max"), 1)
  expect_equal(stat_value(summ, "prop_gt_10"), 0)
  expect_equal(stat_value(summ, "prop_gt_50"), 0)
  expect_equal(stat_value(summ, "ess"), 10)
  expect_equal(stat_value(summ, "ess_fraction"), 1)
})

test_that("Kish ESS matches hand-computed values and is scale invariant", {
  ten <- generics::tidy(check_density_ratios(rep(1, 10)))
  expect_equal(stat_value(ten, "ess"), 10)

  scaled <- generics::tidy(check_density_ratios(rep(2, 4)))
  expect_equal(stat_value(scaled, "ess"), 4)
  expect_equal(stat_value(scaled, "ess_fraction"), 1)

  outlier <- generics::tidy(check_density_ratios(c(1, 1, 1, 1, 100)))
  expect_equal(stat_value(outlier, "ess"), 1.081168, tolerance = 1e-5)
})

test_that("Kish ESS is scale invariant at extreme magnitudes", {
  # Equal weights give a full effective sample size at any overall scale, even
  # where the naive sum of squares would overflow or underflow to give NaN.
  huge <- generics::tidy(check_density_ratios(rep(1e160, 4)))
  expect_equal(stat_value(huge, "ess"), 4)
  expect_equal(stat_value(huge, "ess_fraction"), 1)

  tiny <- generics::tidy(check_density_ratios(rep(1e-200, 4)))
  expect_equal(stat_value(tiny, "ess"), 4)
  expect_equal(stat_value(tiny, "ess_fraction"), 1)
})

# ---- Lognormal generator behavior ----------------------------------------

test_that("the mean ratio stays at 1 for every severity", {
  for (k in c(0.25, 0.5, 1.0)) {
    summ <- generics::tidy(check_density_ratios(
      gen_lognormal_ratios(50000, k = k, seed = 1)
    ))
    expect_equal(stat_value(summ, "mean"), 1, tolerance = 0.03)
  }
})

test_that("quantiles, ESS fraction, and exceedance match lognormal targets at k = 1", {
  summ <- generics::tidy(check_density_ratios(
    gen_lognormal_ratios(50000, k = 1, seed = 1)
  ))

  expect_equal(stat_value(summ, "quantile_50"), 0.6065, tolerance = 0.05)
  expect_equal(stat_value(summ, "quantile_95"), 3.142, tolerance = 0.15)
  expect_equal(stat_value(summ, "quantile_99"), 6.211, tolerance = 0.6)
  expect_equal(stat_value(summ, "ess_fraction"), 0.368, tolerance = 0.03)
  # waldo treats tolerance relatively at this magnitude, so the target is
  # checked absolutely instead.
  expect_lt(abs(stat_value(summ, "prop_gt_10") - 0.0025), 0.002)
})

test_that("the tail summaries respond monotonically to severity", {
  ks <- c(0, 0.5, 1, 1.5, 2)
  summaries <- lapply(ks, function(k) {
    generics::tidy(check_density_ratios(
      gen_lognormal_ratios(50000, k = k, seed = 1)
    ))
  })

  # The 95th percentile of the shift-MTP ratio peaks near k = 1.6 and then
  # declines, so it is not monotone across this range; the 99th percentile, the
  # exceedance proportion, and the maximum are.
  q99 <- vapply(summaries, stat_value, numeric(1), statistic = "quantile_99")
  prop10 <- vapply(summaries, stat_value, numeric(1), statistic = "prop_gt_10")
  maxes <- vapply(summaries, stat_value, numeric(1), statistic = "max")
  ess_frac <- vapply(
    summaries,
    stat_value,
    numeric(1),
    statistic = "ess_fraction"
  )

  expect_true(all(diff(q99) >= 0))
  expect_true(all(diff(prop10) >= 0))
  expect_true(all(diff(maxes) >= 0))
  expect_true(all(diff(ess_frac) < 0))
  expect_lt(ess_frac[length(ess_frac)], 0.05)
})

# ---- Structural signature --------------------------------------------------

test_that("all-zero ratios give a zero Kish effective sample size", {
  # When every weight is zero the Kish denominator collapses, so the effective
  # sample size is defined to be zero rather than NaN.
  summ <- generics::tidy(check_density_ratios(rep(0, 10)))

  expect_equal(stat_value(summ, "ess"), 0)
  expect_equal(stat_value(summ, "ess_fraction"), 0)
  expect_equal(stat_value(summ, "prop_zero"), 1)
})

test_that("a structural violation shows a below-one mean and exact-zero mass", {
  # A leaked-support scenario: 30% of the intervention mass sits outside the
  # observed support, so those ratios are exactly zero and the mean drops below
  # the RN-derivative invariant of 1.
  ratios <- c(rep(0, 300), rep(1, 700))
  summ <- generics::tidy(check_density_ratios(ratios))

  expect_lt(stat_value(summ, "mean"), 1)
  expect_equal(stat_value(summ, "prop_zero"), 0.3)
})

# ---- Cumulative product across time ---------------------------------------

test_that("the cumulative product collapses ESS while per-time ESS stays flat", {
  withr::local_seed(1)
  n <- 50000
  times <- 6
  k <- 0.5
  m <- matrix(
    exp(k * stats::rnorm(n * times) - k^2 / 2),
    nrow = n,
    ncol = times
  )
  summ <- generics::tidy(check_density_ratios(m))

  cum_ess <- vapply(
    seq_len(times),
    function(t) stat_value(summ, "cumulative_ess_fraction", time = t),
    numeric(1)
  )
  per_time_ess <- vapply(
    seq_len(times),
    function(t) stat_value(summ, "ess_fraction", time = t),
    numeric(1)
  )
  cum_mean <- vapply(
    seq_len(times),
    function(t) stat_value(summ, "cumulative_mean", time = t),
    numeric(1)
  )
  cum_max <- vapply(
    seq_len(times),
    function(t) stat_value(summ, "cumulative_max", time = t),
    numeric(1)
  )

  expect_true(all(diff(cum_ess) < 0))
  expect_lt(cum_ess[times], 0.30)
  for (t in seq_len(times)) {
    expect_equal(cum_ess[t], exp(-t * k^2), tolerance = 0.03)
    expect_equal(per_time_ess[t], 0.779, tolerance = 0.03)
    expect_equal(cum_mean[t], 1, tolerance = 0.03)
  }
  # The running maximum of cumulative products is not monotone draw by draw, so
  # only its overall growth across time is asserted.
  expect_gt(cum_max[times], cum_max[1])
})

# ---- Cumulative product overflow ------------------------------------------

test_that("a cumulative product that overflows to Inf aborts with a range error", {
  # Two finite time points of 1e200 for the first row overflow the cumulative
  # product to 1e400, which is Inf in double precision, even though every
  # per-time ratio validate_ratios sees is finite. The overflow poisons the
  # cumulative ess (NaN) and quantiles (Inf) downstream, so the diagnostic must
  # abort at the point where the time-point product leaves double range.
  m <- matrix(c(1e200, 1, 1e200, 1), nrow = 2)
  expect_error(
    check_density_ratios(m),
    class = "positively_range_error"
  )
})

test_that("the cumulative-product overflow message is stable", {
  m <- matrix(c(1e200, 1, 1e200, 1), nrow = 2)
  expect_snapshot_abort(check_density_ratios(m))
})

# ---- lmtp accessor ---------------------------------------------------------

test_that("an lmtp fit is a pass-through to the matrix summaries", {
  m <- matrix(
    gen_lognormal_ratios(300, k = 0.5, seed = 2),
    nrow = 100,
    ncol = 3
  )
  from_matrix <- check_density_ratios(m)
  from_fit <- check_density_ratios(fake_lmtp(m))

  expect_identical(from_fit@results, from_matrix@results)
})

# ---- psw pass-through ------------------------------------------------------

test_that("a psw exposure is a pass-through to the numeric path", {
  skip_if_not_installed("propensity")
  ratios <- gen_lognormal_ratios(200, k = 1)
  from_psw <- check_density_ratios(propensity::psw(ratios))
  from_numeric <- check_density_ratios(ratios)
  expect_identical(from_psw@results, from_numeric@results)
})

# ---- Autoplot contract -----------------------------------------------------

test_that("autoplot() returns a ggplot for point and matrix inputs", {
  res_point <- check_density_ratios(gen_lognormal_ratios(200, k = 1))
  expect_s3_class(ggplot2::autoplot(res_point), "ggplot")

  m <- matrix(gen_lognormal_ratios(300, k = 0.5), nrow = 100, ncol = 3)
  res_matrix <- check_density_ratios(m)
  expect_s3_class(ggplot2::autoplot(res_matrix), "ggplot")
})

test_that("the cumulative view is a ggplot for a matrix and aborts for a point", {
  m <- matrix(gen_lognormal_ratios(300, k = 0.5), nrow = 100, ncol = 3)
  res_matrix <- check_density_ratios(m)
  expect_s3_class(
    ggplot2::autoplot(res_matrix, type = "cumulative"),
    "ggplot"
  )

  res_point <- check_density_ratios(gen_lognormal_ratios(200, k = 1))
  expect_error(
    ggplot2::autoplot(res_point, type = "cumulative"),
    class = "positively_error"
  )
})

test_that("autoplot() rejects an unknown type as a classed error", {
  res <- check_density_ratios(gen_lognormal_ratios(200, k = 1))

  # A view name is chosen from a fixed menu like any other argument, and the
  # call rendering a figure does not make its failure a lesser one.
  expect_error(
    ggplot2::autoplot(res, type = "bogus"),
    class = "positively_args_error"
  )
})

test_that("density-ratio autoplot views render as expected", {
  announce_doppelganger(
    "Density ratios point histogram",
    "Density ratios time boxplot",
    "Density ratios cumulative quantiles"
  )
  res_point <- check_density_ratios(gen_lognormal_ratios(200, k = 1))
  expect_doppelganger(
    "Density ratios point histogram",
    ggplot2::autoplot(res_point, type = "distribution")
  )

  m <- matrix(gen_lognormal_ratios(300, k = 0.5), nrow = 100, ncol = 3)
  res_matrix <- check_density_ratios(m)
  expect_doppelganger(
    "Density ratios time boxplot",
    ggplot2::autoplot(res_matrix, type = "distribution")
  )
  expect_doppelganger(
    "Density ratios cumulative quantiles",
    ggplot2::autoplot(res_matrix, type = "cumulative")
  )
})

test_that("plot() draws the view and returns the result invisibly", {
  local_null_device()
  res <- check_density_ratios(gen_lognormal_ratios(50, k = 1))
  expect_identical(plot(res), res)
})

# ---- Print method ----------------------------------------------------------

test_that("the print method is stable", {
  res <- check_density_ratios(gen_lognormal_ratios(200, k = 1))
  expect_snapshot(print(res))
})

test_that("the print method shows a nonzero proportion exactly zero", {
  # A structural violation leaves exact zeros in the ratios, the case the
  # "Proportion exactly zero" line exists for.
  res <- check_density_ratios(c(0, 0, 1, 1, 2, 5))
  expect_snapshot(print(res))
})

test_that("the multi-time print method is stable", {
  m <- matrix(gen_lognormal_ratios(300, k = 0.5), nrow = 100, ncol = 3)
  res <- check_density_ratios(m)
  expect_snapshot(print(res))
})

test_that("the percentile print labels use correct ordinals", {
  res <- check_density_ratios(
    gen_lognormal_ratios(200, k = 1),
    probs = c(0.01, 0.02, 0.03, 0.21, 0.995)
  )
  expect_snapshot(print(res))
})

test_that("the teens percentiles take a plain 'th' ordinal", {
  # The 11th, 12th, and 13th percentiles are the exception to the ordinal rule
  # and always take "th", not "st", "nd", or "rd".
  res <- check_density_ratios(
    gen_lognormal_ratios(200, k = 1),
    probs = c(0.11, 0.12, 0.13)
  )
  expect_snapshot(print(res))
})

test_that("print omits the quantile block when probs is the maximum alone", {
  # With probs = 1 there are no interior quantiles to display, so the printed
  # summary jumps from the mean straight to the maximum.
  res <- check_density_ratios(gen_lognormal_ratios(200, k = 1), probs = 1)
  expect_snapshot(print(res))
})

# ---- Error snapshots -------------------------------------------------------

test_that("the classed error messages are stable", {
  expect_snapshot_abort(check_density_ratios(c(1, -0.5, 2)))
  expect_snapshot_abort(check_density_ratios(c(1, NA, 2)))
  expect_snapshot_abort(check_density_ratios(numeric(0)))
  expect_snapshot_abort(check_density_ratios(c("a", "b")))
  expect_snapshot_abort(check_density_ratios(tibble::tibble(a = 1:3)))
  expect_snapshot_abort(check_density_ratios(c(1, 1, 1), probs = c(0.5, 1.5)))
  expect_snapshot_abort(check_density_ratios(c(1, 1, 1), probs = c(0.5, 0.5)))
  expect_snapshot_abort(check_density_ratios(
    c(1, 1, 1),
    thresholds = c(-1, 50)
  ))
  expect_snapshot_abort(check_density_ratios(c(1, 1, 1), thresholds = "10"))
  expect_snapshot_abort(check_density_ratios(
    c(1, 1, 1),
    thresholds = c(10, 10)
  ))
  malformed <- structure(list(other = 1), class = "lmtp")
  expect_snapshot_abort(check_density_ratios(malformed))

  res_point <- check_density_ratios(gen_lognormal_ratios(50, k = 1))
  expect_snapshot_abort(ggplot2::autoplot(res_point, type = "cumulative"))
  expect_snapshot_abort(check_density_ratios(
    c(1, 1, 1),
    thresholds = numeric(0)
  ))
  expect_snapshot_abort(check_density_ratios(
    c(1, 1, 1),
    thresholds = c(10, NA)
  ))
  expect_snapshot_abort(check_density_ratios(matrix(
    c("a", "b", "c", "d"),
    nrow = 2
  )))
  expect_snapshot_abort(check_density_ratios(c(1, Inf, 2)))
  expect_snapshot_abort(check_density_ratios(matrix(c(1, 2, Inf, 3), nrow = 2)))
})

# ---- Display methods -------------------------------------------------------

# Two hundred point-treatment ratios, so the reading covers a single time point.
density_ratios_display_result <- function() {
  check_density_ratios(gen_lognormal_ratios(200, k = 1))
}

test_that("diagnostic_label() names the ratio check without naming its class", {
  res <- density_ratios_display_result()
  label <- expect_readable_label(res)

  printed <- printed_text(res)
  expect_no_match(printed, S7::S7_class(res)@name, fixed = TRUE)
  expect_match(printed, label, fixed = TRUE)
})

test_that("diagnostic_headline() reads the ratios it summarizes", {
  res <- density_ratios_display_result()
  headline <- expect_readable_headline(res)
  text <- rendered_text(headline)

  expect_match(text, "\\b200\\b")
})

test_that("the density ratio label and headline are stable", {
  res <- density_ratios_display_result()
  expect_snapshot({
    diagnostic_label(res)
    writeLines(diagnostic_headline(res))
  })
})

test_that("diagnostic_headline() marks the multi-time statistics as cumulative", {
  # A matrix holds one ratio per observation per time point, and the headline
  # statistics are the cumulative-product values at the final one, so neither
  # the count nor the statistics describe the whole set of ratios.
  m <- matrix(gen_lognormal_ratios(300, k = 0.5), nrow = 100, ncol = 3)
  res <- check_density_ratios(m)
  headline <- expect_readable_headline(res)
  text <- rendered_text(headline)

  expect_match(text, "\\b100\\b")
  expect_match(text, "cumulative", fixed = TRUE)
  # The plural follows the number of time points rather than the ratio count.
  expect_match(text, "time points", fixed = TRUE)
})

test_that("the multi-time density ratio headline is stable", {
  m <- matrix(gen_lognormal_ratios(300, k = 0.5), nrow = 100, ncol = 3)
  res <- check_density_ratios(m)
  expect_snapshot(writeLines(diagnostic_headline(res)))
})
