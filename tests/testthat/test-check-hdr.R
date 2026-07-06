# The HDR non-overlap ratio of Bao and Schomaker (2025). For a continuous
# exposure A with covariates L and a support level `mass`, the HDR at profile l
# is A_alpha(l) = { a : f(a | l) >= f_alpha(l) } with mass probability, and the
# non-overlap ratio tau_hat(a) = (1/n) sum_j 1{ a not in A_alpha(l_j) } is the
# fraction of covariate profiles for which the common target a is unsupported.
# Under the default hdr_density_normal() estimator, membership reduces to the
# interval test abs(a - mu_hat(l)) <= z * sigma_hat with z = qnorm((1 + mass) /
# 2), so tau_hat(a) is the fraction of fitted means more than z * sigma_hat from
# a. The fixed results columns are value (a) and nonoverlap (tau_hat), plus
# `time` for the sequential variant. Magnitude claims are anchored to the
# closed-form oracle tau(a) = Phi((a - z*sigma)/s_mu) + 1 - Phi((a + z*sigma)/
# s_mu) with s_mu = |beta| for the linear-Gaussian design, per the sim REPORT.
#
# Key caveat carried from the REPORT: the normal estimator detects mean-shift
# support gaps, not multimodal gaps, and the non-overlap floor is not
# universally zero under strong covariate dependence.

local_quiet <- function(.env = parent.frame()) {
  withr::local_options(positively.quiet = TRUE, .local_envir = .env)
}

# ---- Scenario generators --------------------------------------------------

# Linear-Gaussian design: L ~ N(0, 1), A | L ~ N(beta * L, sigma^2). The default
# estimator recovers mu_hat(l) = beta * l and sigma_hat = sigma, so the whole
# tau curve has a closed form. One seed across beta holds L fixed.
sim_hdr_linear <- function(n, beta = 1, sigma = 1, seed = 1) {
  withr::local_seed(seed)
  l <- stats::rnorm(n)
  a <- stats::rnorm(n, mean = beta * l, sd = sigma)
  tibble::tibble(exposure = a, l = l)
}

# Two-group mixture with a known unsupported fraction at the target a = 0. A
# fraction p sits far from 0 (mean 5) and the rest sit at 0, both with small
# within-group SD, so lm(exposure ~ g) recovers two means and tau_hat(0) = p.
sim_hdr_mixture <- function(n, p = 0.3, seed = 1) {
  withr::local_seed(seed)
  far <- stats::runif(n) < p
  a <- ifelse(far, stats::rnorm(n, 5, 0.3), stats::rnorm(n, 0, 0.3))
  tibble::tibble(exposure = a, g = as.numeric(far))
}

# Bimodal truth with a genuine support hole at a = 0: A = 5 * B + N(0, 0.5^2),
# B in {-1, +1}, modes at -5 and +5. The covariate x carries no information
# about B, so the fitted normal has mean ~0 and sigma ~5 and fills the gap. This
# is the documented misspecification failure mode of the normal estimator.
sim_hdr_bimodal <- function(n, seed = 1) {
  withr::local_seed(seed)
  b <- sample(c(-1, 1), n, replace = TRUE)
  a <- 5 * b + stats::rnorm(n, 0, 0.5)
  tibble::tibble(exposure = a, x = stats::rnorm(n))
}

# Three-time longitudinal design with a mean-shift support gap at t = 2: the
# whole stratum's supported dose moves to mean 5, away from the common target
# a = 0, while t = 1 and t = 3 stay centered at 0. Unlike dgp_longitudinal()'s
# bimodal gap, a mean shift is visible to the normal estimator.
sim_hdr_seq_meanshift <- function(n, seed = 1) {
  withr::local_seed(seed)
  l0 <- stats::rnorm(n)
  a1 <- stats::rnorm(n, mean = 0, sd = 1)
  l1 <- stats::rnorm(n, mean = a1)
  a2 <- stats::rnorm(n, mean = 5, sd = 0.5)
  l2 <- stats::rnorm(n, mean = a2)
  a3 <- stats::rnorm(n, mean = 0, sd = 1)
  tibble::tibble(
    id = seq_len(n),
    l0 = l0,
    a1 = a1,
    l1 = l1,
    a2 = a2,
    l2 = l2,
    a3 = a3
  )
}

# The non-overlap value at a single target, matched exactly on the grid value.
tau_at <- function(res, a) {
  res@results$nonoverlap[res@results$value == a]
}

# ---- Argument validation --------------------------------------------------

test_that("check_hdr() aborts on a binary exposure", {
  local_quiet()
  data <- dgp_good_positivity(n = 100, seed = 1)
  expect_error(
    check_hdr(data, exposure, c(x1, x2)),
    class = "positively_error"
  )
})

test_that("check_hdr() aborts on a categorical exposure", {
  local_quiet()
  data <- dgp_structural_subgroup(n = 200, seed = 1)
  data$exposure <- factor(sample(c("a", "b", "c"), nrow(data), replace = TRUE))
  expect_error(
    check_hdr(data, exposure, x1),
    class = "positively_error"
  )
})

test_that("check_hdr() rejects non-data-frame input", {
  local_quiet()
  expect_error(
    check_hdr(1:10, exposure, x1),
    class = "positively_error"
  )
})

test_that("check_hdr() rejects an empty covariate selection", {
  local_quiet()
  data <- sim_hdr_linear(80, seed = 1)
  expect_error(
    check_hdr(data, exposure, tidyselect::starts_with("zzz")),
    class = "positively_error"
  )
})

test_that("check_hdr() requires a single exposure column", {
  local_quiet()
  data <- sim_hdr_linear(80, seed = 1)
  expect_error(
    check_hdr(data, c(exposure, l), l),
    class = "positively_error"
  )
})

test_that("check_hdr() rejects a mass outside the open unit interval", {
  local_quiet()
  data <- sim_hdr_linear(80, seed = 1)
  expect_error(
    check_hdr(data, exposure, l, mass = 0),
    class = "positively_error"
  )
  expect_error(
    check_hdr(data, exposure, l, mass = 1),
    class = "positively_error"
  )
  expect_error(
    check_hdr(data, exposure, l, mass = 1.5),
    class = "positively_error"
  )
  expect_error(
    check_hdr(data, exposure, l, mass = -0.1),
    class = "positively_error"
  )
})

test_that("check_hdr() rejects a non-scalar mass", {
  local_quiet()
  data <- sim_hdr_linear(80, seed = 1)
  expect_error(
    check_hdr(data, exposure, l, mass = c(0.9, 0.95)),
    class = "positively_error"
  )
})

test_that("check_hdr() aborts on missing exposure or covariate values", {
  local_quiet()
  data <- sim_hdr_linear(80, seed = 1)

  na_exposure <- data
  na_exposure$exposure[1] <- NA
  expect_error(
    check_hdr(na_exposure, exposure, l),
    class = "positively_error"
  )

  na_covariate <- data
  na_covariate$l[1] <- NA
  expect_error(
    check_hdr(na_covariate, exposure, l),
    class = "positively_error"
  )
})

test_that("check_hdr() aborts on fewer than two observations", {
  local_quiet()
  data <- tibble::tibble(exposure = 0.5, l = -0.2)
  expect_error(
    check_hdr(data, exposure, l),
    class = "positively_error"
  )
})

test_that("check_hdr() validates a user-supplied target grid", {
  local_quiet()
  data <- sim_hdr_linear(80, seed = 1)
  expect_error(
    check_hdr(data, exposure, l, values = "a"),
    class = "positively_error"
  )
  expect_error(
    check_hdr(data, exposure, l, values = c(0, NA)),
    class = "positively_error"
  )
  expect_error(
    check_hdr(data, exposure, l, values = numeric(0)),
    class = "positively_error"
  )
})

test_that("check_hdr() rejects a density_estimator that is not an hdr_density", {
  local_quiet()
  data <- sim_hdr_linear(80, seed = 1)
  expect_error(
    check_hdr(data, exposure, l, density_estimator = "not an estimator"),
    class = "positively_error"
  )
  expect_error(
    check_hdr(data, exposure, l, density_estimator = list()),
    class = "positively_error"
  )
})

# ---- Closed-form oracle and curve helpers ---------------------------------

# Closed-form non-overlap for the linear-Gaussian design (REPORT expectation 7):
# tau(a) = Phi((a - z*sigma)/s_mu) + 1 - Phi((a + z*sigma)/s_mu), s_mu = |beta|.
tau_oracle <- function(a, beta, sigma, mass) {
  z <- stats::qnorm((1 + mass) / 2)
  s_mu <- abs(beta)
  stats::pnorm((a - z * sigma) / s_mu) +
    1 -
    stats::pnorm((a + z * sigma) / s_mu)
}

# The non-overlap curve ordered by target value, so monotonicity and oracle
# comparisons read against a sorted grid regardless of row order.
tau_curve <- function(res) {
  res@results$nonoverlap[order(res@results$value)]
}

# A structurally faithful normal-equivalent estimator built through the public
# constructor. Test files do not share helper scope, so this mirrors the one in
# test-hdr-density.R. with_threshold = FALSE leaves the cutoff NULL to exercise
# the numeric-grid fallback through check_hdr().
make_normal_like <- function(with_threshold = FALSE) {
  fit <- function(formula, data) {
    model <- stats::lm(formula, data = data)
    list(model = model, sigma = summary(model)$sigma)
  }
  density <- function(state, a, newdata) {
    mu <- stats::predict(state$model, newdata = newdata)
    stats::dnorm(a, mean = mu, sd = state$sigma)
  }
  if (!with_threshold) {
    return(new_hdr_density(fit = fit, density = density))
  }
  hdr_threshold <- function(state, newdata, mass) {
    z <- stats::qnorm((1 + mass) / 2)
    rep(stats::dnorm(z) / state$sigma, nrow(newdata))
  }
  new_hdr_density(fit = fit, density = density, hdr_threshold = hdr_threshold)
}

# ---- Result class and structure -------------------------------------------

test_that("check_hdr() returns an hdr_result diagnostic", {
  local_quiet()
  data <- sim_hdr_linear(300, beta = 1, seed = 1)
  res <- check_hdr(data, exposure, l)

  expect_true(S7::S7_inherits(res, positivity_diagnostic))
  expect_identical(S7::S7_class(res)@name, "hdr_result")
  expect_identical(res@exposure_type, "continuous")
  expect_identical(res@n, 300L)
})

test_that("hdr_result carries the mass and estimator-label properties", {
  local_quiet()
  data <- sim_hdr_linear(300, beta = 1, seed = 1)
  res <- check_hdr(data, exposure, l, mass = 0.9)

  expect_type(res@mass, "double")
  expect_identical(res@mass, 0.9)
  expect_type(res@density_estimator, "character")
  expect_length(res@density_estimator, 1)
})

test_that("hdr_result has the fixed results columns", {
  local_quiet()
  data <- sim_hdr_linear(300, beta = 1, seed = 1)
  res <- check_hdr(data, exposure, l, values = c(-1, 0, 1))

  expect_s3_class(res@results, "tbl_df")
  expect_setequal(names(res@results), c("value", "nonoverlap"))
  expect_type(res@results$value, "double")
  expect_type(res@results$nonoverlap, "double")
})

test_that("results hold one aggregate row per target value", {
  local_quiet()
  data <- sim_hdr_linear(300, beta = 1, seed = 1)

  res <- check_hdr(data, exposure, l, values = c(-1, 0, 1))
  expect_identical(nrow(res@results), 3L)
  expect_setequal(res@results$value, c(-1, 0, 1))
})

test_that("the default target grid spans 100 points over the observed range", {
  local_quiet()
  data <- sim_hdr_linear(300, beta = 1, seed = 1)
  res <- check_hdr(data, exposure, l)

  expect_identical(nrow(res@results), 100L)
  expect_gte(min(res@results$value), min(data$exposure))
  expect_lte(max(res@results$value), max(data$exposure))
})

test_that("tidy() and glance() follow the shared diagnostic contract", {
  local_quiet()
  data <- sim_hdr_linear(200, beta = 1, seed = 1)
  res <- check_hdr(data, exposure, l, values = c(-1, 0, 1))

  expect_identical(generics::tidy(res), res@results)
  glanced <- generics::glance(res)
  expect_s3_class(glanced, "tbl_df")
  expect_identical(nrow(glanced), 1L)
})

# ---- Range and validity (expectation 1) -----------------------------------

test_that("every non-overlap value lies in the unit interval", {
  local_quiet()
  # Mixed designs: well-supported linear-Gaussian, a planted mixture, and a
  # misspecified bimodal truth all keep tau within [0, 1].
  linear <- check_hdr(sim_hdr_linear(500, beta = 1, seed = 1), exposure, l)
  mixture <- check_hdr(sim_hdr_mixture(500, p = 0.3, seed = 1), exposure, g)
  bimodal <- check_hdr(sim_hdr_bimodal(500, seed = 1), exposure, x)

  for (res in list(linear, mixture, bimodal)) {
    expect_true(all(res@results$nonoverlap >= 0 & res@results$nonoverlap <= 1))
  }
})

# ---- Null / quiet (expectation 2) -----------------------------------------

test_that("the central target is quiet under weak covariate dependence", {
  local_quiet()
  data <- sim_hdr_linear(5000, beta = 0.5, sigma = 1, seed = 1)
  res <- check_hdr(data, exposure, l, mass = 0.95, values = 0)

  expect_lt(tau_at(res, 0), 0.02)
})

# ---- Severe violation (expectation 3) -------------------------------------

test_that("a target beyond all support gives tau exactly one", {
  local_quiet()
  data <- sim_hdr_linear(5000, beta = 1, sigma = 1, seed = 1)
  model <- stats::lm(exposure ~ l, data = data)
  mu <- stats::predict(model)
  sigma <- summary(model)$sigma
  a_extreme <- max(mu) + 5 * sigma

  res <- check_hdr(data, exposure, l, values = a_extreme)
  expect_equal(tau_at(res, a_extreme), 1)
})

# ---- Known-fraction magnitude (expectation 4) -----------------------------

test_that("tau recovers a planted unsupported fraction", {
  local_quiet()
  data <- sim_hdr_mixture(20000, p = 0.3, seed = 1)
  res <- check_hdr(data, exposure, g, mass = 0.95, values = 0)

  expect_equal(tau_at(res, 0), 0.30, tolerance = 0.02)
})

# ---- Monotonicity in the target (expectation 5) ---------------------------

test_that("tau is non-decreasing on a one-sided grid above the center", {
  local_quiet()
  data <- sim_hdr_linear(20000, beta = 1, sigma = 1, seed = 1)
  # The center sits at mean(mu) ~ 0; above it the curve is non-decreasing, with
  # small slack for Monte-Carlo noise in the estimated endpoints.
  res <- check_hdr(data, exposure, l, values = seq(0, 3, by = 0.25))

  expect_true(all(diff(tau_curve(res)) >= -0.01))
})

# ---- Monotonicity in mass (expectation 6) ---------------------------------

test_that("tau is non-increasing as mass widens the region", {
  local_quiet()
  data <- sim_hdr_linear(10000, beta = 1.5, sigma = 1, seed = 1)
  masses <- c(0.5, 0.8, 0.9, 0.95, 0.99)
  taus <- vapply(
    masses,
    function(m) {
      tau_at(check_hdr(data, exposure, l, mass = m, values = 3), 3)
    },
    numeric(1)
  )

  expect_true(all(diff(taus) <= 0.01))
})

# ---- Closed-form oracle agreement (expectation 7) -------------------------

test_that("the whole tau curve matches the closed-form oracle", {
  local_quiet()
  beta <- 1
  sigma <- 1
  mass <- 0.95
  data <- sim_hdr_linear(20000, beta = beta, sigma = sigma, seed = 1)
  grid <- seq(-3, 3, by = 0.5)
  res <- check_hdr(data, exposure, l, mass = mass, values = grid)

  oracle <- tau_oracle(grid, beta, sigma, mass)
  expect_equal(tau_curve(res), oracle, tolerance = 0.03)
})

# ---- Numeric-grid threshold fallback (expectation 8) ----------------------

test_that("the numeric-grid cutoff reproduces the closed-form normal curve", {
  local_quiet()
  # An estimator with hdr_threshold = NULL forces the numeric-grid fallback; its
  # tau curve must match the default estimator's closed-form cutoff. Comparison
  # is on the same target grid over a well-supported range.
  data <- sim_hdr_linear(2000, beta = 1, sigma = 1, seed = 1)
  grid <- seq(-2, 2, by = 0.5)

  closed_form <- check_hdr(data, exposure, l, values = grid)
  numeric_grid <- check_hdr(
    data,
    exposure,
    l,
    values = grid,
    density_estimator = make_normal_like(with_threshold = FALSE)
  )

  expect_equal(
    tau_curve(numeric_grid),
    tau_curve(closed_form),
    tolerance = 0.01
  )
})

# ---- Pluggable-density flow ------------------------------------------------

test_that("a custom estimator with a closed-form threshold round-trips", {
  local_quiet()
  # A user-supplied estimator that mirrors the default geometry must reproduce
  # the default tau curve, confirming the pluggable path carries through
  # check_hdr() end to end.
  data <- sim_hdr_linear(2000, beta = 1, sigma = 1, seed = 1)
  grid <- c(-1.5, -0.5, 0.5, 1.5)

  custom <- check_hdr(
    data,
    exposure,
    l,
    values = grid,
    density_estimator = make_normal_like(with_threshold = TRUE)
  )
  default <- check_hdr(data, exposure, l, values = grid)

  expect_equal(tau_curve(custom), tau_curve(default), tolerance = 1e-6)
})

test_that("the estimator label reflects a user-supplied estimator", {
  local_quiet()
  data <- sim_hdr_linear(200, beta = 1, seed = 1)
  res <- check_hdr(
    data,
    exposure,
    l,
    density_estimator = make_normal_like(with_threshold = TRUE)
  )
  expect_type(res@density_estimator, "character")
  expect_identical(res@density_estimator, "custom")
})

test_that("the default estimator carries the normal label", {
  local_quiet()
  data <- sim_hdr_linear(200, beta = 1, seed = 1)
  res <- check_hdr(data, exposure, l)
  expect_identical(res@density_estimator, "normal")
})

# ---- Bimodal misspecification (expectation 11) ----------------------------

test_that("the normal estimator misses a bimodal support gap", {
  local_quiet()
  # Known failure mode: the fitted normal fills the hole between the modes, so it
  # reports the gap target (a = 0) as supported when the truth is fully
  # unsupported. A flexible density estimator would be needed to recover it.
  data <- sim_hdr_bimodal(20000, seed = 1)
  res <- check_hdr(data, exposure, x, values = 0)

  expect_lt(tau_at(res, 0), 0.05)
})

# ---- Sequential structure --------------------------------------------------

test_that("check_hdr_seq() returns an hdr_result with a time column", {
  local_quiet()
  data <- dgp_longitudinal(n = 300, seed = 5)
  res <- check_hdr_seq(
    data,
    c(a1, a2, a3),
    list(l0, l1, l2),
    values = 0
  )

  expect_true(S7::S7_inherits(res, positivity_diagnostic))
  expect_identical(S7::S7_class(res)@name, "hdr_result")
  expect_identical(res@exposure_type, "continuous")
  expect_setequal(names(res@results), c("time", "value", "nonoverlap"))
  expect_length(sort(unique(res@results$time)), 3)
})

test_that("check_hdr_seq() carries one aggregate row per time and value", {
  local_quiet()
  data <- dgp_longitudinal(n = 300, seed = 5)
  res <- check_hdr_seq(
    data,
    c(a1, a2, a3),
    list(l0, l1, l2),
    values = c(-1, 0, 1)
  )

  expect_identical(nrow(res@results), 3L * 3L)
  expect_true(all(res@results$nonoverlap >= 0 & res@results$nonoverlap <= 1))
})

# ---- Sequential per-time isolation (expectation 10) -----------------------

test_that("check_hdr_seq() isolates a mean-shift gap to its time point", {
  local_quiet()
  data <- sim_hdr_seq_meanshift(10000, seed = 1)
  res <- check_hdr_seq(
    data,
    c(a1, a2, a3),
    list(l0, l1, l2),
    values = 0
  )

  times <- sort(unique(res@results$time))
  tau <- vapply(
    times,
    function(t) {
      res@results$nonoverlap[res@results$time == t & res@results$value == 0]
    },
    numeric(1)
  )

  expect_gt(tau[2], 0.9)
  expect_lt(tau[1], 0.05)
  expect_lt(tau[3], 0.05)
})

# ---- Conditioning-set construction ----------------------------------------

test_that("a finite lag shrinks the conditioning set", {
  local_quiet()
  withr::local_seed(1)
  n <- 2000
  data <- tibble::tibble(
    id = seq_len(n),
    l0 = stats::rnorm(n),
    a1 = stats::rnorm(n),
    l1 = stats::rnorm(n),
    a2 = stats::rnorm(n),
    l2 = stats::rnorm(n)
  )
  # The final exposure depends only on the baseline covariate l0, which enters
  # the t = 3 conditioning set at lag = Inf but not at lag = 0 (only l2 does).
  # A visible l0 gives spread-out fitted means, so tau(0) is high; hidden, the
  # fitted means collapse to the marginal mean and tau(0) is near zero.
  data$a3 <- 3 * data$l0 + stats::rnorm(n, 0, 0.3)

  tau3 <- function(lag) {
    res <- check_hdr_seq(
      data,
      c(a1, a2, a3),
      list(l0, l1, l2),
      values = 0,
      lag = lag
    )
    res@results$nonoverlap[res@results$time == 3 & res@results$value == 0]
  }

  expect_gt(tau3(Inf), 0.5)
  expect_lt(tau3(0), 0.1)
})

test_that(".baseline enters every conditioning set", {
  local_quiet()
  withr::local_seed(1)
  n <- 2000
  b <- stats::rnorm(n)
  data <- tibble::tibble(
    b = b,
    l0 = stats::rnorm(n),
    l1 = stats::rnorm(n),
    l2 = stats::rnorm(n),
    a1 = 3 * b + stats::rnorm(n, 0, 0.3),
    a2 = 3 * b + stats::rnorm(n, 0, 0.3),
    a3 = 3 * b + stats::rnorm(n, 0, 0.3)
  )

  tau_at_0 <- function(res) {
    res@results$nonoverlap[res@results$value == 0]
  }

  # lag = 0 keeps prior exposures out of the conditioning sets, so the baseline
  # covariate b is the only route to the exposure-determining signal.
  with_baseline <- check_hdr_seq(
    data,
    c(a1, a2, a3),
    list(l0, l1, l2),
    .baseline = b,
    values = 0,
    lag = 0
  )
  without_baseline <- check_hdr_seq(
    data,
    c(a1, a2, a3),
    list(l0, l1, l2),
    values = 0,
    lag = 0
  )

  expect_true(all(tau_at_0(with_baseline) > 0.5))
  expect_true(all(tau_at_0(without_baseline) < 0.1))
})

test_that("a length-one covariate list recycles across time points", {
  local_quiet()
  data <- dgp_longitudinal(n = 300, seed = 5)

  recycled <- check_hdr_seq(
    data,
    c(a1, a2, a3),
    list(l1),
    values = c(-1, 0, 1)
  )
  explicit <- check_hdr_seq(
    data,
    c(a1, a2, a3),
    list(l1, l1, l1),
    values = c(-1, 0, 1)
  )

  expect_length(sort(unique(recycled@results$time)), 3)
  expect_equal(recycled@results$nonoverlap, explicit@results$nonoverlap)
})

# ---- Autoplot contract ----------------------------------------------------

test_that("autoplot() returns a ggplot for the point diagnostic", {
  local_quiet()
  data <- sim_hdr_linear(300, beta = 1, seed = 1)
  res <- check_hdr(data, exposure, l)

  expect_s3_class(ggplot2::autoplot(res), "ggplot")
})

test_that("autoplot() returns a ggplot for the sequential diagnostic", {
  local_quiet()
  data <- dgp_longitudinal(n = 300, seed = 5)
  res <- check_hdr_seq(data, c(a1, a2, a3), list(l0, l1, l2), values = 0)

  expect_s3_class(ggplot2::autoplot(res), "ggplot")
})

test_that("plot() draws the HDR view and returns the result invisibly", {
  local_quiet()
  data <- sim_hdr_linear(300, beta = 1, seed = 1)
  res <- check_hdr(data, exposure, l)

  path <- withr::local_tempfile(fileext = ".pdf")
  grDevices::pdf(path)
  withr::defer(grDevices::dev.off())
  expect_identical(plot(res), res)
})

# ---- Snapshots ------------------------------------------------------------

test_that("the point print method is stable", {
  local_quiet()
  data <- sim_hdr_linear(300, beta = 1, seed = 1)
  res <- check_hdr(data, exposure, l, values = c(-1, 0, 1))
  expect_snapshot(print(res))
})

test_that("the sequential print method is stable", {
  local_quiet()
  data <- dgp_longitudinal(n = 300, seed = 5)
  res <- check_hdr_seq(
    data,
    c(a1, a2, a3),
    list(l0, l1, l2),
    values = c(-1, 0, 1)
  )
  expect_snapshot(print(res))
})

test_that("check_hdr() argument validation messages are stable", {
  local_quiet()
  data <- sim_hdr_linear(80, seed = 1)
  binary <- dgp_good_positivity(n = 100, seed = 1)

  expect_snapshot(check_hdr(1:10, exposure, x1), error = TRUE)
  expect_snapshot(check_hdr(binary, exposure, c(x1, x2)), error = TRUE)
  expect_snapshot(
    check_hdr(data, exposure, tidyselect::starts_with("zzz")),
    error = TRUE
  )
  expect_snapshot(check_hdr(data, c(exposure, l), l), error = TRUE)
  expect_snapshot(check_hdr(data, exposure, l, mass = 0), error = TRUE)
  expect_snapshot(check_hdr(data, exposure, l, mass = 1.5), error = TRUE)
  expect_snapshot(
    check_hdr(data, exposure, l, mass = c(0.9, 0.95)),
    error = TRUE
  )
  expect_snapshot(check_hdr(data, exposure, l, values = "a"), error = TRUE)
  expect_snapshot(
    check_hdr(data, exposure, l, values = c(0, NA)),
    error = TRUE
  )
  expect_snapshot(
    check_hdr(data, exposure, l, values = numeric(0)),
    error = TRUE
  )
  expect_snapshot(
    check_hdr(data, exposure, l, density_estimator = "not an estimator"),
    error = TRUE
  )

  na_covariate <- data
  na_covariate$l[1] <- NA
  expect_snapshot(check_hdr(na_covariate, exposure, l), error = TRUE)

  one_row <- tibble::tibble(exposure = 0.5, l = -0.2)
  expect_snapshot(check_hdr(one_row, exposure, l), error = TRUE)
})

test_that("check_hdr_seq() argument validation messages are stable", {
  local_quiet()
  data <- dgp_longitudinal(n = 300, seed = 5)

  expect_snapshot(
    check_hdr_seq(data, c(a1, a2, a3), list(l0, l1)),
    error = TRUE
  )
  expect_snapshot(
    check_hdr_seq(data, c(a1, a2, a3), list(l0, l1, l2), lag = -1),
    error = TRUE
  )
  expect_snapshot(
    check_hdr_seq(data, c(a1, a2, a3), list(l0, l1, l2), lag = 1.5),
    error = TRUE
  )
  expect_snapshot(
    check_hdr_seq(data, c(a1, a2, a3), list(l0, l1, l2), lag = "a"),
    error = TRUE
  )

  binary <- data
  binary$a2 <- rep(c(0, 1), length.out = nrow(binary))
  expect_snapshot(
    check_hdr_seq(binary, c(a1, a2, a3), list(l0, l1, l2)),
    error = TRUE
  )
})
