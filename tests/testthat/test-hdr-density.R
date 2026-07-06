# The pluggable conditional-density contract (api-design section 1.4). A
# density estimator is an internal S7 object of class `hdr_density` carrying
# three function properties:
#   fit(formula, data)             -> fitted state (a list)
#   density(state, a, newdata)     -> f_hat(a | l), a recycled against newdata
#   hdr_threshold(state, newdata, mass) -> f_alpha(l_j), the mass-capturing cutoff
# `hdr_density_normal()` supplies the default: lm(exposure ~ covariates) with a
# Gaussian conditional density of residual SD sigma, and the closed-form cutoff
# dnorm(z) / sigma where z = qnorm((1 + mass) / 2). `new_hdr_density()` is the
# developer constructor; hdr_threshold is optional (NULL requests the numeric
# grid fallback), so the user contract is two functions minimum.

# A structurally faithful normal-equivalent estimator built through the public
# constructor. It mirrors the default estimator's geometry so the pluggable path
# can be checked against the closed-form path. With `with_threshold = FALSE` the
# cutoff is left NULL to request the numeric grid fallback.
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

# ---- new_hdr_density() construction and validation ------------------------

test_that("new_hdr_density() returns an hdr_density object with the three slots", {
  est <- make_normal_like()

  expect_identical(S7::S7_class(est)@name, "hdr_density")
  expect_type(est@fit, "closure")
  expect_type(est@density, "closure")
})

test_that("new_hdr_density() accepts a NULL hdr_threshold", {
  est <- new_hdr_density(
    fit = function(formula, data) list(),
    density = function(state, a, newdata) rep(1, nrow(newdata))
  )
  expect_null(est@hdr_threshold)
})

test_that("new_hdr_density() stores a supplied hdr_threshold", {
  est <- make_normal_like(with_threshold = TRUE)
  expect_type(est@hdr_threshold, "closure")
})

test_that("new_hdr_density() defaults the label to \"custom\"", {
  est <- make_normal_like()
  expect_identical(est@label, "custom")
})

test_that("new_hdr_density() stores a supplied label", {
  est <- new_hdr_density(
    fit = function(formula, data) list(),
    density = function(state, a, newdata) rep(1, nrow(newdata)),
    label = "spline"
  )
  expect_identical(est@label, "spline")
})

test_that("hdr_density_normal() carries the \"normal\" label", {
  expect_identical(hdr_density_normal()@label, "normal")
})

test_that("new_hdr_density() rejects a non-string label", {
  density <- function(state, a, newdata) rep(1, nrow(newdata))
  fit <- function(formula, data) list()

  expect_error(
    new_hdr_density(fit = fit, density = density, label = 1),
    class = "positively_error"
  )
  expect_error(
    new_hdr_density(fit = fit, density = density, label = c("a", "b")),
    class = "positively_error"
  )
  expect_error(
    new_hdr_density(fit = fit, density = density, label = NA_character_),
    class = "positively_error"
  )
})

test_that("new_hdr_density() rejects a non-function fit", {
  expect_error(
    new_hdr_density(
      fit = "not a function",
      density = function(state, a, newdata) rep(1, nrow(newdata))
    ),
    class = "positively_error"
  )
})

test_that("new_hdr_density() rejects a non-function density", {
  expect_error(
    new_hdr_density(
      fit = function(formula, data) list(),
      density = 1
    ),
    class = "positively_error"
  )
})

test_that("new_hdr_density() rejects a non-function, non-NULL hdr_threshold", {
  expect_error(
    new_hdr_density(
      fit = function(formula, data) list(),
      density = function(state, a, newdata) rep(1, nrow(newdata)),
      hdr_threshold = "not a function"
    ),
    class = "positively_error"
  )
})

# ---- hdr_density_normal() structure ---------------------------------------

test_that("hdr_density_normal() returns an hdr_density estimator", {
  est <- hdr_density_normal()

  expect_identical(S7::S7_class(est)@name, "hdr_density")
  expect_type(est@fit, "closure")
  expect_type(est@density, "closure")
  expect_type(est@hdr_threshold, "closure")
})

# ---- Threshold form and membership (expectation 9) ------------------------

test_that("the normal threshold equals dnorm(z) / sigma and is constant", {
  withr::local_seed(1)
  data <- tibble::tibble(
    exposure = stats::rnorm(400, mean = 0.5 * stats::rnorm(400)),
    l = stats::rnorm(400)
  )
  # Refit the same linear-Gaussian model the default estimator uses so the
  # residual SD is available independently.
  data$l <- stats::rnorm(400)
  data$exposure <- stats::rnorm(400, mean = 1.2 * data$l, sd = 0.8)
  model <- stats::lm(exposure ~ l, data = data)
  sigma <- summary(model)$sigma

  est <- hdr_density_normal()
  state <- est@fit(exposure ~ l, data)

  for (mass in c(0.80, 0.90, 0.95, 0.99)) {
    z <- stats::qnorm((1 + mass) / 2)
    threshold <- est@hdr_threshold(state, data, mass)
    # The Gaussian cutoff is homoskedastic, so every observation shares it.
    expect_equal(length(unique(round(threshold, 12))), 1L)
    expect_equal(
      threshold[1],
      stats::dnorm(z) / sigma,
      tolerance = 1e-6
    )
  }
})

test_that("normal HDR membership reduces to abs(a - mu) <= z * sigma", {
  withr::local_seed(2)
  data <- tibble::tibble(l = stats::rnorm(300))
  data$exposure <- stats::rnorm(300, mean = data$l, sd = 1)
  model <- stats::lm(exposure ~ l, data = data)
  mu <- stats::predict(model, newdata = data)
  sigma <- summary(model)$sigma
  mass <- 0.95
  z <- stats::qnorm((1 + mass) / 2)

  est <- hdr_density_normal()
  state <- est@fit(exposure ~ l, data)
  threshold <- est@hdr_threshold(state, data, mass)

  # For several targets, the density-vs-threshold membership test must agree
  # with the interval test on the fitted mean. Targets are held off the exact
  # boundary so a residual-SD denominator choice cannot flip a borderline case.
  for (a in c(-2, -0.5, 0.5, 2)) {
    inside_density <- est@density(state, a, data) >= threshold
    inside_interval <- abs(a - mu) <= z * sigma
    expect_identical(inside_density, inside_interval)
  }
})

test_that("the normal density matches a plain Gaussian around fitted means", {
  withr::local_seed(3)
  data <- tibble::tibble(l = stats::rnorm(200))
  data$exposure <- stats::rnorm(200, mean = 0.7 * data$l, sd = 1.3)
  model <- stats::lm(exposure ~ l, data = data)
  mu <- stats::predict(model, newdata = data)
  sigma <- summary(model)$sigma

  est <- hdr_density_normal()
  state <- est@fit(exposure ~ l, data)

  expect_equal(
    est@density(state, 0.25, data),
    stats::dnorm(0.25, mean = mu, sd = sigma),
    tolerance = 1e-6
  )
})

# ---- Snapshots ------------------------------------------------------------

test_that("new_hdr_density() validation messages are stable", {
  good_density <- function(state, a, newdata) rep(1, nrow(newdata))
  good_fit <- function(formula, data) list()

  expect_snapshot(
    new_hdr_density(fit = "not a function", density = good_density),
    error = TRUE
  )
  expect_snapshot(
    new_hdr_density(fit = good_fit, density = 1),
    error = TRUE
  )
  expect_snapshot(
    new_hdr_density(
      fit = good_fit,
      density = good_density,
      hdr_threshold = "not a function"
    ),
    error = TRUE
  )
  expect_snapshot(
    new_hdr_density(fit = good_fit, density = good_density, label = 1),
    error = TRUE
  )
})
