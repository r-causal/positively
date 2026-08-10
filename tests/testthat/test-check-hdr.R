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
# s_mu), which holds with s_mu = |beta| for the linear-Gaussian design below.
#
# The normal estimator bounds what those claims can be: it detects mean-shift
# support gaps, not multimodal gaps, and its non-overlap floor is not
# universally zero under strong covariate dependence.

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

# Three waves of a continuous dose dispensed on coarse_dose_grid(), the same
# grid dgp_coarse_dose() uses. Each wave's dose is snapped from a prescribing
# score that tracks the current covariate, so no exposure column can hold more
# than eight distinct values and detection reads every one of them as
# categorical, while the central dose stays supported and the grid endpoints do
# not.
sim_hdr_seq_coarse <- function(n = 150, seed = 1) {
  withr::local_seed(seed)
  dose_grid <- coarse_dose_grid()
  step <- dose_grid[2] - dose_grid[1]
  center <- mean(range(dose_grid))
  dispense <- function(score) {
    dose_grid[pmin(pmax(round(score / step), 1L), length(dose_grid))]
  }
  l0 <- stats::rnorm(n)
  a1 <- dispense(center + 2 * l0 + stats::rnorm(n, sd = 5))
  l1 <- stats::rnorm(n, mean = 0.3 * (a1 - center))
  a2 <- dispense(center + 2 * l1 + stats::rnorm(n, sd = 5))
  l2 <- stats::rnorm(n, mean = 0.3 * (a2 - center))
  a3 <- dispense(center + 2 * l2 + stats::rnorm(n, sd = 5))
  tibble::tibble(l0 = l0, a1 = a1, l1 = l1, a2 = a2, l2 = l2, a3 = a3)
}

# A perfectly determined exposure: A = 2 * X with no residual noise at all. The
# exposure varies, so the diagnostic has a grid to span, but the fitted density
# has no width and no profile's HDR reaches any other profile's dose, so every
# target inside the observed range reads one.
sim_hdr_determined <- function(n = 200, seed = 1) {
  withr::local_seed(seed)
  x <- stats::rnorm(n)
  tibble::tibble(exposure = 2 * x, x = x)
}

# Three waves whose second and third exposures are determined by the covariate
# entering their own conditioning sets, while the first carries genuine residual
# noise. Two waves of three read flat at one, so a report that named every wave
# or only the first would be visible.
sim_hdr_seq_determined <- function(n = 200, seed = 1) {
  withr::local_seed(seed)
  l0 <- stats::rnorm(n)
  a1 <- stats::rnorm(n, mean = l0)
  l1 <- stats::rnorm(n)
  l2 <- stats::rnorm(n)
  tibble::tibble(l0 = l0, a1 = a1, l1 = l1, a2 = 2 * l1, l2 = l2, a3 = 3 * l2)
}

# Two well-separated covariate groups with a genuine support gap between them.
# The exposure is far from deterministic, with an R-squared near 0.9 against a
# residual standard deviation of about 1, and the two groups' regions are
# roughly [-4.9, -1.1] and [1.2, 5.0], so a target in the middle is unsupported
# at every profile because the profiles genuinely disagree about the dose.
sim_hdr_separated <- function(n = 400, seed = 1) {
  withr::local_seed(seed)
  g <- rep(c(0, 1), length.out = n)
  tibble::tibble(g = g, exposure = 3 * (2 * g - 1) + stats::rnorm(n, 0, 1))
}

# Three waves of separated groups, so every wave carries a real support gap and
# none is deterministic. The waves' ranges differ, so the pooled target grid
# leaves a wave with few targets of its own to be judged on.
sim_hdr_seq_separated <- function(n = 400, seed = 1) {
  withr::local_seed(seed)
  g <- rep(c(0, 1), length.out = n)
  wave <- function(center, spread) {
    center + spread * (2 * g - 1) + stats::rnorm(n, 0, 1)
  }
  tibble::tibble(
    l0 = g,
    a1 = wave(0, 3),
    l1 = g,
    a2 = wave(6, 3),
    l2 = g,
    a3 = wave(-6, 3)
  )
}

# A three-level factor covariate carried alongside the treatment-contrast
# indicators lm() would expand it into, plus a logical and a character encoding
# of the same grouping. The group shifts are large against the residual spread,
# so a target sitting on one group's dose is unsupported at the other two
# groups' profiles and no equivalence below is read off a flat curve.
sim_hdr_grouped <- function(n = 200, seed = 1) {
  withr::local_seed(seed)
  g <- factor(sample(c("a", "b", "c"), n, replace = TRUE))
  l <- stats::rnorm(n)
  shift <- unname(c(a = 0, b = 5, c = -5)[as.character(g)])
  tibble::tibble(
    exposure = shift + 0.5 * l + stats::rnorm(n, 0, 0.5),
    l = l,
    g = g,
    g_b = as.numeric(g == "b"),
    g_c = as.numeric(g == "c"),
    g_is_b = g == "b",
    g_chr = as.character(g)
  )
}

# Three waves whose doses all shift with one two-level factor, carried alongside
# the single indicator lm() would expand it into. Each wave's two groups sit
# about eight residual standard deviations apart, so either target on the pooled
# grid is unsupported at roughly half the profiles.
sim_hdr_seq_grouped <- function(n = 200, seed = 1) {
  withr::local_seed(seed)
  g <- factor(rep(c("ctl", "trt"), length.out = n))
  shift <- 4 * (g == "trt")
  tibble::tibble(
    g = g,
    g_trt = as.numeric(g == "trt"),
    l1 = stats::rnorm(n),
    l2 = stats::rnorm(n),
    l3 = stats::rnorm(n),
    a1 = shift + stats::rnorm(n, 0, 0.5),
    a2 = shift + stats::rnorm(n, 0, 0.5),
    a3 = shift + stats::rnorm(n, 0, 0.5)
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
  data <- dgp_categorical(n = 200, seed = 1)
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

test_that("a renamed exposure inside the covariate selection is rejected", {
  local_quiet()
  data <- sim_hdr_linear(80, seed = 1)
  # Renaming the exposure column inside the selection hides it from a name-based
  # overlap check, so the guard must compare resolved positions. Otherwise the
  # renamed column conditions the density on the exposure itself and the failure
  # surfaces later as a confusing missing-column error.
  expect_error(
    check_hdr(data, exposure, c(foo = exposure, l)),
    class = "positively_selection_error"
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

test_that("check_hdr() rejects a non-finite target value", {
  local_quiet()
  data <- sim_hdr_linear(80, seed = 1)
  expect_error(
    check_hdr(data, exposure, l, values = c(0, Inf)),
    class = "positively_range_error"
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

# ---- Covariate types ------------------------------------------------------

# lm() fits a factor through its treatment-contrast indicators, so the design
# matrix behind a factor selection and the design matrix behind the same
# variable's hand-encoded indicator columns are the same matrix. Every
# equivalence in this section rests on that identity.

test_that("check_hdr() accepts a factor covariate", {
  local_quiet()
  data <- sim_hdr_grouped(200, seed = 1)
  values <- c(-5, 0, 5)

  factor_fit <- check_hdr(data, exposure, c(l, g), values = values)
  indicator_fit <- check_hdr(data, exposure, c(l, g_b, g_c), values = values)

  expect_equal(generics::tidy(factor_fit), generics::tidy(indicator_fit))
  # Each target belongs to one of the three groups, so most profiles exclude it.
  # A curve pinned at zero would satisfy the equality above without testing it.
  expect_true(all(generics::tidy(factor_fit)$nonoverlap > 0.5))
})

test_that("check_hdr() accepts logical and character covariates", {
  local_quiet()
  data <- sim_hdr_grouped(200, seed = 1)
  values <- c(-5, 0, 5)

  # The character column holds the factor's labels, and lm() reads it as a
  # factor over the same alphabetically ordered levels.
  expect_equal(
    generics::tidy(check_hdr(data, exposure, c(l, g_chr), values = values)),
    generics::tidy(check_hdr(data, exposure, c(l, g), values = values))
  )
  # The logical column is the "b" indicator, so it stands in for that one dummy.
  expect_equal(
    generics::tidy(check_hdr(data, exposure, c(l, g_is_b), values = values)),
    generics::tidy(check_hdr(data, exposure, c(l, g_b), values = values))
  )
})

test_that("check_hdr() rejects a Date covariate", {
  local_quiet()
  data <- sim_hdr_linear(80, seed = 1)
  data$d <- as.Date("2020-01-01") + seq_len(nrow(data))

  expect_snapshot_abort(
    check_hdr(data, exposure, c(l, d)),
    class = "positively_type_error"
  )
})

test_that("check_hdr() rejects missing values in an accepted covariate type", {
  local_quiet()
  data <- sim_hdr_grouped(200, seed = 1)
  data$g[1] <- NA
  data$g_chr[2] <- NA

  expect_error(
    check_hdr(data, exposure, c(l, g)),
    class = "positively_missing_error"
  )
  expect_error(
    check_hdr(data, exposure, c(l, g_chr)),
    class = "positively_missing_error"
  )
})

test_that("an unaccepted covariate type is reported before a missing value", {
  local_quiet()
  data <- sim_hdr_linear(80, seed = 1)
  data$d <- as.Date("2020-01-01") + seq_len(nrow(data))
  data$l[1] <- NA

  # Both gates have grounds to fire on this selection, so the class says which
  # one runs first. A type failure is the more specific report of the two.
  expect_error(
    check_hdr(data, exposure, c(l, d)),
    class = "positively_type_error"
  )
})

test_that("check_hdr_seq() accepts a factor covariate", {
  local_quiet()
  data <- sim_hdr_seq_grouped(200, seed = 1)
  values <- c(0, 4)

  factor_fit <- check_hdr_seq(
    data,
    c(a1, a2, a3),
    list(c(l1, g), c(l2, g), c(l3, g)),
    values = values
  )
  indicator_fit <- check_hdr_seq(
    data,
    c(a1, a2, a3),
    list(c(l1, g_trt), c(l2, g_trt), c(l3, g_trt)),
    values = values
  )

  expect_equal(generics::tidy(factor_fit), generics::tidy(indicator_fit))
  expect_true(all(generics::tidy(factor_fit)$nonoverlap > 0.3))
})

test_that("check_hdr_seq() accepts a factor .baseline", {
  local_quiet()
  data <- sim_hdr_seq_grouped(200, seed = 1)
  values <- c(0, 4)

  # lag = 0 keeps prior exposures out of the conditioning sets, so the baseline
  # is the only route to the grouping that determines the dose.
  factor_fit <- check_hdr_seq(
    data,
    c(a1, a2, a3),
    list(l1, l2, l3),
    .baseline = g,
    values = values,
    lag = 0
  )
  indicator_fit <- check_hdr_seq(
    data,
    c(a1, a2, a3),
    list(l1, l2, l3),
    .baseline = g_trt,
    values = values,
    lag = 0
  )

  expect_equal(generics::tidy(factor_fit), generics::tidy(indicator_fit))
  expect_true(all(generics::tidy(factor_fit)$nonoverlap > 0.3))
})

test_that("check_hdr_seq() rejects a Date covariate or baseline", {
  local_quiet()
  data <- sim_hdr_seq_grouped(200, seed = 1)
  data$d <- as.Date("2020-01-01") + seq_len(nrow(data))

  expect_snapshot_abort(
    check_hdr_seq(data, c(a1, a2, a3), list(c(l1, d), l2, l3)),
    class = "positively_type_error"
  )
  expect_snapshot_abort(
    check_hdr_seq(data, c(a1, a2, a3), list(l1, l2, l3), .baseline = d),
    class = "positively_type_error"
  )

  # The two routes share a message body, so the argument the failure is
  # attributed to is the only thing that separates them. Snapshot expectations
  # skip on CRAN, so this match is the attribution guard there.
  baseline_error <- expect_error(
    check_hdr_seq(data, c(a1, a2, a3), list(l1, l2, l3), .baseline = d),
    class = "positively_type_error"
  )
  expect_match(
    conditionMessage(baseline_error),
    "`.baseline` must select",
    fixed = TRUE
  )
})

test_that("check_hdr_seq() rejects missing values in a covariate or baseline", {
  local_quiet()
  data <- sim_hdr_seq_grouped(200, seed = 1)
  data$g[1] <- NA

  # The covariate sets and the baseline reach the completeness gate by separate
  # routes, so each needs its own call to be exercised.
  expect_error(
    check_hdr_seq(data, c(a1, a2, a3), list(c(l1, g), l2, l3)),
    class = "positively_missing_error"
  )
  expect_error(
    check_hdr_seq(data, c(a1, a2, a3), list(l1, l2, l3), .baseline = g),
    class = "positively_missing_error"
  )
})

# ---- Declared exposure type ------------------------------------------------

test_that("a declared continuous type runs the diagnostic on a coarse dose", {
  local_quiet()
  data <- dgp_coarse_dose(n = 150, seed = 1)
  res <- check_hdr(data, exposure, x1, exposure_type = "continuous")

  expect_true(S7::S7_inherits(res, positivity_diagnostic))
  expect_identical(S7::S7_class(res)@name, "hdr_result")
  expect_identical(res@exposure_type, "continuous")
  expect_identical(res@n, 150L)
  expect_identical(nrow(res@results), 100L)
  expect_true(all(res@results$nonoverlap >= 0 & res@results$nonoverlap <= 1))

  # The dose tracks the covariate, so the middle of the dispensing grid is well
  # supported while both ends are not. Pinning that contrast, rather than the
  # absence of an error, is what keeps the declared type from merely reaching
  # degenerate arithmetic.
  targeted <- check_hdr(
    data,
    exposure,
    x1,
    exposure_type = "continuous",
    values = c(2.5, 11.25, 20)
  )
  expect_lt(tau_at(targeted, 11.25), 0.2)
  expect_gt(tau_at(targeted, 2.5), 0.4)
  expect_gt(tau_at(targeted, 20), 0.4)
})

test_that("the coarse dose a declaration rescues is one detection misreads", {
  local_quiet()
  # Paired with the test above: without this half, raising the unique-value
  # cutoff in is_categorical() would make the fixture detect continuous and leave
  # the declaration test passing while it covered nothing.
  data <- dgp_coarse_dose(n = 150, seed = 1)
  expect_identical(detect_exposure_type(data$exposure), "categorical")

  err <- expect_error(
    check_hdr(data, exposure, x1, exposure_type = "auto"),
    class = "positively_exposure_type_error"
  )
  expect_match(conditionMessage(err), "categorical", fixed = TRUE)
})

test_that("a declared continuous type rejects a non-numeric exposure", {
  local_quiet()
  data <- dgp_coarse_dose(n = 150, seed = 1)
  data$exposure <- factor(rep(c("low", "mid", "high"), length.out = nrow(data)))

  err <- expect_error(
    check_hdr(data, exposure, x1, exposure_type = "continuous"),
    class = "positively_exposure_type_error"
  )
  # The declaration is what fails, so the error names the type the column cannot
  # carry rather than arriving later as a generic non-numeric complaint from
  # validate_numeric_columns().
  expect_match(conditionMessage(err), "continuous", fixed = TRUE)
})

test_that("check_hdr() rejects a type outside its supported menu", {
  local_quiet()
  data <- dgp_coarse_dose(n = 150, seed = 1)
  err <- expect_error(
    check_hdr(data, exposure, x1, exposure_type = "binary"),
    class = "positively_args_error"
  )
  expect_match(
    conditionMessage(err),
    '`exposure_type` must be one of "auto" or "continuous", not "binary".',
    fixed = TRUE
  )
})

test_that("check_hdr() announces nothing on either exposure_type path", {
  withr::local_options(positively.quiet = FALSE)
  data <- sim_hdr_linear(100, beta = 1, seed = 1)

  # check_hdr() supports one exposure type, so a returning call took the only
  # path there is and the announcement would restate what the function already
  # promises. A declared type does not announce, so both paths must be silent.
  expect_no_message(check_hdr(data, exposure, l, values = 0))
  expect_no_message(
    check_hdr(data, exposure, l, values = 0, exposure_type = "continuous")
  )
})

# ---- Degenerate exposures -------------------------------------------------

test_that("check_hdr() aborts on a constant exposure", {
  local_quiet()
  # A constant exposure collapses the documented 100-point grid to one target
  # and leaves the fitted density no width, so the ratio reads zero at the one
  # observed dose. That is not an imprecise answer, it is the opposite of the
  # truth, which is why this aborts rather than warns.
  declared <- tibble::tibble(
    exposure = rep(3, 50),
    l = seq(-1, 1, length.out = 50)
  )
  expect_error(
    check_hdr(declared, exposure, l, exposure_type = "continuous"),
    class = "positively_variance_error"
  )

  # The unique-value heuristic reads a constant column as categorical only once
  # the ratio of one distinct value to the row count falls below the cutoff, so
  # at five rows or fewer the auto path reaches the same arithmetic.
  detected <- tibble::tibble(exposure = rep(3, 4), l = c(-1, 0, 1, 2))
  expect_identical(detect_exposure_type(detected$exposure), "continuous")
  expect_error(
    check_hdr(detected, exposure, l),
    class = "positively_variance_error"
  )
})

test_that("check_hdr_seq() names every constant wave in one abort", {
  local_quiet()
  withr::local_options(cli.width = 200)
  data <- dgp_longitudinal(n = 300, seed = 5)
  data$a1 <- 2
  data$a3 <- 7

  err <- expect_error(
    check_hdr_seq(
      data,
      c(a1, a2, a3),
      list(l0, l1, l2),
      exposure_type = "continuous",
      values = 0
    ),
    class = "positively_variance_error"
  )
  # Both offending waves are named together, so a user does not correct one
  # column, rerun, and be told about the next.
  expect_match(conditionMessage(err), "`a1` and `a3`", fixed = TRUE)
})

test_that("no unclassed condition escapes a perfectly determined exposure", {
  local_quiet()
  # summary.lm() calls the fit "essentially perfect" and warns about it, which
  # would put a bare base R warning in front of a user of a package whose every
  # condition is classed. Nothing the diagnostic raises here may be unclassed.
  data <- sim_hdr_determined(200, seed = 1)

  raised <- list()
  withCallingHandlers(
    check_hdr(data, exposure, x),
    warning = function(cnd) {
      raised[[length(raised) + 1L]] <<- cnd
      invokeRestart("muffleWarning")
    }
  )

  expect_gt(length(raised), 0)
  for (cnd in raised) {
    expect_s3_class(cnd, "positively_warning")
  }
})

test_that("a non-overlap curve flat at one inside the observed range warns", {
  local_quiet()
  # The exposure varies, so a ratio of one is literally correct: a dose fixed by
  # the covariates is supported at no other profile. What the warning reports is
  # that the reading is driven by the width of the fitted density rather than by
  # any comparison between profiles.
  data <- sim_hdr_determined(200, seed = 1)

  expect_warning(
    check_hdr(data, exposure, x, values = c(-1, 0, 1)),
    class = "positively_degenerate_density_warning"
  )
  expect_warning(
    check_hdr(data, exposure, x),
    class = "positively_degenerate_density_warning"
  )
})

test_that("a saturated fit returns a result rather than a bare R condition", {
  local_quiet()
  # With as many parameters as observations the residual standard deviation is
  # NaN, so every ratio is NA. A flat-at-one test written as a comparison over
  # those ratios is itself NA, and branching on it raises an unclassed base R
  # error from inside a function whose whole purpose is to keep those away.
  saturated <- tibble::tibble(
    exposure = c(1, 2, 4),
    x = c(0, 1, 3),
    y = c(1, 0, 2)
  )
  res <- expect_no_warning(
    check_hdr(saturated, exposure, c(x, y), exposure_type = "continuous")
  )
  expect_identical(S7::S7_class(res)@name, "hdr_result")
  expect_true(all(is.na(res@results$nonoverlap)))

  # The second wave conditions on three columns over three rows, so it saturates
  # while the first wave keeps a residual standard deviation and reads normally.
  seq_saturated <- tibble::tibble(
    l0 = c(0, 1, 2),
    a1 = c(0, 3, 1),
    l1 = c(1, 0, 2),
    a2 = c(2, 5, 1)
  )
  seq_res <- expect_no_warning(
    check_hdr_seq(
      seq_saturated,
      c(a1, a2),
      list(l0, l1),
      exposure_type = "continuous"
    )
  )
  expect_identical(S7::S7_class(seq_res)@name, "hdr_result")
  expect_true(anyNA(seq_res@results$nonoverlap))
})

test_that("a genuine support gap between profiles does not warn", {
  local_quiet()
  # The exposure is nowhere near deterministic and the ratio of one in the gap
  # is set precisely by a comparison between profiles, which is the reading the
  # diagnostic exists to report. Whether the warning fires must be a property of
  # the fit, so it cannot depend on which doses the user happened to ask about.
  data <- sim_hdr_separated(400, seed = 1)

  expect_no_warning(check_hdr(data, exposure, g, values = c(-0.5, 0, 0.5)))
  expect_no_warning(check_hdr(data, exposure, g, values = 0))
  expect_no_warning(check_hdr(data, exposure, g))

  # The probe list really does land in the gap, so the test would pass on a
  # criterion that had simply stopped firing.
  probed <- suppressWarnings(
    check_hdr(data, exposure, g, values = c(-0.5, 0, 0.5))
  )
  expect_true(all(probed@results$nonoverlap == 1))
})

test_that("a sequential support gap does not warn on the pooled grid", {
  local_quiet()
  # Each wave is judged against its own observed range while the target grid
  # spans the pooled one, so a wave can be left with very few targets of its own.
  # A criterion reading the shared grid would let that decide the verdict.
  data <- sim_hdr_seq_separated(400, seed = 1)

  expect_no_warning(check_hdr_seq(data, c(a1, a2, a3), list(l0, l1, l2)))
  expect_no_warning(
    check_hdr_seq(data, c(a1, a2, a3), list(l0, l1, l2), values = c(-6, 0, 6))
  )
})

test_that("check_hdr_seq() names every degenerate wave in one warning", {
  local_quiet()
  withr::local_options(cli.width = 200)
  data <- sim_hdr_seq_determined(200, seed = 1)

  cnd <- expect_warning(
    check_hdr_seq(
      data,
      c(a1, a2, a3),
      list(l0, l1, l2),
      values = c(-1, 0, 1)
    ),
    class = "positively_degenerate_density_warning"
  )
  # The first wave carries residual noise, so naming every wave would be as
  # wrong as naming only the first.
  expect_match(conditionMessage(cnd), "2 and 3", fixed = TRUE)
})

# ---- Closed-form oracle and curve helpers ---------------------------------

# Closed-form non-overlap for the linear-Gaussian design:
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

test_that("tidy() returns the results tibble", {
  local_quiet()
  data <- sim_hdr_linear(200, beta = 1, seed = 1)
  res <- check_hdr(data, exposure, l, values = c(-1, 0, 1))

  expect_identical(generics::tidy(res), res@results)
})

test_that("glance() reports the HDR settings and the non-overlap range", {
  local_quiet()
  data <- sim_hdr_linear(200, beta = 1, seed = 1)
  values <- c(-1, 0, 1)
  res <- check_hdr(data, exposure, l, values = values)
  glanced <- generics::glance(res)

  expect_s3_class(glanced, "tbl_df")
  expect_identical(nrow(glanced), 1L)
  expect_setequal(
    names(glanced),
    c(
      "n",
      "mass",
      "density_estimator",
      "n_values",
      "nonoverlap_min",
      "nonoverlap_max"
    )
  )

  expect_identical(glanced$n, 200L)
  expect_identical(glanced$mass, 0.95)
  expect_identical(glanced$density_estimator, "normal")
  expect_identical(glanced$n_values, 3L)

  # Under the normal estimator tau(a) is the fraction of fitted means further
  # than z * sigma_hat from a, so the whole curve follows from lm() alone.
  model <- stats::lm(exposure ~ l, data = data)
  radius <- stats::qnorm((1 + 0.95) / 2) * stats::sigma(model)
  fitted_means <- stats::fitted(model)
  tau <- vapply(
    values,
    function(value) mean(abs(value - fitted_means) > radius),
    numeric(1)
  )
  expect_equal(glanced$nonoverlap_min, min(tau))
  expect_equal(glanced$nonoverlap_max, max(tau))
})

test_that("glance() adds the time-point count on a sequential result", {
  local_quiet()
  data <- dgp_longitudinal(n = 200, seed = 5)
  res <- check_hdr_seq(data, c(a1, a2, a3), list(l0, l1, l2), values = c(0, 3))
  glanced <- generics::glance(res)

  expect_identical(nrow(glanced), 1L)
  expect_setequal(
    names(glanced),
    c(
      "n",
      "n_times",
      "mass",
      "density_estimator",
      "n_values",
      "nonoverlap_min",
      "nonoverlap_max"
    )
  )

  # Three exposure columns and two targets were supplied.
  expect_identical(glanced$n_times, 3L)
  expect_identical(glanced$n_values, 2L)

  # The range spans every wave. dgp_longitudinal() gives its three waves
  # different support, the second of them bimodal, so the extremes come from
  # different time points and a summary that read one wave would report the
  # wrong range.
  results <- generics::tidy(res)
  expect_equal(glanced$nonoverlap_min, min(results$nonoverlap))
  expect_equal(glanced$nonoverlap_max, max(results$nonoverlap))
  expect_false(identical(
    results$time[which.min(results$nonoverlap)],
    results$time[which.max(results$nonoverlap)]
  ))
})

# ---- Range and validity ---------------------------------------------------

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

# ---- Null / quiet ---------------------------------------------------------

test_that("the central target is quiet under weak covariate dependence", {
  local_quiet()
  data <- sim_hdr_linear(5000, beta = 0.5, sigma = 1, seed = 1)
  res <- check_hdr(data, exposure, l, mass = 0.95, values = 0)

  expect_lt(tau_at(res, 0), 0.02)
})

# ---- Severe violation -----------------------------------------------------

test_that("a target beyond all support gives tau exactly one", {
  local_quiet()
  data <- sim_hdr_linear(5000, beta = 1, sigma = 1, seed = 1)
  model <- stats::lm(exposure ~ l, data = data)
  mu <- stats::predict(model)
  sigma <- summary(model)$sigma
  a_extreme <- max(mu) + 5 * sigma

  # Every target here reads one, but the target sits beyond all observed
  # support, so it is a deliberate probe at an unobserved dose and is entitled
  # to that ratio. Restricting the flat-at-one reading to targets inside the
  # observed range is what keeps this case quiet.
  res <- expect_no_warning(check_hdr(data, exposure, l, values = a_extreme))
  expect_equal(tau_at(res, a_extreme), 1)
})

# ---- Known-fraction magnitude ---------------------------------------------

test_that("tau recovers a planted unsupported fraction", {
  local_quiet()
  data <- sim_hdr_mixture(20000, p = 0.3, seed = 1)
  res <- check_hdr(data, exposure, g, mass = 0.95, values = 0)

  expect_equal(tau_at(res, 0), 0.30, tolerance = 0.02)
})

# ---- Monotonicity in the target -------------------------------------------

test_that("tau is non-decreasing on a one-sided grid above the center", {
  local_quiet()
  data <- sim_hdr_linear(20000, beta = 1, sigma = 1, seed = 1)
  # The center sits at mean(mu) ~ 0; above it the curve is non-decreasing, with
  # small slack for Monte-Carlo noise in the estimated endpoints.
  res <- check_hdr(data, exposure, l, values = seq(0, 3, by = 0.25))

  expect_true(all(diff(tau_curve(res)) >= -0.01))
})

# ---- Monotonicity in mass -------------------------------------------------

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

# ---- Closed-form oracle agreement -----------------------------------------

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

# ---- Numeric-grid threshold fallback --------------------------------------

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

# ---- Bimodal misspecification ---------------------------------------------

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

# ---- Sequential declared exposure type -------------------------------------

test_that("a declared continuous type runs the sequential diagnostic", {
  local_quiet()
  data <- sim_hdr_seq_coarse(n = 150, seed = 1)

  res <- check_hdr_seq(
    data,
    c(a1, a2, a3),
    list(l0, l1, l2),
    exposure_type = "continuous",
    values = c(2.5, 11.25, 20)
  )

  expect_identical(S7::S7_class(res)@name, "hdr_result")
  expect_identical(res@exposure_type, "continuous")
  expect_identical(nrow(res@results), 3L * 3L)
  expect_length(sort(unique(res@results$time)), 3)
  expect_true(all(res@results$nonoverlap >= 0 & res@results$nonoverlap <= 1))

  # The central dose is supported at every wave and the grid endpoints are not,
  # so the sequential curve carries signal rather than a flat zero.
  central <- res@results$nonoverlap[res@results$value == 11.25]
  ends <- res@results$nonoverlap[res@results$value != 11.25]
  expect_true(all(central < 0.2))
  expect_true(all(ends > 0.2))
})

test_that("the sequential doses a declaration rescues are ones detection misreads", {
  local_quiet()
  # A width wide enough that no bullet wraps, so the assertions below read the
  # message rather than the console geometry they happen to run under.
  withr::local_options(cli.width = 200)
  # Paired with the test above, as in the point diagnostics: without this half,
  # raising the unique-value cutoff in is_categorical() would make the waves
  # detect continuous and leave the declaration test passing while it covered
  # nothing. Asserting the auto path by condition class rather than by snapshot
  # also keeps the sequential gate standing when the snapshot is re-recorded.
  data <- sim_hdr_seq_coarse(n = 150, seed = 1)
  for (name in c("a1", "a2", "a3")) {
    expect_identical(detect_exposure_type(data[[name]]), "categorical")
  }

  err <- expect_error(
    check_hdr_seq(
      data,
      c(a1, a2, a3),
      list(l0, l1, l2),
      exposure_type = "auto",
      values = c(2.5, 11.25, 20)
    ),
    class = "positively_exposure_type_error"
  )
  # Every wave is coarse, so every wave offends. Naming them together is what
  # keeps a user from correcting one column, rerunning, and being told about the
  # next one.
  message <- conditionMessage(err)
  expect_match(
    message,
    "`a1`, `a2`, and `a3` were detected as \"categorical\".",
    fixed = TRUE
  )
  expect_match(
    message,
    "If `a1`, `a2`, and `a3` are continuous, set `exposure_type = \"continuous\"`.",
    fixed = TRUE
  )
})

test_that("a declared continuous type names every factor exposure it cannot carry", {
  local_quiet()
  withr::local_options(cli.width = 200)
  data <- sim_hdr_seq_coarse(n = 150, seed = 1)
  data$a2 <- factor(rep(c("low", "high"), length.out = nrow(data)))
  data$a3 <- factor(rep(c("low", "high"), length.out = nrow(data)))

  err <- expect_error(
    check_hdr_seq(
      data,
      c(a1, a2, a3),
      list(l0, l1, l2),
      exposure_type = "continuous",
      values = 0
    ),
    class = "positively_exposure_type_error"
  )
  # Every offending column is named, so a user with several exposures knows
  # which ones the declared type cannot describe.
  expect_match(
    conditionMessage(err),
    "`a2` and `a3` are <factor>.",
    fixed = TRUE
  )
})

test_that("sequential offenders read differently are named apart", {
  local_quiet()
  withr::local_options(cli.width = 200)
  data <- dgp_longitudinal(n = 300, seed = 5)
  data$a2 <- rep(c(0, 1), length.out = nrow(data))
  data$a3 <- rep(seq_len(8), length.out = nrow(data))
  expect_identical(detect_exposure_type(data$a2), "binary")
  expect_identical(detect_exposure_type(data$a3), "categorical")

  err <- expect_error(
    check_hdr_seq(data, c(a1, a2, a3), list(l0, l1, l2), values = 0),
    class = "positively_exposure_type_error"
  )
  message <- conditionMessage(err)

  # One bullet per distinct reading, so a user is told what each column was read
  # as rather than being handed one verdict covering both.
  readings <- gregexpr("detected as", message, fixed = TRUE)[[1]]
  expect_length(readings, 2L)
  expect_match(message, "`a2` was detected as \"binary\".", fixed = TRUE)
  expect_match(message, "`a3` was detected as \"categorical\".", fixed = TRUE)

  # The hint is single, because the escape hatch is the same for either reading.
  expect_match(
    message,
    "If `a2` and `a3` are continuous, set `exposure_type = \"continuous\"`.",
    fixed = TRUE
  )
})

test_that("a detection failure is reported before a structural one", {
  local_quiet()
  withr::local_options(cli.width = 200)
  # A Date column is read as continuous and then cannot carry it; a two-valued
  # column is read as binary, which the diagnostic does not support at all. The
  # structural offender comes first in selection order, so position cannot be
  # what decides which of the two failures is reported.
  data <- dgp_longitudinal(n = 300, seed = 5)
  data$a2 <- as.Date("2020-01-01") + seq_len(nrow(data))
  data$a3 <- rep(c(0, 1), length.out = nrow(data))
  expect_identical(detect_exposure_type(data$a2), "continuous")
  expect_identical(detect_exposure_type(data$a3), "binary")

  err <- expect_error(
    check_hdr_seq(data, c(a1, a2, a3), list(l0, l1, l2), values = 0),
    class = "positively_exposure_type_error"
  )
  message <- conditionMessage(err)
  expect_match(message, "`a3` was detected as \"binary\".", fixed = TRUE)
  expect_no_match(message, "a2", fixed = TRUE)
})

test_that("a column name containing braces survives the detection message", {
  local_quiet()
  withr::local_options(cli.width = 200)
  # Column names are user data. They reach the message as values rather than as
  # template text, so a name cli would otherwise read as an expression stays a
  # name, and the package's own condition classes survive instead of being
  # replaced by an internal cli failure.
  data <- dgp_longitudinal(n = 300, seed = 5)
  names(data)[match(c("a2", "a3"), names(data))] <- c("a{b}", "a{1}")
  data[["a{b}"]] <- rep(c(0, 1), length.out = nrow(data))
  data[["a{1}"]] <- rep(c(0, 1), length.out = nrow(data))

  err <- expect_error(
    check_hdr_seq(data, c(a1, `a{b}`, `a{1}`), list(l0, l1, l2), values = 0),
    class = "positively_exposure_type_error"
  )
  expect_match(
    conditionMessage(err),
    "`a{b}` and `a{1}` were detected as \"binary\".",
    fixed = TRUE
  )
})

test_that("a column name containing braces survives the structural message", {
  local_quiet()
  withr::local_options(cli.width = 200)
  data <- dgp_longitudinal(n = 300, seed = 5)
  names(data)[match(c("a2", "a3"), names(data))] <- c("a{b}", "a{1}")
  data[["a{b}"]] <- factor(rep(c("low", "high"), length.out = nrow(data)))
  data[["a{1}"]] <- factor(rep(c("low", "high"), length.out = nrow(data)))

  err <- expect_error(
    check_hdr_seq(
      data,
      c(a1, `a{b}`, `a{1}`),
      list(l0, l1, l2),
      exposure_type = "continuous",
      values = 0
    ),
    class = "positively_exposure_type_error"
  )
  expect_match(
    conditionMessage(err),
    "`a{b}` and `a{1}` are <factor>.",
    fixed = TRUE
  )
})

test_that("structural offenders of different classes are named apart", {
  local_quiet()
  withr::local_options(cli.width = 200)
  # The structural bullets group by the column's class vector, so a run whose
  # offenders all share one class would leave a collapsed grouping green.
  data <- dgp_longitudinal(n = 300, seed = 5)
  data$a1 <- factor(rep(c("low", "high"), length.out = nrow(data)))
  data$a2 <- rep(c("low", "high"), length.out = nrow(data))
  data$a3 <- as.Date("2020-01-01") + seq_len(nrow(data))

  err <- expect_error(
    check_hdr_seq(
      data,
      c(a1, a2, a3),
      list(l0, l1, l2),
      exposure_type = "continuous",
      values = 0
    ),
    class = "positively_exposure_type_error"
  )
  message <- conditionMessage(err)
  expect_match(message, "`a1` is <factor>.", fixed = TRUE)
  expect_match(message, "`a2` is <character>.", fixed = TRUE)
  expect_match(message, "`a3` is <Date>.", fixed = TRUE)
})

test_that("check_hdr_seq() rejects a type outside its supported menu", {
  local_quiet()
  data <- sim_hdr_seq_coarse(n = 150, seed = 1)

  err <- expect_error(
    check_hdr_seq(
      data,
      c(a1, a2, a3),
      list(l0, l1, l2),
      exposure_type = "binary",
      values = 0
    ),
    class = "positively_args_error"
  )
  expect_match(
    conditionMessage(err),
    '`exposure_type` must be one of "auto" or "continuous", not "binary".',
    fixed = TRUE
  )
})

test_that("check_hdr_seq() resolves every exposure without announcing", {
  withr::local_options(positively.quiet = FALSE)
  data <- dgp_longitudinal(n = 300, seed = 5)

  # Detection is run once per exposure column, and its announcement names
  # `.exposure`, which is not an argument of this function. Resolving the
  # sequential types silently is what keeps one call from emitting a misnamed
  # message per time point.
  expect_no_message(
    check_hdr_seq(data, c(a1, a2, a3), list(l0, l1, l2), values = 0)
  )
})

# ---- Sequential per-time isolation ----------------------------------------

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

# ---- Exposure in the conditioning set -------------------------------------

test_that("check_hdr() rejects covariates overlapping the exposure", {
  local_quiet()
  data <- sim_hdr_linear(80, seed = 1)
  expect_error(
    check_hdr(data, exposure, c(exposure, l)),
    class = "positively_selection_error"
  )
  expect_error(
    check_hdr(data, exposure, tidyselect::everything()),
    class = "positively_selection_error"
  )
})

test_that("check_hdr_seq() rejects an exposure in its own conditioning set", {
  local_quiet()
  data <- dgp_longitudinal(n = 300, seed = 5)
  # a1 is the time-1 exposure, so listing it in the time-1 covariate set puts the
  # response into its own model.
  expect_error(
    check_hdr_seq(data, c(a1, a2), list(c(l0, a1), l1), values = 0),
    class = "positively_selection_error"
  )
  # A baseline exposure enters every conditioning set, including its own time.
  expect_error(
    check_hdr_seq(
      data,
      c(a1, a2, a3),
      list(l0, l1, l2),
      .baseline = a2,
      values = 0
    ),
    class = "positively_selection_error"
  )
})

test_that("check_hdr_seq() keeps a prior exposure as a valid conditioning column", {
  local_quiet()
  data <- dgp_longitudinal(n = 300, seed = 5)
  # At the full lag the time-1 exposure a1 legitimately conditions the time-2
  # model, so the run succeeds with every time point represented.
  res <- check_hdr_seq(data, c(a1, a2, a3), list(l0, l1, l2), values = 0)
  expect_identical(S7::S7_class(res)@name, "hdr_result")
  expect_length(sort(unique(res@results$time)), 3)
})

# ---- Empty conditioning sets ----------------------------------------------

test_that("check_hdr_seq() rejects an empty conditioning set", {
  local_quiet()
  data <- dgp_longitudinal(n = 300, seed = 5)
  expect_error(
    check_hdr_seq(
      data,
      c(a1, a2, a3),
      list(tidyselect::starts_with("zzz")),
      values = 0
    ),
    class = "positively_empty_error"
  )
})

test_that("an empty per-time selection is allowed when .baseline covers it", {
  local_quiet()
  data <- dgp_longitudinal(n = 300, seed = 5)
  # The per-time covariate selection is empty, but the baseline covariate fills
  # every conditioning set, so all time points are still returned.
  res <- check_hdr_seq(
    data,
    c(a1, a2, a3),
    list(tidyselect::starts_with("zzz")),
    .baseline = l0,
    values = 0
  )
  expect_length(sort(unique(res@results$time)), 3)
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

test_that("the point view marks each result with a point", {
  local_quiet()
  # A single target value per group draws no polyline, so a point layer is what
  # keeps the result visible.
  data <- sim_hdr_linear(300, beta = 1, seed = 1)
  res <- check_hdr(data, exposure, l, values = 0)

  built <- ggplot2::ggplot_build(ggplot2::autoplot(res))
  geoms <- vapply(
    built$plot$layers,
    function(ly) class(ly$geom)[1],
    character(1)
  )
  point_index <- which(geoms == "GeomPoint")
  expect_length(point_index, 1)
  expect_identical(nrow(built$data[[point_index]]), nrow(res@results))
})

test_that("the sequential view marks each result with a point", {
  local_quiet()
  data <- dgp_longitudinal(n = 300, seed = 5)
  res <- check_hdr_seq(data, c(a1, a2, a3), list(l0, l1, l2), values = 0)

  built <- ggplot2::ggplot_build(ggplot2::autoplot(res))
  geoms <- vapply(
    built$plot$layers,
    function(ly) class(ly$geom)[1],
    character(1)
  )
  point_index <- which(geoms == "GeomPoint")
  expect_length(point_index, 1)
  expect_identical(nrow(built$data[[point_index]]), nrow(res@results))
})

test_that("plot() draws the HDR view and returns the result invisibly", {
  local_quiet()
  local_null_device()
  data <- sim_hdr_linear(300, beta = 1, seed = 1)
  res <- check_hdr(data, exposure, l)

  expect_identical(plot(res), res)
})

test_that("HDR autoplot views render as expected", {
  local_quiet()
  announce_doppelganger(
    "HDR non-overlap point curve",
    "HDR non-overlap sequential curves"
  )
  data <- sim_hdr_linear(300, beta = 1, seed = 1)
  res_point <- check_hdr(data, exposure, l)
  expect_doppelganger(
    "HDR non-overlap point curve",
    ggplot2::autoplot(res_point)
  )

  long <- dgp_longitudinal(n = 300, seed = 5)
  res_seq <- check_hdr_seq(
    long,
    c(a1, a2, a3),
    list(l0, l1, l2),
    values = seq(-2, 2, length.out = 25)
  )
  expect_doppelganger(
    "HDR non-overlap sequential curves",
    ggplot2::autoplot(res_seq)
  )
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

  expect_snapshot_abort(check_hdr(1:10, exposure, x1))
  expect_snapshot_abort(check_hdr(binary, exposure, c(x1, x2)))
  expect_snapshot_abort(check_hdr(
    data,
    exposure,
    tidyselect::starts_with("zzz")
  ))
  expect_snapshot_abort(check_hdr(data, c(exposure, l), l))
  expect_snapshot_abort(check_hdr(data, exposure, l, mass = 0))
  expect_snapshot_abort(check_hdr(data, exposure, l, mass = 1.5))
  expect_snapshot_abort(check_hdr(data, exposure, l, mass = c(0.9, 0.95)))
  expect_snapshot_abort(check_hdr(data, exposure, l, values = "a"))
  expect_snapshot_abort(check_hdr(data, exposure, l, values = c(0, NA)))
  expect_snapshot_abort(check_hdr(data, exposure, l, values = numeric(0)))
  expect_snapshot_abort(check_hdr(data, exposure, l, values = c(0, Inf)))
  expect_snapshot_abort(check_hdr(
    data,
    exposure,
    l,
    density_estimator = "not an estimator"
  ))

  na_covariate <- data
  na_covariate$l[1] <- NA
  expect_snapshot_abort(check_hdr(na_covariate, exposure, l))

  one_row <- tibble::tibble(exposure = 0.5, l = -0.2)
  expect_snapshot_abort(check_hdr(one_row, exposure, l))

  expect_snapshot_abort(check_hdr(data, exposure, c(exposure, l)))

  constant <- tibble::tibble(
    exposure = rep(3, 50),
    l = seq(-1, 1, length.out = 50)
  )
  expect_snapshot_abort(check_hdr(
    constant,
    exposure,
    l,
    exposure_type = "continuous"
  ))
})

test_that("the flat-at-one non-overlap messages are stable", {
  local_quiet()
  withr::local_options(warn = 0)

  point <- sim_hdr_determined(200, seed = 1)
  expect_snapshot(res <- check_hdr(point, exposure, x, values = c(-1, 0, 1)))

  sequential <- sim_hdr_seq_determined(200, seed = 1)
  expect_snapshot(
    res <- check_hdr_seq(
      sequential,
      c(a1, a2, a3),
      list(l0, l1, l2),
      values = c(-1, 0, 1)
    )
  )
})

test_that("check_hdr_seq() argument validation messages are stable", {
  local_quiet()
  data <- dgp_longitudinal(n = 300, seed = 5)

  expect_snapshot_abort(check_hdr_seq(data, c(a1, a2, a3), list(l0, l1)))
  expect_snapshot_abort(check_hdr_seq(
    data,
    c(a1, a2, a3),
    list(l0, l1, l2),
    lag = -1
  ))
  expect_snapshot_abort(check_hdr_seq(
    data,
    c(a1, a2, a3),
    list(l0, l1, l2),
    lag = 1.5
  ))
  expect_snapshot_abort(check_hdr_seq(
    data,
    c(a1, a2, a3),
    list(l0, l1, l2),
    lag = "a"
  ))

  binary <- data
  binary$a2 <- rep(c(0, 1), length.out = nrow(binary))
  expect_snapshot_abort(check_hdr_seq(binary, c(a1, a2, a3), list(l0, l1, l2)))

  binary_pair <- binary
  binary_pair$a3 <- rep(c(0, 1), length.out = nrow(binary_pair))
  expect_snapshot_abort(check_hdr_seq(
    binary_pair,
    c(a1, a2, a3),
    list(l0, l1, l2)
  ))

  factors <- data
  factors$a2 <- factor(rep(c("low", "high"), length.out = nrow(factors)))
  factors$a3 <- factor(rep(c("low", "high"), length.out = nrow(factors)))
  expect_snapshot_abort(check_hdr_seq(
    factors,
    c(a1, a2, a3),
    list(l0, l1, l2),
    exposure_type = "continuous"
  ))

  mixed <- data
  mixed$a1 <- factor(rep(c("low", "high"), length.out = nrow(mixed)))
  mixed$a2 <- rep(c("low", "high"), length.out = nrow(mixed))
  mixed$a3 <- as.Date("2020-01-01") + seq_len(nrow(mixed))
  expect_snapshot_abort(check_hdr_seq(
    mixed,
    c(a1, a2, a3),
    list(l0, l1, l2),
    exposure_type = "continuous"
  ))

  ranked <- data
  ranked$a1 <- factor(
    rep(c("low", "high"), length.out = nrow(ranked)),
    ordered = TRUE
  )
  ranked$a2 <- factor(rep(c("low", "high"), length.out = nrow(ranked)))
  expect_snapshot_abort(check_hdr_seq(
    ranked,
    c(a1, a2, a3),
    list(l0, l1, l2),
    exposure_type = "continuous"
  ))

  expect_snapshot_abort(check_hdr_seq(
    data,
    c(a1, a2),
    list(c(l0, a1), l1),
    values = 0
  ))
  expect_snapshot_abort(check_hdr_seq(
    data,
    c(a1, a2, a3),
    list(l0, l1, l2),
    .baseline = a2,
    values = 0
  ))
  expect_snapshot_abort(check_hdr_seq(
    data,
    c(a1, a2, a3),
    list(tidyselect::starts_with("zzz")),
    values = 0
  ))

  constant <- data
  constant$a1 <- 2
  constant$a3 <- 7
  expect_snapshot_abort(check_hdr_seq(
    constant,
    c(a1, a2, a3),
    list(l0, l1, l2),
    exposure_type = "continuous",
    values = 0
  ))
})

# ---- Display methods -------------------------------------------------------

# Three intervention targets under the normal density estimator, whose choice
# decides how the non-overlap curve was computed.
hdr_display_result <- function() {
  check_hdr(
    sim_hdr_linear(300, beta = 1, seed = 1),
    exposure,
    l,
    values = c(-1, 0, 1)
  )
}

test_that("diagnostic_label() names the HDR check without naming its class", {
  local_quiet()
  res <- hdr_display_result()
  label <- expect_readable_label(res)

  printed <- printed_text(res)
  expect_no_match(printed, S7::S7_class(res)@name, fixed = TRUE)
  expect_match(printed, label, fixed = TRUE)
})

test_that("diagnostic_headline() reads the estimator and the target count", {
  local_quiet()
  res <- hdr_display_result()
  headline <- expect_readable_headline(res)
  text <- rendered_text(headline)

  expect_match(text, "normal", fixed = TRUE)
  expect_match(text, "\\b3\\b")
})

test_that("the HDR label and headline are stable", {
  local_quiet()
  res <- hdr_display_result()
  expect_snapshot({
    diagnostic_label(res)
    writeLines(diagnostic_headline(res))
  })
})
