# Seeded data-generating processes with known positivity properties.
#
# These simulators are test infrastructure. Each documents the structure it
# plants so tests can assert that a diagnostic recovers it. Every helper takes a
# seed and uses withr::local_seed() so draws are deterministic and the global
# RNG state is left untouched.

# Good positivity: the true propensity score is bounded away from 0 and 1, so
# every covariate stratum has both exposure levels represented.
dgp_good_positivity <- function(n = 500, seed = 1) {
  withr::local_seed(seed)
  x1 <- stats::rnorm(n)
  x2 <- stats::rnorm(n)
  # Squeeze the propensity into (0.2, 0.8): plogis() lands in (0, 1) and the
  # affine transform bounds it.
  ps <- 0.2 + 0.6 * stats::plogis(0.5 * x1 - 0.5 * x2)
  exposure <- stats::rbinom(n, 1L, ps)
  tibble::tibble(exposure = exposure, x1 = x1, x2 = x2, ps = ps)
}

# Practical positivity violation (Petersen style): a steep propensity model
# pushes the tails of the covariate distribution to near-deterministic
# treatment, so some fitted scores fall below 0.05 and some exceed 0.95.
dgp_practical_violation <- function(n = 500, seed = 2) {
  withr::local_seed(seed)
  x1 <- stats::rnorm(n)
  x2 <- stats::rnorm(n)
  ps <- stats::plogis(3 * x1)
  exposure <- stats::rbinom(n, 1L, ps)
  tibble::tibble(exposure = exposure, x1 = x1, x2 = x2, ps = ps)
}

# Structural violation: a planted subgroup rule where treatment is impossible.
# No subject with x1 > 2 and x2 == "b" is ever treated (true propensity 0).
dgp_structural_subgroup <- function(n = 1000, seed = 3) {
  withr::local_seed(seed)
  x1 <- stats::rnorm(n)
  x2 <- factor(sample(c("a", "b"), n, replace = TRUE))
  ps <- stats::plogis(0.5 * x1)
  in_gap <- x1 > 2 & x2 == "b"
  ps[in_gap] <- 0
  exposure <- stats::rbinom(n, 1L, ps)
  tibble::tibble(exposure = exposure, x1 = x1, x2 = x2, ps = ps)
}

# Continuous exposure with a hard support gap: the exposure is drawn on
# [0, 2] or [4, 6], so the interval (2, 4) is never observed.
dgp_continuous_support_gap <- function(n = 500, seed = 4) {
  withr::local_seed(seed)
  low <- stats::runif(n, 0, 2)
  high <- stats::runif(n, 4, 6)
  in_low <- stats::runif(n) < 0.5
  exposure <- ifelse(in_low, low, high)
  x1 <- 0.5 * exposure + stats::rnorm(n)
  tibble::tibble(exposure = exposure, x1 = x1)
}

# Longitudinal wide-format data over three time points with a per-time
# violation: the time-2 exposure has a hard support gap on (2, 4), while the
# time-1 and time-3 exposures are unconstrained.
dgp_longitudinal <- function(n = 500, seed = 5) {
  withr::local_seed(seed)
  l0 <- stats::rnorm(n)
  a1 <- stats::rnorm(n, mean = l0)
  l1 <- stats::rnorm(n, mean = a1)
  low <- stats::runif(n, 0, 2)
  high <- stats::runif(n, 4, 6)
  in_low <- stats::runif(n) < 0.5
  a2 <- ifelse(in_low, low, high)
  l2 <- stats::rnorm(n, mean = a2)
  a3 <- stats::rnorm(n, mean = l2)
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

# Longitudinal wide-format data over three time points with monotone binary
# exposures, built so the sPoRT follower and risk-set mechanics are visible.
# Once treated, a subject stays treated, so the risk set of subjects still
# following the regime shrinks across time. Two deterministic regions are
# planted at t = 2 among the followers: treatment never initiates where the
# current covariate l1 is above 1, and never initiates where the lagged
# covariate l0 is above 1. Each region is defined by a single covariate, so the
# reading rule surfaces both. The l0 region enters the t = 2 conditioning set
# only when the lag window reaches back one time point, so it is flagged under
# the full history and absent under a lag of zero.
dgp_longitudinal_binary <- function(n = 2000, seed = 6) {
  withr::local_seed(seed)
  l0 <- stats::rnorm(n)
  a1 <- stats::rbinom(n, 1L, 0.3)
  l1 <- stats::rnorm(n, mean = 0.3 * a1)

  a2 <- a1
  followers2 <- a1 == 0
  p2 <- rep(0.5, n)
  p2[l0 > 1] <- 0
  p2[l1 > 1] <- 0
  a2[followers2] <- stats::rbinom(sum(followers2), 1L, p2[followers2])

  l2 <- stats::rnorm(n, mean = 0.3 * a2)

  a3 <- a2
  followers3 <- a2 == 0
  a3[followers3] <- stats::rbinom(sum(followers3), 1L, 0.4)

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

# Longitudinal wide-format data over three time points with monotone binary
# exposures and per-wave censoring events c1, c2, c3. Each is coded 1 when the
# subject is censored at that wave among those still uncensored going in, and 0
# otherwise; a subject censored at an earlier wave carries 0 at every later wave
# because it has already left the risk set. It plants a censoring violation at
# t = 2: among the
# subjects uncensored through t = 1, censoring is near-certain where the baseline
# covariate l0 is above 1, while it is light everywhere else and at the first and
# last time points. Because censoring is analyzed on the uncensored-so-far risk
# set, that set shrinks across time as earlier censoring accrues. The exposure
# side carries its own planted region at t = 2 where treatment never initiates
# for l1 above 1, kept distinct from the l0 censoring region so exposure and
# censoring subgroups are separable in the tidy output.
dgp_longitudinal_binary_censoring <- function(n = 2000, seed = 7) {
  withr::local_seed(seed)
  l0 <- stats::rnorm(n)
  a1 <- stats::rbinom(n, 1L, 0.3)
  c1 <- stats::rbinom(n, 1L, 0.08)

  l1 <- stats::rnorm(n, mean = 0.3 * a1)

  a2 <- a1
  followers2 <- a1 == 0
  p2 <- rep(0.5, n)
  p2[l1 > 1] <- 0
  a2[followers2] <- stats::rbinom(sum(followers2), 1L, p2[followers2])

  uncensored1 <- c1 == 0
  q2 <- rep(0.1, n)
  q2[l0 > 1] <- 0.98
  c2 <- rep(0L, n)
  c2[uncensored1] <- stats::rbinom(sum(uncensored1), 1L, q2[uncensored1])

  l2 <- stats::rnorm(n, mean = 0.3 * a2)

  a3 <- a2
  followers3 <- a2 == 0
  a3[followers3] <- stats::rbinom(sum(followers3), 1L, 0.4)

  uncensored2 <- c1 == 0 & c2 == 0
  c3 <- rep(0L, n)
  c3[uncensored2] <- stats::rbinom(sum(uncensored2), 1L, 0.1)

  tibble::tibble(
    id = seq_len(n),
    l0 = l0,
    a1 = a1,
    c1 = c1,
    l1 = l1,
    a2 = a2,
    c2 = c2,
    l2 = l2,
    a3 = a3,
    c3 = c3
  )
}
