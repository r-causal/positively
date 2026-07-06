# Builds the two lazy-loaded datasets that ship with positively. Run this
# script interactively to regenerate data/pos_violations.rda and
# data/pos_violations_long.rda. A fixed seed per dataset keeps the output
# deterministic and byte-identical across runs.

# ---- pos_violations --------------------------------------------------------

set.seed(2024L)

n_point <- 1000L
x1 <- stats::rnorm(n_point)
x2 <- stats::rnorm(n_point)
region <- factor(sample(c("a", "b"), n_point, replace = TRUE))

# A steep coefficient on x1 drives the tails toward near-deterministic
# treatment, the practical near-violation.
ps <- stats::plogis(3 * x1 - 0.5 * x2)

# The structural subgroup: region "b" with x2 above 1 can never be exposed.
structural <- region == "b" & x2 > 1
ps[structural] <- 0

exposure <- stats::rbinom(n_point, 1L, ps)

pos_violations <- tibble::tibble(
  exposure = exposure,
  x1 = x1,
  x2 = x2,
  region = region
)

# ---- pos_violations_long ---------------------------------------------------

set.seed(2025L)

n_long <- 500L
l0 <- stats::rnorm(n_long)
a1 <- stats::rnorm(n_long, mean = l0)
l1 <- stats::rnorm(n_long, mean = a1)

# The time-2 exposure is drawn on [0, 2] or [4, 6], leaving a support gap on
# the open interval (2, 4) that a1 and a3 populate.
low <- stats::runif(n_long, 0, 2)
high <- stats::runif(n_long, 4, 6)
in_low <- stats::runif(n_long) < 0.5
a2 <- ifelse(in_low, low, high)

l2 <- stats::rnorm(n_long, mean = a2)
a3 <- stats::rnorm(n_long, mean = l2)

pos_violations_long <- tibble::tibble(
  id = seq_len(n_long),
  l0 = l0,
  a1 = a1,
  l1 = l1,
  a2 = a2,
  l2 = l2,
  a3 = a3
)

usethis::use_data(
  pos_violations,
  pos_violations_long,
  overwrite = TRUE
)
