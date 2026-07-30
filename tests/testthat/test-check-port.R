# PoRT reads a fixed rule off a regression tree: a terminal subgroup is flagged
# when its exposure prevalence falls below beta or above 1 - beta AND its size is
# at least alpha of the sample. The reading rule is deterministic given the tree;
# only the tree carries stochasticity, so magnitude claims use tolerances while
# the hand-built anchor uses exact fields. Every assertion about which subgroups
# were found is written against the flagged-row subset of @results, so it holds
# whether @results reports only flagged subgroups or every terminal leaf with its
# status.

# Rows of @results that met the PoRT reading rule.
flagged_rows <- function(res) {
  res@results[res@results$flagged, , drop = FALSE]
}

# All numbers appearing in a rule description, used only where a test asserts a
# recovered split threshold within tolerance.
rule_numbers <- function(x) {
  as.numeric(unlist(regmatches(x, gregexpr("[0-9]+\\.?[0-9]*", x))))
}

# A directly constructed port_result, used to exercise the class validator and
# the print and plot branches that depend on property shapes a fitted run does
# not readily produce.
make_port_result <- function(alpha = 0.05, beta = 0.05, gamma = 2) {
  port_result(
    results = empty_port_tibble(sequential = FALSE),
    exposure = "exposure",
    exposure_type = "binary",
    n = 1L,
    params = list(),
    call = quote(check_port()),
    trees = list(),
    alpha = alpha,
    beta = beta,
    gamma = gamma
  )
}

# ---- Local scenario generators --------------------------------------------
# The shared helpers in helper-dgp.R fix a single exposure-covariate dependence.
# PoRT tests need to dial subgroup size, subgroup prevalence, and the joint
# structure independently, so these seeded generators live here, mirroring the
# per-file generator idiom used by the hat-values and HDR tests.

# Deterministic reading-rule anchor (handcalc.json). One binary predictor g.
# g == 1: 300 rows, exactly 3 treated (prevalence 0.010, size 0.300).
# g == 0: 700 rows, exactly 350 treated (prevalence 0.500). No randomness.
port_anchor_data <- function() {
  g <- c(rep(1L, 300), rep(0L, 700))
  exposure <- c(
    rep(1L, 3),
    rep(0L, 297), # g == 1: 3 of 300 treated
    rep(1L, 350),
    rep(0L, 350) # g == 0: 350 of 700 treated
  )
  tibble::tibble(exposure = exposure, g = g)
}

# Single-covariate structural violation: never treated when x3 > 70. The planted
# threshold is 70; recovered thresholds cluster just above it.
sim_port_structural <- function(n = 2000, seed = 1) {
  withr::local_seed(seed)
  x3 <- stats::runif(n, 0, 100)
  x1 <- stats::rnorm(n)
  exposure <- stats::rbinom(n, 1L, 0.5)
  exposure[x3 > 70] <- 0L
  tibble::tibble(exposure = exposure, x3 = x3, x1 = x1)
}

# A binary covariate z marks a never-treated subgroup covering `viol_frac` of the
# sample; elsewhere the prevalence is 0.5. Used for the alpha size gate.
sim_port_sized <- function(n = 3000, viol_frac = 0.06, seed = 1) {
  withr::local_seed(seed)
  z <- as.integer(stats::runif(n) < viol_frac)
  ps <- ifelse(z == 1L, 0, 0.5)
  exposure <- stats::rbinom(n, 1L, ps)
  tibble::tibble(exposure = exposure, z = z, x1 = stats::rnorm(n))
}

# A binary covariate b marks a quarter of the sample in which the exposure
# prevalence is a tunable near-beta value; elsewhere the prevalence is 0.5. Used
# for the beta bound, where 0.051 straddles beta = 0.05.
sim_port_prev <- function(n = 3000, prev = 0.03, frac = 0.25, seed = 1) {
  withr::local_seed(seed)
  b <- as.integer(stats::runif(n) < frac)
  ps <- ifelse(b == 1L, prev, 0.5)
  exposure <- stats::rbinom(n, 1L, ps)
  tibble::tibble(exposure = exposure, b = b, x1 = stats::rnorm(n))
}

# Joint violation: never treated when x1 > 0 AND x2 == "b". Neither predictor
# alone isolates a violating subgroup, so gamma must reach 2 to find it.
sim_port_joint <- function(n = 3000, seed = 1) {
  withr::local_seed(seed)
  x1 <- stats::rnorm(n)
  x2 <- factor(sample(c("a", "b"), n, replace = TRUE))
  in_gap <- x1 > 0 & x2 == "b"
  ps <- ifelse(in_gap, 0, 0.5)
  exposure <- stats::rbinom(n, 1L, ps)
  tibble::tibble(exposure = exposure, x1 = x1, x2 = x2)
}

# Three-level exposure where the "high" level is never assigned in a subgroup
# (x2 == "b"), exercising the per-level (one-versus-rest) categorical check.
sim_port_categorical <- function(n = 3000, seed = 1) {
  withr::local_seed(seed)
  x2 <- factor(sample(c("a", "b"), n, replace = TRUE))
  draw_level <- function(is_b) {
    probs <- if (is_b) {
      c(low = 0.5, medium = 0.5, high = 0)
    } else {
      c(low = 1 / 3, medium = 1 / 3, high = 1 / 3)
    }
    sample(names(probs), 1, prob = probs)
  }
  level <- vapply(x2 == "b", draw_level, character(1))
  tibble::tibble(
    exposure = factor(level, levels = c("low", "medium", "high")),
    x2 = x2
  )
}

# A binary exposure with one never-treated subgroup and one always-treated
# subgroup, so both tails of the reading rule fire.
sim_port_extremes <- function(n = 3000, seed = 1) {
  withr::local_seed(seed)
  w <- factor(sample(c("a", "b", "c"), n, replace = TRUE))
  ps <- c(a = 0.5, b = 0, c = 1)[as.character(w)]
  exposure <- stats::rbinom(n, 1L, ps)
  tibble::tibble(exposure = exposure, w = w)
}

# Continuous predictor carrying a hard threshold violation at 0.70, used to show
# that a coarse quantile categorization of the predictor dilutes and misses it.
sim_port_granularity <- function(n = 3000, seed = 1) {
  withr::local_seed(seed)
  x3 <- stats::runif(n)
  exposure <- stats::rbinom(n, 1L, 0.5)
  exposure[x3 > 0.70] <- 0L
  tibble::tibble(exposure = exposure, x3 = x3)
}

# ---- Argument validation ---------------------------------------------------

test_that("check_port() rejects non-data-frame input", {
  local_quiet()
  expect_error(
    check_port(1:10, exposure, x1),
    class = "positively_error"
  )
})

test_that("check_port() rejects an empty covariate selection", {
  local_quiet()
  data <- dgp_good_positivity(n = 200, seed = 1)
  expect_error(
    check_port(data, exposure, tidyselect::starts_with("zzz")),
    class = "positively_error"
  )
})

test_that("check_port() requires a single exposure column", {
  local_quiet()
  data <- dgp_good_positivity(n = 200, seed = 1)
  expect_error(
    check_port(data, c(exposure, x1), x2),
    class = "positively_error"
  )
})

test_that("check_port() aborts on a single-level exposure", {
  local_quiet()
  data <- dgp_good_positivity(n = 200, seed = 1)
  data$exposure <- 1L
  expect_error(
    check_port(data, exposure, c(x1, x2)),
    class = "positively_error"
  )
})

test_that("declaring binary on a three-level exposure aborts rather than dropping levels", {
  local_quiet()
  # port_level_responses() reduces a binary exposure to the indicator of its
  # highest level alone, so an unchecked binary declaration would examine one
  # level of three and report a third of the rows the same column yields when
  # read as categorical, with no sign that two levels were never looked at.
  data <- sim_port_categorical(n = 600, seed = 1)

  categorical <- check_port(data, exposure, x2, exposure_type = "categorical")
  expect_setequal(
    unique(categorical@results$exposure_level),
    c("low", "medium", "high")
  )

  expect_error(
    check_port(data, exposure, x2, exposure_type = "binary"),
    class = "positively_exposure_type_error"
  )
})

test_that("check_port() rejects alpha outside the unit interval", {
  local_quiet()
  data <- dgp_good_positivity(n = 200, seed = 1)
  expect_error(
    check_port(data, exposure, c(x1, x2), alpha = 0),
    class = "positively_error"
  )
  expect_error(
    check_port(data, exposure, c(x1, x2), alpha = 1.5),
    class = "positively_error"
  )
  expect_error(
    check_port(data, exposure, c(x1, x2), alpha = c(0.05, 0.1)),
    class = "positively_error"
  )
})

test_that("check_port() rejects a numeric beta outside its range", {
  local_quiet()
  data <- dgp_good_positivity(n = 200, seed = 1)
  expect_error(
    check_port(data, exposure, c(x1, x2), beta = 0),
    class = "positively_error"
  )
  expect_error(
    check_port(data, exposure, c(x1, x2), beta = -0.1),
    class = "positively_error"
  )
  expect_error(
    check_port(data, exposure, c(x1, x2), beta = 1.5),
    class = "positively_error"
  )
})

test_that("check_port() accepts beta = \"gruber\" but rejects other strings", {
  local_quiet()
  data <- dgp_good_positivity(n = 400, seed = 1)
  # "gruber" is the one permitted string; it must not raise.
  expect_no_error(check_port(data, exposure, c(x1, x2), beta = "gruber"))
  expect_error(
    check_port(data, exposure, c(x1, x2), beta = "banana"),
    class = "positively_error"
  )
})

test_that("check_port() rejects a gamma that is not a positive whole number", {
  local_quiet()
  data <- dgp_good_positivity(n = 200, seed = 1)
  expect_error(
    check_port(data, exposure, c(x1, x2), gamma = 0),
    class = "positively_error"
  )
  expect_error(
    check_port(data, exposure, c(x1, x2), gamma = 2.5),
    class = "positively_error"
  )
  expect_error(
    check_port(data, exposure, c(x1, x2), gamma = -1),
    class = "positively_error"
  )
})

test_that("check_port() rejects a bad n_bins for continuous exposures", {
  local_quiet()
  data <- dgp_continuous_support_gap(n = 300, seed = 4)
  expect_error(
    check_port(data, exposure, x1, n_bins = 0),
    class = "positively_error"
  )
  # A single bin makes every observation share one level, flagging every
  # subgroup, so n_bins must be at least two.
  expect_error(
    check_port(data, exposure, x1, n_bins = 1),
    class = "positively_error"
  )
  expect_error(
    check_port(data, exposure, x1, n_bins = 2.5),
    class = "positively_error"
  )
})

test_that("check_port() rejects non-numeric breaks", {
  local_quiet()
  data <- dgp_continuous_support_gap(n = 300, seed = 4)
  expect_error(
    check_port(data, exposure, x1, breaks = "a"),
    class = "positively_error"
  )
})

# ---- Class validator ------------------------------------------------------

test_that("port_result rejects a non-scalar alpha", {
  expect_error(make_port_result(alpha = c(0.05, 0.1)), regexp = "@alpha")
})

test_that("port_result rejects an empty beta", {
  expect_error(make_port_result(beta = double(0)), regexp = "@beta")
})

test_that("port_result rejects a non-scalar gamma", {
  expect_error(make_port_result(gamma = c(1, 2)), regexp = "@gamma")
})

# ---- Assembly helpers -----------------------------------------------------

test_that("empty_port_tibble carries the fixed schema", {
  point <- empty_port_tibble(sequential = FALSE)
  expect_s3_class(point, "tbl_df")
  expect_identical(nrow(point), 0L)
  expect_setequal(
    names(point),
    c(
      "subgroup",
      "description",
      "exposure_level",
      "n",
      "proportion",
      "prevalence",
      "flagged"
    )
  )
  expect_type(point$flagged, "logical")

  sequential <- empty_port_tibble(sequential = TRUE)
  expect_identical(names(sequential)[[1]], "time")
  expect_type(sequential$time, "integer")
})

test_that("assemble_port_results returns the empty schema when no rows survive", {
  from_nothing <- assemble_port_results(list(), sequential = FALSE)
  expect_identical(nrow(from_nothing), 0L)
  expect_true("flagged" %in% names(from_nothing))

  from_empties <- assemble_port_results(
    list(empty_port_tibble(sequential = TRUE)),
    sequential = TRUE
  )
  expect_identical(nrow(from_empties), 0L)
  expect_identical(names(from_empties)[[1]], "time")
})

# ---- Result class and properties ------------------------------------------

test_that("check_port() returns a port_result diagnostic", {
  local_quiet()
  data <- dgp_good_positivity(n = 500, seed = 1)
  res <- check_port(data, exposure, c(x1, x2))

  expect_true(S7::S7_inherits(res, positivity_diagnostic))
  expect_identical(S7::S7_class(res)@name, "port_result")
  expect_identical(res@exposure_type, "binary")
  expect_identical(res@n, 500L)
})

test_that("port_result carries the alpha, beta, and gamma properties", {
  local_quiet()
  data <- dgp_good_positivity(n = 500, seed = 1)
  res <- check_port(
    data,
    exposure,
    c(x1, x2),
    alpha = 0.1,
    beta = 0.08,
    gamma = 3
  )

  expect_identical(res@alpha, 0.1)
  expect_identical(res@beta, 0.08)
  expect_identical(res@gamma, 3)
})

test_that("port_result@trees is a list of rpart objects", {
  local_quiet()
  data <- dgp_good_positivity(n = 500, seed = 1)
  res <- check_port(data, exposure, c(x1, x2))

  expect_type(res@trees, "list")
  expect_gte(length(res@trees), 1)
  expect_true(all(vapply(
    res@trees,
    function(tree) inherits(tree, "rpart"),
    logical(1)
  )))
})

test_that("port_result has the fixed point results columns", {
  local_quiet()
  data <- dgp_good_positivity(n = 500, seed = 1)
  res <- check_port(data, exposure, c(x1, x2))

  expect_s3_class(res@results, "tbl_df")
  expect_setequal(
    names(res@results),
    c(
      "subgroup",
      "description",
      "exposure_level",
      "n",
      "proportion",
      "prevalence",
      "flagged"
    )
  )
  expect_type(res@results$flagged, "logical")
})

test_that("tidy() returns the results tibble and glance() a one-row summary", {
  local_quiet()
  data <- dgp_good_positivity(n = 500, seed = 1)
  res <- check_port(data, exposure, c(x1, x2))

  expect_identical(generics::tidy(res), res@results)
  glanced <- generics::glance(res)
  expect_s3_class(glanced, "tbl_df")
  expect_identical(nrow(glanced), 1L)
})

# ---- Deterministic reading-rule anchor ------------------------------------

test_that("the reading rule flags exactly the hand-built subgroup", {
  local_quiet()
  data <- port_anchor_data()
  res <- check_port(data, exposure, g, alpha = 0.05, beta = 0.05)

  flagged <- flagged_rows(res)
  expect_identical(nrow(flagged), 1L)
  expect_identical(flagged$n, 300L)
  expect_equal(flagged$proportion, 0.30, tolerance = 1e-6)
  expect_equal(flagged$prevalence, 0.01, tolerance = 1e-6)
  expect_true(flagged$flagged)
  expect_match(flagged$description, "g")
})

test_that("the alpha boundary gates the anchor subgroup out", {
  local_quiet()
  data <- port_anchor_data()
  # The anchor subgroup has size 0.30, so alpha = 0.31 excludes it.
  res <- check_port(data, exposure, g, alpha = 0.31, beta = 0.05)
  expect_identical(nrow(flagged_rows(res)), 0L)
})

test_that("the beta boundary gates the anchor subgroup out", {
  local_quiet()
  data <- port_anchor_data()
  # The anchor subgroup has prevalence 0.01, so beta = 0.005 no longer flags it.
  res <- check_port(data, exposure, g, alpha = 0.05, beta = 0.005)
  expect_identical(nrow(flagged_rows(res)), 0L)
})

# ---- Structural recovery --------------------------------------------------

test_that("check_port() recovers a planted single-covariate violation", {
  local_quiet()
  data <- sim_port_structural(n = 2000, seed = 1)
  res <- check_port(data, exposure, c(x3, x1))

  flagged <- flagged_rows(res)
  on_x3 <- grepl("x3", flagged$description)
  expect_true(any(on_x3))
  expect_true(all(flagged$prevalence[on_x3] < res@beta))
})

test_that("the recovered split threshold sits near the planted value", {
  local_quiet()
  data <- sim_port_structural(n = 2000, seed = 1)
  res <- check_port(data, exposure, c(x3, x1))

  flagged <- flagged_rows(res)
  on_x3 <- flagged[grepl("x3", flagged$description), , drop = FALSE]
  # The planted threshold is 70; recovered splits cluster just above it. If the
  # description exposes the split, at least one number should sit within 1 of 70.
  numbers <- unlist(lapply(on_x3$description, rule_numbers))
  near_planted <- any(abs(numbers - 70) <= 1)
  expect_true(near_planted)
})

# ---- Quiet on good positivity ---------------------------------------------

test_that("check_port() stays quiet under good positivity at n = 1000", {
  local_quiet()
  data <- dgp_good_positivity(n = 1000, seed = 1)
  res <- check_port(data, exposure, c(x1, x2), alpha = 0.05, beta = 0.05)
  expect_equal(sum(res@results$flagged), 0)
})

test_that("average flags stay low across seeds under good positivity", {
  local_quiet()
  counts <- vapply(
    1:5,
    function(s) {
      res <- check_port(
        dgp_good_positivity(n = 1000, seed = s),
        exposure,
        c(x1, x2)
      )
      sum(res@results$flagged)
    },
    numeric(1)
  )
  expect_lte(mean(counts), 0.1)
})

# ---- Direction-only n-scaling ---------------------------------------------

test_that("spurious flags are more common at small n than large n", {
  local_quiet()
  mean_flags <- function(n) {
    mean(vapply(
      1:20,
      function(s) {
        res <- check_port(
          dgp_good_positivity(n = n, seed = s),
          exposure,
          c(x1, x2)
        )
        sum(res@results$flagged)
      },
      numeric(1)
    ))
  }
  expect_gt(mean_flags(150), mean_flags(3000))
})

# ---- alpha size gate ------------------------------------------------------

test_that("a subgroup near 6 percent is gated by alpha", {
  local_quiet()
  data <- sim_port_sized(n = 3000, viol_frac = 0.06, seed = 1)

  at_05 <- check_port(data, exposure, c(z, x1), alpha = 0.05)
  at_10 <- check_port(data, exposure, c(z, x1), alpha = 0.10)

  expect_true(any(grepl("z", flagged_rows(at_05)$description)))
  expect_false(any(grepl("z", flagged_rows(at_10)$description)))
})

test_that("flag counts are non-increasing across an increasing alpha grid", {
  local_quiet()
  data <- sim_port_sized(n = 3000, viol_frac = 0.06, seed = 1)
  alphas <- c(0.01, 0.05, 0.10, 0.20)
  counts <- vapply(
    alphas,
    function(a) {
      sum(check_port(data, exposure, c(z, x1), alpha = a)@results$flagged)
    },
    numeric(1)
  )
  expect_true(all(diff(counts) <= 0))
})

# ---- beta bound -----------------------------------------------------------

test_that("a near-beta subgroup is flagged reliably at the looser bound", {
  local_quiet()
  # Prevalence 0.051 straddles beta = 0.05, so it is flagged at beta = 0.10 but
  # only a coin flip at 0.05. Assert the reliable side.
  data <- sim_port_prev(n = 4000, prev = 0.051, seed = 1)
  at_10 <- check_port(data, exposure, c(b, x1), beta = 0.10)
  expect_true(any(grepl("b", flagged_rows(at_10)$description)))
})

test_that("detection is non-decreasing across an increasing beta grid", {
  local_quiet()
  # Descendant suppression means the raw flag count is not monotone in beta, so
  # the monotone signal is detection (any flag), not the count.
  data <- sim_port_prev(n = 4000, prev = 0.051, seed = 1)
  betas <- c(0.02, 0.05, 0.10, 0.20)
  detected <- vapply(
    betas,
    function(bt) {
      any(check_port(data, exposure, c(b, x1), beta = bt)@results$flagged)
    },
    logical(1)
  )
  expect_true(all(diff(as.integer(detected)) >= 0))
})

# ---- gamma unlocks joint violations ---------------------------------------

test_that("a joint violation needs gamma = 2 to be found", {
  local_quiet()
  data <- sim_port_joint(n = 3000, seed = 1)

  uses_both <- function(res) {
    desc <- flagged_rows(res)$description
    any(grepl("x1", desc) & grepl("x2", desc))
  }

  expect_false(uses_both(check_port(data, exposure, c(x1, x2), gamma = 1)))
  expect_true(uses_both(check_port(data, exposure, c(x1, x2), gamma = 2)))
})

# ---- Gruber default beta --------------------------------------------------

test_that("beta = \"gruber\" resolves to 5 / (sqrt(n) * log(n))", {
  local_quiet()
  data <- dgp_good_positivity(n = 1000, seed = 1)
  res <- check_port(data, exposure, c(x1, x2), beta = "gruber")

  expected <- 5 / (sqrt(1000) * log(1000))
  expect_equal(res@beta, expected, tolerance = 1e-8)
  # 0.0229 is what that formula evaluates to at n = 1000. Pinning the number
  # separately means a change to the expression above cannot pass by moving
  # both sides of the equality together.
  expect_equal(res@beta, 0.0229, tolerance = 1e-3)
})

test_that("the Gruber bound scales with n", {
  local_quiet()
  res4000 <- check_port(
    dgp_good_positivity(n = 4000, seed = 1),
    exposure,
    c(x1, x2),
    beta = "gruber"
  )
  expect_equal(res4000@beta, 5 / (sqrt(4000) * log(4000)), tolerance = 1e-8)
})

# ---- Degenerate resolved Gruber bound -------------------------------------

# A small sample pushes the Gruber bound to or above 0.5, where the reading rule
# flags every subgroup because a prevalence below beta and above 1 - beta no
# longer bracket a gap. A resolved bound of 0.5 or more is a hard error, just as
# a numeric beta of 0.5 or more is rejected up front.
port_gruber_data <- function(n, zeros) {
  withr::with_seed(
    3,
    tibble::tibble(
      exposure = c(rep(0L, zeros), rep(1L, n - zeros)),
      x1 = stats::rnorm(n)
    )
  )
}

test_that("a resolved Gruber bound of 0.5 or more is rejected", {
  local_quiet()
  # n = 12 resolves the Gruber bound to about 0.581, at or above the 0.5 limit.
  data <- port_gruber_data(n = 12, zeros = 7)
  expect_error(
    check_port(data, exposure, x1, beta = "gruber"),
    class = "positively_range_error"
  )
})

test_that("the degenerate Gruber bound message pins the bound and sample size", {
  local_quiet()
  data <- port_gruber_data(n = 12, zeros = 7)
  expect_snapshot_abort(check_port(data, exposure, x1, beta = "gruber"))
})

test_that("a Gruber bound below 0.5 still runs", {
  local_quiet()
  # n = 15 resolves the Gruber bound to about 0.477, below the 0.5 limit.
  data <- port_gruber_data(n = 15, zeros = 9)
  res <- check_port(data, exposure, x1, beta = "gruber")
  expect_identical(S7::S7_class(res)@name, "port_result")
  expect_equal(res@beta, 5 / (sqrt(15) * log(15)), tolerance = 1e-8)
})

# ---- Practical and structural flagged identically -------------------------

test_that("a low-prevalence subgroup is flagged like an empty one", {
  local_quiet()
  empty <- sim_port_prev(n = 3000, prev = 0.00, seed = 1)
  practical <- sim_port_prev(n = 3000, prev = 0.03, seed = 1)

  expect_true(any(grepl(
    "b",
    flagged_rows(check_port(empty, exposure, c(b, x1)))$description
  )))
  expect_true(any(grepl(
    "b",
    flagged_rows(check_port(practical, exposure, c(b, x1)))$description
  )))
})

test_that("a subgroup above beta is not flagged", {
  local_quiet()
  # Prevalence 0.12 is above beta = 0.05 and below 1 - beta, so it never flags.
  data <- sim_port_prev(n = 3000, prev = 0.12, seed = 1)
  res <- check_port(data, exposure, c(b, x1), beta = 0.05)
  expect_false(any(grepl("b", flagged_rows(res)$description)))
})

# ---- Categorical exposure per level ---------------------------------------

test_that("a never-assigned level in a subgroup is flagged at that level", {
  local_quiet()
  data <- sim_port_categorical(n = 3000, seed = 1)
  res <- check_port(data, exposure, x2)

  expect_identical(res@exposure_type, "categorical")
  flagged <- flagged_rows(res)
  high_flag <- flagged[flagged$exposure_level == "high", , drop = FALSE]
  expect_gte(nrow(high_flag), 1)
  expect_true(all(high_flag$prevalence < res@beta))
  expect_true(any(grepl("x2", high_flag$description)))
})

test_that("both extremes of a binary exposure are flagged", {
  local_quiet()
  data <- sim_port_extremes(n = 3000, seed = 1)
  res <- check_port(data, exposure, w)

  flagged <- flagged_rows(res)
  expect_true(any(flagged$prevalence < res@beta))
  expect_true(any(flagged$prevalence > 1 - res@beta))
})

# ---- Threshold granularity failure mode -----------------------------------

test_that("a continuous predictor recovers a threshold a coarse binning misses", {
  local_quiet()
  data <- sim_port_granularity(n = 3000, seed = 1)

  # The exposure is binary; n_bins applies to a continuous exposure, so here the
  # granularity caveat is exercised on the predictor by pre-binning it coarsely.
  fine <- check_port(data, exposure, x3)
  expect_true(any(grepl("x3", flagged_rows(fine)$description)))

  coarse <- data
  # A two-way split at 0.5 puts treated (0.5 to 0.7) and never-treated (above
  # 0.7) in one category, diluting its prevalence to about 0.2 so the reading
  # rule can no longer isolate the violation.
  coarse$x3 <- cut(data$x3, breaks = c(-Inf, 0.5, Inf), labels = FALSE)
  coarse_res <- check_port(coarse, exposure, x3)
  expect_false(any(grepl("x3", flagged_rows(coarse_res)$description)))
})

# ---- Continuous binning ----------------------------------------------------

test_that("explicit breaks bin a continuous exposure", {
  local_quiet()
  data <- dgp_continuous_support_gap(n = 500, seed = 4)
  # The exposure occupies [0, 2] and [4, 6]; explicit cut points override the
  # quantile binning, and the empty (2, 4) bin simply contributes no level.
  res <- check_port(
    data,
    exposure,
    x1,
    breaks = c(0, 2, 4, 6),
    exposure_type = "continuous"
  )
  expect_identical(res@exposure_type, "continuous")
  expect_gte(nrow(res@results), 1L)

  # The lowest exposure bin "[0,2]" carries a covariate subgroup on x1 whose
  # exposure prevalence falls below beta, so it flags with an x1 description.
  low_flag <- flagged_rows(res)
  low_flag <- low_flag[low_flag$exposure_level == "[0,2]", , drop = FALSE]
  expect_gte(nrow(low_flag), 1L)
  expect_true(any(grepl("x1", low_flag$description)))
})

test_that("quantile bins recover the continuous support gap", {
  local_quiet()
  # Default quantile binning cuts the exposure into three levels. Because
  # x1 = 0.5 * exposure + noise, high-x1 subjects concentrate in the highest
  # exposure bin and are nearly absent from the lowest, so the lowest bin flags a
  # subgroup carrying a high x1 lower bound.
  data <- dgp_continuous_support_gap(n = 1000, seed = 4)
  res <- check_port(data, exposure, x1, exposure_type = "continuous")

  flagged <- flagged_rows(res)
  expect_identical(nrow(flagged), 2L)
  # Every flagged subgroup sits below beta.
  expect_true(all(flagged$prevalence < res@beta))

  # An interval label like "[0.00102,1.27]" whose bounds both fall in [0, 2] is a
  # low exposure bin. At least one flagged low bin carries a high x1 lower bound.
  bounds <- lapply(
    flagged$exposure_level,
    function(level) {
      as.numeric(regmatches(
        level,
        gregexpr("-?[0-9]+\\.?[0-9]*", level)
      )[[1]])
    }
  )
  low_bin <- vapply(bounds, function(b) all(b >= 0 & b <= 2), logical(1))
  high_x1 <- grepl("x1>=", flagged$description, fixed = TRUE)
  expect_true(any(low_bin & high_x1))
})

# ---- Autoplot contract -----------------------------------------------------

test_that("autoplot() returns a ggplot", {
  local_quiet()
  data <- sim_port_structural(n = 500, seed = 1)
  res <- check_port(data, exposure, c(x3, x1))
  expect_s3_class(ggplot2::autoplot(res), "ggplot")
})

test_that("PoRT autoplot renders as expected", {
  local_quiet()
  data <- sim_port_structural(n = 500, seed = 1)
  res <- check_port(data, exposure, c(x3, x1))
  expect_doppelganger(
    "PoRT subgroup prevalence bars",
    ggplot2::autoplot(res)
  )
})

test_that("duplicate subgroups draw a bar of their own prevalence, not a sum", {
  local_quiet()
  res <- check_port(port_anchor_data(), exposure, g)
  # Duplicate every reported row: the plot must still draw each bar to its own
  # prevalence rather than stacking duplicates to twice the length.
  duplicated_res <- res
  duplicated_res@results <- vctrs::vec_rbind(res@results, res@results)

  built <- ggplot2::ggplot_build(ggplot2::autoplot(duplicated_res))
  bars <- built$data[[1]]
  expect_equal(max(bars$xmax), max(res@results$prevalence), tolerance = 1e-8)
})

test_that("bar width is proportional to subgroup size", {
  local_quiet()
  res <- check_port(port_anchor_data(), exposure, g)
  plot_data <- port_plot_data(res)
  expect_equal(
    plot_data$width,
    plot_data$n / max(plot_data$n),
    tolerance = 1e-8
  )
  # The anchor's two subgroups differ in size, so their bars differ in width.
  expect_gt(length(unique(plot_data$width)), 1)
})

test_that("port_plot_data returns the display columns for an empty result", {
  # With no reported subgroups the plotting frame keeps the label, width, and
  # flagged columns but carries no rows.
  plot_data <- port_plot_data(make_port_result())
  expect_identical(nrow(plot_data), 0L)
  expect_true(all(c("label", "width", "flagged") %in% names(plot_data)))
  expect_identical(levels(plot_data$flagged), c("FALSE", "TRUE"))
})

test_that("plot() draws the PoRT view and returns the result invisibly", {
  local_quiet()
  local_null_device()
  data <- sim_port_structural(n = 500, seed = 1)
  res <- check_port(data, exposure, c(x3, x1))

  expect_identical(plot(res), res)
})

test_that("flagged_only draws only the flagged subgroups", {
  local_quiet()
  data <- sim_port_structural(n = 500, seed = 1)
  res <- check_port(data, exposure, c(x3, x1))

  full <- port_plot_data(res)
  flagged_data <- port_plot_data(res, flagged_only = TRUE)
  expect_gt(nrow(flagged_data), 0)
  expect_lt(nrow(flagged_data), nrow(full))
  expect_identical(flagged_data$label, full$label[full$flagged == "TRUE"])
  expect_true(all(flagged_data$flagged == "TRUE"))

  built <- ggplot2::ggplot_build(ggplot2::autoplot(res, flagged_only = TRUE))
  expect_identical(nrow(built$data[[1]]), nrow(flagged_data))
})

test_that("flagged_only recomputes widths relative to the shown subgroups", {
  local_quiet()
  data <- sim_port_structural(n = 500, seed = 1)
  res <- check_port(data, exposure, c(x3, x1))

  flagged_data <- port_plot_data(res, flagged_only = TRUE)
  expect_equal(
    flagged_data$width,
    flagged_data$n / max(flagged_data$n),
    tolerance = 1e-8
  )
  expect_equal(max(flagged_data$width), 1, tolerance = 1e-8)
})

test_that("flagged_only suppresses the fill legend", {
  local_quiet()
  data <- sim_port_structural(n = 500, seed = 1)
  res <- check_port(data, exposure, c(x3, x1))

  flagged_plot <- ggplot2::autoplot(res, flagged_only = TRUE)
  expect_identical(flagged_plot$scales$get_scales("fill")$guide, "none")

  default_plot <- ggplot2::autoplot(res)
  expect_identical(default_plot$scales$get_scales("fill")$guide, "legend")
})

test_that("flagged_only with nothing flagged draws the empty panel", {
  local_quiet()
  data <- sim_port_structural(n = 500, seed = 1)
  res <- check_port(data, exposure, c(x3, x1))
  results <- res@results
  results$flagged <- FALSE
  res@results <- results

  expect_identical(nrow(port_plot_data(res, flagged_only = TRUE)), 0L)

  built <- NULL
  expect_no_error(
    built <- ggplot2::ggplot_build(ggplot2::autoplot(res, flagged_only = TRUE))
  )
  # The reference lines survive the empty bar frame.
  expect_identical(nrow(built$data[[2]]), 2L)
})

test_that("autoplot() rejects a non-flag flagged_only as a classed error", {
  local_quiet()
  data <- sim_port_structural(n = 500, seed = 1)
  res <- check_port(data, exposure, c(x3, x1))

  expect_error(
    ggplot2::autoplot(res, flagged_only = "yes"),
    class = "positively_type_error"
  )
  expect_error(
    ggplot2::autoplot(res, flagged_only = NA),
    class = "positively_type_error"
  )
  expect_error(
    ggplot2::autoplot(res, flagged_only = c(TRUE, FALSE)),
    class = "positively_type_error"
  )
  expect_snapshot_abort(
    ggplot2::autoplot(res, flagged_only = "yes"),
    class = "positively_type_error"
  )
})

test_that("PoRT flagged-only autoplot renders as expected", {
  local_quiet()
  data <- sim_port_structural(n = 500, seed = 1)
  res <- check_port(data, exposure, c(x3, x1))
  expect_doppelganger(
    "PoRT flagged subgroups only",
    ggplot2::autoplot(res, flagged_only = TRUE)
  )
})

test_that("plot() forwards flagged_only", {
  local_quiet()
  local_null_device()
  data <- sim_port_structural(n = 500, seed = 1)
  res <- check_port(data, exposure, c(x3, x1))

  expect_identical(plot(res, flagged_only = TRUE), res)
  expect_error(
    plot(res, flagged_only = "yes"),
    class = "positively_type_error"
  )
})

# ---- rpart controls and passthrough ---------------------------------------

test_that("the fitted trees pin the PoRT rpart controls", {
  local_quiet()
  data <- sim_port_structural(n = 500, seed = 1)
  res <- check_port(data, exposure, c(x3, x1))
  control <- res@trees[[1]]$control
  expect_identical(control$minsplit, 20L)
  expect_identical(control$minbucket, 6L)
  expect_identical(control$maxdepth, 30L)
  expect_identical(control$cp, 0)
})

test_that("the dots reach rpart.control", {
  local_quiet()
  data <- sim_port_structural(n = 500, seed = 1)
  default <- check_port(data, exposure, c(x3, x1))
  pruned <- check_port(data, exposure, c(x3, x1), cp = 0.5)
  expect_identical(pruned@trees[[1]]$control$cp, 0.5)
  # A heavy complexity penalty changes the fitted trees.
  expect_false(identical(default@trees, pruned@trees))
})

test_that("check_port() rejects a dot that is not an rpart.control option", {
  local_quiet()
  data <- sim_port_structural(n = 500, seed = 1)
  expect_snapshot_abort(
    check_port(data, exposure, c(x3, x1), alpa = 0.4),
    class = "positively_args_error"
  )
})

test_that("correlated covariates surface a flagged subgroup for each near-copy", {
  local_quiet()
  withr::local_seed(1)
  n <- 3000
  xa <- stats::runif(n, 0, 100)
  xb <- xa + stats::rnorm(n, sd = 0.01)
  exposure <- stats::rbinom(n, 1L, 0.5)
  exposure[xa > 70] <- 0L
  data <- tibble::tibble(exposure = exposure, xa = xa, xb = xb)

  res <- check_port(data, exposure, c(xa, xb))
  flagged <- flagged_rows(res)
  # Duplicates are surfaced, not de-duplicated: both near-copies flag.
  expect_true(any(grepl("xa", flagged$description)))
  expect_true(any(grepl("xb", flagged$description)))
})

# ---- Description simplification -------------------------------------------

test_that("descriptions collapse chained conditions to the binding bounds", {
  local_quiet()
  # A deep tree on a single continuous covariate stacks many splits on one
  # variable; the description must keep at most one lower and one upper bound.
  withr::local_seed(1)
  n <- 2000
  x3 <- stats::runif(n, 0, 100)
  exposure <- stats::rbinom(n, 1L, 0.5)
  exposure[x3 > 70] <- 0L
  data <- tibble::tibble(exposure = exposure, x3 = x3)

  res <- check_port(data, exposure, x3)
  counts <- vapply(
    res@results$description,
    function(desc) {
      pieces <- strsplit(desc, " & ", fixed = TRUE)[[1]]
      c(
        lower = sum(grepl(">", pieces)),
        upper = sum(grepl("<", pieces))
      )
    },
    numeric(2)
  )
  expect_true(all(counts["lower", ] <= 1))
  expect_true(all(counts["upper", ] <= 1))
})

test_that("numeric beta of 0.5 or more is rejected", {
  local_quiet()
  data <- dgp_good_positivity(n = 200, seed = 1)
  expect_error(
    check_port(data, exposure, c(x1, x2), beta = 0.5),
    class = "positively_error"
  )
  expect_error(
    check_port(data, exposure, c(x1, x2), beta = 0.6),
    class = "positively_error"
  )
})

# ---- Snapshots ------------------------------------------------------------

test_that("the point print method is stable", {
  local_quiet()
  data <- port_anchor_data()
  res <- check_port(data, exposure, g, alpha = 0.05, beta = 0.05)
  expect_snapshot(print(res))
})

test_that("the print method reports an unresolved beta", {
  # When every per-time beta is non-finite, as arises if all sequential time
  # points are skipped, the threshold prints as not resolved.
  res <- make_port_result(beta = NA_real_)
  expect_snapshot(print(res))
})

test_that("check_port() argument validation messages are stable", {
  local_quiet()
  data <- dgp_good_positivity(n = 200, seed = 1)
  single_level <- data
  single_level$exposure <- 1L
  three_level <- data
  three_level$exposure <- rep(c("a", "b", "c"), length.out = nrow(data))

  expect_snapshot_abort(check_port(1:10, exposure, x1))
  expect_snapshot_abort(check_port(
    data,
    exposure,
    tidyselect::starts_with("zzz")
  ))
  expect_snapshot_abort(check_port(data, c(exposure, x1), x2))
  expect_snapshot_abort(check_port(single_level, exposure, c(x1, x2)))
  expect_snapshot_abort(check_port(data, exposure, c(x1, x2), alpha = 0))
  expect_snapshot_abort(check_port(data, exposure, c(x1, x2), alpha = 1.5))
  expect_snapshot_abort(check_port(
    data,
    exposure,
    c(x1, x2),
    alpha = c(0.05, 0.1)
  ))
  expect_snapshot_abort(check_port(data, exposure, c(x1, x2), beta = 0))
  expect_snapshot_abort(check_port(data, exposure, c(x1, x2), beta = 1.5))
  expect_snapshot_abort(check_port(data, exposure, c(x1, x2), beta = 0.6))
  expect_snapshot_abort(check_port(data, exposure, c(x1, x2), beta = "banana"))
  expect_snapshot_abort(check_port(data, exposure, c(x1, x2), gamma = 0))
  expect_snapshot_abort(check_port(data, exposure, c(x1, x2), gamma = 2.5))
  expect_snapshot_abort(check_port(data, exposure, c(x1, x2), n_bins = 0))
  expect_snapshot_abort(check_port(data, exposure, c(x1, x2), breaks = "a"))
  expect_snapshot_abort(check_port(
    three_level,
    exposure,
    c(x1, x2),
    exposure_type = "binary"
  ))
})

# ---- Degenerate breaks and quantile collapse ------------------------------

# The Pair 2 fixtures share one seeded continuous exposure on [0, 6].
port_binning_data <- function() {
  withr::with_seed(
    4,
    tibble::tibble(exposure = stats::runif(300, 0, 6), x1 = stats::rnorm(300))
  )
}

# A zero-inflated exposure: 280 of 300 observations at zero, so the requested
# quantile cut points coincide and three bins collapse to one.
port_zero_inflated_data <- function() {
  withr::with_seed(
    4,
    tibble::tibble(
      exposure = c(rep(0, 280), stats::runif(20, 1, 6)),
      x1 = stats::rnorm(300)
    )
  )
}

test_that("a single break point is rejected rather than read as an interval count", {
  local_quiet()
  # A length-one breaks vector defines no interval; cut() would otherwise read
  # it as a count of equally spaced intervals.
  expect_error(
    check_port(
      port_binning_data(),
      exposure,
      x1,
      breaks = 4,
      exposure_type = "continuous"
    ),
    class = "positively_range_error"
  )
})

test_that("breaks defining a single interval are rejected", {
  local_quiet()
  # Two cut points bound one interval, so every observation shares one exposure
  # level and the reading rule has nothing to contrast.
  expect_error(
    check_port(
      port_binning_data(),
      exposure,
      x1,
      breaks = c(0, 6),
      exposure_type = "continuous"
    ),
    class = "positively_range_error"
  )
})

test_that("missing values in breaks are rejected", {
  local_quiet()
  expect_error(
    check_port(
      port_binning_data(),
      exposure,
      x1,
      breaks = c(0, NA, 6),
      exposure_type = "continuous"
    ),
    class = "positively_missing_error"
  )
})

test_that("a zero-inflated exposure whose quantile bins collapse is rejected", {
  local_quiet()
  expect_error(
    check_port(
      port_zero_inflated_data(),
      exposure,
      x1,
      n_bins = 3,
      exposure_type = "continuous"
    ),
    class = "positively_binning_error"
  )
})

test_that("well-separated explicit breaks yield one exposure level per interval", {
  local_quiet()
  res <- check_port(
    port_binning_data(),
    exposure,
    x1,
    breaks = c(0, 2, 4, 6),
    exposure_type = "continuous"
  )
  expect_identical(res@exposure_type, "continuous")
  expect_length(unique(res@results$exposure_level), 3L)
})

# ---- Missing covariates ---------------------------------------------------

test_that("check_port() aborts on a missing covariate value", {
  local_quiet()
  data <- dgp_good_positivity(n = 300, seed = 1)
  data$x1[1] <- NA
  # rpart would silently drop or surrogate the missing row, so the reported
  # subgroup sizes would no longer sum to the sample.
  expect_error(
    check_port(data, exposure, c(x1, x2)),
    class = "positively_missing_error"
  )
})

# ---- Missing exposure ------------------------------------------------------

# A 0/1 exposure with a large block of missing values. rpart would otherwise
# treat the exposure as a third category and report deflated prevalences, so a
# missing exposure value must abort.
port_missing_exposure_data <- function() {
  withr::with_seed(1, {
    exposure <- stats::rbinom(1000, 1L, 0.5)
    exposure[seq_len(400)] <- NA
    tibble::tibble(exposure = exposure, g = stats::rnorm(1000))
  })
}

test_that("check_port() aborts on a missing exposure value", {
  local_quiet()
  expect_error(
    check_port(port_missing_exposure_data(), exposure, g),
    class = "positively_missing_error"
  )
})

test_that("check_port() missing-exposure message is stable", {
  local_quiet()
  expect_snapshot_abort(check_port(port_missing_exposure_data(), exposure, g))
})

# ---- Response-name collision ----------------------------------------------

# A covariate named .port_response collides with the internal response column
# rpart is handed, degenerating the fit; the covariate must be read like any
# other.
test_that("a covariate named .port_response is handled like any other", {
  local_quiet()
  data <- withr::with_seed(1, {
    v <- stats::runif(800, 0, 2)
    exposure <- stats::rbinom(800, 1L, 0.5)
    exposure[v > 1] <- 0L
    data.frame(exposure = exposure, .port_response = v, check.names = FALSE)
  })
  res <- check_port(data, exposure, tidyselect::all_of(".port_response"))

  flagged <- flagged_rows(res)
  expect_true(any(grepl(".port_response", flagged$description, fixed = TRUE)))
})

# ---- Degenerate binning snapshots -----------------------------------------

test_that("check_port() degenerate binning messages are stable", {
  local_quiet()
  expect_snapshot_abort(check_port(
    port_binning_data(),
    exposure,
    x1,
    breaks = 4,
    exposure_type = "continuous"
  ))
  expect_snapshot_abort(check_port(
    port_binning_data(),
    exposure,
    x1,
    breaks = c(0, 6),
    exposure_type = "continuous"
  ))
  expect_snapshot_abort(check_port(
    port_binning_data(),
    exposure,
    x1,
    breaks = c(0, NA, 6),
    exposure_type = "continuous"
  ))
  expect_snapshot_abort(check_port(
    port_zero_inflated_data(),
    exposure,
    x1,
    n_bins = 3,
    exposure_type = "continuous"
  ))
})

test_that("check_port() missing-covariate message is stable", {
  local_quiet()
  data <- dgp_good_positivity(n = 300, seed = 1)
  data$x1[1] <- NA
  expect_snapshot_abort(check_port(data, exposure, c(x1, x2)))
})
