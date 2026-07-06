# sPoRT applies the PoRT reading rule per time point along a treatment regime.
# Under the stratified strategy, at each time t the risk set is the subjects
# still following the rule (not yet initiated under a monotone rule), the
# exposure is the time-t treatment, and the conditioning set is baseline plus the
# current confounders. Subsetting to followers is what exposes a violation that
# the already-initiated would otherwise mask, so these tests contrast the
# stratified fit against a naive point fit on the full sample.

local_quiet <- function(.env = parent.frame()) {
  withr::local_options(positively.quiet = TRUE, .local_envir = .env)
}

flagged_rows <- function(res) {
  res@results[res@results$flagged, , drop = FALSE]
}

# Monotone binary treatment over three time points with a planted sequential
# violation at t = 2 (REPORT claim 6). Region C (baseline c1 > 1) never
# initiates at t = 2 among subjects still untreated after t = 1. Because the
# treatment is monotone, already-initiated subjects carry a2 == 1 forward and
# inflate the region-C prevalence on the full sample, masking the violation
# unless the risk set is subset to followers.
sim_port_seq <- function(n = 3000, seed = 1) {
  withr::local_seed(seed)
  c1 <- stats::rnorm(n)
  x <- stats::rnorm(n)
  region_c <- c1 > 1

  a1 <- stats::rbinom(n, 1L, 0.3)

  a2 <- a1
  followers2 <- a1 == 0
  p2 <- ifelse(region_c, 0, 0.4)
  a2[followers2] <- stats::rbinom(sum(followers2), 1L, p2[followers2])

  a3 <- a2
  followers3 <- a2 == 0
  a3[followers3] <- stats::rbinom(sum(followers3), 1L, 0.4)

  tibble::tibble(id = seq_len(n), c1 = c1, x = x, a1 = a1, a2 = a2, a3 = a3)
}

# ---- Argument validation ---------------------------------------------------

test_that("check_port_seq() rejects non-data-frame input", {
  local_quiet()
  expect_error(
    check_port_seq(1:10, c(a1, a2), list(c1)),
    class = "positively_error"
  )
})

test_that("check_port_seq() rejects a covariate list of the wrong length", {
  local_quiet()
  data <- sim_port_seq(n = 300, seed = 1)
  # Three exposures but a two-element covariate list is neither matched nor a
  # single recyclable element.
  expect_error(
    check_port_seq(data, c(a1, a2, a3), list(c1, x)),
    class = "positively_error"
  )
})

test_that("check_port_seq() rejects an empty exposure selection", {
  local_quiet()
  data <- sim_port_seq(n = 300, seed = 1)
  expect_error(
    check_port_seq(data, tidyselect::starts_with("zzz"), list(c1)),
    class = "positively_error"
  )
})

test_that("check_port_seq() validates the lag", {
  local_quiet()
  data <- sim_port_seq(n = 300, seed = 1)
  expect_error(
    check_port_seq(data, c(a1, a2, a3), list(c1), lag = -1),
    class = "positively_error"
  )
  expect_error(
    check_port_seq(data, c(a1, a2, a3), list(c1), lag = 1.5),
    class = "positively_error"
  )
  expect_error(
    check_port_seq(data, c(a1, a2, a3), list(c1), lag = "a"),
    class = "positively_error"
  )
})

test_that("check_port_seq() validates alpha, beta, and gamma", {
  local_quiet()
  data <- sim_port_seq(n = 300, seed = 1)
  expect_error(
    check_port_seq(data, c(a1, a2, a3), list(c1), alpha = 1.5),
    class = "positively_error"
  )
  expect_error(
    check_port_seq(data, c(a1, a2, a3), list(c1), beta = 1.5),
    class = "positively_error"
  )
  expect_error(
    check_port_seq(data, c(a1, a2, a3), list(c1), gamma = 0),
    class = "positively_error"
  )
})

test_that("check_port_seq() rejects an unknown strategy", {
  local_quiet()
  data <- sim_port_seq(n = 300, seed = 1)
  expect_error(
    check_port_seq(data, c(a1, a2, a3), list(c1), strategy = "bogus")
  )
})

# ---- Result class and structure -------------------------------------------

test_that("check_port_seq() returns a port_result with a time column", {
  local_quiet()
  data <- sim_port_seq(n = 1000, seed = 1)
  res <- check_port_seq(data, c(a1, a2, a3), list(c1))

  expect_true(S7::S7_inherits(res, positivity_diagnostic))
  expect_identical(S7::S7_class(res)@name, "port_result")
  expect_identical(res@exposure_type, "binary")
  expect_setequal(
    names(res@results),
    c(
      "time",
      "subgroup",
      "description",
      "exposure_level",
      "n",
      "proportion",
      "prevalence",
      "flagged"
    )
  )
})

test_that("check_port_seq() carries one rpart tree per time point", {
  local_quiet()
  data <- sim_port_seq(n = 1000, seed = 1)
  res <- check_port_seq(data, c(a1, a2, a3), list(c1))

  expect_type(res@trees, "list")
  expect_gte(length(res@trees), 3)
  expect_true(all(vapply(
    res@trees,
    function(tree) inherits(tree, "rpart"),
    logical(1)
  )))
})

test_that("a length-one covariate list recycles across time points", {
  local_quiet()
  data <- sim_port_seq(n = 1000, seed = 1)
  recycled <- check_port_seq(data, c(a1, a2, a3), list(c1))
  explicit <- check_port_seq(data, c(a1, a2, a3), list(c1, c1, c1))

  expect_equal(recycled@results, explicit@results)
})

test_that("tidy() and glance() work on the sequential result", {
  local_quiet()
  data <- sim_port_seq(n = 1000, seed = 1)
  res <- check_port_seq(data, c(a1, a2, a3), list(c1))

  expect_identical(generics::tidy(res), res@results)
  expect_identical(nrow(generics::glance(res)), 1L)
})

test_that("check_port_seq() runs on a continuous exposure by categorizing it", {
  local_quiet()
  # dgp_longitudinal() carries continuous time-varying exposures. sPoRT
  # categorizes each into n_bins before applying the reading rule, so the run
  # should return a well-formed per-time result.
  data <- dgp_longitudinal(n = 400, seed = 5)
  res <- check_port_seq(data, c(a1, a2, a3), list(l0, l1, l2), n_bins = 3)

  expect_identical(res@exposure_type, "continuous")
  expect_length(sort(unique(res@results$time)), 3)
})

# ---- sPoRT localization (expectation 10) ----------------------------------

test_that("flags localize to the time point of the planted violation", {
  local_quiet()
  data <- sim_port_seq(n = 3000, seed = 1)
  res <- check_port_seq(data, c(a1, a2, a3), list(c1))

  flagged <- flagged_rows(res)
  expect_true(any(flagged$time == 2))
  expect_false(any(flagged$time %in% c(1, 3)))
  # The planted region is defined by the baseline covariate c1.
  expect_true(any(grepl("c1", flagged$description[flagged$time == 2])))
})

# ---- History subsetting guard (expectation 12) ----------------------------

test_that("stratified subsetting reveals a violation that naive pooling masks", {
  local_quiet()
  data <- sim_port_seq(n = 3000, seed = 1)

  # The stratified strategy subsets to followers (a1 == 0) at t = 2 and flags
  # region C, whose treated prevalence among followers is zero.
  seq_res <- check_port_seq(data, c(a1, a2, a3), list(c1))
  seq_flagged <- flagged_rows(seq_res)
  expect_true(any(seq_flagged$time == 2 & grepl("c1", seq_flagged$description)))

  # A naive point fit on the full sample sees the already-initiated carrying
  # a2 == 1 forward, which inflates the region-C prevalence above beta, so the
  # violation is not flagged.
  naive <- check_port(data, a2, c1)
  expect_false(any(grepl("c1", flagged_rows(naive)$description)))
})

# ---- Risk-set shrinkage and per-time Gruber beta (expectation 11) ---------

test_that("the per-time Gruber beta grows as the monotone risk set shrinks", {
  local_quiet()
  data <- sim_port_seq(n = 3000, seed = 1)
  res <- check_port_seq(data, c(a1, a2, a3), list(c1), beta = "gruber")

  # Risk-set sizes under the monotone rule: everyone at t = 1, the untreated
  # after t = 1 at t = 2, the untreated after t = 2 at t = 3.
  n_t <- c(nrow(data), sum(data$a1 == 0), sum(data$a2 == 0))
  expect_true(all(diff(n_t) < 0))

  expected_beta <- 5 / (sqrt(n_t) * log(n_t))
  # @beta holds the resolved per-time thresholds, ordered by time.
  expect_equal(res@beta, expected_beta, tolerance = 1e-8)
  expect_true(all(diff(res@beta) > 0))
})

# ---- Autoplot contract -----------------------------------------------------

test_that("autoplot() returns a ggplot faceted over time", {
  local_quiet()
  data <- sim_port_seq(n = 1000, seed = 1)
  res <- check_port_seq(data, c(a1, a2, a3), list(c1))
  expect_s3_class(ggplot2::autoplot(res), "ggplot")
})

test_that("sPoRT autoplot renders faceted over time", {
  local_quiet()
  data <- sim_port_seq(n = 1000, seed = 1)
  res <- check_port_seq(data, c(a1, a2, a3), list(c1))
  expect_doppelganger(
    "sPoRT subgroups faceted by time",
    ggplot2::autoplot(res)
  )
})

test_that("each facet shows its own resolved beta reference lines", {
  local_quiet()
  data <- sim_port_seq(n = 3000, seed = 1)
  res <- check_port_seq(data, c(a1, a2, a3), list(c1), beta = "gruber")

  built <- ggplot2::ggplot_build(ggplot2::autoplot(res))
  vlines <- built$data[[2]]
  per_panel <- split(vlines$xintercept, vlines$PANEL)
  for (t in seq_along(res@beta)) {
    panel <- as.character(t)
    if (!is.null(per_panel[[panel]])) {
      expected <- sort(c(res@beta[[t]], 1 - res@beta[[t]]))
      expect_equal(sort(per_panel[[panel]]), expected, tolerance = 1e-8)
    }
  }
  # The upper reference line moves across time, so the panels are not identical.
  expect_gt(length(unique(vlines$xintercept)), 2)
})

# ---- Degenerate risk sets (finding 5) -------------------------------------

# Collect and muffle the risk-set warnings, since more than one wave may skip.
collect_risk_set_warnings <- function(expr) {
  seen <- character(0)
  value <- withCallingHandlers(
    expr,
    positively_risk_set_warning = function(w) {
      seen <<- c(seen, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  list(value = value, warnings = seen)
}

test_that("an all-treated prior wave is skipped with a classed warning", {
  local_quiet()
  data <- sim_port_seq(n = 500, seed = 1)
  # Everyone is treated at time 1, so the time-2 follower set is empty.
  data$a1 <- 1L
  data$a2 <- 1L

  caught <- collect_risk_set_warnings(
    check_port_seq(data, c(a1, a2, a3), list(c1))
  )
  expect_gte(length(caught$warnings), 1)
  expect_false(any(caught$value@results$time == 2))
  expect_true(is.na(caught$value@beta[[2]]))
})

test_that("a single-follower risk set is skipped under the Gruber bound", {
  local_quiet()
  data <- sim_port_seq(n = 500, seed = 1)
  # Leave exactly one untreated subject after time 1.
  data$a1 <- 1L
  data$a1[1] <- 0L

  caught <- collect_risk_set_warnings(
    check_port_seq(data, c(a1, a2, a3), list(c1), beta = "gruber")
  )
  expect_gte(length(caught$warnings), 1)
  expect_true(is.na(caught$value@beta[[2]]))
  expect_false(any(caught$value@results$time == 2))
})

test_that("the untreated level is resolved globally, not per wave", {
  local_quiet()
  data <- sim_port_seq(n = 1000, seed = 1)
  # If every survivor at time 1 were treated, a per-wave min() would flip the
  # follower set. A global untreated level keeps followers as the untreated.
  global <- check_port_seq(data, c(a1, a2, a3), list(c1))
  n_followers_t2 <- sum(data$a1 == 0)
  expect_gt(n_followers_t2, 0)
  expect_true(any(global@results$time == 2))
})

# ---- Snapshots ------------------------------------------------------------

test_that("the sequential print method is stable", {
  local_quiet()
  data <- sim_port_seq(n = 1000, seed = 1)
  res <- check_port_seq(data, c(a1, a2, a3), list(c1), beta = "gruber")
  expect_snapshot(print(res))
})

test_that("check_port_seq() argument validation messages are stable", {
  local_quiet()
  data <- sim_port_seq(n = 300, seed = 1)

  expect_snapshot(check_port_seq(1:10, c(a1, a2), list(c1)), error = TRUE)
  expect_snapshot(
    check_port_seq(data, tidyselect::starts_with("zzz"), list(c1)),
    error = TRUE
  )
  expect_snapshot(
    check_port_seq(data, c(a1, a2, a3), list(c1, x)),
    error = TRUE
  )
  expect_snapshot(
    check_port_seq(data, c(a1, a2, a3), list(c1), lag = -1),
    error = TRUE
  )
  expect_snapshot(
    check_port_seq(data, c(a1, a2, a3), list(c1), lag = 1.5),
    error = TRUE
  )
  expect_snapshot(
    check_port_seq(data, c(a1, a2, a3), list(c1), lag = "a"),
    error = TRUE
  )
  expect_snapshot(
    check_port_seq(data, c(a1, a2, a3), list(c1), alpha = 1.5),
    error = TRUE
  )
  expect_snapshot(
    check_port_seq(data, c(a1, a2, a3), list(c1), beta = 1.5),
    error = TRUE
  )
  expect_snapshot(
    check_port_seq(data, c(a1, a2, a3), list(c1), beta = 0.6),
    error = TRUE
  )
  expect_snapshot(
    check_port_seq(data, c(a1, a2, a3), list(c1), gamma = 0),
    error = TRUE
  )
  expect_snapshot(
    check_port_seq(data, c(a1, a2, a3), list(c1), strategy = "pooled"),
    error = TRUE
  )
  expect_snapshot(
    check_port_seq(data, c(a1, a2, a3), list(c1), .censoring = c1),
    error = TRUE
  )
})
