# sPoRT applies the PoRT reading rule per time point along a treatment regime.
# Under the stratified strategy, at each time t the risk set is the subjects
# still following the rule (not yet initiated under a monotone rule), the
# exposure is the time-t treatment, and the conditioning set is baseline plus the
# current confounders. Subsetting to followers is what exposes a violation that
# the already-initiated would otherwise mask, so these tests contrast the
# stratified fit against a naive point fit on the full sample.

low_support_rows <- function(res) {
  res@results[res@results$low_support, , drop = FALSE]
}

# Monotone binary treatment over three time points with a planted sequential
# violation at t = 2. Region C (baseline c1 > 1) never initiates at t = 2
# among subjects still untreated after t = 1. Because the treatment is
# monotone, already-initiated subjects carry a2 == 1 forward and inflate the
# region-C prevalence on the full sample, masking the violation unless the
# risk set is subset to followers.
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

test_that("check_port_seq() rejects a dot that is not an rpart.control option", {
  local_quiet()
  data <- sim_port_seq(n = 300, seed = 1)
  expect_error(
    check_port_seq(data, c(a1, a2, a3), list(c1), alpa = 0.4),
    class = "positively_args_error"
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

test_that("check_port_seq() rejects an empty conditioning set", {
  local_quiet()
  data <- sim_port_seq(n = 300, seed = 1)
  expect_error(
    check_port_seq(data, c(a1, a2), list(tidyselect::starts_with("zzz"))),
    class = "positively_empty_error"
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
    check_port_seq(data, c(a1, a2, a3), list(c1), strategy = "bogus"),
    class = "positively_args_error"
  )
})

test_that("the pooled strategy is rejected ahead of an empty conditioning set", {
  local_quiet()
  data <- sim_port_seq(n = 300, seed = 1)
  # An empty covariate selection leaves the first time point with nothing to
  # condition on, but the pooled strategy is the more fundamental obstruction:
  # it is not implemented at all, so it must be reported first rather than the
  # downstream empty-conditioning error.
  expect_error(
    check_port_seq(
      data,
      c(a1, a2),
      list(tidyselect::starts_with("zzz")),
      strategy = "pooled"
    ),
    class = "positively_strategy_error"
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
      "low_support"
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

# ---- Declared exposure type ------------------------------------------------

# Every reported row carries proportion = n / n_t against its own time point's
# risk set, so the risk-set size is recoverable exactly from the tidy output
# whatever shape the trees took.
risk_set_sizes <- function(res) {
  sizes <- res@results$n / res@results$proportion
  as.integer(round(tapply(sizes, res@results$time, max)))
}

test_that("the detected type flips with the wave count on one dose column", {
  local_quiet()
  data <- dgp_longitudinal_dose()

  # One wave of doses clears the unique-value cutoff, so detection reads the
  # column as continuous and the reading rule runs on quantile bins.
  one_wave <- check_port_seq(data, a1, list(l0), n_bins = 3, cp = 0.01)
  expect_identical(one_wave@exposure_type, "continuous")
  expect_length(unique(one_wave@results$exposure_level), 3L)

  # Adding a wave changes nothing about any dose column, yet detection reads the
  # exposure pooled across waves, so the same doses now fall below the cutoff.
  # The type follows the wave count, and with it the binning path is abandoned:
  # one response per distinct dose replaces one per bin, and n_bins goes unread.
  two_waves <- check_port_seq(
    data,
    c(a1, a2),
    list(l0, l1),
    n_bins = 3,
    cp = 0.01
  )
  expect_identical(two_waves@exposure_type, "categorical")
  expect_length(unique(two_waves@results$exposure_level), 30L)
})

test_that("the announced type names the argument check_port_seq() takes", {
  withr::local_options(positively.quiet = FALSE)
  data <- sim_port_seq(n = 3000, seed = 1)

  # Detection runs once on the pooled exposure, so the announcement belongs to
  # `.exposures`. Naming `.exposure` would point a user at an argument this
  # function does not have, which is why the announced name is threaded from the
  # caller rather than fixed in the detector.
  expect_message(
    check_port_seq(data, c(a1, a2), list(c1)),
    "Treating `.exposures` as binary",
    fixed = TRUE
  )
})

test_that("a declared continuous type bins a dose the pooled reading misreads", {
  local_quiet()
  data <- dgp_longitudinal_dose()
  res <- check_port_seq(
    data,
    c(a1, a2, a3),
    list(l0, l1, l2),
    n_bins = 3,
    exposure_type = "continuous",
    cp = 0.01
  )

  expect_identical(res@exposure_type, "continuous")

  # One exposure level per bin at every wave is what proves n_bins is read. A
  # run that merely succeeded would say nothing: the categorical branch also
  # succeeds, and silently reports one level per distinct dose instead.
  per_time <- tapply(
    res@results$exposure_level,
    res@results$time,
    function(level) length(unique(level))
  )
  expect_identical(as.integer(per_time), rep(3L, 3))

  # The waves share a marginal dose distribution, so their quantile cut points
  # and bin labels coincide, and the labels are intervals rather than doses.
  levels_seen <- unique(res@results$exposure_level)
  expect_length(levels_seen, 3L)
  expect_true(all(startsWith(levels_seen, "[") | startsWith(levels_seen, "(")))

  # n_bins now moves the result, which is the whole point: under the pooled
  # detection it was inert.
  five <- check_port_seq(
    data,
    c(a1, a2, a3),
    list(l0, l1, l2),
    n_bins = 5,
    exposure_type = "continuous",
    cp = 0.01
  )
  expect_length(unique(five@results$exposure_level), 5L)
})

test_that("explicit breaks bin a declared continuous sequential dose", {
  local_quiet()
  data <- dgp_longitudinal_dose()
  res <- check_port_seq(
    data,
    c(a1, a2, a3),
    list(l0, l1, l2),
    breaks = c(0, 25, 50, 75),
    exposure_type = "continuous",
    cp = 0.01
  )

  # The cut points are the user's, not the wave's quantiles, so the labels pin
  # that breaks reached the binning rather than being dropped with n_bins.
  expect_setequal(
    unique(res@results$exposure_level),
    c("[0,25]", "(25,50]", "(50,75]")
  )
})

test_that("a declared binary type rejects a three-level pooled exposure", {
  local_quiet()
  data <- sim_port_seq(n = 300, seed = 1)
  # A third level at time 2 leaves the pooled exposure with three values, so no
  # binary regime exists for the follower restriction to follow.
  data$a2[seq_len(40)] <- 2L
  expect_length(unique(c(data$a1, data$a2)), 3L)

  err <- expect_error(
    check_port_seq(data, c(a1, a2), list(c1), exposure_type = "binary"),
    class = "positively_exposure_type_error"
  )
  expect_match(
    conditionMessage(err),
    "`.exposures` has 3 distinct values",
    fixed = TRUE
  )

  # Detection reads the same pooled column as categorical and runs, so the abort
  # comes from the declaration alone.
  detected <- check_port_seq(data, c(a1, a2), list(c1), cp = 0.01)
  expect_identical(detected@exposure_type, "categorical")
})

test_that("a declared continuous type rejects a factor exposure", {
  local_quiet()
  data <- sim_port_seq(n = 300, seed = 1)
  data$a1 <- factor(data$a1, labels = c("no", "yes"))
  data$a2 <- factor(data$a2, labels = c("no", "yes"))

  err <- expect_error(
    check_port_seq(data, c(a1, a2), list(c1), exposure_type = "continuous"),
    class = "positively_exposure_type_error"
  )
  # A factor pools to a factor, and its level codes are not doses, so the
  # declaration fails on the column rather than binning the codes.
  expect_match(conditionMessage(err), "continuous", fixed = TRUE)
  expect_match(conditionMessage(err), "`.exposures` is <factor>", fixed = TRUE)
})

test_that("check_port_seq() offers every exposure type it can compute", {
  local_quiet()
  data <- sim_port_seq(n = 300, seed = 1)

  err <- expect_error(
    check_port_seq(data, c(a1, a2), list(c1), exposure_type = "ordinal"),
    class = "positively_args_error"
  )
  # The sequential reader computes all three types, so the menu is the full one.
  # cli wraps it across lines at the console width, so the comparison collapses
  # whitespace.
  expect_match(
    gsub("\\s+", " ", conditionMessage(err)),
    paste(
      '`exposure_type` must be one of "auto", "binary", "categorical", or',
      '"continuous", not "ordinal".'
    ),
    fixed = TRUE
  )
})

test_that("a declared non-binary type drops the monotone follower restriction", {
  local_quiet()
  data <- dgp_longitudinal_binary()
  followers <- as.integer(c(
    nrow(data),
    sum(data$a1 == 0),
    sum(data$a2 == 0)
  ))

  detected <- check_port_seq(data, c(a1, a2, a3), list(l0, l1, l2), cp = 0.01)
  expect_identical(detected@exposure_type, "binary")
  expect_identical(risk_set_sizes(detected), followers)

  # Declaring the same 0/1 column categorical takes the non-binary branch of the
  # risk set, which applies no follower restriction, so every subject stays in
  # every wave.
  declared <- check_port_seq(
    data,
    c(a1, a2, a3),
    list(l0, l1, l2),
    exposure_type = "categorical",
    cp = 0.01
  )
  expect_identical(declared@exposure_type, "categorical")
  expect_identical(risk_set_sizes(declared), rep(nrow(data), 3L))

  # The restriction is what exposes the planted time-2 violation, so dropping it
  # loses the flag: the already-initiated carry their treatment forward and the
  # trees flag that carry-forward instead.
  detected_low_support <- low_support_rows(detected)
  declared_low_support <- low_support_rows(declared)
  expect_true(any(
    detected_low_support$time == 2 &
      grepl("l1", detected_low_support$description)
  ))
  expect_false(any(
    declared_low_support$time == 2 &
      grepl("l1", declared_low_support$description)
  ))
  expect_true(any(
    declared_low_support$time == 2 &
      grepl("a1", declared_low_support$description)
  ))
})

test_that("a declared binary type keeps the follower restriction", {
  local_quiet()
  data <- dgp_longitudinal_binary()

  detected <- check_port_seq(data, c(a1, a2, a3), list(l0, l1, l2), cp = 0.01)
  declared <- check_port_seq(
    data,
    c(a1, a2, a3),
    list(l0, l1, l2),
    exposure_type = "binary",
    cp = 0.01
  )

  # A pooled exposure with exactly two values is the one case detection cannot
  # get wrong, so declaring binary unlocks no run detection would refuse. What
  # it does carry is the follower restriction, which every non-binary type
  # drops.
  expect_identical(declared@exposure_type, "binary")
  expect_equal(declared@results, detected@results)
  expect_identical(
    risk_set_sizes(declared),
    as.integer(c(nrow(data), sum(data$a1 == 0), sum(data$a2 == 0)))
  )
})

test_that("a degenerate pooled exposure is rejected before the declared type", {
  local_quiet()
  data <- sim_port_seq(n = 300, seed = 1)
  data$a1 <- 0L
  data$a2 <- 0L

  # A single distinct value cannot carry a binary type either, but the
  # degenerate exposure is the more specific complaint and is reported first.
  expect_error(
    check_port_seq(data, c(a1, a2), list(c1), exposure_type = "binary"),
    class = "positively_exposure_error"
  )
})

# ---- sPoRT localization ---------------------------------------------------

test_that("flags localize to the time point of the planted violation", {
  local_quiet()
  data <- sim_port_seq(n = 3000, seed = 1)
  res <- check_port_seq(data, c(a1, a2, a3), list(c1))

  low_support <- low_support_rows(res)
  expect_true(any(low_support$time == 2))
  expect_false(any(low_support$time %in% c(1, 3)))
  # The planted region is defined by the baseline covariate c1.
  expect_true(any(grepl("c1", low_support$description[low_support$time == 2])))
})

# ---- History subsetting guard ---------------------------------------------

test_that("stratified subsetting reveals a violation that naive pooling masks", {
  local_quiet()
  data <- sim_port_seq(n = 3000, seed = 1)

  # The stratified strategy subsets to followers (a1 == 0) at t = 2 and flags
  # region C, whose treated prevalence among followers is zero.
  seq_res <- check_port_seq(data, c(a1, a2, a3), list(c1))
  seq_low_support <- low_support_rows(seq_res)
  expect_true(any(
    seq_low_support$time == 2 & grepl("c1", seq_low_support$description)
  ))

  # A naive point fit on the full sample sees the already-initiated carrying
  # a2 == 1 forward, which inflates the region-C prevalence above beta, so the
  # violation does not read as low support.
  naive <- check_port(data, a2, c1)
  expect_false(any(grepl("c1", low_support_rows(naive)$description)))
})

# ---- Risk-set shrinkage and per-time Gruber beta --------------------------

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

test_that("low_support_only keeps a facet for a time point it empties", {
  local_quiet()
  data <- sim_port_seq(n = 1000, seed = 1)
  res <- check_port_seq(data, c(a1, a2, a3), list(c1))

  # Pin the shape under test: time 2 carries no low-support rows, time 1 carries
  # at least one.
  results <- res@results
  results$low_support[results$time == 2] <- FALSE
  results$low_support[which(results$time == 1)[1]] <- TRUE
  res@results <- results

  built <- ggplot2::ggplot_build(ggplot2::autoplot(
    res,
    low_support_only = TRUE
  ))
  layout <- built$layout$layout
  expect_setequal(layout$time, 1:3)

  bars <- built$data[[1]]
  bar_times <- layout$time[match(bars$PANEL, layout$PANEL)]
  expect_true(any(bar_times == 1))
  expect_false(any(bar_times == 2))

  # The time-2 facet still shows its reference lines.
  vlines <- built$data[[2]]
  vline_times <- layout$time[match(vlines$PANEL, layout$PANEL)]
  expect_true(any(vline_times == 2))
})

# ---- Censoring autoplot -----------------------------------------------------

test_that("the censoring run stores resolved per-time censoring thresholds", {
  local_quiet()
  data <- dgp_longitudinal_binary_censoring()
  res <- check_port_seq(
    data,
    c(a1, a2, a3),
    list(l0, l1, l2),
    .censoring = c(c1, c2, c3),
    beta = "gruber",
    cp = 0.01
  )

  # The censoring thresholds resolve from the censoring risk-set sizes, the
  # subjects uncensored through the previous time point, which differ from the
  # exposure risk-set sizes.
  risk_sets <- c(
    nrow(data),
    sum(data$c1 == 0),
    sum(data$c1 == 0 & data$c2 == 0)
  )
  expected <- 5 / (sqrt(risk_sets) * log(risk_sets))
  expect_equal(res@censoring_beta, expected, tolerance = 1e-8)
  expect_false(isTRUE(all.equal(res@censoring_beta, res@beta)))
})

test_that("autoplot on a censoring run facets by type", {
  local_quiet()
  data <- dgp_longitudinal_binary_censoring()
  res <- check_port_seq(
    data,
    c(a1, a2, a3),
    list(l0, l1, l2),
    .censoring = c(c1, c2, c3),
    cp = 0.01
  )

  built <- ggplot2::ggplot_build(ggplot2::autoplot(res))
  # The facet layout carries the type column, so exposure and censoring
  # subgroups occupy distinct panels.
  expect_true("type" %in% names(built$layout$layout))
  expect_setequal(unique(built$layout$layout$type), c("exposure", "censoring"))
})

test_that("censoring panels draw their own resolved reference lines", {
  local_quiet()
  data <- dgp_longitudinal_binary_censoring()
  res <- check_port_seq(
    data,
    c(a1, a2, a3),
    list(l0, l1, l2),
    .censoring = c(c1, c2, c3),
    beta = "gruber",
    cp = 0.01
  )

  built <- ggplot2::ggplot_build(ggplot2::autoplot(res))
  layout <- built$layout$layout
  vlines <- built$data[[2]]

  for (i in seq_len(nrow(layout))) {
    panel <- layout$PANEL[[i]]
    time <- layout$time[[i]]
    beta_t <- if (layout$type[[i]] == "censoring") {
      res@censoring_beta[[time]]
    } else {
      res@beta[[time]]
    }
    got <- sort(vlines$xintercept[vlines$PANEL == panel])
    if (length(got) > 0 && is.finite(beta_t)) {
      expect_equal(got, sort(c(beta_t, 1 - beta_t)), tolerance = 1e-8)
    }
  }
})

test_that("sPoRT censoring autoplot renders by type", {
  local_quiet()
  data <- dgp_longitudinal_binary_censoring()
  res <- check_port_seq(
    data,
    c(a1, a2, a3),
    list(l0, l1, l2),
    .censoring = c(c1, c2, c3),
    cp = 0.01
  )
  expect_doppelganger(
    "sPoRT censoring subgroups by type",
    ggplot2::autoplot(res)
  )
})

test_that("the censoring run carries a family-neutral prevalence axis label", {
  local_quiet()
  data <- dgp_longitudinal_binary_censoring()

  censoring <- check_port_seq(
    data,
    c(a1, a2, a3),
    list(l0, l1, l2),
    .censoring = c(c1, c2, c3),
    cp = 0.01
  )
  # The censoring run pools exposure and censoring subgroups into one plot, so
  # the axis label must not claim exposure prevalence for every bar.
  expect_identical(
    ggplot2::autoplot(censoring)$labels$x,
    "Prevalence in subgroup"
  )

  exposure_only <- check_port_seq(
    data,
    c(a1, a2, a3),
    list(l0, l1, l2),
    cp = 0.01
  )
  expect_identical(
    ggplot2::autoplot(exposure_only)$labels$x,
    "Exposure prevalence in subgroup"
  )
})

test_that("low_support_only keeps both type facets on a censoring run", {
  local_quiet()
  data <- dgp_longitudinal_binary_censoring()
  res <- check_port_seq(
    data,
    c(a1, a2, a3),
    list(l0, l1, l2),
    .censoring = c(c1, c2, c3),
    cp = 0.01
  )

  built <- ggplot2::ggplot_build(ggplot2::autoplot(
    res,
    low_support_only = TRUE
  ))
  layout <- built$layout$layout
  expect_setequal(unique(layout$type), c("exposure", "censoring"))

  # Every drawn bar is one of the deduplicated low-support rows.
  low_support_data <- port_plot_data(res, low_support_only = TRUE)
  expect_identical(nrow(built$data[[1]]), nrow(low_support_data))
  expect_true(all(low_support_data$low_support == "TRUE"))
})

# ---- Degenerate risk sets -------------------------------------------------

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

test_that("autoplot builds when every time point is skipped", {
  local_quiet()
  # A tiny risk set puts every wave's resolved Gruber bound at or above 0.5, so
  # each wave is skipped and @results is the empty sequential tibble. Faceting an
  # empty frame aborts the build, so the empty result must render unfaceted.
  make_data <- function() {
    withr::local_seed(1)
    n <- 12
    d <- tibble::tibble(id = seq_len(n), c1 = rnorm(n), a1 = rbinom(n, 1L, 0.3))
    followers2 <- d$a1 == 0
    d$a2 <- d$a1
    d$a2[followers2] <- rbinom(sum(followers2), 1L, 0.4)
    followers3 <- d$a2 == 0
    d$a3 <- d$a2
    d$a3[followers3] <- rbinom(sum(followers3), 1L, 0.4)
    d
  }
  data <- make_data()

  caught <- collect_risk_set_warnings(
    check_port_seq(data, c(a1, a2, a3), list(c1), beta = "gruber")
  )
  res <- caught$value
  expect_identical(nrow(res@results), 0L)
  expect_true(all(is.na(res@beta)))

  expect_no_error(ggplot2::ggplot_build(ggplot2::autoplot(res)))
  built <- ggplot2::ggplot_build(ggplot2::autoplot(res))
  expect_identical(nrow(built$layout$layout), 1L)
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

# ---- Binary longitudinal fixture ------------------------------------------

test_that("the stratified risk set shrinks over time in the tidy output", {
  local_quiet()
  data <- dgp_longitudinal_binary()
  # Grow shallow trees so each quiet wave reports its risk set as one node,
  # rather than fragmenting into noise leaves under the default cp of zero.
  res <- check_port_seq(data, c(a1, a2, a3), list(l0, l1, l2), cp = 0.01)

  # Under the monotone rule the follower set is everyone at t = 1, the untreated
  # after t = 1 at t = 2, and the untreated after t = 2 at t = 3.
  risk_sets <- c(nrow(data), sum(data$a1 == 0), sum(data$a2 == 0))
  expect_true(all(diff(risk_sets) < 0))

  # The largest reported subgroup at each time is the whole risk set, so the
  # per-time size in the tidy output tracks the shrinking follower set.
  per_time_n <- as.integer(tapply(res@results$n, res@results$time, max))
  expect_identical(per_time_n, as.integer(risk_sets))
  expect_true(all(diff(per_time_n) < 0))
})

test_that("the planted time-2 subgroup reads as low support", {
  local_quiet()
  data <- dgp_longitudinal_binary()
  res <- check_port_seq(data, c(a1, a2, a3), list(l0, l1, l2))

  low_support <- low_support_rows(res)
  # Among the t = 2 followers, treatment is deterministic in the l1 region, so
  # the reading rule surfaces it.
  expect_true(any(low_support$time == 2 & grepl("l1", low_support$description)))
  # The quiet first and last waves carry no flag.
  expect_false(any(low_support$time %in% c(1, 3)))
})

test_that("the lag argument changes the conditioning set", {
  local_quiet()
  data <- dgp_longitudinal_binary()

  # With the full history the t = 2 conditioning set holds the lagged covariate
  # l0, and its deterministic region reads as low support.
  full <- check_port_seq(data, c(a1, a2, a3), list(l0, l1, l2))
  full_t2 <- full@results$description[full@results$time == 2]
  expect_true(any(grepl("l0", full_t2)))
  expect_true(any(grepl("l0", low_support_rows(full)$description)))

  # With lag zero the t = 2 conditioning set drops l0, so it appears in no
  # subgroup, and only the l1 region still reads as low support.
  lagged <- check_port_seq(data, c(a1, a2, a3), list(l0, l1, l2), lag = 0)
  lagged_t2 <- lagged@results$description[lagged@results$time == 2]
  expect_false(any(grepl("l0", lagged_t2)))
  lagged_low_support <- low_support_rows(lagged)
  expect_false(any(grepl("l0", lagged_low_support$description)))
  expect_true(any(
    lagged_low_support$time == 2 & grepl("l1", lagged_low_support$description)
  ))
})

# ---- Censoring -------------------------------------------------------------

test_that("check_port_seq() rejects a censoring selection of the wrong length", {
  local_quiet()
  data <- dgp_longitudinal_binary_censoring(n = 400, seed = 7)
  expect_error(
    check_port_seq(
      data,
      c(a1, a2, a3),
      list(l0, l1, l2),
      .censoring = c(c1, c2)
    ),
    class = "positively_error"
  )
})

test_that("check_port_seq() rejects non-binary censoring indicators", {
  local_quiet()
  data <- dgp_longitudinal_binary_censoring(n = 400, seed = 7)
  # The l0, l1, l2 columns are continuous, so they are not valid 0/1 indicators.
  expect_error(
    check_port_seq(
      data,
      c(a1, a2, a3),
      list(l0, l1, l2),
      .censoring = c(l0, l1, l2)
    ),
    class = "positively_error"
  )
})

test_that("a censoring result distinguishes censoring rows from exposure rows", {
  local_quiet()
  data <- dgp_longitudinal_binary_censoring()
  res <- check_port_seq(
    data,
    c(a1, a2, a3),
    list(l0, l1, l2),
    .censoring = c(c1, c2, c3),
    cp = 0.01
  )
  expect_true("type" %in% names(res@results))
  expect_setequal(unique(res@results$type), c("exposure", "censoring"))
})

test_that("the planted censoring subgroup has low support at its time", {
  local_quiet()
  data <- dgp_longitudinal_binary_censoring()
  res <- check_port_seq(
    data,
    c(a1, a2, a3),
    list(l0, l1, l2),
    .censoring = c(c1, c2, c3)
  )

  censoring_low_support <- res@results[
    res@results$type == "censoring" & res@results$low_support,
    ,
    drop = FALSE
  ]
  # Censoring is near-certain where l0 > 1 among the uncensored-so-far at t = 2.
  expect_true(any(
    censoring_low_support$time == 2 &
      grepl("l0", censoring_low_support$description)
  ))
})

test_that("the censoring risk set shrinks as prior censoring accrues", {
  local_quiet()
  data <- dgp_longitudinal_binary_censoring()
  res <- check_port_seq(
    data,
    c(a1, a2, a3),
    list(l0, l1, l2),
    .censoring = c(c1, c2, c3),
    cp = 0.01
  )

  # The censoring tree at time t runs on the subjects uncensored through t - 1.
  risk_sets <- c(
    nrow(data),
    sum(data$c1 == 0),
    sum(data$c1 == 0 & data$c2 == 0)
  )
  expect_true(all(diff(risk_sets) < 0))

  censoring <- res@results[res@results$type == "censoring", , drop = FALSE]
  per_time_n <- as.integer(tapply(censoring$n, censoring$time, max))
  expect_identical(per_time_n, as.integer(risk_sets))
  expect_true(all(diff(per_time_n) < 0))
})

test_that("supplying .censoring shrinks the exposure risk sets", {
  local_quiet()
  data <- dgp_longitudinal_binary_censoring()

  without <- check_port_seq(data, c(a1, a2, a3), list(l0, l1, l2), cp = 0.01)
  with_censoring <- check_port_seq(
    data,
    c(a1, a2, a3),
    list(l0, l1, l2),
    .censoring = c(c1, c2, c3),
    cp = 0.01
  )

  # A run without censoring carries no type column and reports exposure rows.
  expect_false("type" %in% names(without@results))

  # The largest reported subgroup at each time is the whole risk set, so the
  # per-time size tracks the exposure risk set.
  max_n <- function(rows) as.integer(tapply(rows$n, rows$time, max))
  exposure_rows <- with_censoring@results[
    with_censoring@results$type == "exposure",
    ,
    drop = FALSE
  ]

  # Without censoring the exposure risk set is the regime followers alone.
  followers <- c(nrow(data), sum(data$a1 == 0), sum(data$a2 == 0))
  expect_identical(max_n(without@results), followers)

  # With censoring the exposure risk set is the followers who are also
  # uncensored through the previous time point, so it is strictly smaller from
  # time 2 on and unchanged at time 1, which has no prior censoring.
  followers_uncensored <- c(
    nrow(data),
    sum(data$a1 == 0 & data$c1 == 0),
    sum(data$a2 == 0 & data$c1 == 0 & data$c2 == 0)
  )
  expect_identical(max_n(exposure_rows), followers_uncensored)
  expect_identical(followers_uncensored[[1]], followers[[1]])
  expect_true(all(followers_uncensored[-1] < followers[-1]))
})

test_that("pooled exposure type detection skips NA exposures", {
  local_quiet()
  data <- dgp_longitudinal_binary_censoring(n = 1000, seed = 7)
  # Censored subjects carry no observed exposure at later time points, so the
  # pooled exposure column holds 0, 1, and NA.
  data$a2[data$c1 == 1] <- NA
  data$a3[data$c1 == 1 | data$c2 == 1] <- NA

  res <- check_port_seq(
    data,
    c(a1, a2, a3),
    list(l0, l1, l2),
    .censoring = c(c1, c2, c3),
    cp = 0.01
  )

  # Skipping the NA values keeps the pooled exposure binary rather than
  # resolving it to categorical on the three distinct values 0, 1, and NA.
  expect_identical(res@exposure_type, "binary")
})

# ---- Factor-coded binary exposures ----------------------------------------

# check_port() accepts a factor exposure by resolving its levels with sort(), so
# the sequential reader must resolve the untreated (reference) level the same
# way rather than with min(), which is undefined for unordered factors.

test_that("check_port_seq() reads factor-coded binary exposures like integers", {
  local_quiet()
  data <- sim_port_seq(n = 3000, seed = 1)
  integer_res <- check_port_seq(data, c(a1, a2, a3), list(c1))

  factored <- data
  factored$a1 <- factor(factored$a1, labels = c("no", "yes"))
  factored$a2 <- factor(factored$a2, labels = c("no", "yes"))
  factored$a3 <- factor(factored$a3, labels = c("no", "yes"))
  factor_res <- check_port_seq(factored, c(a1, a2, a3), list(c1))

  expect_identical(factor_res@exposure_type, "binary")

  # The planted time-2 violation and its localization survive the factor coding.
  factor_low_support <- low_support_rows(factor_res)
  expect_true(any(factor_low_support$time == 2))
  expect_false(any(factor_low_support$time %in% c(1, 3)))

  # The readings fall on the same time points and subgroups as the integer run.
  low_support_keys <- function(res) {
    low_support <- low_support_rows(res)
    sort(paste(low_support$time, low_support$description))
  }
  expect_identical(
    low_support_keys(factor_res),
    low_support_keys(integer_res)
  )

  # The per-time follower risk set, read as the largest reported subgroup at
  # each time, is identical to the integer-coded run.
  max_n <- function(res) {
    as.integer(tapply(res@results$n, res@results$time, max))
  }
  expect_identical(max_n(factor_res), max_n(integer_res))
})

test_that("the untreated level is the first factor level, not the sorted label", {
  local_quiet()
  data <- sim_port_seq(n = 3000, seed = 1)

  # Reverse the alphabetical order of the labels relative to the level order:
  # value 0 maps to "untreated" (the reference level) and value 1 to "treated",
  # yet "treated" sorts before "untreated". A reader that resolved the untreated
  # level by sorted label order would flip the follower set to the treated.
  labelled <- data
  labelled$a1 <- factor(labelled$a1, labels = c("untreated", "treated"))
  labelled$a2 <- factor(labelled$a2, labels = c("untreated", "treated"))
  labelled$a3 <- factor(labelled$a3, labels = c("untreated", "treated"))

  res <- check_port_seq(labelled, c(a1, a2, a3), list(c1))

  # The time-2 follower risk set is the subjects untreated after time 1, read as
  # the largest reported subgroup at that time. Resolving the untreated level to
  # the alphabetically first label would instead size it from the treated.
  followers_t2 <- sum(data$a1 == 0)
  treated_t1 <- sum(data$a1 == 1)
  expect_false(followers_t2 == treated_t1)

  time2_max_n <- max(res@results$n[res@results$time == 2])
  expect_identical(as.integer(time2_max_n), as.integer(followers_t2))
})

# ---- Snapshots ------------------------------------------------------------

test_that("the sequential print method is stable", {
  local_quiet()
  data <- sim_port_seq(n = 1000, seed = 1)
  res <- check_port_seq(data, c(a1, a2, a3), list(c1), beta = "gruber")
  expect_snapshot(print(res))
})

test_that("the sequential print is stable on the binary longitudinal fixture", {
  local_quiet()
  data <- dgp_longitudinal_binary()
  res <- check_port_seq(data, c(a1, a2, a3), list(l0, l1, l2), cp = 0.01)
  expect_snapshot(print(res))
})

test_that("the sequential print reports censoring subgroups distinctly", {
  local_quiet()
  data <- dgp_longitudinal_binary_censoring()
  res <- check_port_seq(
    data,
    c(a1, a2, a3),
    list(l0, l1, l2),
    .censoring = c(c1, c2, c3),
    cp = 0.01
  )
  expect_snapshot(print(res))
})

test_that("check_port_seq() censoring validation messages are stable", {
  local_quiet()
  data <- dgp_longitudinal_binary_censoring(n = 400, seed = 7)

  expect_snapshot_abort(check_port_seq(
    data,
    c(a1, a2, a3),
    list(l0, l1, l2),
    .censoring = c(c1, c2)
  ))
  expect_snapshot_abort(check_port_seq(
    data,
    c(a1, a2, a3),
    list(l0, l1, l2),
    .censoring = c(l0, l1, l2)
  ))
  expect_snapshot_abort(check_port_seq(
    data,
    c(a1, a2, a3),
    list(l0, l1, l2),
    .censoring = c(c1, c2, c3),
    strategy = "pooled"
  ))
})

test_that("check_port_seq() argument validation messages are stable", {
  local_quiet()
  data <- sim_port_seq(n = 300, seed = 1)

  expect_snapshot_abort(check_port_seq(1:10, c(a1, a2), list(c1)))
  expect_snapshot_abort(check_port_seq(
    data,
    tidyselect::starts_with("zzz"),
    list(c1)
  ))
  expect_snapshot_abort(check_port_seq(data, c(a1, a2, a3), list(c1, x)))
  expect_snapshot_abort(check_port_seq(
    data,
    c(a1, a2),
    list(tidyselect::starts_with("zzz"))
  ))
  expect_snapshot_abort(check_port_seq(data, c(a1, a2, a3), list(c1), lag = -1))
  expect_snapshot_abort(check_port_seq(
    data,
    c(a1, a2, a3),
    list(c1),
    lag = 1.5
  ))
  expect_snapshot_abort(check_port_seq(
    data,
    c(a1, a2, a3),
    list(c1),
    lag = "a"
  ))
  expect_snapshot_abort(check_port_seq(
    data,
    c(a1, a2, a3),
    list(c1),
    alpha = 1.5
  ))
  expect_snapshot_abort(check_port_seq(
    data,
    c(a1, a2, a3),
    list(c1),
    beta = 1.5
  ))
  expect_snapshot_abort(check_port_seq(
    data,
    c(a1, a2, a3),
    list(c1),
    beta = 0.6
  ))
  expect_snapshot_abort(check_port_seq(
    data,
    c(a1, a2, a3),
    list(c1),
    gamma = 0
  ))
  expect_snapshot_abort(check_port_seq(
    data,
    c(a1, a2, a3),
    list(c1),
    strategy = "pooled"
  ))
})

# ---- Missing covariates ----------------------------------------------------

test_that("check_port_seq() aborts on a missing covariate in a risk set", {
  local_quiet()
  data <- sim_port_seq(n = 300, seed = 1)
  data$c1[1] <- NA
  # The missing value sits in the time-1 risk set, where every subject is a
  # follower, so rpart would silently drop or surrogate it and the reported
  # subgroup sizes would no longer sum to the risk set.
  expect_error(
    check_port_seq(data, c(a1, a2, a3), list(c1)),
    class = "positively_missing_error"
  )
})

test_that("a censored subject may carry NA in later covariates without aborting", {
  local_quiet()
  data <- dgp_longitudinal_binary_censoring()
  # A subject censored at time 1 has no observed time-2 or time-3 covariates, so
  # its later covariate values are missing. Because it leaves every later risk
  # set, those missing values never enter a tree and the run proceeds.
  data$l1[data$c1 == 1] <- NA
  data$l2[data$c1 == 1 | data$c2 == 1] <- NA

  res <- check_port_seq(
    data,
    c(a1, a2, a3),
    list(l0, l1, l2),
    .censoring = c(c1, c2, c3),
    cp = 0.01
  )
  expect_true(S7::S7_inherits(res, positivity_diagnostic))
  expect_true("type" %in% names(res@results))
})

test_that("check_port_seq() missing-covariate message names the column and time", {
  local_quiet()
  data <- sim_port_seq(n = 300, seed = 1)
  data$c1[1] <- NA
  expect_snapshot_abort(check_port_seq(data, c(a1, a2, a3), list(c1)))
})

# ---- Missing time-t exposure and censoring status -------------------------

# A missing time-t exposure marks a subject with no observed treatment at that
# time, distinct from a missing conditioning covariate. It must leave the
# exposure risk set for that time point rather than enter the tree as an NA
# response, so the risk-set size and its resolved Gruber threshold count only the
# subjects with an observed exposure.
test_that("a missing time-t exposure leaves that time's exposure risk set", {
  local_quiet()
  data <- withr::with_seed(2, {
    n <- 800
    c1 <- stats::rnorm(n)
    a1 <- stats::rbinom(n, 1L, 0.3)
    a2 <- a1
    f <- a1 == 0
    a2[f] <- stats::rbinom(sum(f), 1L, 0.4)
    drop <- f & stats::rbinom(n, 1L, 0.5) == 1
    a2[drop] <- NA
    tibble::tibble(id = seq_len(n), c1 = c1, a1 = a1, a2 = a2)
  })
  # The conditioning covariate c1 is complete, so the missingness is confined to
  # the time-2 exposure column and never trips the risk-set completeness check.
  expect_false(anyNA(data$c1))
  k <- sum(data$a1 == 0 & !is.na(data$a2))

  res <- check_port_seq(
    data,
    c(a1, a2),
    list(c1),
    beta = "gruber",
    cp = 0.01
  )

  time2 <- res@results[res@results$time == 2, , drop = FALSE]
  # The largest reported subgroup at time 2 is the whole risk set, so its size is
  # the count of followers with an observed exposure and its proportion is one.
  expect_identical(as.integer(max(time2$n)), as.integer(k))
  expect_equal(max(time2$proportion), 1)
  # The Gruber bound resolves from the observed-exposure risk-set size.
  expect_equal(res@beta[[2]], 5 / (sqrt(k) * log(k)), tolerance = 1e-8)
})

test_that("a missing time-t censoring status leaves that censoring risk set", {
  local_quiet()
  data <- dgp_longitudinal_binary_censoring()
  # Set the time-2 censoring status to unknown for a slice of the subjects
  # uncensored through time 1, the time-2 censoring risk set.
  in_risk <- which(data$c1 == 0)
  data$c2[in_risk[seq_len(300)]] <- NA
  k2 <- sum(data$c1 == 0 & !is.na(data$c2))
  t3_risk <- sum(data$c1 == 0 & data$c2 == 0, na.rm = TRUE)

  res <- check_port_seq(
    data,
    c(a1, a2, a3),
    list(l0, l1, l2),
    .censoring = c(c1, c2, c3),
    beta = "gruber",
    cp = 0.01
  )

  # The time-2 censoring threshold resolves from the observed-status count, not
  # from the full uncensored-through-time-1 cohort.
  expect_equal(
    res@censoring_beta[[2]],
    5 / (sqrt(k2) * log(k2)),
    tolerance = 1e-8
  )

  # A subject with an unknown time-2 status is absent from the time-3 censoring
  # risk set, so the largest reported time-3 censoring subgroup is the
  # observed-uncensored cohort.
  censoring <- res@results[res@results$type == "censoring", , drop = FALSE]
  time3 <- censoring[censoring$time == 3, , drop = FALSE]
  expect_identical(as.integer(max(time3$n)), as.integer(t3_risk))
})

# ---- Degenerate pooled exposure -------------------------------------------

test_that("an identically constant pooled exposure is rejected", {
  local_quiet()
  data <- sim_port_seq(n = 300, seed = 1)
  data$a1 <- 0L
  data$a2 <- 0L
  # A single distinct value across every time point leaves the reading rule with
  # no exposed level to contrast against the unexposed.
  expect_error(
    check_port_seq(data, c(a1, a2), list(c1)),
    class = "positively_error"
  )
})

test_that("an all-missing pooled exposure is rejected", {
  local_quiet()
  data <- sim_port_seq(n = 300, seed = 1)
  data$a1 <- NA_integer_
  data$a2 <- NA_integer_
  expect_error(
    check_port_seq(data, c(a1, a2), list(c1)),
    class = "positively_error"
  )
})

test_that("check_port_seq() degenerate pooled-exposure messages are stable", {
  local_quiet()
  data <- sim_port_seq(n = 300, seed = 1)
  constant <- data
  constant$a1 <- 0L
  constant$a2 <- 0L
  all_missing <- data
  all_missing$a1 <- NA_integer_
  all_missing$a2 <- NA_integer_

  expect_snapshot_abort(check_port_seq(constant, c(a1, a2), list(c1)))
  expect_snapshot_abort(check_port_seq(all_missing, c(a1, a2), list(c1)))
})
