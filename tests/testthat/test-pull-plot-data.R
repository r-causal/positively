# `pull_plot_data()` hands back exactly the frame one autoplot view draws, so
# nearly every block here is a drift test: the accessor's frame against what
# ggplot2 reports the corresponding layers drawing. A frame that stopped
# matching its view would still look reasonable on its own, which is what these
# comparisons are for.

# ---- Scenario generators --------------------------------------------------

# Gaussian exposure with an independent Gaussian covariate. Both intervention
# values sit inside the support, so every EDP is positive and finite.
sim_pull_gaussian <- function(n, seed = 1) {
  withr::local_seed(seed)
  tibble::tibble(exposure = stats::rnorm(n), x1 = stats::rnorm(n))
}

# Three strata under exact categorical matching: two are treated at different
# rates and get distinct finite ideal weights, and the third is never treated,
# so its outcome support at a* = 1 is zero and its ideal weight is infinite. The
# scatter view draws the finite and the infinite rows in separate layers.
sim_pull_mixed_support <- function() {
  data.frame(
    exposure = c(
      rep(c(0L, 1L), c(4L, 4L)),
      rep(c(0L, 1L), c(4L, 2L)),
      rep(0L, 5L)
    ),
    s = factor(rep(c("a", "b", "c"), c(8L, 6L, 5L)))
  )
}

# ---- View menus, read off the methods -------------------------------------

# The menus are read from the methods rather than restated, so a block cannot
# pass against a view the plot no longer offers. The S7 methods on `autoplot()`
# register as S3 methods on ggplot2's generic; `pull_plot_data()` is an S7
# generic and its methods are reached through S7.
autoplot_method <- function(object) {
  utils::getS3method("autoplot", class(object)[[1]])
}

autoplot_types <- function(object) {
  eval(formals(autoplot_method(object))$type)
}

pull_plot_data_types <- function(class) {
  eval(formals(S7::method(pull_plot_data, class))$type)
}

deparsed_call <- function(condition) {
  paste(deparse(conditionCall(condition)), collapse = " ")
}

# ---- Reading a built plot -------------------------------------------------

# Built layers are read by what draws them rather than by the order they were
# added in: a histogram or a bar carries a count, a reference line an intercept,
# points a shape, a box its outliers, a ribbon a ymax.
layers_drawing <- function(built, aesthetic) {
  Filter(function(layer) aesthetic %in% names(layer), built$data)
}

layer_drawing <- function(built, aesthetic) {
  layers_drawing(built, aesthetic)[[1]]
}

# Every strip label a faceted plot renders, with its labeller applied. The
# labels are read off the built layout, one facet dimension at a time, so what
# a block compares against is the text the figure draws rather than the lookup
# behind it.
rendered_strip_labels <- function(plot) {
  built <- ggplot2::ggplot_build(plot)
  layout <- built$layout$layout
  # ggplot2 spells its own layout bookkeeping columns in upper case, so dropping
  # those leaves the variables the faceting was built on, whatever they are
  # named.
  keys <- layout[!grepl("^[A-Z_]+$", names(layout))]
  labeller <- match.fun(built$layout$facet_params$labeller)
  unlist(lapply(names(keys), function(key) {
    as.character(unlist(labeller(unique(keys[key]))))
  }))
}

# ---- The generic ----------------------------------------------------------

test_that("pull_plot_data() is an exported generic with an edp_result method", {
  expect_true("pull_plot_data" %in% getNamespaceExports("positively"))
  expect_s3_class(pull_plot_data, "S7_generic")
  expect_identical(names(formals(pull_plot_data))[[1]], "x")

  data <- sim_pull_gaussian(60)
  res <- check_edp(
    data,
    exposure,
    x1,
    values = c(0, 1),
    exposure_type = "continuous"
  )

  expect_true(is.function(S7::method(pull_plot_data, edp_result)))
  expect_s3_class(pull_plot_data(res), "tbl_df")
})

test_that("pull_plot_data() defaults to the frame the default view draws", {
  data <- sim_pull_gaussian(60)
  res <- check_edp(
    data,
    exposure,
    x1,
    values = c(0, 1),
    exposure_type = "continuous"
  )

  # Whichever view heads the menu is the default of both, so the two defaults
  # are compared rather than named.
  default_type <- autoplot_types(res)[[1]]
  expect_identical(
    pull_plot_data(res),
    pull_plot_data(res, type = default_type)
  )

  # The default view draws a single frame at the plot level, so what the plot
  # holds is what it draws.
  expect_identical(pull_plot_data(res), ggplot2::autoplot(res)$data)
})

test_that("the type menu mirrors autoplot's and refuses an unknown view", {
  data <- sim_pull_gaussian(60)
  res <- check_edp(
    data,
    exposure,
    x1,
    values = c(0, 1),
    exposure_type = "continuous"
  )

  # The accessor exists to hand back what a view draws, so it offers exactly the
  # views autoplot() offers, in the same order.
  expect_identical(
    autoplot_types(res),
    c("boxplot", "histogram", "ecdf", "density", "scatter")
  )
  expect_identical(pull_plot_data_types(edp_result), autoplot_types(res))

  # A view name is chosen from a fixed menu like any other argument, and asking
  # for the frame rather than the figure does not make the refusal a lesser one.
  # The reader is shown the call they made rather than the method dispatch
  # reached.
  menu_error <- rlang::catch_cnd(pull_plot_data(res, type = "bogus"))
  expect_s3_class(menu_error, "positively_args_error")
  expect_match(deparsed_call(menu_error), "pull_plot_data", fixed = TRUE)
})

# ---- Frames against what the views draw -----------------------------------

test_that("the boxplot frame is the hand-computed five-number summary", {
  # Five observations whose exposure and covariate both run 0 to 4, with both
  # half-distances at 1. Every kernel factor is then
  # 0.5 ^ ((a_j - a*) ^ 2 + (x_j - x_i) ^ 2), a negative power of two, so each
  # EDP is a sum of powers of two and is exact in double precision.
  data <- tibble::tibble(exposure = c(0, 1, 2, 3, 4), x1 = c(0, 1, 2, 3, 4))
  res <- check_edp(
    data,
    exposure,
    x1,
    bw_exposure = 1,
    bw_covariates = 1,
    values = 0,
    exposure_type = "continuous"
  )
  edp <- c(
    2^0 + 2^-2 + 2^-8 + 2^-18 + 2^-32,
    2^0 + 2^-5 + 2^-13 + 2^-25,
    2^-2 + 2^-3 + 2^-10 + 2^-20,
    2^-4 + 2^-8 + 2^-17,
    2^-8 + 2^-9 + 2^-15
  )
  expect_identical(res@results$edp, edp)

  frame <- pull_plot_data(res, type = "boxplot")

  expect_identical(
    names(frame),
    c(
      "intervention",
      "n",
      "ymin",
      "lower",
      "middle",
      "upper",
      "ymax",
      "outliers"
    )
  )
  expect_identical(nrow(frame), 1L)
  expect_s3_class(frame$intervention, "factor")
  expect_identical(levels(frame$intervention), "0")
  expect_type(frame$n, "integer")
  expect_identical(frame$n, 5L)
  for (stat in c("ymin", "lower", "middle", "upper", "ymax")) {
    expect_type(frame[[stat]], "double")
  }
  expect_type(frame$outliers, "list")
  expect_type(frame$outliers[[1]], "double")

  # Type 7 reads a quantile at h = (n - 1)p + 1, which for these five
  # observations is 1.2, 2, 3, 4, and 4.8: the quartiles and the median are
  # order statistics outright, and each whisker interpolates a fifth of the way
  # in from the end. The interpolated ends are compared rather than pinned bit
  # for bit, because a fifth is not a binary fraction.
  sorted <- sort(edp)
  expect_equal(
    frame$ymin,
    sorted[[1]] + 0.2 * (sorted[[2]] - sorted[[1]]),
    tolerance = 1e-12
  )
  expect_identical(frame$lower, sorted[[2]])
  expect_identical(frame$middle, sorted[[3]])
  expect_identical(frame$upper, sorted[[4]])
  expect_equal(
    frame$ymax,
    sorted[[4]] + 0.8 * (sorted[[5]] - sorted[[4]]),
    tolerance = 1e-12
  )

  # The smallest and the largest observation are the two lying strictly outside
  # the whiskers.
  expect_identical(sort(frame$outliers[[1]]), sorted[c(1L, 5L)])
})

test_that("the boxplot frame is one row per intervention in grid order", {
  data <- sim_pull_gaussian(60)
  res <- check_edp(
    data,
    exposure,
    x1,
    values = list(down = function(.x) .x - 1, 0),
    exposure_type = "continuous"
  )

  frame <- pull_plot_data(res, type = "boxplot")

  # The levels follow the order the interventions were given rather than the
  # order their labels sort in, and so do the rows.
  expect_identical(
    levels(frame$intervention),
    unique(names(res@params$values))
  )
  expect_identical(levels(frame$intervention), c("down", "0"))
  expect_identical(as.character(frame$intervention), c("down", "0"))
  expect_identical(frame$n, c(60L, 60L))

  plot <- ggplot2::autoplot(res, type = "boxplot")
  expect_identical(plot$data, frame)

  built <- ggplot2::ggplot_build(plot)
  box <- built$data[[1]]
  # One box per row of the frame, at the heights the frame gives rather than at
  # the ones a boxplot would compute for itself.
  expect_identical(nrow(box), nrow(frame))
  expect_identical(as.numeric(box$x), c(1, 2))
  for (stat in c("ymin", "lower", "middle", "upper", "ymax")) {
    expect_equal(box[[stat]], frame[[stat]])
  }

  # The point layer draws the outliers list-column unchopped: every value, with
  # its own box, and nothing besides.
  points <- built$data[[2]]
  expect_gt(nrow(points), 0L)
  expect_identical(nrow(points), length(unlist(frame$outliers)))
  expect_equal(
    unname(lapply(split(points$y, as.numeric(points$x)), sort)),
    unname(lapply(frame$outliers, sort))
  )
})

test_that("the estimator boxplot frame covers the two EDP measures", {
  data <- sim_pull_gaussian(60)
  res <- check_edp(
    data,
    exposure,
    x1,
    variant = "estimator",
    values = list(down = function(.x) .x - 1, 0),
    exposure_type = "continuous"
  )

  frame <- pull_plot_data(res, type = "boxplot")

  expect_identical(
    names(frame),
    c(
      "intervention",
      "measure",
      "n",
      "ymin",
      "lower",
      "middle",
      "upper",
      "ymax",
      "outliers"
    )
  )
  expect_s3_class(frame$measure, "factor")
  expect_s3_class(frame$intervention, "factor")
  # The ideal weight shares no scale with the two EDP measures, so this view
  # leaves it to the scatter view and to tidy().
  expect_identical(levels(frame$measure), c("edp_outcome", "edp_treatment"))
  expect_identical(levels(frame$intervention), c("down", "0"))

  # One row per measure and intervention, and every pairing is present.
  expect_identical(nrow(frame), 4L)
  expect_setequal(
    paste(frame$intervention, frame$measure),
    c(
      "down edp_outcome",
      "0 edp_outcome",
      "down edp_treatment",
      "0 edp_treatment"
    )
  )
  expect_identical(frame$n, rep(60L, 4L))
  expect_type(frame$outliers, "list")

  # Each row summarizes its own measure within its own intervention, which the
  # median pins independently of the quantile call.
  for (row in seq_len(nrow(frame))) {
    measure <- as.character(frame$measure[[row]])
    label <- as.character(frame$intervention[[row]])
    values <- res@results[[measure]][res@results$intervention == label]
    expect_identical(length(values), frame$n[[row]])
    expect_equal(frame$middle[[row]], stats::median(values))
  }
})

test_that("the density frame is the observation frame and needs no ggridges", {
  data <- sim_pull_gaussian(60)
  res <- check_edp(
    data,
    exposure,
    x1,
    values = c(0, 1),
    exposure_type = "continuous"
  )

  frame <- pull_plot_data(res, type = "density")
  expect_identical(frame, pull_plot_data(res, type = "histogram"))

  # The frame is the numbers rather than the figure, so the package that draws
  # the ridgelines is not needed to read them.
  local_mocked_bindings(is_installed = function(...) FALSE, .package = "rlang")
  expect_identical(pull_plot_data(res, type = "density"), frame)
})

test_that("the histogram frame is what the histogram view draws", {
  data <- sim_pull_gaussian(60)
  res <- check_edp(
    data,
    exposure,
    x1,
    values = c(0, 1),
    exposure_type = "continuous"
  )

  frame <- pull_plot_data(res, type = "histogram")

  # The observation-grain results, keyed on the intervention label in grid
  # order. A function intervention varies `value` from row to row, so the label
  # is the only key that separates the interventions. Pinning the whole frame
  # rather than a column at a time leaves no room for a column to be dropped,
  # reordered, or quietly recomputed.
  expected <- res@results
  expected$intervention <- factor(
    expected$intervention,
    levels = names(res@params$values)
  )
  expect_identical(frame, expected)

  plot <- ggplot2::autoplot(res, type = "histogram")
  expect_identical(plot$data, frame)

  built <- ggplot2::ggplot_build(plot)
  layout <- built$layout$layout
  layer <- built$data[[1]]
  # One panel per intervention, in level order, and every row of the frame lands
  # in the bins of its own panel.
  expect_identical(
    as.character(layout$intervention),
    levels(frame$intervention)
  )
  panel_counts <- as.numeric(tapply(layer$count, layer$PANEL, sum))
  expect_identical(panel_counts, as.numeric(table(frame$intervention)))
})

test_that("the ecdf frame is what the ecdf view draws", {
  data <- sim_pull_gaussian(60)
  res <- check_edp(
    data,
    exposure,
    x1,
    values = c(0, 1),
    exposure_type = "continuous"
  )

  frame <- pull_plot_data(res, type = "ecdf")

  expect_s3_class(frame$intervention, "factor")
  expect_identical(levels(frame$intervention), names(res@params$values))
  expect_identical(nrow(frame), nrow(res@results))
  # Anchored to the results rather than only to the curve drawn below, so a
  # frame that agreed with its own view but not with the diagnostic would
  # still be caught.
  expect_identical(frame$edp, res@results$edp)

  plot <- ggplot2::autoplot(res, type = "ecdf")
  expect_identical(plot$data, frame)

  built <- ggplot2::ggplot_build(plot)
  layer <- built$data[[1]]
  # One curve per intervention, each stepping at exactly that intervention's EDP
  # values. `stat_ecdf()` pads a curve's ends with infinities.
  expect_length(unique(layer$group), nlevels(frame$intervention))
  drawn <- lapply(split(layer$x, layer$group), function(x) {
    sort(x[is.finite(x)])
  })
  expected <- lapply(
    split(frame$edp, frame$intervention),
    function(edp) sort(unique(edp))
  )
  expect_equal(unname(drawn), unname(expected))
})

test_that("the scatter frame is the full estimator results and covers both layers", {
  local_quiet()
  data <- sim_pull_mixed_support()
  res <- check_edp(
    data,
    exposure,
    s,
    variant = "estimator",
    values = 1,
    categorical_similarity = 0
  )

  frame <- pull_plot_data(res, type = "scatter")

  # The view splits on `ideal_weight` itself, so the whole results tibble goes
  # back and the split the plot draws stays derivable from it.
  expect_identical(frame, res@results)
  n_finite <- sum(is.finite(frame$ideal_weight))
  n_infinite <- sum(is.infinite(frame$ideal_weight))
  expect_gt(n_finite, 0L)
  expect_gt(n_infinite, 0L)

  built <- ggplot2::ggplot_build(ggplot2::autoplot(res, type = "scatter"))
  layer_rows <- vapply(built$data, nrow, integer(1))
  expect_length(built$data, 2)
  expect_setequal(layer_rows, c(n_finite, n_infinite))
  expect_identical(sum(layer_rows), nrow(frame))

  # The two layers together draw the frame's rows and nothing else.
  drawn <- do.call(
    rbind,
    lapply(built$data, function(layer) layer[c("x", "y")])
  )
  expect_equal(sort(drawn$x), sort(frame$edp_treatment))
  expect_equal(sort(drawn$y), sort(frame$edp_outcome))
})

test_that("a function intervention keys the frame by its label", {
  data <- sim_pull_gaussian(60)
  res <- check_edp(
    data,
    exposure,
    x1,
    values = list(down = function(.x) .x - 1, 0),
    exposure_type = "continuous"
  )

  frame <- pull_plot_data(res, type = "histogram")

  expect_identical(levels(frame$intervention), c("down", "0"))
  expect_identical(as.integer(table(frame$intervention)), rep(nrow(data), 2L))
  # The shift gives nearly every observation its own intervened exposure, so
  # `value` cannot separate the interventions and the label must.
  shifted <- frame$value[frame$intervention == "down"]
  expect_gt(length(unique(shifted)), 1L)
  expect_identical(unique(frame$value[frame$intervention == "0"]), 0)
})

test_that("interventions sharing a label share a level", {
  data <- sim_pull_gaussian(60)
  # The vector path labels each intervention by its own value, so a repeated
  # value gives two interventions one label.
  res <- check_edp(
    data,
    exposure,
    x1,
    values = c(1, 1),
    exposure_type = "continuous"
  )

  frame <- pull_plot_data(res, type = "histogram")

  # The shared label is one level, and neither intervention's rows are dropped
  # or turned into missing levels along the way.
  expect_identical(levels(frame$intervention), "1")
  expect_identical(nrow(frame), 2L * nrow(data))
  expect_false(anyNA(frame$intervention))

  for (type in c("histogram", "ecdf", "density")) {
    expect_identical(pull_plot_data(res, type = type), frame)
  }

  # Both observation views draw the merged rows together, one panel and one
  # curve rather than two of either.
  histogram <- ggplot2::ggplot_build(ggplot2::autoplot(res, type = "histogram"))
  expect_identical(nrow(histogram$layout$layout), 1L)
  expect_identical(sum(histogram$data[[1]]$count), as.numeric(nrow(frame)))

  ecdf <- ggplot2::ggplot_build(ggplot2::autoplot(res, type = "ecdf"))
  expect_length(unique(ecdf$data[[1]]$group), 1L)
})

# ---- Gate parity ----------------------------------------------------------

test_that("pull_plot_data() and autoplot() refuse the scatter view alike", {
  local_quiet()
  data <- sim_pull_gaussian(60)
  plain <- check_edp(
    data,
    exposure,
    x1,
    values = 0,
    exposure_type = "continuous"
  )

  pull_error <- rlang::catch_cnd(pull_plot_data(plain, type = "scatter"))
  plot_error <- rlang::catch_cnd(ggplot2::autoplot(plain, type = "scatter"))

  # One gate behind two entry points: the reader is told the same thing either
  # way and is shown the call they made.
  expect_identical(conditionMessage(pull_error), conditionMessage(plot_error))
  expect_match(deparsed_call(pull_error), "pull_plot_data", fixed = TRUE)
  expect_match(deparsed_call(plot_error), "autoplot", fixed = TRUE)

  expect_snapshot_abort(
    pull_plot_data(plain, type = "scatter"),
    class = "positively_variant_error"
  )
  expect_snapshot_abort(
    ggplot2::autoplot(plain, type = "scatter"),
    class = "positively_variant_error"
  )
})

# ---- Density ratios -------------------------------------------------------

# The exact Gaussian shift-MTP density ratio, log r = k Z - k ^ 2 / 2. A vector
# of these is a point treatment; a matrix is a time-varying treatment, one
# column per time point.
sim_pull_ratios <- function(n, k = 0.5, seed = 1) {
  withr::local_seed(seed)
  exp(k * stats::rnorm(n) - k^2 / 2)
}

test_that("the density-ratio menu mirrors autoplot's", {
  res <- check_density_ratios(sim_pull_ratios(40))

  expect_identical(autoplot_types(res), c("distribution", "cumulative"))
  expect_identical(
    pull_plot_data_types(density_ratios_result),
    autoplot_types(res)
  )
  expect_identical(
    pull_plot_data(res),
    pull_plot_data(res, type = "distribution")
  )
})

test_that("the distribution frame unchops the ratios by time", {
  m <- matrix(sim_pull_ratios(60), nrow = 20, ncol = 3)
  res <- check_density_ratios(m)

  frame <- pull_plot_data(res, type = "distribution")

  # One row per ratio per time point, the grain the boxplot view summarizes and
  # the histogram view bins. `@ratios` is the only place the raw ratios live, so
  # nothing else can hand a reader the numbers behind either figure.
  expect_identical(names(frame), c("time", "ratio"))
  expect_type(frame$time, "integer")
  expect_type(frame$ratio, "double")
  expect_identical(nrow(frame), 60L)
  expect_identical(frame$time, rep(1:3, each = 20L))
  expect_identical(frame$ratio, unlist(res@ratios))

  plot <- ggplot2::autoplot(res, type = "distribution")
  expect_identical(plot$data, frame)

  built <- ggplot2::ggplot_build(plot)
  box <- layer_drawing(built, "outliers")
  # One box per time point, at the positions the times take, each spanning that
  # time's ratios with its outliers included.
  expect_identical(nrow(box), 3L)
  expect_identical(as.numeric(box$x), c(1, 2, 3))
  by_time <- split(frame$ratio, frame$time)
  expect_equal(box$ymin_final, unname(vapply(by_time, min, numeric(1))))
  expect_equal(box$ymax_final, unname(vapply(by_time, max, numeric(1))))
  expect_equal(box$middle, unname(vapply(by_time, stats::median, numeric(1))))
})

test_that("the distribution frame is what the point histogram draws", {
  res <- check_density_ratios(sim_pull_ratios(50, k = 1))

  frame <- pull_plot_data(res, type = "distribution")

  # A point treatment is one time point, and the frame keeps the time column
  # either way, so one view has one shape whichever figure it draws.
  expect_identical(names(frame), c("time", "ratio"))
  expect_identical(frame$time, rep(1L, 50L))
  expect_identical(frame$ratio, res@ratios[[1]])

  plot <- ggplot2::autoplot(res, type = "distribution")
  expect_identical(plot$data, frame)

  built <- ggplot2::ggplot_build(plot)
  bins <- layer_drawing(built, "count")
  # Every ratio is binned once and the bins span the frame, so the histogram
  # accounts for the frame's rows and for nothing besides.
  expect_identical(sum(bins$count), as.numeric(nrow(frame)))
  expect_gte(min(frame$ratio), min(bins$xmin))
  expect_lte(max(frame$ratio), max(bins$xmax))
})

test_that("the cumulative frame is the series the view draws", {
  m <- matrix(sim_pull_ratios(60), nrow = 20, ncol = 3)
  res <- check_density_ratios(m)

  frame <- pull_plot_data(res, type = "cumulative")

  expect_identical(names(frame), c("time", "statistic", "value", "series"))
  expect_type(frame$time, "integer")
  expect_type(frame$statistic, "character")
  expect_type(frame$value, "double")
  expect_s3_class(frame$series, "factor")
  # The series read as text and are ordered by the quantile behind them, with
  # the maximum last rather than where its label happens to sort.
  expect_identical(
    levels(frame$series),
    c(
      "50th percentile",
      "90th percentile",
      "95th percentile",
      "99th percentile",
      "Maximum"
    )
  )
  # Only the cumulative rows, and among those only the quantiles and the
  # maximum: a mean or a proportion shares no scale with them.
  expect_true(all(startsWith(frame$statistic, "cumulative_")))
  expect_identical(nrow(frame), 15L)

  plot <- ggplot2::autoplot(res, type = "cumulative")
  expect_identical(plot$data, frame)

  built <- ggplot2::ggplot_build(plot)
  points <- layer_drawing(built, "shape")
  # One point per row and one group per series, each at the time and the value
  # the frame gives it.
  expect_identical(nrow(points), nrow(frame))
  expect_length(unique(points$group), nlevels(frame$series))
  expect_equal(sort(points$x), sort(as.numeric(frame$time)))
  expect_equal(
    unname(lapply(split(points$y, points$group), sort)),
    unname(lapply(split(frame$value, frame$series), sort))
  )
})

test_that("pull_plot_data() and autoplot() refuse the cumulative view alike", {
  res <- check_density_ratios(sim_pull_ratios(40, k = 1))

  pull_error <- rlang::catch_cnd(pull_plot_data(res, type = "cumulative"))
  plot_error <- rlang::catch_cnd(ggplot2::autoplot(res, type = "cumulative"))

  # One gate behind two entry points: the reader is told the same thing either
  # way and is shown the call they made.
  expect_identical(conditionMessage(pull_error), conditionMessage(plot_error))
  expect_match(deparsed_call(pull_error), "pull_plot_data", fixed = TRUE)
  expect_match(deparsed_call(plot_error), "autoplot", fixed = TRUE)

  expect_snapshot_abort(
    pull_plot_data(res, type = "cumulative"),
    class = "positively_type_error"
  )
  expect_snapshot_abort(
    ggplot2::autoplot(res, type = "cumulative"),
    class = "positively_type_error"
  )
})

# ---- ETA.Bias -------------------------------------------------------------

# A steep propensity model pushes fitted scores toward zero and one, the
# practical violation ETA.Bias is aimed at. The outcome is linear in the
# exposure, so the estimand holds one term at every truncation level.
sim_pull_eta <- function(n = 300, seed = 1) {
  withr::local_seed(seed)
  x1 <- stats::rnorm(n)
  x2 <- stats::rnorm(n)
  a <- stats::rbinom(n, 1L, stats::plogis(2 * (x1 + x2)))
  tibble::tibble(a = a, y = a + x1 + x2 + stats::rnorm(n), x1 = x1, x2 = x2)
}

# The same mechanism with a continuous exposure, which caps a stabilized weight
# rather than bounding a probability and so reports no lower bound at any level.
sim_pull_eta_continuous <- function(n = 300, seed = 1) {
  withr::local_seed(seed)
  x1 <- stats::rnorm(n)
  x2 <- stats::rnorm(n)
  a <- x1 + x2 + stats::rnorm(n)
  tibble::tibble(a = a, y = a + x1 + x2 + stats::rnorm(n), x1 = x1, x2 = x2)
}

# The bootstrap draw is seeded so the estimates a frame carries are the ones the
# figure drawn beside it holds.
fit_pull_eta <- function(data = sim_pull_eta(), ..., seed = 2024) {
  withr::local_seed(seed)
  check_eta_bias(data, a, y, c(x1, x2), ...)
}

test_that("the eta-bias menu mirrors autoplot's", {
  local_quiet()
  res <- fit_pull_eta(n_boot = 10)

  expect_identical(autoplot_types(res), c("bootstrap", "sweep"))
  expect_identical(pull_plot_data_types(eta_bias_result), autoplot_types(res))
  expect_identical(pull_plot_data(res), pull_plot_data(res, type = "bootstrap"))
})

test_that("the bootstrap frame unchops the bootstrap estimates", {
  local_quiet()
  res <- fit_pull_eta(truncation_grid = c(0, 0.05, 0.1), n_boot = 20)

  frame <- pull_plot_data(res, type = "bootstrap")

  expect_identical(
    names(frame),
    c("term", "level", "label", "estimate", "truth")
  )
  expect_s3_class(frame$term, "factor")
  expect_s3_class(frame$level, "factor")
  expect_s3_class(frame$label, "factor")
  expect_type(frame$estimate, "double")
  expect_type(frame$truth, "double")

  # One row per draw per reported reading. `@boot_estimates` is the only place
  # the draws live, so nothing else hands a reader the distributions the facets
  # show.
  expect_identical(nrow(frame), sum(lengths(res@boot_estimates)))
  expect_identical(frame$estimate, unlist(res@boot_estimates))
  expect_identical(levels(frame$term), unique(res@results$term))
  # Both keys are factors so the panels follow the sweep and the order the terms
  # were built in. Left as text the level index would sort lexically and put a
  # tenth level ahead of a second one.
  expect_identical(levels(frame$level), c("1", "2", "3"))
  expect_identical(
    levels(frame$label),
    c("lower = 0", "lower = 0.05", "lower = 0.1")
  )
  # The truth is a property of the term alone, so every row of a term carries
  # the one value the estimator in it was aimed at.
  expect_identical(frame$truth, unname(res@truth[as.character(frame$term)]))

  plot <- ggplot2::autoplot(res, type = "bootstrap")
  expect_identical(plot$data, frame)

  built <- ggplot2::ggplot_build(plot)

  # The strips read the labels the frame carries, so a reader rebuilding the
  # figure keys the facets on a column rather than on a lookup they cannot see.
  # The term strip rides along in the same list, so the labels are counted
  # rather than compared to the whole of it.
  strips <- rendered_strip_labels(plot)
  expect_true(all(levels(frame$label) %in% strips))
  expect_identical(
    sum(strips %in% levels(frame$label)),
    nlevels(frame$label)
  )

  bins <- layer_drawing(built, "count")
  # Every draw is binned once, inside the panel of its own level.
  expect_identical(sum(bins$count), as.numeric(nrow(frame)))
  expect_identical(
    unname(vapply(split(bins$count, bins$PANEL), sum, numeric(1))),
    as.numeric(table(frame$level))
  )

  # One truth per term, so the lines draw the frame's distinct truths and no
  # target a panel's estimator was never aimed at.
  lines <- layer_drawing(built, "xintercept")
  expect_setequal(lines$xintercept, unique(frame$truth))
})

test_that("the sweep frame is the band the sweep view draws", {
  local_quiet()
  res <- fit_pull_eta(truncation_grid = c(0, 0.05, 0.1), n_boot = 20)

  frame <- pull_plot_data(res, type = "sweep")

  expect_identical(
    names(frame),
    c("term", "truncation", "bias", "lower", "upper")
  )
  expect_s3_class(frame$term, "factor")
  expect_identical(levels(frame$term), unique(res@results$term))
  expect_type(frame$truncation, "double")
  expect_identical(nrow(frame), nrow(res@results))
  expect_identical(frame$bias, res@results$bias)
  # The band is two Monte Carlo standard errors either side of the bias, which
  # is the interval the figure fills rather than one derived from it.
  expect_identical(frame$lower, res@results$bias - 2 * res@results$mc_se)
  expect_identical(frame$upper, res@results$bias + 2 * res@results$mc_se)
  # A discrete exposure truncates a fitted probability, so its sweep is drawn
  # against the lower bound it applied.
  expect_identical(frame$truncation, res@results$truncation_lower)

  plot <- ggplot2::autoplot(res, type = "sweep")
  expect_identical(plot$data, frame)

  built <- ggplot2::ggplot_build(plot)
  points <- layer_drawing(built, "shape")
  expect_equal(points$x, frame$truncation)
  expect_equal(points$y, frame$bias)

  band <- layer_drawing(built, "ymax")
  expect_equal(band$x, frame$truncation)
  expect_equal(band$ymin, frame$lower)
  expect_equal(band$ymax, frame$upper)
})

test_that("a continuous sweep frame is drawn against the quantile levels", {
  local_quiet()
  grid <- c(0, 0.05, 0.1)
  res <- fit_pull_eta(
    data = sim_pull_eta_continuous(),
    exposure_type = "continuous",
    truncation_grid = grid,
    n_boot = 20
  )

  frame <- pull_plot_data(res, type = "sweep")

  # A continuous run caps a stabilized weight rather than bounding a
  # probability, so it reports no lower bound at any level and the quantile
  # levels the caller swept are what the sweep can be drawn against.
  expect_true(all(is.na(res@results$truncation_lower)))
  expect_identical(frame$truncation, grid)
  expect_identical(frame$bias, res@results$bias)

  built <- ggplot2::ggplot_build(ggplot2::autoplot(res, type = "sweep"))
  points <- layer_drawing(built, "shape")
  expect_equal(points$x, frame$truncation)
  expect_equal(points$y, frame$bias)
})

test_that("pull_plot_data() and autoplot() refuse the sweep view alike", {
  local_quiet()
  res <- fit_pull_eta(n_boot = 10)

  pull_error <- rlang::catch_cnd(pull_plot_data(res, type = "sweep"))
  plot_error <- rlang::catch_cnd(ggplot2::autoplot(res, type = "sweep"))

  expect_identical(conditionMessage(pull_error), conditionMessage(plot_error))
  expect_match(deparsed_call(pull_error), "pull_plot_data", fixed = TRUE)
  expect_match(deparsed_call(plot_error), "autoplot", fixed = TRUE)

  expect_snapshot_abort(
    pull_plot_data(res, type = "sweep"),
    class = "positively_sweep_absent_error"
  )
  expect_snapshot_abort(
    ggplot2::autoplot(res, type = "sweep"),
    class = "positively_sweep_absent_error"
  )
})

# ---- Hat values -----------------------------------------------------------

# A dose that tracks its one covariate, so leverage rises at the ends of the
# exposure range and the profile has structure to show.
sim_pull_hat <- function(n = 150, seed = 1) {
  withr::local_seed(seed)
  x1 <- stats::rnorm(n)
  tibble::tibble(dose = stats::rnorm(n, mean = x1), x1 = x1)
}

test_that("the hat-value menu mirrors autoplot's", {
  local_quiet()
  res <- check_hat_values(sim_pull_hat(), dose, x1, null_reps = 5)

  expect_identical(autoplot_types(res), c("null", "profile"))
  expect_identical(pull_plot_data_types(hat_values_result), autoplot_types(res))
  expect_identical(pull_plot_data(res), pull_plot_data(res, type = "null"))
})

test_that("the null frame carries the null draws and the lines beside them", {
  local_quiet()
  # The null draw is seeded so the quantile the frame carries is the one the
  # figure marks.
  res <- withr::with_seed(
    2024,
    check_hat_values(sim_pull_hat(), dose, x1, null_reps = 40)
  )

  frame <- pull_plot_data(res, type = "null")

  expect_identical(
    names(frame),
    c("phi", "null_quantile", "phi_hat", "conf_level")
  )
  for (column in names(frame)) {
    expect_type(frame[[column]], "double")
  }
  # `@null_dist` is the only place the null replicates live. The three scalars
  # the figure marks ride as repeated columns, because a view returns one tibble
  # whatever the figure does with it.
  expect_identical(nrow(frame), 40L)
  expect_identical(frame$phi, res@null_dist)
  expect_identical(frame$null_quantile, rep(res@null_quantile, 40L))
  expect_identical(frame$phi_hat, rep(res@phi_hat, 40L))
  expect_identical(frame$conf_level, rep(res@params$conf_level, 40L))

  plot <- ggplot2::autoplot(res, type = "null")
  expect_identical(plot$data, frame)

  built <- ggplot2::ggplot_build(plot)

  bins <- layer_drawing(built, "count")
  expect_identical(sum(bins$count), as.numeric(nrow(frame)))
  expect_gte(min(frame$phi), min(bins$xmin))
  expect_lte(max(frame$phi), max(bins$xmax))

  # Two reference lines, read by what draws them rather than by the order they
  # were added in: the observed reading is the solid coloured one and the null
  # quantile the dashed one.
  lines <- layers_drawing(built, "xintercept")
  expect_length(lines, 2)
  observed <- Filter(function(line) all(line$colour == "firebrick"), lines)
  dashed <- Filter(function(line) all(line$linetype == "dashed"), lines)
  expect_length(observed, 1)
  expect_length(dashed, 1)
  expect_identical(observed[[1]]$xintercept, unique(frame$phi_hat))
  expect_identical(dashed[[1]]$xintercept, unique(frame$null_quantile))
})

test_that("the profile frame is the fraction the profile view draws", {
  local_quiet()
  res <- check_hat_values(sim_pull_hat(), dose, x1, null_reps = 5)

  frame <- pull_plot_data(res, type = "profile")

  expect_identical(names(frame), c("prob", "fraction"))
  expect_type(frame$prob, "double")
  expect_type(frame$fraction, "double")
  # One row per exposure percentile the candidates were built at, holding the
  # share of that percentile's candidates that read as high leverage. The
  # results carry one row per candidate, so the aggregation is the whole of what
  # the view adds and it belongs in the frame.
  expect_equal(frame$prob, sort(unique(res@results$prob)))
  expect_equal(
    frame$fraction,
    as.double(tapply(res@results$high_leverage, res@results$prob, mean))
  )

  plot <- ggplot2::autoplot(res, type = "profile")
  expect_identical(plot$data, frame)

  built <- ggplot2::ggplot_build(plot)
  points <- layer_drawing(built, "shape")
  line <- layer_drawing(built, "flipped_aes")
  expect_equal(points$x, frame$prob)
  expect_equal(points$y, frame$fraction)
  expect_equal(line$x, frame$prob)
  expect_equal(line$y, frame$fraction)
})

# ---- PoRT -----------------------------------------------------------------

# A deterministic reading-rule anchor. One binary covariate g: 300 rows at
# g == 1 with 3 treated, a prevalence of 0.010, and 700 rows at g == 0 with 350
# treated, a prevalence of 0.500. Nothing is drawn at random, so the reported
# subgroups are the same two on every run.
sim_pull_port <- function() {
  tibble::tibble(
    exposure = c(rep(1L, 3), rep(0L, 297), rep(1L, 350), rep(0L, 350)),
    g = c(rep(1L, 300), rep(0L, 700))
  )
}

# A censoring result assembled directly, so the exposure and the censoring
# family carry thresholds that differ from each other and from wave to wave.
# What is under test is which pair of thresholds a row was judged against, and a
# run whose thresholds all resolve to one number cannot show it.
make_pull_port_seq <- function() {
  results <- tibble::tibble(
    time = c(1L, 1L, 2L, 2L),
    type = c("exposure", "censoring", "exposure", "censoring"),
    subgroup = rep("g", 4L),
    description = c("g>=1", "g>=2", "g>=3", "g>=4"),
    exposure_level = rep("1", 4L),
    n = c(10L, 20L, 30L, 40L),
    proportion = c(0.1, 0.2, 0.3, 0.4),
    prevalence = c(0.01, 0.02, 0.03, 0.04),
    low_support = c(TRUE, TRUE, FALSE, TRUE)
  )
  port_result(
    results = results,
    exposure = c("a1", "a2"),
    exposure_type = "binary",
    n = 100L,
    params = list(),
    call = quote(check_port_seq()),
    trees = list(),
    alpha = 0.05,
    beta = c(0.05, 0.1),
    censoring_beta = c(0.2, 0.3),
    gamma = 2
  )
}

test_that("the PoRT frame is the plotting frame and the thresholds behind it", {
  local_quiet()
  res <- check_port(sim_pull_port(), exposure, g)

  # One view, so no view to name. The accessor takes the argument the plot takes
  # instead, and in the same position.
  expect_identical(
    names(formals(S7::method(pull_plot_data, port_result))),
    c("x", "low_support_only", "...")
  )
  expect_identical(
    names(formals(autoplot_method(res)))[[2]],
    "low_support_only"
  )

  frame <- pull_plot_data(res)

  expect_identical(
    names(frame),
    c(
      "subgroup",
      "description",
      "exposure_level",
      "n",
      "proportion",
      "prevalence",
      "low_support",
      "label",
      "width",
      "beta_lower",
      "beta_upper"
    )
  )
  expect_s3_class(frame$low_support, "factor")
  expect_identical(levels(frame$low_support), c("FALSE", "TRUE"))
  expect_type(frame$beta_lower, "double")
  expect_type(frame$beta_upper, "double")

  # Duplicate rows, which arise when the same leaf is reported by more than one
  # tree, are collapsed to one bar, so the frame is the reported subgroups
  # rather than the results.
  expect_identical(nrow(frame), 2L)
  expect_equal(frame$prevalence, c(0.01, 0.5))
  expect_identical(as.character(frame$low_support), c("TRUE", "FALSE"))

  # Every row carries the pair of thresholds its prevalence was read against,
  # which is what a reader needs to see why a subgroup reads as low support.
  expect_identical(frame$beta_lower, rep(res@beta[[1]], nrow(frame)))
  expect_identical(frame$beta_upper, rep(1 - res@beta[[1]], nrow(frame)))

  plot <- ggplot2::autoplot(res)
  expect_identical(plot$data, frame)

  built <- ggplot2::ggplot_build(plot)
  bars <- layer_drawing(built, "xmax")
  expect_identical(nrow(bars), nrow(frame))
  drawn <- bars[order(bars$xmax), ]
  expected <- frame[order(frame$prevalence), ]
  expect_equal(drawn$xmax, expected$prevalence)
  # The bars run horizontally, so the width proportional to subgroup size is the
  # extent each bar takes in y.
  expect_equal(drawn$ymax - drawn$ymin, expected$width)

  lines <- layer_drawing(built, "xintercept")
  expect_setequal(
    lines$xintercept,
    c(unique(frame$beta_lower), unique(frame$beta_upper))
  )
})

test_that("each PoRT row carries the thresholds its own facet is drawn with", {
  res <- make_pull_port_seq()

  frame <- pull_plot_data(res)
  reference <- port_beta_reference(res)

  expect_true(all(c("time", "type") %in% names(frame)))
  expect_type(frame$beta_lower, "double")
  expect_type(frame$beta_upper, "double")
  # A censoring run judges the two families on different risk sets, so a row's
  # thresholds follow its family as well as its wave.
  expect_identical(frame$beta_lower, c(0.05, 0.2, 0.1, 0.3))
  expect_identical(frame$beta_upper, 1 - c(0.05, 0.2, 0.1, 0.3))

  # autoplot() draws its reference lines from a frame built per facet. The two
  # have to name the same pair wherever a facet has a bar in it, or a row would
  # be read against one pair and drawn against another.
  for (row in seq_len(nrow(frame))) {
    facet <- reference[
      reference$time == frame$time[[row]] &
        reference$type == frame$type[[row]],
    ]
    expect_setequal(
      facet$xintercept,
      c(frame$beta_lower[[row]], frame$beta_upper[[row]])
    )
  }
})

test_that("low_support_only filters the frame and the bars alike", {
  local_quiet()
  res <- check_port(sim_pull_port(), exposure, g)

  full <- pull_plot_data(res)
  filtered <- pull_plot_data(res, low_support_only = TRUE)

  expect_identical(names(filtered), names(full))
  expect_identical(nrow(filtered), 1L)
  expect_identical(as.character(filtered$low_support), "TRUE")
  # The widths are read against the largest subgroup still shown rather than
  # against one the filter dropped.
  expect_equal(max(filtered$width), 1)

  plot <- ggplot2::autoplot(res, low_support_only = TRUE)
  expect_identical(plot$data, filtered)

  built <- ggplot2::ggplot_build(plot)
  expect_identical(nrow(layer_drawing(built, "xmax")), nrow(filtered))
})

test_that("pull_plot_data() and autoplot() refuse a non-flag alike", {
  local_quiet()
  res <- check_port(sim_pull_port(), exposure, g)

  pull_error <- rlang::catch_cnd(pull_plot_data(res, low_support_only = "yes"))
  plot_error <- rlang::catch_cnd(
    ggplot2::autoplot(res, low_support_only = "yes")
  )

  # The argument the two entry points share is validated once behind both, so
  # the reader is told the same thing either way and is shown the call they
  # made.
  expect_identical(conditionMessage(pull_error), conditionMessage(plot_error))
  expect_match(deparsed_call(pull_error), "pull_plot_data", fixed = TRUE)
  expect_match(deparsed_call(plot_error), "autoplot", fixed = TRUE)

  expect_snapshot_abort(
    pull_plot_data(res, low_support_only = "yes"),
    class = "positively_type_error"
  )
  expect_snapshot_abort(
    ggplot2::autoplot(res, low_support_only = "yes"),
    class = "positively_type_error"
  )
})

# ---- HDR ------------------------------------------------------------------

sim_pull_hdr <- function(n = 150, seed = 1) {
  withr::local_seed(seed)
  l <- stats::rnorm(n)
  tibble::tibble(exposure = stats::rnorm(n, mean = l), l = l)
}

test_that("the point HDR frame is the results the curve draws", {
  local_quiet()
  res <- check_hdr(sim_pull_hdr(), exposure, l, values = c(-1, 0, 1))

  # One view, so there is no view to name and neither signature offers a menu.
  expect_null(formals(S7::method(pull_plot_data, hdr_result))$type)
  expect_null(formals(autoplot_method(res))$type)

  frame <- pull_plot_data(res)

  # A point result has nothing to key the curves on, so the frame is the results
  # as they stand.
  expect_identical(frame, res@results)
  expect_identical(names(frame), c("value", "nonoverlap"))

  plot <- ggplot2::autoplot(res)
  expect_identical(plot$data, frame)

  built <- ggplot2::ggplot_build(plot)
  points <- layer_drawing(built, "shape")
  expect_equal(points$x, frame$value)
  expect_equal(points$y, frame$nonoverlap)
})

test_that("the sequential HDR frame keys the curves on the time point", {
  local_quiet()
  data <- dgp_longitudinal(n = 200, seed = 5)
  res <- check_hdr_seq(
    data,
    c(a1, a2, a3),
    list(l0, l1, l2),
    values = c(-1, 0, 1)
  )

  frame <- pull_plot_data(res)

  # The wave separates the curves, so it is a factor here where the results
  # carry it as the integer index. Left numeric it would colour the curves off a
  # continuous scale, which reads a wave as a quantity.
  expect_identical(names(frame), names(res@results))
  expect_s3_class(frame$time, "factor")
  expect_identical(levels(frame$time), c("1", "2", "3"))
  expect_identical(as.integer(as.character(frame$time)), res@results$time)
  expect_identical(frame$value, res@results$value)
  expect_identical(frame$nonoverlap, res@results$nonoverlap)

  plot <- ggplot2::autoplot(res)
  expect_identical(plot$data, frame)

  built <- ggplot2::ggplot_build(plot)
  points <- layer_drawing(built, "shape")
  # One curve per wave, together drawing every row of the frame.
  expect_length(unique(points$group), nlevels(frame$time))
  expect_identical(nrow(points), nrow(frame))
  expect_equal(sort(points$y), sort(frame$nonoverlap))
})

# ---- Extrapolation --------------------------------------------------------

# Two exposure groups drawn from one Gaussian covariate cloud, so both groups
# cover each other and the hull test has a mix of inside and outside rows to
# report.
sim_pull_extrapolation <- function(n = 60, p = 3, seed = 1) {
  withr::local_seed(seed)
  covariates <- matrix(stats::rnorm(n * p), nrow = n, ncol = p)
  colnames(covariates) <- paste0("x", seq_len(p))
  out <- tibble::as_tibble(as.data.frame(covariates))
  n0 <- n %/% 2L
  out$exposure <- c(rep(0L, n0), rep(1L, n - n0))
  out[c("exposure", paste0("x", seq_len(p)))]
}

test_that("the extrapolation menu mirrors autoplot's", {
  local_quiet()
  res <- check_extrapolation(
    sim_pull_extrapolation(),
    exposure,
    tidyselect::starts_with("x"),
    hull = FALSE
  )

  expect_identical(autoplot_types(res), c("distribution", "hull"))
  expect_identical(
    pull_plot_data_types(extrapolation_result),
    autoplot_types(res)
  )
  expect_identical(
    pull_plot_data(res),
    pull_plot_data(res, type = "distribution")
  )
})

test_that("the distribution frame groups the results by exposure", {
  local_quiet()
  res <- check_extrapolation(
    sim_pull_extrapolation(),
    exposure,
    tidyselect::starts_with("x"),
    hull = FALSE
  )

  frame <- pull_plot_data(res, type = "distribution")

  # The exposure panels the view rather than measuring anything, so it is a
  # factor here where the results carry the raw level.
  expect_identical(names(frame), names(res@results))
  expect_s3_class(frame$exposure, "factor")
  expect_identical(levels(frame$exposure), c("0", "1"))
  expect_identical(frame$frac_nearby, res@results$frac_nearby)

  plot <- ggplot2::autoplot(res, type = "distribution")
  expect_identical(plot$data, frame)

  built <- ggplot2::ggplot_build(plot)
  # One panel per exposure group, in level order, with every row of the frame
  # binned inside the panel of its own group.
  expect_identical(
    as.character(built$layout$layout$exposure),
    levels(frame$exposure)
  )
  bins <- layer_drawing(built, "count")
  expect_identical(
    unname(vapply(split(bins$count, bins$PANEL), sum, numeric(1))),
    as.numeric(table(frame$exposure))
  )
})

test_that("the hull frame is the membership the bars draw", {
  local_quiet()
  skip_if_not_installed("lpSolve")
  res <- check_extrapolation(
    sim_pull_extrapolation(),
    exposure,
    tidyselect::starts_with("x"),
    hull = TRUE
  )

  frame <- pull_plot_data(res, type = "hull")

  # The view stacks the bars on a reading of `in_hull` rather than on the
  # logical itself, so the reading is a column of the frame and reads the way
  # the legend does.
  expect_identical(names(frame), c(names(res@results), "membership"))
  expect_s3_class(frame$exposure, "factor")
  expect_s3_class(frame$membership, "factor")
  expect_identical(levels(frame$membership), c("Inside hull", "Outside hull"))
  expect_identical(frame$membership == "Inside hull", res@results$in_hull)

  plot <- ggplot2::autoplot(res, type = "hull")
  expect_identical(plot$data, frame)

  built <- ggplot2::ggplot_build(plot)
  bars <- layer_drawing(built, "count")
  # One bar per exposure group per membership, counting the frame's rows and no
  # others.
  expect_identical(nrow(bars), 4L)
  expect_identical(sum(bars$count), as.numeric(nrow(frame)))
  expect_identical(
    unname(vapply(split(bars$count, bars$x), sum, numeric(1))),
    as.numeric(table(frame$exposure))
  )
  expect_setequal(
    bars$count,
    as.numeric(table(frame$exposure, frame$membership))
  )

  # The key names the memberships the frame carries, in the frame's own order.
  guide <- ggplot2::get_guide_data(plot, "fill")
  expect_identical(as.character(guide$.label), levels(frame$membership))
})

test_that("pull_plot_data() and autoplot() refuse the hull view alike", {
  local_quiet()
  res <- check_extrapolation(
    sim_pull_extrapolation(),
    exposure,
    tidyselect::starts_with("x"),
    hull = FALSE
  )

  pull_error <- rlang::catch_cnd(pull_plot_data(res, type = "hull"))
  plot_error <- rlang::catch_cnd(ggplot2::autoplot(res, type = "hull"))

  expect_identical(conditionMessage(pull_error), conditionMessage(plot_error))
  expect_match(deparsed_call(pull_error), "pull_plot_data", fixed = TRUE)
  expect_match(deparsed_call(plot_error), "autoplot", fixed = TRUE)

  expect_snapshot_abort(
    pull_plot_data(res, type = "hull"),
    class = "positively_hull_absent_error"
  )
  expect_snapshot_abort(
    ggplot2::autoplot(res, type = "hull"),
    class = "positively_hull_absent_error"
  )
})

# ---- Positivity check -----------------------------------------------------

sim_pull_check <- function(n = 200, seed = 1) {
  withr::local_seed(seed)
  x1 <- stats::rnorm(n)
  tibble::tibble(
    exposure = stats::rbinom(n, 1L, stats::plogis(3 * x1)),
    x1 = x1
  )
}

test_that("a positivity check delegates to the diagnostic it is named", {
  local_quiet()
  check <- check_positivity(
    sim_pull_check(),
    exposure,
    x1,
    diagnostics = c("port", "edp")
  )

  # A container holds children rather than a frame of its own, so naming one is
  # the whole of what this does: the frame is the child's, unchanged.
  expect_identical(pull_plot_data(check, "port"), pull_plot_data(check$port))
  expect_identical(pull_plot_data(check, "edp"), pull_plot_data(check$edp))

  # A view name reaches the child through the dots, so a child's menu is
  # reachable without extracting it first.
  expect_identical(
    pull_plot_data(check, "edp", type = "ecdf"),
    pull_plot_data(check$edp, type = "ecdf")
  )
})

test_that("a positivity check refuses to guess which diagnostic to pull", {
  local_quiet()
  check <- check_positivity(
    sim_pull_check(),
    exposure,
    x1,
    diagnostics = c("port", "edp")
  )

  # A container holds one frame per child and this returns one tibble, so
  # without a name there is nothing to return. autoplot() has a panel to fall
  # back on; there is no such thing for a frame.
  missing_error <- rlang::catch_cnd(pull_plot_data(check))
  expect_s3_class(missing_error, "positively_diagnostic_error")
  for (name in names(check)) {
    expect_match(conditionMessage(missing_error), name, fixed = TRUE)
  }

  # autoplot() spells the panel as a `NULL` diagnostic, so a reader carrying the
  # argument over lands on the same refusal rather than on an internal failure.
  null_error <- rlang::catch_cnd(pull_plot_data(check, NULL))
  expect_s3_class(null_error, "positively_diagnostic_error")
  expect_identical(
    conditionMessage(null_error),
    conditionMessage(missing_error)
  )

  # extract_named_check() is what `$` and `[[` refuse an unknown name with, so
  # every way of naming a child spells the refusal the same, and the reader is
  # shown the call they made.
  unknown_error <- rlang::catch_cnd(pull_plot_data(check, "bogus"))
  extract_error <- rlang::catch_cnd(check[["bogus"]])
  expect_s3_class(unknown_error, "positively_diagnostic_error")
  expect_identical(
    conditionMessage(unknown_error),
    conditionMessage(extract_error)
  )
  expect_match(deparsed_call(unknown_error), "pull_plot_data", fixed = TRUE)

  expect_snapshot_abort(
    pull_plot_data(check),
    class = "positively_diagnostic_error"
  )
  expect_snapshot_abort(
    pull_plot_data(check, "bogus"),
    class = "positively_diagnostic_error"
  )
})
