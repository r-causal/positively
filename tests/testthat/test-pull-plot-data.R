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
autoplot_types <- function(object) {
  method <- utils::getS3method("autoplot", class(object)[[1]])
  eval(formals(method)$type)
}

pull_plot_data_types <- function(class) {
  eval(formals(S7::method(pull_plot_data, class))$type)
}

deparsed_call <- function(condition) {
  paste(deparse(conditionCall(condition)), collapse = " ")
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

  frame <- pull_plot_data(res)

  # The shared label is one level, and neither intervention's rows are dropped
  # or turned into missing levels along the way.
  expect_identical(levels(frame$intervention), "1")
  expect_identical(nrow(frame), 2L * nrow(data))
  expect_false(anyNA(frame$intervention))

  for (type in c("histogram", "ecdf")) {
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
