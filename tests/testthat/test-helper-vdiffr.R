# `expect_doppelganger()` has to name a figure's snapshot file before it is
# allowed to skip, and vdiffr derives that name with the unexported
# str_standardise(), so `doppelganger_file()` reimplements the rule rather than
# calling it. The two have to agree exactly: announcing a name that does not
# match the file on disk leaves the real baseline looking unused, and the
# cleanup pass at the end of the run deletes it.

test_that("doppelganger_file() reproduces the vdiffr slug rule", {
  skip_if_not_installed("vdiffr")
  titles <- c(
    # A sample of the titles this suite pins, not the whole set.
    "Extrapolation frac nearby distribution",
    "EDP estimator scatter with infinite weights",
    "sPoRT subgroups faceted by time",
    "HDR non-overlap point curve",
    "Density ratios cumulative quantiles",
    # Shapes the live titles never reach.
    "Consecutive   separators -- collapse",
    "2024 leads with a digit",
    "Trailing punctuation!!!"
  )

  # Reaching the unexported function is the point of the test, so it comes
  # through getFromNamespace() rather than `:::`.
  str_standardise <- utils::getFromNamespace("str_standardise", "vdiffr")
  upstream <- vapply(
    titles,
    str_standardise,
    character(1),
    sep = "-",
    USE.NAMES = FALSE
  )
  expect_identical(
    vapply(titles, doppelganger_file, character(1), USE.NAMES = FALSE),
    paste0(upstream, ".svg")
  )
})
