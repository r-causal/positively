# A thin wrapper over vdiffr::expect_doppelganger() that skips when vdiffr is
# not installed, so machines without it, including CRAN, skip the visual
# comparison rather than failing.
expect_doppelganger <- function(title, fig, ...) {
  testthat::skip_if_not_installed("vdiffr")
  vdiffr::expect_doppelganger(title, fig, ...)
}
