# A thin wrapper over vdiffr::expect_doppelganger() that skips when vdiffr is
# not installed, so machines without it, including CRAN, skip the visual
# comparison rather than failing.
#
# testthat keeps the `.md` snapshot entries of a skipped test, but a file
# snapshot survives only when the test announces it. An unannounced `.svg`
# looks unused to the cleanup pass at the end of the run and is deleted, so
# every skip has to come after the announcement. The skip inside the wrapper
# abandons the rest of its block, which strands any later figure in the same
# block: a block that renders more than one figure, or that skips on some other
# condition, must open with an `announce_doppelganger()` call naming every
# figure it pins.

# Reproduces vdiffr:::str_standardise(), the unexported function that turns a
# doppelganger title into its snapshot file name.
doppelganger_file <- function(title) {
  slug <- gsub("[^a-z0-9]", "-", tolower(title))
  slug <- gsub("--+", "-", slug)
  slug <- gsub("^-|-$", "", slug)
  paste0(slug, ".svg")
}

announce_doppelganger <- function(...) {
  for (title in c(...)) {
    testthat::announce_snapshot_file(name = doppelganger_file(title))
  }
  invisible(NULL)
}

expect_doppelganger <- function(title, fig, ...) {
  announce_doppelganger(title)
  testthat::skip_if_not_installed("vdiffr")
  vdiffr::expect_doppelganger(title, fig, ...)
}
