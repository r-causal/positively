# The contract every result class's display methods hold to. The container
# report is assembled from `diagnostic_label()` and `diagnostic_headline()`, so
# both are checked against fitted objects: what a reader sees is the subject,
# not how the package spells anything internally.

# Collapse rendered output to one string, so that a match is not defeated by
# wherever cli wrapped a line, and drop the big mark, so that a count matches
# whatever grouping the renderer applied to it.
rendered_text <- function(lines) {
  gsub(",", "", paste(lines, collapse = " "), fixed = TRUE)
}

printed_lines <- function(x) {
  utils::capture.output(print(x))
}

printed_text <- function(x) {
  rendered_text(printed_lines(x))
}

# The report names each section for the diagnostic that produced it, which is
# the name `names()` lists and `$` extracts by. Positions are read off the whole
# rendered block rather than line by line, because cli wraps a long line at a
# space, and the name is matched on word boundaries so that "port" is not found
# inside "support" or "reported".
section_positions <- function(container) {
  text <- printed_text(container)
  vapply(
    names(container),
    function(name) regexpr(paste0("\\b", name, "\\b"), text)[[1]],
    integer(1)
  )
}

# The diagnostic whose section the report opens with, read off the positions so
# that a block asserting on the leader has already guarded them: an unfound name
# reports -1, which would otherwise win `which.min()`.
leading_section <- function(positions) {
  names(positions)[[which.min(positions)]]
}

# The label is what a diagnostic is called wherever it is named for a reader, so
# it is prose and not an identifier. That rules out the S7 class name, which is
# internal and undocumented, and it rules out the punctuation the package spells
# symbols with: no underscore and no period anywhere in the label.
expect_readable_label <- function(x) {
  label <- diagnostic_label(x)
  testthat::expect_type(label, "character")
  testthat::expect_length(label, 1L)
  testthat::expect_false(is.na(label))
  testthat::expect_gt(nchar(label), 0L)
  testthat::expect_false(identical(label, S7::S7_class(x)@name))
  testthat::expect_no_match(label, "[_.]")
  invisible(label)
}

# The headline is the reading printed under a diagnostic's section of the
# container report. Each element is one line of that section, so none carries a
# line break of its own or falls back to the internal class name, and the
# reading states a number: a section that named its diagnostic and stopped would
# have reported nothing.
expect_readable_headline <- function(x) {
  headline <- diagnostic_headline(x)
  testthat::expect_type(headline, "character")
  testthat::expect_gte(length(headline), 1L)
  testthat::expect_false(anyNA(headline))
  testthat::expect_true(all(nchar(headline) > 0L))
  testthat::expect_no_match(headline, "\n", fixed = TRUE)
  rendered <- rendered_text(headline)
  testthat::expect_no_match(rendered, S7::S7_class(x)@name, fixed = TRUE)
  testthat::expect_match(rendered, "[0-9]")
  invisible(headline)
}
