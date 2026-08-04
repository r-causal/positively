# The abstract positivity_diagnostic parent cannot be instantiated, so a
# throwaway concrete subclass lets these specs exercise the inherited
# properties, validator, and shared tidy()/glance()/print wiring.
make_test_diagnostic <- function(
  results = tibble::tibble(.id = 1:3, value = c(0.1, 0.2, 0.3)),
  exposure = "a",
  exposure_type = "binary",
  n = 3L,
  params = list(alpha = 0.05),
  call = quote(check_test())
) {
  test_diagnostic <- S7::new_class(
    "test_diagnostic",
    parent = positivity_diagnostic
  )
  test_diagnostic(
    results = results,
    exposure = exposure,
    exposure_type = exposure_type,
    n = n,
    params = params,
    call = call
  )
}

# positivity_check carries its own metadata alongside a named list of children,
# so a hand-built container needs all six values. These defaults stand in for
# the ones a block is not exercising; every block overrides what it asserts on.
make_test_container <- function(
  checks = list(test = make_test_diagnostic()),
  exposure = "a",
  exposure_type = "binary",
  covariates = c("x1", "x2"),
  n = 3L,
  call = quote(check_positivity())
) {
  positivity_check(
    checks = checks,
    exposure = exposure,
    exposure_type = exposure_type,
    covariates = covariates,
    n = n,
    call = call
  )
}

test_that("positivity_diagnostic is an S7 class", {
  expect_true(S7::S7_inherits(make_test_diagnostic(), positivity_diagnostic))
})

test_that("positivity_diagnostic is abstract and cannot be instantiated", {
  expect_error(
    positivity_diagnostic(
      results = tibble::tibble(),
      exposure = "a",
      exposure_type = "binary",
      n = 0L,
      params = list(),
      call = quote(f())
    ),
    regexp = "abstract"
  )
})

test_that("positivity_diagnostic declares the fixed property set", {
  expect_setequal(
    names(positivity_diagnostic@properties),
    c("results", "exposure", "exposure_type", "n", "params", "call")
  )
})

test_that("the parent validator accepts a valid exposure_type", {
  obj <- make_test_diagnostic(exposure_type = "continuous")
  expect_identical(obj@exposure_type, "continuous")
})

test_that("the parent validator rejects an unknown exposure_type", {
  expect_error(
    make_test_diagnostic(exposure_type = "ordinal"),
    regexp = "exposure_type"
  )
})

test_that("tidy() returns the results tibble of a diagnostic", {
  results <- tibble::tibble(.id = 1:2, value = c(0.4, 0.6))
  obj <- make_test_diagnostic(results = results)
  expect_identical(generics::tidy(obj), results)
})

test_that("glance() returns a one-row tibble for a diagnostic", {
  glanced <- generics::glance(make_test_diagnostic())
  expect_s3_class(glanced, "tbl_df")
  expect_identical(nrow(glanced), 1L)
})

test_that("the inherited glance() carries the sample size and nothing else", {
  # A subclass that adds no statistics of its own has only the sample size to
  # report. The exposure, its type, and the results row count are metadata, so
  # the default must not smuggle them in as statistics.
  glanced <- generics::glance(make_test_diagnostic(n = 42L))

  expect_setequal(names(glanced), "n")
  expect_identical(glanced$n, 42L)
})

test_that("the default summary() reports a diagnostic's own statistics", {
  # The inherited glance() reports only the sample size, and the sample size is
  # metadata the container states once rather than a statistic a diagnostic
  # computed, so a subclass that adds nothing of its own summarizes to no rows.
  # The column set and its types still have to hold, because the container
  # stacks these rows across its children.
  summarized <- summary(make_test_diagnostic(n = 42L))

  expect_s3_class(summarized, "tbl_df")
  expect_identical(names(summarized), c("statistic", "value", "threshold"))
  expect_identical(nrow(summarized), 0L)
  expect_type(summarized$statistic, "character")
  expect_type(summarized$value, "double")
  expect_type(summarized$threshold, "double")
})

test_that("summary() on a container with no children keeps the column set", {
  summarized <- summary(make_test_container(checks = list()))

  expect_s3_class(summarized, "tbl_df")
  expect_identical(
    names(summarized),
    c("diagnostic", "statistic", "value", "threshold")
  )
  expect_identical(nrow(summarized), 0L)
  expect_type(summarized$diagnostic, "character")
  expect_type(summarized$statistic, "character")
  expect_type(summarized$value, "double")
  expect_type(summarized$threshold, "double")
})

test_that("the parent validator rejects a non-scalar exposure_type", {
  expect_error(
    make_test_diagnostic(exposure_type = c("binary", "continuous")),
    regexp = "exposure_type"
  )
})

test_that("the parent validator rejects an empty exposure_type", {
  expect_error(
    make_test_diagnostic(exposure_type = character(0)),
    regexp = "exposure_type"
  )
})

test_that("the parent validator rejects a non-scalar n", {
  expect_error(
    make_test_diagnostic(n = c(1L, 2L)),
    regexp = "@n"
  )
})

test_that("the parent validator rejects a missing n", {
  expect_error(
    make_test_diagnostic(n = NA_integer_),
    regexp = "@n"
  )
})

test_that("the parent validator rejects an empty exposure", {
  expect_error(
    make_test_diagnostic(exposure = character(0)),
    regexp = "@exposure"
  )
})

# ---- The container's property set and metadata ----------------------------

test_that("positivity_check declares its property set", {
  expect_setequal(
    names(positivity_check@properties),
    c("checks", "exposure", "exposure_type", "covariates", "n", "call")
  )
})

test_that("positivity_check stores validated diagnostic children under their names", {
  obj <- make_test_diagnostic()
  container <- make_test_container(checks = list(test = obj))
  expect_length(container@checks, 1)
  expect_identical(names(container@checks), "test")
  expect_identical(container@checks[["test"]], obj)
})

test_that("positivity_check carries the metadata it was constructed with", {
  originating_call <- quote(check_positivity(df, dose, c(x1, x2, region)))
  container <- make_test_container(
    exposure = "dose",
    exposure_type = "continuous",
    covariates = c("x1", "x2", "region"),
    n = 500L,
    call = originating_call
  )
  expect_identical(container@exposure, "dose")
  expect_identical(container@exposure_type, "continuous")
  expect_identical(container@covariates, c("x1", "x2", "region"))
  expect_identical(container@n, 500L)
  expect_identical(container@call, originating_call)
})

# ---- The container validator ----------------------------------------------

test_that("positivity_check rejects non-diagnostic children", {
  expect_error(
    make_test_container(checks = list(bad = 1L)),
    regexp = "positivity_diagnostic"
  )
})

test_that("positivity_check rejects an unnamed checks list", {
  obj <- make_test_diagnostic()
  expect_error(
    make_test_container(checks = list(obj, obj)),
    regexp = "named"
  )
})

test_that("positivity_check rejects an empty or missing diagnostic name", {
  obj <- make_test_diagnostic()

  blank <- list(obj, obj)
  names(blank) <- c("edp", "")
  expect_error(make_test_container(checks = blank), regexp = "named")

  missing_name <- list(obj, obj)
  names(missing_name) <- c("edp", NA_character_)
  expect_error(make_test_container(checks = missing_name), regexp = "named")
})

test_that("positivity_check rejects duplicate diagnostic names", {
  obj <- make_test_diagnostic()
  repeated <- list(obj, obj)
  names(repeated) <- c("port", "port")
  expect_error(make_test_container(checks = repeated), regexp = "repeat")
})

test_that("the container validator rejects an unknown exposure_type", {
  expect_error(
    make_test_container(exposure_type = "ordinal"),
    regexp = "@exposure_type"
  )
})

test_that("the container validator rejects a non-scalar exposure_type", {
  expect_error(
    make_test_container(exposure_type = c("binary", "continuous")),
    regexp = "@exposure_type"
  )
})

test_that("the container validator rejects an empty exposure", {
  expect_error(
    make_test_container(exposure = character(0)),
    regexp = "@exposure"
  )
})

test_that("the container validator rejects empty covariates", {
  expect_error(
    make_test_container(covariates = character(0)),
    regexp = "@covariates"
  )
})

test_that("the container validator rejects a non-scalar n", {
  expect_error(
    make_test_container(n = c(1L, 2L)),
    regexp = "@n"
  )
})

test_that("the container validator rejects a missing n", {
  expect_error(
    make_test_container(n = NA_integer_),
    regexp = "@n"
  )
})

# ---- Naming and extracting children ---------------------------------------

test_that("names() returns the names of the checks list", {
  res <- make_test_container(
    checks = list(edp = make_test_diagnostic(), port = make_test_diagnostic())
  )
  expect_identical(names(res), c("edp", "port"))
  expect_identical(names(res), names(res@checks))
})

test_that("names() on a container with no children is character(0)", {
  expect_identical(names(make_test_container(checks = list())), character(0))
})

test_that("`[[` extracts a diagnostic by name and by whole-number position", {
  edp <- make_test_diagnostic()
  port <- make_test_diagnostic()
  res <- make_test_container(checks = list(edp = edp, port = port))
  expect_identical(res[["port"]], port)
  expect_identical(res[[2]], port)
  expect_identical(res[[2]], res[["port"]])
})

test_that("`[[` rejects invalid non-character indices", {
  res <- make_test_container(
    checks = list(edp = make_test_diagnostic(), port = make_test_diagnostic())
  )
  expect_error(res[[NA]], class = "positively_error")
  expect_error(res[[TRUE]], class = "positively_error")
  expect_error(res[[1.9]], class = "positively_error")
})

test_that("`[[` rejects a non-scalar index", {
  res <- make_test_container(
    checks = list(edp = make_test_diagnostic(), port = make_test_diagnostic())
  )
  expect_error(res[[c("edp", "port")]], class = "positively_diagnostic_error")
})

test_that("`[[` rejects an unknown diagnostic name", {
  res <- make_test_container(
    checks = list(edp = make_test_diagnostic(), port = make_test_diagnostic())
  )
  expect_error(res[["nonexistent"]], class = "positively_diagnostic_error")
})

test_that("`[[` rejects an out-of-bounds position", {
  res <- make_test_container(
    checks = list(edp = make_test_diagnostic(), port = make_test_diagnostic())
  )
  expect_error(res[[3]], class = "positively_bounds_error")
})

test_that("`$` extracts a diagnostic by name", {
  edp <- make_test_diagnostic()
  port <- make_test_diagnostic()
  res <- make_test_container(checks = list(edp = edp, port = port))
  expect_identical(res$port, port)
  expect_identical(res$port, res[["port"]])
})

test_that("`$` rejects an unknown diagnostic name rather than returning NULL", {
  res <- make_test_container(
    checks = list(edp = make_test_diagnostic(), port = make_test_diagnostic())
  )
  expect_error(res$nonexistent, class = "positively_diagnostic_error")
})

test_that("`$` does not partially match a diagnostic name", {
  res <- make_test_container(
    checks = list(edp = make_test_diagnostic(), port = make_test_diagnostic())
  )
  expect_error(res$por, class = "positively_diagnostic_error")
})

# Completion after `check$` runs through .DollarNames.default, which matches the
# supplied pattern against `names(x)`. The container defines no .DollarNames
# method of its own, so what completion offers depends entirely on the
# container's `names()` method. Completion passes a `^`-anchored regular
# expression rather than a bare prefix, so the pattern is matched as a regexp.
test_that(".DollarNames offers the diagnostic names and honours a pattern", {
  res <- make_test_container(
    checks = list(edp = make_test_diagnostic(), port = make_test_diagnostic())
  )
  expect_identical(utils::.DollarNames(res, ""), c("edp", "port"))
  expect_identical(utils::.DollarNames(res, "^p"), "port")
  expect_identical(utils::.DollarNames(res, "^z"), character(0))
})

# ---- Printing and pinned messages -----------------------------------------

test_that("printing a positivity_check produces sectioned output", {
  container <- make_test_container(checks = list(test = make_test_diagnostic()))
  expect_output(print(container), "test")
})

test_that("the invalid-index message is stable", {
  res <- make_test_container(
    checks = list(edp = make_test_diagnostic(), port = make_test_diagnostic())
  )
  expect_snapshot_abort(res[[1.9]])
})

test_that("the unknown-diagnostic message is stable", {
  res <- make_test_container(
    checks = list(edp = make_test_diagnostic(), port = make_test_diagnostic())
  )
  expect_snapshot_abort(res[["nonexistent"]])
  expect_snapshot_abort(res$nonexistent)
})

test_that("printing a diagnostic is stable", {
  expect_snapshot(print(make_test_diagnostic()))
})

test_that("printing a container with two children is stable", {
  container <- make_test_container(
    checks = list(edp = make_test_diagnostic(), port = make_test_diagnostic())
  )
  expect_snapshot(print(container))
})

# ---- The container report --------------------------------------------------

# The report states the exposure, its type, the sample size, and the covariates
# from the container's own metadata. A block that counts how often the report
# states one of them needs tokens nothing else in the output carries, so the
# shared generators' columns are renamed: counting "exposure" would otherwise
# count a header label as well as the column that label names.
report_data <- function(data) {
  names(data)[1:3] <- c("statin", "bmi", "ldl")
  data
}

# One child holding low-support rows. Good positivity leaves PoRT with no
# subgroup over its reading rule, while two of the 400 units sit further than
# one geometric variability from the nearest unit of the opposite exposure, so
# extrapolation holds the only rows and is the last diagnostic requested.
one_low_support_container <- function() {
  check_positivity(
    report_data(dgp_good_positivity(n = 400, seed = 1)),
    statin,
    c(bmi, ldl)
  )
}

# A continuous run holding only EDP and the leverage check, so the section with
# something to report is the one requested second. The leverage check reports a
# single statistic rather than a set of rows, and the null draw is seeded because
# the quantile phi-hat is compared against is drawn inside the run.
scalar_finding_container <- function() {
  withr::local_seed(2024)
  check_positivity(
    dgp_continuous_support_gap(n = 200, seed = 4),
    exposure,
    x1,
    exposure_type = "continuous",
    diagnostics = c("edp", "hat_values"),
    args = list(hat_values = list(null_reps = 25))
  )
}

# Two children holding low-support rows, where the one requested later holds
# far more of them.
two_low_support_container <- function() {
  check_positivity(
    report_data(dgp_practical_violation(n = 317, seed = 2)),
    statin,
    c(bmi, ldl)
  )
}

# Rebuild a container with its children in a different order and nothing else
# changed, so that a block can vary the order the diagnostics were requested in.
reorder_checks <- function(container, order) {
  positivity_check(
    checks = container@checks[order],
    exposure = container@exposure,
    exposure_type = container@exposure_type,
    covariates = container@covariates,
    n = container@n,
    call = container@call
  )
}

# The block the report ends on. The footer is one sentence, which wraps at the
# width testthat pins, so the search covers the block rather than its last
# physical line. cli separates that block from the section above it with an
# empty line.
trailing_block <- function(lines) {
  filled <- which(nzchar(trimws(lines)))
  lines <- lines[seq_len(max(filled))]
  breaks <- which(!nzchar(trimws(lines)))
  start <- if (length(breaks) > 0) max(breaks) + 1L else 1L
  rendered_text(utils::tail(lines, length(lines) - start + 1L))
}

# How many times a rendered report states a value.
count_stated <- function(text, value) {
  found <- gregexpr(value, text, fixed = TRUE)[[1]]
  if (found[[1]] == -1L) 0L else length(found)
}

test_that("the report states the run's metadata once, not once per section", {
  local_quiet()
  container <- check_positivity(
    report_data(dgp_practical_violation(n = 317, seed = 2)),
    statin,
    c(bmi, ldl),
    diagnostics = c("edp", "port")
  )
  text <- printed_text(container)

  # The container owns all four values, so the report states them for the run
  # rather than letting every child repeat what it was handed.
  expect_identical(count_stated(text, "statin"), 1L)
  expect_identical(count_stated(text, "binary"), 1L)
  expect_identical(count_stated(text, "317"), 1L)
  expect_identical(count_stated(text, "bmi"), 1L)
  expect_identical(count_stated(text, "ldl"), 1L)
})

test_that("a section with a finding to report leads", {
  local_quiet()
  container <- one_low_support_container()

  # Extrapolation is the only child with a finding here: two of its units sit
  # beyond one geometric variability of an opposite unit. PoRT found no subgroup
  # over its reading rule and EDP declares no finding of any kind, so the order
  # follows what each child reports rather than what its results contain.
  expect_identical(sum(container$extrapolation@results$low_support), 2L)
  expect_identical(sum(container$port@results$low_support), 0L)
  expect_gt(nrow(sniff_violations(container$extrapolation)), 0L)
  expect_identical(nrow(sniff_violations(container$port)), 0L)
  expect_identical(nrow(sniff_violations(container$edp)), 0L)

  positions <- section_positions(container)
  expect_true(all(positions > 0L))
  expect_lt(positions[["extrapolation"]], positions[["edp"]])
  expect_lt(positions[["extrapolation"]], positions[["port"]])
  # The two holding none keep the order they were requested in.
  expect_lt(positions[["edp"]], positions[["port"]])
})

test_that("the requested order does not decide which section leads", {
  local_quiet()
  container <- one_low_support_container()
  reversed <- reorder_checks(container, rev(names(container)))

  requested <- section_positions(container)
  after_reversal <- section_positions(reversed)
  expect_true(all(requested > 0L))
  expect_true(all(after_reversal > 0L))

  # Extrapolation is the only child with a finding to report, so it leads
  # whether it was requested last or first.
  expect_identical(leading_section(after_reversal), leading_section(requested))
  expect_identical(leading_section(requested), "extrapolation")

  # The sections with nothing to report follow the requested order, which the
  # reversal changed.
  expect_lt(after_reversal[["port"]], after_reversal[["edp"]])
})

test_that("sections with findings keep the requested order", {
  local_quiet()
  container <- two_low_support_container()

  # Extrapolation reports far more than PoRT and was requested after it. Having
  # a finding is the whole of what the order records, not how much was found, so
  # between two sections that report the requested order stands.
  expect_gt(
    nrow(sniff_violations(container$extrapolation)),
    0L
  )
  expect_gt(
    sum(container$extrapolation@results$low_support),
    sum(container$port@results$low_support)
  )

  positions <- section_positions(container)
  expect_true(all(positions > 0L))
  expect_lt(positions[["port"]], positions[["extrapolation"]])
  expect_lt(positions[["extrapolation"]], positions[["edp"]])
})

test_that("a diagnostic with nothing to report does not lead", {
  local_quiet()
  # A diagnostic subclassed outside the package reports no findings, so it takes
  # no section ahead of one that does, however its results are shaped and
  # whatever it was requested before.
  container <- one_low_support_container()
  mixed <- positivity_check(
    checks = c(list(plain = make_test_diagnostic()), container@checks),
    exposure = container@exposure,
    exposure_type = container@exposure_type,
    covariates = container@covariates,
    n = container@n,
    call = container@call
  )

  positions <- section_positions(mixed)
  expect_true(all(positions > 0L))
  expect_identical(leading_section(positions), "extrapolation")
  expect_gt(positions[["plain"]], positions[["extrapolation"]])
})

test_that("a section whose finding is a single statistic leads", {
  local_quiet()
  container <- scalar_finding_container()

  # The leverage check reports phi-hat against its own permutation null and
  # declares no low_support column, so a report reading only that column would
  # leave it behind EDP, which reports nothing and was requested first.
  expect_true(container$hat_values@exceeds_null)
  expect_false("low_support" %in% names(container$hat_values@results))
  expect_gt(nrow(sniff_violations(container$hat_values)), 0L)
  expect_identical(nrow(sniff_violations(container$edp)), 0L)

  positions <- section_positions(container)
  expect_true(all(positions > 0L))
  expect_lt(positions[["hat_values"]], positions[["edp"]])
})

test_that("no report line runs past the console width", {
  local_quiet()

  # Every line the report emits is wrapped to the width it is rendered at,
  # including the footer, which cli leaves unwrapped unless asked. The footer
  # names the section that leads, so its length follows that name: the first
  # container leads with `extrapolation`, which overruns the default width on its
  # own, and the second leads with the shorter `port`, which only overruns a
  # narrow one.
  cases <- list(
    list(width = 80, container = one_low_support_container()),
    list(width = 45, container = two_low_support_container())
  )
  for (case in cases) {
    lines <- withr::with_options(
      list(cli.width = case$width, width = case$width),
      printed_lines(case$container)
    )
    expect_lte(max(nchar(lines)), case$width)
  }
})

test_that("the report states the rule behind a count, not its parameters", {
  local_quiet()
  container <- two_low_support_container()
  text <- printed_text(container)

  # PoRT reports a subgroup whose exposure prevalence falls outside
  # [beta, 1 - beta] and that covers at least alpha of the sample. A reader
  # shown alpha, beta, and gamma has to reassemble that rule, so the report
  # states the rule and leaves the parameters on the child.
  expect_no_match(text, "alpha", fixed = TRUE)
  expect_no_match(text, "beta", fixed = TRUE)
  expect_no_match(text, "gamma", fixed = TRUE)
  expect_match(text, "0[.]05")
  expect_match(text, "0[.]95")
})

test_that("the report closes by naming what to do next", {
  local_quiet()
  container <- one_low_support_container()

  footer <- trailing_block(printed_lines(container))

  # The last thing the report says is what the reader can call from here, so it
  # names a call rather than restating a finding.
  next_steps <- c("$", "tidy(", "glance(", "summary(", "plot(")
  named <- vapply(
    next_steps,
    function(step) grepl(step, footer, fixed = TRUE),
    logical(1)
  )
  expect_true(any(named))
})

test_that("the report is stable", {
  local_quiet()
  skip_if_not_installed("lpSolve")
  container <- two_low_support_container()
  expect_snapshot(print(container))
})

test_that("a diagnostic that defines neither display method still has both", {
  # Every section of the report is built from these two generics, so a
  # diagnostic subclassed outside the package falls back rather than fails.
  obj <- make_test_diagnostic()

  label <- diagnostic_label(obj)
  expect_type(label, "character")
  expect_length(label, 1L)
  expect_false(is.na(label))

  headline <- diagnostic_headline(obj)
  expect_type(headline, "character")
  expect_gte(length(headline), 1L)
  expect_false(anyNA(headline))
})
