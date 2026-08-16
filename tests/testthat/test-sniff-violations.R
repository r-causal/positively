# sniff_violations() reports what a run found, one row per finding, in a schema
# every diagnostic shares. Participation is declared rather than inferred: a
# class contributes rows by defining a method, and the inherited default reports
# none whatever its results happen to hold.

# The columns every result carries, in order.
sniff_columns <- function() {
  c("diagnostic", "scope", "label", "n", "statistic", "value", "threshold")
}

# Assertions every result satisfies, whichever object or scope produced it. The
# value and the threshold are the two the container used to coerce to character,
# so they are pinned as numeric wherever they appear.
expect_sniff_schema <- function(found) {
  testthat::expect_s3_class(found, "tbl_df")
  testthat::expect_identical(names(found), sniff_columns())
  testthat::expect_type(found$diagnostic, "character")
  testthat::expect_type(found$scope, "character")
  testthat::expect_type(found$label, "character")
  # A count column counts, and every method agrees on that so the rows of one
  # container are not coerced to a wider type by a single method's choice.
  testthat::expect_type(found$n, "integer")
  testthat::expect_type(found$statistic, "character")
  testthat::expect_type(found$value, "double")
  testthat::expect_type(found$threshold, "double")
  invisible(found)
}

# ---- Local fixtures --------------------------------------------------------

# A throwaway concrete subclass carrying every column the package reads a flag
# from anywhere. The default method reports nothing from any of them, so a
# results tibble that a column-reading implementation would take for a finding
# is what exercises it.
make_unclaimed_diagnostic <- function() {
  unclaimed_diagnostic <- S7::new_class(
    "unclaimed_diagnostic",
    parent = positivity_diagnostic
  )
  unclaimed_diagnostic(
    results = tibble::tibble(
      .id = 1:3,
      low_support = c(TRUE, TRUE, FALSE),
      in_hull = c(FALSE, FALSE, TRUE),
      high_leverage = c(TRUE, FALSE, TRUE),
      prop_zero = c(0.4, 0.2, 0)
    ),
    exposure = "a",
    exposure_type = "binary",
    n = 3L,
    params = list(),
    call = quote(check_unclaimed())
  )
}

# A steep propensity model puts the tails of x1 near deterministic treatment. At
# n = 317 the reading rule reports 45 subgroups, 4 of which meet it at the
# default prevalence threshold.
port_finding <- function(...) {
  check_port(
    dgp_practical_violation(n = 317, seed = 2),
    exposure,
    c(x1, x2),
    ...
  )
}

# Good positivity leaves every reported subgroup inside the reading rule.
port_no_finding <- function() {
  check_port(dgp_good_positivity(n = 400, seed = 1), exposure, c(x1, x2))
}

# A sequential run with censoring under the Gruber bound. The censoring risk set
# is not the follower set, so the two sides resolve different per-wave
# thresholds and a row carrying the wrong one is visible.
port_censoring_finding <- function() {
  check_port_seq(
    dgp_longitudinal_binary_censoring(n = 400, seed = 7),
    c(a1, a2, a3),
    list(l0, l1, l2),
    .censoring = c(c1, c2, c3),
    beta = "gruber"
  )
}

# Two of 400 units sit beyond one geometric variability of the nearest opposite
# unit while 19 fall outside the opposite hull, so both readings fire and they
# disagree about which units they name.
extrapolation_finding <- function(hull = TRUE) {
  check_extrapolation(
    dgp_good_positivity(n = 400, seed = 1),
    exposure,
    c(x1, x2),
    hull = hull
  )
}

# Every one of the 200 units has an opposite-exposure unit within one geometric
# variability, so the Gower reading is silent, while 8 of them fall outside the
# opposite hull. Under `hull = FALSE` the hull reading has nothing to say
# either, which is the run on which the diagnostic found nothing at all.
extrapolation_hull_only <- function(hull = TRUE) {
  check_extrapolation(
    dgp_structural_subgroup(n = 200, seed = 3),
    exposure,
    c(x1, x2),
    hull = hull
  )
}

# The null draw is seeded because the quantile phi-hat is compared against is
# drawn inside the run.
hat_values_finding <- function() {
  withr::local_seed(2024)
  check_hat_values(
    dgp_continuous_support_gap(n = 200, seed = 4),
    exposure,
    x1,
    null_reps = 25,
    exposure_type = "continuous"
  )
}

# Constructed directly rather than fitted: a leverage profile that stays inside
# its null needs an exposure drawn independently of the covariates, which the
# shared generators do not provide. The fields agree with each other because
# phi-hat is the share of candidates over the cutoff, so 2 of 100 carrying high
# leverage is a phi-hat of 0.02, below the 0.05 null quantile.
make_quiet_hat_values <- function() {
  hat_values_result(
    results = tibble::tibble(
      .id = rep(1:20, each = 5),
      prob = rep(seq(0.1, 0.9, length.out = 5), times = 20),
      value = rep(seq(-2, 2, length.out = 5), times = 20),
      hat_value = c(rep(0.01, 98), 0.4, 0.5),
      high_leverage = c(rep(FALSE, 98), TRUE, TRUE)
    ),
    exposure = "dose",
    exposure_type = "continuous",
    n = 20L,
    params = list(
      probs = seq(0.1, 0.9, length.out = 5),
      threshold = 2,
      null_reps = 25,
      null_method = "permutation",
      conf_level = 0.95
    ),
    call = quote(check_hat_values()),
    phi_hat = 0.02,
    null_dist = rep(0.05, 25),
    null_quantile = 0.05,
    exceeds_null = FALSE,
    p = 3L
  )
}

# Constructed directly rather than fitted: ETA bias needs an outcome column, and
# neither the shared generators nor the shipped datasets carry one.
make_eta_bias_result <- function() {
  eta_bias_result(
    results = tibble::tibble(
      term = "1 - 0",
      truncation_lower = 0,
      truncation_upper = 1,
      bias = 0.25,
      mc_se = 0.02,
      boot_mean = 0.75
    ),
    exposure = "a",
    exposure_type = "binary",
    n = 300L,
    params = list(n_boot = 50),
    call = quote(check_eta_bias()),
    estimator = "ipw",
    truth = c("1 - 0" = 0.5),
    boot_estimates = list(c(0.7, 0.8))
  )
}

edp_no_finding <- function() {
  check_edp(
    dgp_continuous_support_gap(n = 200, seed = 4),
    exposure,
    x1,
    values = c(0, 1, 3),
    exposure_type = "continuous"
  )
}

hdr_no_finding <- function() {
  check_hdr(
    dgp_continuous_support_gap(n = 200, seed = 4),
    exposure,
    x1,
    values = c(0, 1, 3)
  )
}

# The binary default set, so the container holds a diagnostic contributing
# subgroup rows, one contributing overall rows, and one contributing none.
binary_container <- function() {
  check_positivity(
    dgp_practical_violation(n = 317, seed = 2),
    exposure,
    c(x1, x2)
  )
}

# ---- The schema ------------------------------------------------------------

test_that("sniff_violations() is part of the public surface", {
  expect_true("sniff_violations" %in% getNamespaceExports("positively"))
})

test_that("the row schema is fixed and its numbers stay numbers", {
  local_quiet()
  found <- sniff_violations(binary_container())

  expect_sniff_schema(found)
  expect_gt(nrow(found), 0L)
})

test_that("every scope is one of the three kinds", {
  local_quiet()
  found <- sniff_violations(
    binary_container(),
    scope = c("subgroup", "overall", "unit")
  )

  expect_true(all(found$scope %in% c("subgroup", "overall", "unit")))
})

# ---- Scope -----------------------------------------------------------------

test_that("unit rows are never in the default", {
  local_quiet()
  container <- binary_container()

  # A subgroup with low support is a covariate region a reader can act on, while
  # a unit with low support is one row of a distribution whose aggregate is the
  # statistic, so the two are not stacked into one answer by default.
  default <- sniff_violations(container)
  units <- sniff_violations(container, scope = "unit")

  expect_false("unit" %in% default$scope)
  expect_setequal(units$scope, "unit")
  expect_gt(nrow(units), nrow(default))
})

test_that("asking for every scope returns the union", {
  local_quiet()
  container <- binary_container()

  default <- sniff_violations(container)
  units <- sniff_violations(container, scope = "unit")
  everything <- sniff_violations(
    container,
    scope = c("subgroup", "overall", "unit")
  )

  expect_identical(nrow(everything), nrow(default) + nrow(units))
  expect_setequal(everything$scope, union(default$scope, units$scope))
})

test_that("asking for one scope on its own returns that scope", {
  local_quiet()
  # A single-element request is its own path through the argument, so it is
  # asked for directly rather than only as part of a pair.
  container <- binary_container()
  found <- sniff_violations(container, scope = "overall")

  expect_sniff_schema(found)
  expect_setequal(found$scope, "overall")
  expect_setequal(found$diagnostic, "extrapolation")
})

test_that("asking for a scope no diagnostic supplies reports no rows", {
  local_quiet()
  # The leverage check reports at overall and unit scope and has no subgroups to
  # report, which is an empty answer rather than a failure.
  found <- sniff_violations(hat_values_finding(), scope = "subgroup")

  expect_sniff_schema(found)
  expect_identical(nrow(found), 0L)
})

test_that("an unrecognized scope is rejected", {
  local_quiet()
  res <- port_finding()
  expect_snapshot_abort(
    sniff_violations(res, scope = "region"),
    class = "positively_args_error"
  )
})

# ---- Declared participation ------------------------------------------------

test_that("the default reports nothing and infers nothing from a flag column", {
  # The default never reads low_support, in_hull, high_leverage, or a
  # proportion. A diagnostic joins by declaring what it reports, so one that
  # declares nothing reports nothing even where a column would have been read.
  found <- sniff_violations(make_unclaimed_diagnostic())

  expect_sniff_schema(found)
  expect_identical(nrow(found), 0L)
})

test_that("an empty result keeps its columns and can be filtered", {
  found <- sniff_violations(make_unclaimed_diagnostic())

  # An empty answer is the common case, so it behaves like a populated one
  # rather than needing a guard at every call site.
  filtered <- dplyr::filter(found, value > 0)
  expect_sniff_schema(filtered)
  expect_identical(nrow(filtered), 0L)
})

test_that("EDP reports nothing", {
  local_quiet()
  res <- edp_no_finding()
  found <- sniff_violations(res)

  # The run measured something at every intervention value, so the silence is
  # the class declining to report rather than the fixture having nothing in it.
  expect_gt(nrow(res@results), 0L)
  expect_true(all(res@results$edp > 0))

  expect_sniff_schema(found)
  expect_identical(nrow(found), 0L)
})

test_that("HDR non-overlap reports nothing", {
  local_quiet()
  res <- hdr_no_finding()
  found <- sniff_violations(res)

  # Non-overlap runs above zero at some target here, which a column-reading
  # implementation could take for a finding, so the silence is declared.
  expect_gt(nrow(res@results), 0L)
  expect_gt(max(res@results$nonoverlap), 0)

  expect_sniff_schema(found)
  expect_identical(nrow(found), 0L)
})

test_that("ETA bias reports nothing", {
  found <- sniff_violations(make_eta_bias_result())

  expect_sniff_schema(found)
  expect_identical(nrow(found), 0L)
})

test_that("density ratios report nothing", {
  found <- sniff_violations(check_density_ratios(c(0, 0, 1, 2, 5)))

  expect_sniff_schema(found)
  expect_identical(nrow(found), 0L)
})

# ---- PoRT ------------------------------------------------------------------

test_that("PoRT contributes one row per low-support subgroup", {
  local_quiet()
  res <- port_finding()
  low_support <- res@results[res@results$low_support, ]
  found <- sniff_violations(res)

  expect_sniff_schema(found)
  expect_identical(nrow(found), 4L)
  expect_identical(nrow(found), nrow(low_support))
  expect_setequal(found$scope, "subgroup")
  expect_setequal(found$statistic, "prevalence")

  # A result read straight off one diagnostic names that diagnostic, which is
  # the same name a container would key it under, rather than waiting for a
  # container to fill the column in.
  expect_setequal(found$diagnostic, "port")

  # Each row carries the subgroup's own reading, so a reader can go from a row
  # back to the covariate region it names.
  expect_identical(found$label, low_support$description)
  expect_identical(found$n, low_support$n)
  expect_equal(found$value, low_support$prevalence)
})

test_that("PoRT reports nothing when no subgroup meets the reading rule", {
  local_quiet()
  res <- port_no_finding()

  # Subgroups were reported; none of them met the rule, which is the difference
  # between a diagnostic that found nothing and one that ran nothing.
  expect_gt(nrow(res@results), 0L)
  expect_identical(sum(res@results$low_support), 0L)
  expect_identical(nrow(sniff_violations(res)), 0L)
})

test_that("PoRT contributes no unit rows", {
  local_quiet()
  found <- sniff_violations(port_finding(), scope = "unit")

  expect_sniff_schema(found)
  expect_identical(nrow(found), 0L)
})

test_that("the threshold is the object's own resolved cut", {
  local_quiet()
  # A prevalence cut the caller moved off the default, so a hardcoded 0.05 would
  # show.
  res <- port_finding(beta = 0.1)
  found <- sniff_violations(res)

  expect_identical(res@beta, 0.1)
  expect_setequal(found$threshold, 0.1)
})

test_that("the threshold agrees with statistic_thresholds()", {
  local_quiet()
  res <- port_finding()
  found <- sniff_violations(res)

  # The package now states a cut in two places, and a reader comparing a row
  # against the summary must not be told two different numbers.
  expect_setequal(
    found$threshold,
    unname(statistic_thresholds(res)[["n_low_support"]])
  )
})

test_that("a sequential run contributes its censoring subgroups", {
  local_quiet()
  res <- port_censoring_finding()
  low_support <- res@results[res@results$low_support, ]
  is_censoring <- low_support$type == "censoring"
  found <- sniff_violations(res)

  # A censoring subgroup with low support is a finding in the same sense as an
  # exposure one, so both sides are reported.
  expect_gt(sum(is_censoring), 0L)
  expect_gt(sum(!is_censoring), 0L)
  expect_identical(nrow(found), nrow(low_support))

  # Each side carries the cut its own risk set resolved. The two vectors differ
  # here, so a row carrying the other side's threshold would be visible. The
  # comparison is order-free because two sources merge into one answer.
  expect_false(identical(res@beta, res@censoring_beta))
  expected <- c(
    res@beta[low_support$time[!is_censoring]],
    res@censoring_beta[low_support$time[is_censoring]]
  )
  expect_equal(sort(found$threshold), sort(expected))
})

# ---- Hat values ------------------------------------------------------------

test_that("hat values contribute one overall row when the null is exceeded", {
  local_quiet()
  res <- hat_values_finding()
  found <- sniff_violations(res)

  expect_true(res@exceeds_null)
  expect_sniff_schema(found)
  expect_identical(nrow(found), 1L)
  expect_identical(found$diagnostic, "hat_values")
  expect_identical(found$scope, "overall")
  expect_identical(found$statistic, "phi_hat")
  expect_equal(found$value, res@phi_hat)

  # The comparison is against the run's own permutation null, which is also what
  # the summary reports as this statistic's cut.
  expect_equal(found$threshold, res@null_quantile)
  expect_equal(
    found$threshold,
    unname(statistic_thresholds(res)[["phi_hat"]])
  )

  # The reading is a property of the whole run rather than of a set of rows.
  expect_identical(found$n, NA_integer_)
})

test_that("hat values contribute no overall row inside the null", {
  res <- make_quiet_hat_values()
  found <- sniff_violations(res)

  expect_false(res@exceeds_null)
  expect_sniff_schema(found)
  expect_identical(nrow(found), 0L)
})

test_that("hat values contribute their high-leverage units at unit scope", {
  local_quiet()
  res <- hat_values_finding()
  found <- sniff_violations(res, scope = "unit")

  # The column is a conservative indicator feeding phi-hat rather than a
  # calibrated per-unit rule, so the rows are labeled for what they are.
  expect_sniff_schema(found)
  expect_identical(nrow(found), sum(res@results$high_leverage))
  expect_setequal(found$diagnostic, "hat_values")
  expect_setequal(found$scope, "unit")
  expect_setequal(found$label, "high_leverage")

  # A unit's cut is stated in hat-value units, which is what the summary pairs
  # with the count of those units.
  expect_equal(
    unique(found$threshold),
    unname(statistic_thresholds(res)[["n_high_leverage"]])
  )
})

# ---- Extrapolation ---------------------------------------------------------

test_that("extrapolation contributes both of its readings", {
  local_quiet()
  res <- extrapolation_finding()
  found <- sniff_violations(res)

  # The two readings disagree about which units they name, and the contrast is
  # the informative part, so neither stands in for the other.
  expect_sniff_schema(found)
  expect_identical(nrow(found), 2L)
  expect_setequal(found$diagnostic, "extrapolation")
  expect_setequal(found$scope, "overall")
  expect_setequal(found$statistic, c("prop_supported", "prop_in_hull"))
  expect_false(identical(found$label[[1]], found$label[[2]]))

  supported <- found[found$statistic == "prop_supported", ]
  hull <- found[found$statistic == "prop_in_hull", ]
  expect_equal(supported$value, mean(!res@results$low_support))
  expect_equal(hull$value, mean(res@results$in_hull))
  expect_identical(supported$n, sum(res@results$low_support))
  expect_identical(hull$n, sum(!res@results$in_hull))

  # prop_supported counts the units whose nearest opposite unit lies further
  # than one geometric variability away, so the cut behind it is that distance,
  # and it is the same number the summary pairs with the same statistic.
  expect_identical(supported$threshold, res@geometric_variability)
  expect_identical(
    supported$threshold,
    unname(statistic_thresholds(res)[["prop_supported"]])
  )

  # Hull membership is a containment test with no numeric cut behind it.
  expect_identical(hull$threshold, NA_real_)
})

test_that("extrapolation omits the hull reading when the hull did not run", {
  local_quiet()
  res <- extrapolation_finding(hull = FALSE)
  found <- sniff_violations(res)

  # in_hull is missing throughout when the test is skipped, so the reading is
  # absent rather than reported as a row of missing values.
  expect_false(res@hull_run)
  expect_true(all(is.na(res@results$in_hull)))
  expect_identical(nrow(found), 1L)
  expect_identical(found$statistic, "prop_supported")
  expect_false(anyNA(found$value))
})

test_that("extrapolation reports only the reading that fires", {
  local_quiet()
  res <- extrapolation_hull_only()
  found <- sniff_violations(res)

  # Every unit has an opposite-exposure unit within one geometric variability,
  # so that reading found nothing and says nothing, while the hull reading found
  # 8 units and says so.
  expect_identical(sum(res@results$low_support), 0L)
  expect_gt(sum(!res@results$in_hull), 0L)

  expect_identical(nrow(found), 1L)
  expect_identical(found$statistic, "prop_in_hull")
  expect_identical(found$n, sum(!res@results$in_hull))
  expect_identical(found$threshold, NA_real_)
})

test_that("extrapolation reports nothing when neither reading fires", {
  local_quiet()
  res <- extrapolation_hull_only(hull = FALSE)
  found <- sniff_violations(res)

  # A proportion is a reading rather than a finding, so a run where every unit
  # is supported and the hull test did not run reports no rows at all. Reporting
  # the proportion anyway would make the accessor non-empty for a dataset with
  # nothing wrong with it.
  expect_identical(sum(res@results$low_support), 0L)
  expect_false(res@hull_run)

  expect_sniff_schema(found)
  expect_identical(nrow(found), 0L)
})

test_that("a unit contributes a row only for the reading that fires", {
  local_quiet()
  res <- extrapolation_finding()
  found <- sniff_violations(res, scope = "unit")

  # A row per unit per reading would double the answer while saying nothing
  # about the units that fired neither reading.
  expect_sniff_schema(found)
  expect_identical(
    nrow(found),
    sum(res@results$low_support) + sum(!res@results$in_hull)
  )
  expect_setequal(found$scope, "unit")
  expect_setequal(found$diagnostic, "extrapolation")

  gower_rows <- found[found$statistic == "gower_min", ]
  hull_rows <- found[found$statistic == "in_hull", ]
  expect_identical(nrow(gower_rows), sum(res@results$low_support))
  expect_identical(nrow(hull_rows), sum(!res@results$in_hull))

  # A Gower row states the unit's own distance to the nearest opposite unit.
  # Hull membership has no per-unit quantity to state, and zero would be the
  # wrong stand-in: a Gower distance of zero in this same column is a unit
  # sitting on top of an opposite-exposure unit, which is the opposite reading.
  expect_equal(
    sort(gower_rows$value),
    sort(res@results$gower_min[res@results$low_support])
  )
  expect_true(all(is.na(hull_rows$value)))

  # A Gower row states the distance it was compared against; a hull row has no
  # numeric cut to state.
  expect_setequal(gower_rows$threshold, res@geometric_variability)
  expect_true(all(is.na(hull_rows$threshold)))
})

test_that("with no hull test a unit contributes only the Gower reading", {
  local_quiet()
  res <- extrapolation_finding(hull = FALSE)
  found <- sniff_violations(res, scope = "unit")

  expect_sniff_schema(found)
  expect_identical(nrow(found), sum(res@results$low_support))
  expect_setequal(found$scope, "unit")
})

# ---- The container ---------------------------------------------------------

test_that("the container names each row for the child it came from", {
  local_quiet()
  container <- binary_container()
  found <- sniff_violations(container)

  expect_sniff_schema(found)
  expect_true(all(found$diagnostic %in% names(container)))

  # The container's rows are its children's rows, so the two cannot drift.
  for (name in names(container)) {
    expect_identical(
      nrow(found[found$diagnostic == name, ]),
      nrow(sniff_violations(container[[name]]))
    )
  }
})

test_that("the container names a child by the name it holds it under", {
  local_quiet()
  # The name comes from the container rather than from the child's own idea of
  # what it is, so a child stored under another name reports under that name.
  child <- port_finding()
  container <- positivity_check(
    checks = list(subgroups = child),
    exposure = "exposure",
    exposure_type = "binary",
    covariates = c("x1", "x2"),
    n = 317L,
    call = quote(check_positivity())
  )
  found <- sniff_violations(container)

  expect_gt(nrow(found), 0L)
  expect_setequal(found$diagnostic, "subgroups")
})

test_that("a container with no children reports no rows", {
  container <- positivity_check(
    checks = list(),
    exposure = "exposure",
    exposure_type = "binary",
    covariates = c("x1", "x2"),
    n = 317L,
    call = quote(check_positivity())
  )
  found <- sniff_violations(container)

  expect_sniff_schema(found)
  expect_identical(nrow(found), 0L)
})

# ---- What the report's ordering follows -------------------------------------

# A container assembled from chosen children, so a block can fix which shapes it
# holds and the order they were requested in. Only the section order is read
# from it, so its own metadata stands in rather than describing any one child.
make_ordering_container <- function(checks) {
  positivity_check(
    checks = checks,
    exposure = "exposure",
    exposure_type = "binary",
    covariates = c("x1", "x2"),
    n = 317L,
    call = quote(check_positivity())
  )
}

# The names of a container's sections, in the order the report printed them.
printed_section_order <- function(container) {
  positions <- section_positions(container)
  testthat::expect_true(all(positions > 0L))
  names(sort(positions))
}

test_that("every shape that reports rows leads every shape that reports none", {
  local_quiet()
  # One child of every result class, plus one subclassed outside the package,
  # requested in an order that interleaves the two groups. PoRT, the leverage
  # check, and extrapolation each have something to report; the rest do not.
  container <- make_ordering_container(list(
    edp = edp_no_finding(),
    port = port_finding(),
    hat_values = hat_values_finding(),
    hdr = hdr_no_finding(),
    extrapolation = extrapolation_finding(),
    eta_bias = make_eta_bias_result(),
    density_ratios = check_density_ratios(c(0, 0, 1, 2, 5)),
    unclaimed = make_unclaimed_diagnostic()
  ))

  requested <- names(container)
  reports <- vapply(
    requested,
    function(name) nrow(sniff_violations(container[[name]])) > 0L,
    logical(1)
  )
  expect_setequal(requested[reports], c("port", "hat_values", "extrapolation"))

  # Sections that report come first and the rest follow, and within each group
  # the order is the one the diagnostics were requested in.
  expect_identical(
    printed_section_order(container),
    c(requested[reports], requested[!reports])
  )
})

test_that("a foreign subclass carrying a support column reports nothing and does not lead", {
  local_quiet()
  # The default reads no column, so a diagnostic defined outside the package
  # neither contributes a row nor takes a section ahead of one that does, even
  # holding a results tibble that a column-reading implementation would flag and
  # even having been requested first.
  unclaimed <- make_unclaimed_diagnostic()
  container <- make_ordering_container(list(
    unclaimed = unclaimed,
    port = port_finding()
  ))

  expect_identical(nrow(sniff_violations(unclaimed)), 0L)
  expect_identical(printed_section_order(container), c("port", "unclaimed"))
})

test_that("a run with nothing to report keeps the order it was requested in", {
  local_quiet()
  # Extrapolation found nothing here, so it stays where it was asked for rather
  # than taking a section on the strength of having a reading to state.
  container <- make_ordering_container(list(
    edp = edp_no_finding(),
    extrapolation = extrapolation_hull_only(hull = FALSE)
  ))

  expect_identical(printed_section_order(container), c("edp", "extrapolation"))
})
