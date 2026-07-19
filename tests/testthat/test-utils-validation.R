# The validate_* helpers share one contract, mirroring halfmoon's validation
# catalogue: a classed positively_error on failure, the input returned
# invisibly on success.

test_that("validate_data_frame() accepts a data frame", {
  data <- tibble::tibble(a = 1:3)
  expect_identical(validate_data_frame(data), data)
})

test_that("validate_data_frame() rejects non-data-frame input", {
  expect_error(
    validate_data_frame(1:10),
    class = "positively_error"
  )
})

test_that("validate_column_selection() returns a non-empty selection", {
  selection <- c(x1 = 1L, x2 = 2L)
  expect_identical(validate_column_selection(selection), selection)
})

test_that("validate_column_selection() rejects an empty selection", {
  expect_error(
    validate_column_selection(integer(0)),
    class = "positively_error"
  )
})

test_that("validate_prob() accepts values in the unit interval", {
  expect_identical(validate_prob(0.05), 0.05)
  expect_identical(validate_prob(c(0, 0.5, 1)), c(0, 0.5, 1))
})

test_that("validate_prob() rejects values outside [0, 1]", {
  expect_error(validate_prob(1.5), class = "positively_error")
  expect_error(validate_prob(-0.1), class = "positively_error")
})

test_that("validate_prob() rejects non-numeric input", {
  expect_error(validate_prob("a"), class = "positively_error")
})

test_that("validate_* failures produce stable messages", {
  expect_snapshot(validate_data_frame(1:10), error = TRUE)
  expect_snapshot(validate_column_selection(integer(0)), error = TRUE)
  expect_snapshot(validate_prob(1.5), error = TRUE)
  expect_snapshot(validate_prob("a"), error = TRUE)
})

# A misspelled column must surface as a classed positively_selection_error from
# every selecting argument, so the package speaks with one voice instead of
# leaking a raw vctrs subscript condition. check_eta_bias() already classes its
# failures and stands as the reference case.

test_that("an unknown column raises a classed selection error from every point-treatment check", {
  local_quiet()
  df <- dgp_good_positivity(n = 80, seed = 1)
  continuous <- dgp_continuous_support_gap(n = 80, seed = 1)

  selections <- list(
    edp_exposure = function() check_edp(df, nope, x1),
    edp_covariate = function() check_edp(df, exposure, nope),
    port_exposure = function() check_port(df, nope, x1),
    port_covariate = function() check_port(df, exposure, nope),
    hat_values_exposure = function() check_hat_values(continuous, nope, x1),
    hat_values_covariate = function() {
      check_hat_values(continuous, exposure, nope)
    },
    hdr_exposure = function() check_hdr(continuous, nope, x1),
    hdr_covariate = function() check_hdr(continuous, exposure, nope),
    extrapolation_exposure = function() check_extrapolation(df, nope, x1),
    extrapolation_covariate = function() {
      check_extrapolation(df, exposure, nope)
    },
    positivity_exposure = function() check_positivity(df, nope, x1),
    positivity_covariate = function() check_positivity(df, exposure, nope),
    # The reference case: check_eta_bias() already classes its selection failure.
    eta_bias_exposure = function() {
      check_eta_bias(df, nope, ps, c(x1, x2), n_boot = 5)
    }
  )

  for (nm in names(selections)) {
    expect_error(
      selections[[nm]](),
      class = "positively_selection_error",
      info = nm
    )
  }
})

test_that("an unknown column raises a classed selection error from every sequential selecting argument", {
  local_quiet()
  binary <- dgp_longitudinal_binary(n = 200, seed = 6)
  continuous <- dgp_longitudinal(n = 100, seed = 5)
  censoring <- dgp_longitudinal_binary_censoring(n = 200, seed = 7)

  selections <- list(
    port_seq_exposure = function() {
      check_port_seq(binary, c(nope, a2, a3), list(l0, l1, l2))
    },
    port_seq_covariate = function() {
      check_port_seq(binary, c(a1, a2, a3), list(nope, l1, l2))
    },
    port_seq_baseline = function() {
      check_port_seq(binary, c(a1, a2, a3), list(l0, l1, l2), .baseline = nope)
    },
    port_seq_censoring = function() {
      check_port_seq(
        censoring,
        c(a1, a2, a3),
        list(l0, l1, l2),
        .censoring = c(c1, nope, c3)
      )
    },
    hdr_seq_exposure = function() {
      check_hdr_seq(continuous, c(nope, a2), list(l0, l1), values = 0)
    },
    hdr_seq_covariate = function() {
      check_hdr_seq(continuous, c(a1, a2), list(nope, l1), values = 0)
    }
  )

  for (nm in names(selections)) {
    expect_error(
      selections[[nm]](),
      class = "positively_selection_error",
      info = nm
    )
  }
})
