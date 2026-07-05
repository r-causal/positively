# The validate_* helpers named in the design doc are validate_data_frame(),
# validate_column_selection(), and validate_prob(). They mirror halfmoon's
# validation catalogue: a classed positively_error on failure, the input
# returned invisibly on success.

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
