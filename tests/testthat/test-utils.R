test_that("abort() signals a classed positively_error condition", {
  expect_error(
    abort("Something went wrong."),
    class = "positively_error"
  )
})

test_that("abort() layers a caller-supplied subclass over positively_error", {
  cnd <- rlang::catch_cnd(
    abort("Bad input.", error_class = "positively_demo_error")
  )
  expect_s3_class(cnd, "positively_demo_error")
  expect_s3_class(cnd, "positively_error")
})

test_that("warn() signals a classed positively_warning condition", {
  expect_warning(
    warn("Heads up."),
    class = "positively_warning"
  )
})

test_that("warn() layers a caller-supplied subclass over positively_warning", {
  cnd <- rlang::catch_cnd(
    warn("Careful.", warning_class = "positively_demo_warning")
  )
  expect_s3_class(cnd, "positively_demo_warning")
  expect_s3_class(cnd, "positively_warning")
})

test_that("alert_info() emits an informational message by default", {
  withr::local_options(positively.quiet = FALSE)
  expect_message(alert_info("informational text"), "informational text")
})

test_that("alert_info() is suppressed when positively.quiet is TRUE", {
  withr::local_options(positively.quiet = TRUE)
  expect_silent(alert_info("should not be shown"))
})
