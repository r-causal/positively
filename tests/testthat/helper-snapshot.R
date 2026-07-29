# A wrapper over testthat::expect_snapshot() for expressions that must abort,
# pairing the pinned message with an assertion on the condition's class.
#
# `expect_snapshot(expr, error = TRUE)` pins the rendered message and nothing
# else, so a condition class can drift whenever the wording survives the change
# and the suite will not notice. Asserting the class first closes that gap. An
# error whose class does not match is not caught by `expect_error()`, so it
# propagates: the block ends at the offending line, no later snapshot in that
# block is compared, and the recorded entries are left as they were.
#
# The expression is evaluated twice, once for the class assertion and once for
# the snapshot, so it must not carry side effects that change the second run.
#
# `class` defaults to the package's base condition class, which every wrapper in
# `R/utils.R` attaches. Matching is by inheritance and follows `parent` chains,
# so a foreign error wrapping a package error still matches. A site whose
# subclass is worth pinning on its own passes that subclass instead.
#
# The expression is spliced into the `expect_snapshot()` call while that call is
# being built, because `expect_snapshot()` itself stopped processing `!!` in
# testthat 3.0.3: `expect_snapshot(!!code, error = TRUE)` would record the
# literal `` `!!`(code) `` as the snapshot source.
expect_snapshot_abort <- function(expr, class = "positively_error") {
  quo <- rlang::enquo0(expr)
  code <- rlang::quo_get_expr(quo)
  env <- rlang::quo_get_env(quo)

  eval(rlang::expr(testthat::expect_error(!!code, class = !!class)), env)
  eval(rlang::expr(testthat::expect_snapshot(!!code, error = TRUE)), env)
}
