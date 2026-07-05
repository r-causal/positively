test_that("the package namespace loads and attaches", {
  expect_true("positively" %in% loadedNamespaces())
  expect_true("package:positively" %in% search())
})

test_that(".onLoad registers S7 methods without error", {
  # .onLoad() calls S7::methods_register(), which is idempotent, so invoking
  # it again on the loaded namespace is a genuine check that it runs cleanly.
  expect_no_error(
    .onLoad(
      libname = dirname(find.package("positively")),
      pkgname = "positively"
    )
  )
})
