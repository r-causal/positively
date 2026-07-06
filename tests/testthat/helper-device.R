# Route base-graphics output to a throwaway device so tests that print or draw a
# plot leave no Rplots.pdf in the package root. The device closes when the
# calling frame exits.
local_null_device <- function(.env = parent.frame()) {
  withr::local_pdf(nullfile(), .local_envir = .env)
}
