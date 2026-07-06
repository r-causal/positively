local_quiet <- function(.env = parent.frame()) {
  withr::local_options(positively.quiet = TRUE, .local_envir = .env)
}
