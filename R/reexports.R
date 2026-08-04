# Re-export the broom-ecosystem generics so that tidy() and glance() are
# available with `library(positively)` alone, without qualifying them as
# `generics::tidy()`. positively registers S7 methods on both generics for its
# result classes and the container.

#' @importFrom generics tidy
#' @export
generics::tidy

#' @importFrom generics glance
#' @export
generics::glance

#' @importFrom ggplot2 autoplot
#' @export
ggplot2::autoplot
