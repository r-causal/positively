# plot() calls autoplot() and prints, so both entry points draw the same views.

#' Plot a hat-value diagnostic
#'
#' Draws one of two views of a [check_hat_values()] result. `type = "null"`
#' shows the null distribution of \eqn{\hat{\phi}} as a histogram with the
#' observed \eqn{\hat{\phi}} and the `conf_level` null quantile marked.
#' `type = "profile"` shows the proportion of high-leverage candidates at each
#' exposure percentile, the leverage profile across the exposure range.
#'
#' @param object A `hat_values_result` from [check_hat_values()].
#' @param type One of `"null"` or `"profile"`.
#' @param ... Not used.
#'
#' @return A [ggplot2::ggplot] object.
#'
#' @examples
#' set.seed(1)
#' n <- 200
#' x1 <- rnorm(n)
#' dose <- rnorm(n, mean = x1)
#' df <- data.frame(dose = dose, x1 = x1)
#' # null_reps is small here to keep the example fast.
#' result <- check_hat_values(df, dose, x1, null_reps = 50)
#'
#' ggplot2::autoplot(result, type = "null")
#' ggplot2::autoplot(result, type = "profile")
#'
#' @name autoplot.hat_values_result
#' @aliases plot.hat_values_result
NULL

method(autoplot, hat_values_result) <- function(
  object,
  type = c("null", "profile"),
  ...
) {
  type <- rlang::arg_match(type)
  if (type == "null") {
    autoplot_hat_values_null(object)
  } else {
    autoplot_hat_values_profile(object)
  }
}

#' The null-distribution view of a hat-value diagnostic
#'
#' @param object A `hat_values_result`.
#'
#' @return A [ggplot2::ggplot] object.
#' @keywords internal
#' @noRd
autoplot_hat_values_null <- function(object) {
  null_data <- tibble::tibble(phi = object@null_dist)
  conf_level <- object@params$conf_level
  # Fewer bins for small null_reps keeps the histogram from reading as ragged;
  # capped at 30 so large null distributions are not oversmoothed.
  bins <- min(30L, max(1L, ceiling(sqrt(length(object@null_dist)))))
  ggplot2::ggplot(null_data, ggplot2::aes(x = .data$phi)) +
    ggplot2::geom_histogram(bins = bins, fill = "grey70", color = "white") +
    ggplot2::geom_vline(
      xintercept = object@null_quantile,
      linetype = "dashed"
    ) +
    ggplot2::geom_vline(
      xintercept = object@phi_hat,
      color = "firebrick",
      linewidth = 1
    ) +
    ggplot2::labs(
      x = expression(hat(phi)),
      y = "Null replicates",
      title = "Observed leverage summary against its null",
      subtitle = paste0(
        "Dashed line: ",
        conf_level,
        " null quantile. Solid line: observed."
      )
    )
}

#' The leverage-profile view of a hat-value diagnostic
#'
#' @param object A `hat_values_result`.
#'
#' @return A [ggplot2::ggplot] object.
#' @keywords internal
#' @noRd
autoplot_hat_values_profile <- function(object) {
  results <- object@results
  fraction <- tapply(results$high_leverage, results$prob, mean)
  profile <- tibble::tibble(
    prob = as.double(names(fraction)),
    fraction = as.double(fraction)
  )
  ggplot2::ggplot(profile, ggplot2::aes(x = .data$prob, y = .data$fraction)) +
    ggplot2::geom_line() +
    ggplot2::geom_point() +
    ggplot2::labs(
      x = "Exposure percentile",
      y = "Proportion of high-leverage candidates",
      title = "Leverage profile across the exposure range"
    )
}

method(plot, hat_values_result) <- function(x, y, ...) {
  print(autoplot(x, ...))
  invisible(x)
}
