# plot() calls autoplot() and prints, so both entry points draw the same views.

#' Plot an ETA.Bias diagnostic
#'
#' Draws one of two views of a [check_eta_bias()] result. `type = "bootstrap"`
#' shows the distribution of bootstrap estimates as a histogram with the target
#' `truth` marked, one facet per truncation level. `type = "sweep"` shows
#' ETA.Bias with a two-Monte-Carlo-standard-error band against the truncation
#' lower bound, the bias-variance tradeoff of weight truncation.
#'
#' @param object An `eta_bias_result` from [check_eta_bias()].
#' @param type One of `"bootstrap"` or `"sweep"`.
#' @param ... Not used.
#'
#' @return A [ggplot2::ggplot] object.
#'
#' @examples
#' set.seed(1)
#' n <- 300
#' x1 <- rnorm(n)
#' x2 <- rnorm(n)
#' a <- rbinom(n, 1, plogis(2 * (x1 + x2)))
#' y <- a + x1 + x2 + rnorm(n)
#' df <- data.frame(a = a, y = y, x1 = x1, x2 = x2)
#' # n_boot is small here to keep the example fast.
#' result <- check_eta_bias(df, a, y, c(x1, x2), n_boot = 25)
#' ggplot2::autoplot(result, type = "bootstrap")
#'
#' swept <- check_eta_bias(
#'   df,
#'   a,
#'   y,
#'   c(x1, x2),
#'   truncation_grid = c(0, 0.05, 0.1),
#'   n_boot = 25
#' )
#' ggplot2::autoplot(swept, type = "sweep")
#'
#' @name autoplot.eta_bias_result
#' @aliases plot.eta_bias_result
NULL

method(autoplot, eta_bias_result) <- function(
  object,
  type = c("bootstrap", "sweep"),
  ...
) {
  type <- rlang::arg_match(type)
  if (type == "bootstrap") {
    autoplot_eta_bias_bootstrap(object)
  } else {
    autoplot_eta_bias_sweep(object, call = rlang::current_call())
  }
}

#' The bootstrap-distribution view of an ETA.Bias diagnostic
#'
#' @param object An `eta_bias_result`.
#'
#' @return A [ggplot2::ggplot] object.
#' @keywords internal
#' @noRd
autoplot_eta_bias_bootstrap <- function(object) {
  results <- object@results
  boot_data <- do.call(
    rbind,
    lapply(seq_len(nrow(results)), function(level) {
      tibble::tibble(
        truncation_lower = results$truncation_lower[[level]],
        estimate = object@boot_estimates[[level]]
      )
    })
  )
  # Facet on a factor ordered by the sweep so numeric lower bounds, including
  # ones that print in scientific notation, order by value rather than
  # lexicographically.
  boot_data$facet <- factor(
    boot_data$truncation_lower,
    levels = unique(results$truncation_lower)
  )
  ggplot2::ggplot(boot_data, ggplot2::aes(x = .data$estimate)) +
    ggplot2::geom_histogram(bins = 30, fill = "grey70", color = "white") +
    ggplot2::geom_vline(xintercept = object@truth, linetype = "dashed") +
    ggplot2::facet_wrap(
      ggplot2::vars(.data$facet),
      labeller = ggplot2::as_labeller(function(value) paste0("lower = ", value))
    ) +
    ggplot2::labs(
      x = "Bootstrap estimate",
      y = "Bootstrap draws",
      title = paste0(
        "Bootstrap estimates for the ",
        object@estimator,
        " estimator"
      ),
      subtitle = "Dashed line: target truth"
    )
}

#' The truncation-sweep view of an ETA.Bias diagnostic
#'
#' @param object An `eta_bias_result`.
#' @param call The calling environment, used to build the error's call.
#'
#' @return A [ggplot2::ggplot] object.
#' @keywords internal
#' @noRd
autoplot_eta_bias_sweep <- function(object, call = rlang::caller_env()) {
  results <- object@results
  if (nrow(results) == 1) {
    abort(
      c(
        "The sweep view needs a truncation sweep of more than one level.",
        i = "Rerun {.fn check_eta_bias} with a {.arg truncation_grid} of multiple lower bounds."
      ),
      error_class = "positively_sweep_absent_error",
      call = call
    )
  }
  sweep_data <- tibble::tibble(
    truncation_lower = results$truncation_lower,
    bias = results$bias,
    lower = results$bias - 2 * results$mc_se,
    upper = results$bias + 2 * results$mc_se
  )
  ggplot2::ggplot(
    sweep_data,
    ggplot2::aes(x = .data$truncation_lower, y = .data$bias)
  ) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed") +
    ggplot2::geom_ribbon(
      ggplot2::aes(ymin = .data$lower, ymax = .data$upper),
      fill = "grey80",
      alpha = 0.5
    ) +
    ggplot2::geom_line() +
    ggplot2::geom_point() +
    ggplot2::labs(
      x = "Truncation lower bound",
      y = "ETA.Bias",
      title = "ETA.Bias across the truncation sweep",
      subtitle = "Band: plus or minus two Monte Carlo standard errors"
    )
}

method(plot, eta_bias_result) <- function(x, y, ...) {
  print(autoplot(x, ...))
  invisible(x)
}
