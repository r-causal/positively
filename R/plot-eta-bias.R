# plot() calls autoplot() and prints, so both entry points draw the same views.

#' Plot an ETA.Bias diagnostic
#'
#' Draws one of two views of a [check_eta_bias()] result. `type = "bootstrap"`
#' shows the distribution of bootstrap estimates as a histogram, one facet per
#' estimand term and truncation level, each facet marking the `truth` its term
#' is aimed at. `type = "sweep"` shows ETA.Bias with a
#' two-Monte-Carlo-standard-error band across the truncation sweep, the
#' bias-variance tradeoff of weight truncation, with one line per estimand term.
#'
#' An estimand of one term, read at one truncation level, draws a single facet.
#' A discrete exposure truncates a fitted probability, so its sweep is drawn
#' against the truncation lower bound. A continuous exposure caps a stabilized
#' weight instead and reports no probability bound at any level, so its sweep is
#' drawn against the quantile levels those caps were read at.
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
#' autoplot(result, type = "bootstrap")
#'
#' swept <- check_eta_bias(
#'   df,
#'   a,
#'   y,
#'   c(x1, x2),
#'   truncation_grid = c(0, 0.05, 0.1),
#'   n_boot = 25
#' )
#' autoplot(swept, type = "sweep")
#'
#' # A continuous exposure caps its weights at a quantile, so the sweep reads
#' # the quantile levels rather than a probability bound.
#' df$a <- x1 + x2 + rnorm(n)
#' df$y <- df$a + x1 + x2 + rnorm(n)
#' continuous <- check_eta_bias(
#'   df,
#'   a,
#'   y,
#'   c(x1, x2),
#'   truncation_grid = c(0, 0.05, 0.1),
#'   n_boot = 25
#' )
#' autoplot(continuous, type = "sweep")
#'
#' @name autoplot.eta_bias_result
#' @aliases plot.eta_bias_result
NULL

method(autoplot, eta_bias_result) <- function(
  object,
  type = c("bootstrap", "sweep"),
  ...
) {
  type <- resolve_arg_match(rlang::arg_match(type))
  if (type == "bootstrap") {
    autoplot_eta_bias_bootstrap(object)
  } else {
    autoplot_eta_bias_sweep(object, call = rlang::current_call())
  }
}

#' Describe the truncation levels a result was read over
#'
#' The `values` are what a sweep is drawn against, the `labels` are what a facet
#' strip reads, and `axis` names the quantity on the sweep's x axis. A discrete
#' run reads all three off the truncation lower bound it applied. A continuous
#' run caps a stabilized weight rather than bounding a probability, so it reads
#' them off the quantile levels its caps came from.
#'
#' @param object An `eta_bias_result`.
#'
#' @return A list of `values`, `labels`, and `axis`.
#' @keywords internal
#' @noRd
eta_bias_levels <- function(object) {
  results <- object@results
  n_levels <- eta_bias_n_levels(results, length(object@truth))
  if (object@exposure_type == "continuous") {
    # The caps themselves are weight values read off the observed fit, and the
    # cap of an untruncated level is infinite, so the quantile levels the caller
    # swept are what the sweep can be drawn against.
    values <- object@params$truncation_grid %||% 0
    return(list(
      values = values,
      labels = paste0("cap quantile = ", values),
      axis = "Weight cap quantile"
    ))
  }
  values <- results$truncation_lower[seq_len(n_levels)]
  list(
    values = values,
    labels = paste0("lower = ", values),
    axis = "Truncation lower bound"
  )
}

#' Index the truncation level each row of a results tibble was read at
#'
#' The rows run level fastest and term slowest, so each term's block of rows
#' walks the sweep once. The position is what identifies a level, because a
#' continuous run reports a missing bound at every one of them.
#'
#' @param n_levels,n_terms The truncation level and estimand term counts.
#'
#' @return An integer vector, one entry per row of the results tibble.
#' @keywords internal
#' @noRd
eta_bias_level_index <- function(n_levels, n_terms) {
  rep(seq_len(n_levels), times = n_terms)
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
  term_levels <- unique(results$term)
  n_terms <- length(object@truth)
  n_levels <- eta_bias_n_levels(results, n_terms)
  sweep_levels <- eta_bias_levels(object)
  level_index <- eta_bias_level_index(n_levels, n_terms)

  boot_data <- do.call(
    rbind,
    lapply(seq_len(nrow(results)), function(cell) {
      tibble::tibble(
        term = results$term[[cell]],
        level = level_index[[cell]],
        estimate = object@boot_estimates[[cell]]
      )
    })
  )
  # Both keys are factors so the panels follow the sweep and the order the terms
  # were built in. Left as they are, the level index would sort lexically and
  # put a tenth level ahead of a second one.
  boot_data$term <- factor(boot_data$term, levels = term_levels)
  boot_data$level <- factor(boot_data$level, levels = seq_len(n_levels))

  # The truth is a property of the term alone, so a term's facets carry the one
  # line the estimator in them was aimed at.
  truth_data <- tibble::tibble(
    term = factor(names(object@truth), levels = term_levels),
    truth = unname(object@truth)
  )

  ggplot2::ggplot(boot_data, ggplot2::aes(x = .data$estimate)) +
    ggplot2::geom_histogram(bins = 30, fill = "grey70", color = "white") +
    ggplot2::geom_vline(
      ggplot2::aes(xintercept = .data$truth),
      data = truth_data,
      linetype = "dashed",
      inherit.aes = FALSE
    ) +
    eta_bias_bootstrap_facets(n_terms, sweep_levels$labels) +
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

#' Cross the estimand terms with the truncation levels
#'
#' One term read at one level has no dimension to cross, so it keeps the single
#' facet naming the level it was read at. Anything wider crosses the two, which
#' keeps a level's terms in separate panels: pooling them would put the bootstrap
#' distributions of two different estimands in one histogram.
#'
#' @param n_terms The number of estimand terms.
#' @param labels The strip label for each truncation level.
#'
#' @return A [ggplot2::ggplot] facet specification.
#' @keywords internal
#' @noRd
eta_bias_bootstrap_facets <- function(n_terms, labels) {
  strip_labeller <- ggplot2::labeller(
    level = ggplot2::as_labeller(function(value) labels[as.integer(value)])
  )
  if (n_terms == 1 && length(labels) == 1) {
    return(ggplot2::facet_wrap(
      ggplot2::vars(.data$level),
      labeller = strip_labeller
    ))
  }
  ggplot2::facet_grid(
    rows = ggplot2::vars(.data$term),
    cols = ggplot2::vars(.data$level),
    labeller = strip_labeller
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
  n_terms <- length(object@truth)
  n_levels <- eta_bias_n_levels(results, n_terms)
  if (n_levels == 1) {
    # An estimand of several terms fills several rows at one truncation level,
    # so the row count reads as a sweep that never ran and the level count is
    # what decides there is nothing to draw.
    entry <- if (object@exposure_type == "continuous") {
      "quantile levels"
    } else {
      "lower bounds"
    }
    abort(
      c(
        "The sweep view needs a truncation sweep of more than one level.",
        i = "Rerun {.fn check_eta_bias} with a {.arg truncation_grid} of multiple {entry}."
      ),
      error_class = "positively_sweep_absent_error",
      call = call
    )
  }
  sweep_levels <- eta_bias_levels(object)
  sweep_data <- tibble::tibble(
    term = factor(results$term, levels = unique(results$term)),
    truncation = sweep_levels$values[eta_bias_level_index(n_levels, n_terms)],
    bias = results$bias,
    lower = results$bias - 2 * results$mc_se,
    upper = results$bias + 2 * results$mc_se
  )
  # One term needs no key to tell terms apart, so it keeps the neutral grey band
  # and the plain line. The colour and fill labels belong to those mappings;
  # naming an aesthetic the plot never maps draws an unknown-label message.
  several_terms <- n_terms > 1
  by_term <- if (several_terms) ggplot2::aes(colour = .data$term) else NULL
  band <- if (several_terms) {
    ggplot2::geom_ribbon(
      ggplot2::aes(ymin = .data$lower, ymax = .data$upper, fill = .data$term),
      alpha = 0.5
    )
  } else {
    ggplot2::geom_ribbon(
      ggplot2::aes(ymin = .data$lower, ymax = .data$upper),
      fill = "grey80",
      alpha = 0.5
    )
  }
  term_label <- if (several_terms) "Term"

  ggplot2::ggplot(
    sweep_data,
    ggplot2::aes(x = .data$truncation, y = .data$bias)
  ) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed") +
    band +
    ggplot2::geom_line(by_term) +
    ggplot2::geom_point(by_term) +
    ggplot2::labs(
      x = sweep_levels$axis,
      y = "ETA.Bias",
      colour = term_label,
      fill = term_label,
      title = "ETA.Bias across the truncation sweep",
      subtitle = "Band: plus or minus two Monte Carlo standard errors"
    )
}

method(plot, eta_bias_result) <- function(x, y, ...) {
  print(autoplot(x, ...))
  invisible(x)
}
