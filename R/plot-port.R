# plot() calls autoplot() and prints, so both entry points draw the same view.

#' Plot a PoRT subgroup diagnostic
#'
#' Draws the reported subgroups as a horizontal bar chart of exposure
#' prevalence, with dashed reference lines at `beta` and `1 - beta`. Bar width
#' is proportional to subgroup size, and flagged subgroups are highlighted. A
#' sequential result is faceted by time point, and each facet carries its own
#' resolved thresholds. A censoring run is faceted by both time point and row
#' type, so the exposure and censoring subgroups are separated, the censoring
#' facets plot censoring prevalence, and each facet's reference lines are drawn
#' from the thresholds that family was judged against.
#'
#' @param object A `port_result` from [check_port()] or [check_port_seq()].
#' @param ... Not used.
#'
#' @return A [ggplot2::ggplot] object.
#'
#' @examples
#' set.seed(1)
#' n <- 400
#' x1 <- rnorm(n)
#' exposure <- rbinom(n, 1, ifelse(x1 > 1, 0, 0.5))
#' df <- data.frame(exposure = exposure, x1 = x1)
#' result <- check_port(df, exposure, x1)
#'
#' ggplot2::autoplot(result)
#'
#' @name autoplot.port_result
#' @aliases plot.port_result
NULL

method(autoplot, port_result) <- function(object, ...) {
  plot_data <- port_plot_data(object)
  reference <- port_beta_reference(object)
  sequential <- "time" %in% names(object@results)
  has_type <- "type" %in% names(object@results)

  # A censoring run pools exposure and censoring subgroups into one plot, so the
  # axis cannot claim exposure prevalence for every bar.
  x_label <- if (has_type) {
    "Prevalence in subgroup"
  } else {
    "Exposure prevalence in subgroup"
  }

  plot <- ggplot2::ggplot(
    plot_data,
    ggplot2::aes(
      x = .data$prevalence,
      y = .data$label,
      fill = .data$flagged
    )
  ) +
    ggplot2::geom_col(ggplot2::aes(width = .data$width)) +
    ggplot2::geom_vline(
      data = reference,
      mapping = ggplot2::aes(xintercept = .data$xintercept),
      linetype = "dashed"
    ) +
    ggplot2::labs(
      x = x_label,
      y = "Subgroup",
      fill = "Flagged",
      title = "PoRT subgroups against the prevalence thresholds"
    )

  # The manual fill scale has no levels to match when there are no bars, and a
  # discrete scale with an all-empty column warns; the empty panel needs no fill.
  if (nrow(plot_data) > 0) {
    plot <- plot +
      ggplot2::scale_fill_manual(
        values = c("FALSE" = "grey70", "TRUE" = "#B2182B"),
        drop = FALSE
      )
  }

  has_facets <- nrow(plot_data) > 0 || nrow(reference) > 0
  if (has_type && has_facets) {
    plot <- plot +
      ggplot2::facet_wrap(
        ggplot2::vars(.data$time, .data$type),
        ncol = 2,
        labeller = ggplot2::label_both,
        scales = "free_y"
      )
  } else if (sequential && has_facets) {
    plot <- plot +
      ggplot2::facet_wrap(
        ggplot2::vars(.data$time),
        ncol = 1,
        labeller = ggplot2::label_both,
        scales = "free_y"
      )
  }

  plot
}

#' Assemble the plotting frame for a PoRT result
#'
#' Duplicate rows, which arise when the same leaf is reported by more than one
#' tree in the combination search, are collapsed to one bar so that prevalences
#' are drawn rather than summed. The duplicates remain in `@results`. Bar width
#' is set proportional to subgroup size.
#'
#' @param object A `port_result`.
#'
#' @return A data frame of the unique reported subgroups with a display `label`,
#'   a `width` column, and a factor `flagged` column.
#' @keywords internal
#' @noRd
port_plot_data <- function(object) {
  plot_data <- object@results
  key_cols <- intersect(
    c(
      "time",
      "type",
      "description",
      "exposure_level",
      "prevalence",
      "n",
      "flagged"
    ),
    names(plot_data)
  )
  plot_data <- plot_data[!duplicated(plot_data[key_cols]), , drop = FALSE]

  if (nrow(plot_data) == 0) {
    plot_data$label <- character(0)
    plot_data$width <- double(0)
    plot_data$flagged <- factor(character(0), levels = c("FALSE", "TRUE"))
    return(plot_data)
  }

  plot_data$label <- paste0(
    plot_data$description,
    " [",
    plot_data$exposure_level,
    "]"
  )
  plot_data$width <- plot_data$n / max(plot_data$n)
  plot_data$flagged <- factor(
    as.character(plot_data$flagged),
    levels = c("FALSE", "TRUE")
  )
  plot_data
}

#' The per-time beta reference-line frame
#'
#' Builds one row per reference line so that each facet shows the thresholds
#' actually applied within it rather than a shared pair. Skipped time points,
#' whose resolved threshold is `NA`, contribute no lines. On a censoring run the
#' exposure and censoring facets are each drawn against their own resolved
#' thresholds, since the two families are judged on different risk-set sizes.
#'
#' @param object A `port_result`.
#'
#' @return A tibble with an `xintercept` column, plus a `time` column for a
#'   sequential result and a `type` column for a censoring result.
#' @keywords internal
#' @noRd
port_beta_reference <- function(object) {
  results <- object@results
  if (!("time" %in% names(results))) {
    beta <- object@beta
    return(tibble::tibble(xintercept = c(beta[[1]], 1 - beta[[1]])))
  }
  if (!("type" %in% names(results))) {
    return(beta_reference_frame(object@beta))
  }
  exposure_reference <- beta_reference_frame(object@beta)
  exposure_reference$type <- rep("exposure", nrow(exposure_reference))
  censoring_reference <- beta_reference_frame(object@censoring_beta)
  censoring_reference$type <- rep("censoring", nrow(censoring_reference))
  vctrs::vec_rbind(exposure_reference, censoring_reference)
}

#' One reference-line row pair per finite per-time threshold
#'
#' @param beta A per-time threshold vector, with `NA` for skipped time points.
#'
#' @return A tibble with `time` and `xintercept` columns.
#' @keywords internal
#' @noRd
beta_reference_frame <- function(beta) {
  times <- seq_along(beta)
  keep <- is.finite(beta)
  tibble::tibble(
    time = rep(times[keep], each = 2),
    xintercept = as.numeric(rbind(beta[keep], 1 - beta[keep]))
  )
}

method(plot, port_result) <- function(x, y, ...) {
  print(autoplot(x, ...))
  invisible(x)
}
