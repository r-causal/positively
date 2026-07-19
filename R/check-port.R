# S7 methods defined here are wired up at load time by S7::methods_register().

#' The PoRT subgroup diagnostic result class
#'
#' `port_result` is the S7 class returned by [check_port()] and
#' [check_port_seq()]. It extends [positivity_diagnostic] with the fitted
#' regression trees and the resolved reading-rule parameters. It is created
#' internally and is not constructed directly by users.
#'
#' @usage NULL
#' @keywords internal
#' @noRd
port_result <- new_class(
  "port_result",
  parent = positivity_diagnostic,
  properties = list(
    trees = class_list,
    alpha = class_double,
    beta = class_double,
    censoring_beta = class_double,
    gamma = class_double
  ),
  validator = function(self) {
    if (length(self@alpha) != 1) {
      "@alpha must be a single value"
    } else if (length(self@beta) < 1) {
      "@beta must hold at least one threshold"
    } else if (length(self@gamma) != 1) {
      "@gamma must be a single value"
    }
  }
)

#' Locate positivity-violating subgroups with the PoRT algorithm
#'
#' `check_port()` runs the Positivity Regression Tree (PoRT) algorithm of
#' Danelian et al. (2023). It grows regression trees on the covariates with the
#' exposure as the response, then reads a fixed rule off each tree to flag
#' subgroups in which one exposure level is rare or absent. It reports a
#' diagnostic only and assigns no severity grade.
#'
#' @details
#' PoRT searches for subgroups, defined by covariate splits, in which positivity
#' is threatened because an exposure level is (near) unobserved. For a binary
#' exposure the response is the indicator of the higher exposure level; for a
#' categorical exposure the algorithm runs once per level with a one-versus-rest
#' indicator; for a continuous exposure the exposure is first categorized into
#' `n_bins` quantile bins (or the bins implied by `breaks`) and each bin is
#' treated as a level. Covariates enter the trees unchanged: `rpart` splits
#' numeric covariates at data-driven thresholds and groups factor levels, so no
#' covariate categorization is performed here.
#'
#' The regression trees follow the control values quoted by Danelian et al.
#' (section 2.5): a minimum of 20 observations to attempt a split, at least six
#' observations in a leaf, a maximum depth of 30, and no cost-complexity pruning
#' (`cp = 0`), so that the largest number of divisions is retained. The paper
#' states a minimum leaf size of six; `rpart`'s own default of
#' `round(minsplit / 3)` is seven, so the value is pinned to six to match the
#' paper. This affects only the smallest leaves and is immaterial at realistic
#' sample sizes, where the `alpha` size gate dominates.
#'
#' The reading rule (Danelian et al., section 2.3) flags a node when its
#' exposure prevalence is extreme, below `beta` or above `1 - beta`, and its
#' size is at least `alpha` of the sample. Once a node is flagged its
#' descendants are suppressed, so the flag is reported at the coarsest subgroup
#' that isolates the violation. Because the rule is deterministic given a tree,
#' only the tree carries stochasticity.
#'
#' The combination search controls how many covariates may jointly define a
#' subgroup (Danelian et al., section 2.4). Step one grows one tree per single
#' covariate; covariates that produce a flagged subgroup are removed, and step
#' `k` grows one tree per `k`-covariate combination of the remaining covariates,
#' up to `gamma`. Raising `gamma` uncovers joint violations that no single
#' covariate isolates, at the cost of more trees.
#'
#' Setting `beta = "gruber"` resolves the threshold to
#' \eqn{5 / (\sqrt{n} \, \ln n)}, the sample-size-adaptive bound used by sPoRT.
#' It is tighter than `0.05` for samples larger than roughly 350 and looser
#' below that, so a small sample is held to a less stringent prevalence
#' threshold.
#'
#' @param .data A data frame.
#' @param .exposure The exposure column, selected with data-masking.
#' @param .covariates The covariate columns, selected with tidyselect. Numeric
#'   and factor covariates are both accepted; the trees split them directly.
#' @param alpha The minimum subgroup size for the reading rule, a fraction of
#'   the sample size strictly between 0 and 1. Defaults to `0.05`.
#' @param beta The prevalence threshold for the reading rule. A single number
#'   strictly between 0 and 0.5 (`0.05` by default), or the string `"gruber"`
#'   for the sample-size-adaptive bound \eqn{5 / (\sqrt{n} \, \ln n)}. The upper
#'   limit is 0.5 because the rule flags a prevalence below `beta` or above
#'   `1 - beta`.
#' @param gamma The maximum number of covariates that may jointly define a
#'   subgroup, a whole number of at least one. Defaults to `2`.
#' @param n_bins For a continuous exposure, the number of quantile bins into
#'   which the exposure is categorized. A whole number of at least two, `3` by
#'   default; a single bin would place every observation at one exposure level
#'   and flag every subgroup. Ignored for binary and categorical exposures.
#' @param breaks For a continuous exposure, at least three distinct numeric cut
#'   points that override `n_bins`. `NULL` (the default) uses quantile bins.
#' @param exposure_type One of `"auto"` (detect from the data, the default),
#'   `"binary"`, `"categorical"`, or `"continuous"`.
#' @param ... Passed to [rpart::rpart.control()], for example `cp`.
#'
#' @return A `port_result` object, an S7 subclass of [positivity_diagnostic].
#'   Its `@results` tibble has one row per reported subgroup with columns
#'   `subgroup` (the covariates defining the rule), `description` (the
#'   human-readable rule), `exposure_level` (the exposure level examined), `n`
#'   (the subgroup size), `proportion` (its share of the sample), `prevalence`
#'   (the exposure prevalence in the subgroup), and `flagged` (the logical
#'   reading-rule outcome). It also carries the `@trees`, `@alpha`, `@beta`, and
#'   `@gamma` properties.
#'
#' @references
#' Danelian G, Foucher Y, Léger M, Le Borgne F, Chatton A (2023). Identification
#' of in-sample positivity violations using regression trees: the PoRT
#' algorithm. \doi{10.1515/jci-2022-0032}
#'
#' Chatton A, Schomaker M, Luque-Fernandez MA, Platt RW, Schnitzer ME (2025).
#' Is checking for sequential positivity violations getting you down? Try
#' sPoRT!
#'
#' @examples
#' set.seed(1)
#' n <- 400
#' x1 <- rnorm(n)
#' x2 <- factor(sample(c("a", "b"), n, replace = TRUE))
#' # No subject with x1 > 1 and x2 == "b" is ever treated.
#' ps <- ifelse(x1 > 1 & x2 == "b", 0, 0.5)
#' exposure <- rbinom(n, 1, ps)
#' df <- data.frame(exposure = exposure, x1 = x1, x2 = x2)
#'
#' result <- check_port(df, exposure, c(x1, x2))
#' result
#'
#' @export
check_port <- function(
  .data,
  .exposure,
  .covariates,
  alpha = 0.05,
  beta = 0.05,
  gamma = 2,
  n_bins = 3,
  breaks = NULL,
  exposure_type = c("auto", "binary", "categorical", "continuous"),
  ...
) {
  validate_data_frame(.data)
  validate_probability(alpha, "alpha")
  beta_spec <- validate_beta(beta)
  validate_count(gamma, "gamma")
  validate_count(n_bins, "n_bins", min = 2)
  validate_breaks(breaks)

  exposure_name <- select_single_exposure(rlang::enquo(.exposure), .data)
  exposure_vec <- .data[[exposure_name]]
  validate_two_level_exposure(exposure_vec)

  resolved_type <- match_exposure_type(exposure_type, exposure_vec)

  covariate_pos <- tidyselect::eval_select(rlang::enquo(.covariates), .data)
  validate_column_selection(covariate_pos, ".covariates")
  covariate_names <- names(covariate_pos)
  validate_complete_columns(.data, covariate_names, ".covariates")

  n <- nrow(.data)
  beta_value <- resolve_beta_scalar(beta_spec, n)
  if (
    identical(beta_spec, "gruber") &&
      (!is.finite(beta_value) || beta_value >= 0.5)
  ) {
    abort(
      c(
        "The {.val gruber} prevalence threshold resolved to {.val {beta_value}} at {.arg n} = {n}.",
        i = "The reading rule needs a threshold below 0.5, but the sample-size-adaptive bound reaches 0.5 or more at small sample sizes.",
        i = "Supply a numeric {.arg beta} below 0.5 for a sample this small."
      ),
      error_class = "positively_range_error"
    )
  }

  level_pairs <- port_level_responses(
    exposure_vec,
    resolved_type,
    n_bins = n_bins,
    breaks = breaks
  )

  search <- run_port(
    covariate_data = .data[covariate_names],
    predictors = covariate_names,
    level_pairs = level_pairs,
    alpha = alpha,
    beta_value = beta_value,
    gamma = gamma,
    n_total = n,
    control = port_control(...)
  )

  results <- assemble_port_results(search$rows, sequential = FALSE)

  port_result(
    results = results,
    exposure = exposure_name,
    exposure_type = resolved_type,
    n = as.integer(n),
    params = list(
      alpha = alpha,
      beta = beta,
      gamma = gamma,
      n_bins = n_bins,
      breaks = breaks
    ),
    call = rlang::current_call(),
    trees = search$trees,
    alpha = as.double(alpha),
    beta = beta_value,
    gamma = as.double(gamma)
  )
}

# ---- Parameter resolution and validation ----------------------------------

#' Validate the `beta` argument of PoRT
#'
#' `beta` is either a single number strictly between 0 and 0.5 or the string
#' `"gruber"`, which selects the sample-size-adaptive bound. The upper limit is
#' 0.5 because the reading rule flags a prevalence below `beta` or above
#' `1 - beta`; at `beta = 0.5` the two thresholds meet and every subgroup is
#' flagged.
#'
#' @param beta The supplied `beta`.
#' @param call The calling environment, used to build the error's call.
#'
#' @return The string `"gruber"` or the validated numeric `beta`.
#' @keywords internal
#' @noRd
validate_beta <- function(beta, call = rlang::caller_env()) {
  if (is.character(beta)) {
    if (length(beta) != 1 || !identical(beta, "gruber")) {
      abort(
        c(
          '{.arg beta} must be a number in (0, 0.5) or {.val gruber}.',
          x = "Received {.val {beta}}."
        ),
        error_class = "positively_beta_error",
        call = call
      )
    }
    return("gruber")
  }
  validate_probability(beta, "beta", call = call)
  if (beta >= 0.5) {
    abort(
      c(
        "{.arg beta} must be below {.val {0.5}}.",
        i = "The reading rule flags prevalence below {.arg beta} or above {.code 1 - beta}, so a value of 0.5 or more flags every subgroup.",
        x = "Found {.val {beta}}."
      ),
      error_class = "positively_range_error",
      call = call
    )
  }
  as.double(beta)
}

#' Validate the `breaks` argument of PoRT
#'
#' @param breaks The supplied `breaks`, either `NULL` or numeric cut points.
#' @param call The calling environment, used to build the error's call.
#'
#' @return `breaks`, invisibly, when it is `NULL` or numeric.
#' @keywords internal
#' @noRd
validate_breaks <- function(breaks, call = rlang::caller_env()) {
  if (is.null(breaks)) {
    return(invisible(breaks))
  }
  if (!is.numeric(breaks)) {
    breaks_class <- class(breaks)[1]
    abort(
      "{.arg breaks} must be numeric or {.code NULL}, not a {.cls {breaks_class}}.",
      error_class = "positively_type_error",
      call = call
    )
  }
  if (anyNA(breaks)) {
    abort(
      "{.arg breaks} must not contain missing values.",
      error_class = "positively_missing_error",
      call = call
    )
  }
  n_distinct <- length(unique(breaks))
  if (n_distinct < 3) {
    abort(
      c(
        "{.arg breaks} must contain at least three distinct cut points.",
        i = "The cut points bound the exposure bins, so three points define the minimum of two bins; a single bin would place every observation at one exposure level and flag every subgroup.",
        x = "Found {n_distinct} distinct cut point{?s}."
      ),
      error_class = "positively_range_error",
      call = call
    )
  }
  invisible(breaks)
}

#' Validate that an exposure has at least two distinct values
#'
#' @param exposure The exposure vector.
#' @param arg The argument name used in error messages.
#' @param call The calling environment, used to build the error's call.
#'
#' @return `exposure`, invisibly, when it has at least two distinct values.
#' @keywords internal
#' @noRd
validate_two_level_exposure <- function(
  exposure,
  arg = ".exposure",
  call = rlang::caller_env()
) {
  n_levels <- length(unique(exposure[!is.na(exposure)]))
  if (n_levels < 2) {
    abort(
      c(
        "{.arg {arg}} must have at least two distinct values.",
        x = "Found {n_levels}."
      ),
      error_class = "positively_exposure_error",
      call = call
    )
  }
  invisible(exposure)
}

#' The Gruber sample-size-adaptive prevalence bound
#'
#' @param n A sample size.
#'
#' @return The bound \eqn{5 / (\sqrt{n} \, \ln n)}.
#' @keywords internal
#' @noRd
gruber_beta <- function(n) {
  5 / (sqrt(n) * log(n))
}

#' Resolve the `beta` specification to a numeric threshold
#'
#' @param beta_spec Either `"gruber"` or a numeric threshold.
#' @param n The sample size, used only for the Gruber bound.
#'
#' @return A single numeric threshold.
#' @keywords internal
#' @noRd
resolve_beta_scalar <- function(beta_spec, n) {
  if (identical(beta_spec, "gruber")) {
    gruber_beta(n)
  } else {
    as.double(beta_spec)
  }
}

#' Build the `rpart` control following the PoRT paper
#'
#' @param ... Overrides passed to [rpart::rpart.control()].
#'
#' @return An `rpart.control` list.
#' @keywords internal
#' @noRd
port_control <- function(...) {
  defaults <- list(
    minsplit = 20L,
    minbucket = 6L,
    maxdepth = 30L,
    cp = 0,
    xval = 0L
  )
  do.call(rpart::rpart.control, utils::modifyList(defaults, list(...)))
}

# ---- Exposure-level responses ---------------------------------------------

#' Build the one-versus-rest response vectors PoRT examines
#'
#' A binary exposure yields a single response, the indicator of its higher
#' level. A categorical exposure yields one response per level. A continuous
#' exposure is first categorized into bins, then handled like a categorical
#' exposure with one response per bin.
#'
#' @param exposure The exposure vector.
#' @param exposure_type The resolved exposure type.
#' @param n_bins The number of quantile bins for a continuous exposure.
#' @param breaks Explicit cut points for a continuous exposure, or `NULL`.
#'
#' @return A list of `list(level, response)` pairs, where `level` is a label and
#'   `response` is an integer 0/1 indicator vector.
#' @keywords internal
#' @noRd
port_level_responses <- function(exposure, exposure_type, n_bins, breaks) {
  if (exposure_type == "binary") {
    present <- sort(unique(exposure[!is.na(exposure)]))
    positive <- present[length(present)]
    return(list(list(
      level = as.character(positive),
      response = as.integer(exposure == positive)
    )))
  }

  membership <- if (exposure_type == "categorical") {
    exposure
  } else {
    bin_continuous_exposure(exposure, n_bins = n_bins, breaks = breaks)
  }

  levels_present <- level_labels(membership)
  lapply(levels_present, function(level) {
    list(
      level = level,
      response = as.integer(as.character(membership) == level)
    )
  })
}

#' The distinct level labels of a membership vector, in a stable order
#'
#' @param membership A factor or vector of level memberships.
#'
#' @return A character vector of level labels.
#' @keywords internal
#' @noRd
level_labels <- function(membership) {
  if (is.factor(membership)) {
    present <- levels(membership)[tabulate(membership, nlevels(membership)) > 0]
    return(as.character(present))
  }
  as.character(sort(unique(membership[!is.na(membership)])))
}

#' Categorize a continuous exposure into bins
#'
#' @param exposure A numeric exposure vector.
#' @param n_bins The number of quantile bins.
#' @param breaks Explicit numeric cut points, or `NULL` for quantile bins.
#'
#' @return A factor of bin memberships.
#' @keywords internal
#' @noRd
bin_continuous_exposure <- function(exposure, n_bins, breaks) {
  cut_points <- if (!is.null(breaks)) {
    sort(unique(breaks))
  } else {
    probs <- seq(0, 1, length.out = n_bins + 1)
    unique(stats::quantile(exposure, probs = probs, na.rm = TRUE, type = 7))
  }
  if (length(cut_points) < 3) {
    abort(
      c(
        "The exposure quantiles define fewer than two bins.",
        i = "Ties in the exposure collapse the quantile cut points, so the requested bins reduce to one.",
        i = "Supply explicit {.arg breaks} that separate the tied values."
      ),
      error_class = "positively_binning_error",
      call = NULL
    )
  }
  cut(exposure, breaks = cut_points, include.lowest = TRUE)
}

# ---- Combination search ---------------------------------------------------

#' Run the PoRT combination search across every exposure level
#'
#' @param covariate_data A data frame of covariate columns.
#' @param predictors The covariate column names.
#' @param level_pairs The `list(level, response)` pairs from
#'   `port_level_responses()`.
#' @param alpha The reading-rule size gate.
#' @param beta_value The resolved numeric prevalence threshold.
#' @param gamma The combination-search depth.
#' @param n_total The sample size the size gate is measured against.
#' @param control The `rpart` control list.
#'
#' @return A list with `rows` (a list of per-subgroup tibbles) and `trees` (a
#'   list of the fitted `rpart` objects).
#' @keywords internal
#' @noRd
run_port <- function(
  covariate_data,
  predictors,
  level_pairs,
  alpha,
  beta_value,
  gamma,
  n_total,
  control
) {
  rows <- list()
  trees <- list()
  for (pair in level_pairs) {
    search <- port_search(
      response = pair$response,
      covariate_data = covariate_data,
      predictors = predictors,
      alpha = alpha,
      beta_value = beta_value,
      gamma = gamma,
      n_total = n_total,
      control = control
    )
    for (subgroup_rows in search$rows) {
      if (nrow(subgroup_rows) > 0) {
        subgroup_rows$exposure_level <- pair$level
        rows[[length(rows) + 1]] <- subgroup_rows
      }
    }
    trees <- c(trees, search$trees)
  }
  list(rows = rows, trees = trees)
}

#' The PoRT combination search for one exposure-level response
#'
#' Grows single-covariate trees, removes covariates that produce a flagged
#' subgroup, then grows trees over `k`-covariate combinations of the survivors
#' up to `gamma`.
#'
#' @inheritParams run_port
#' @param response A single 0/1 indicator response vector.
#'
#' @return A list with `rows` (per-tree tibbles) and `trees` (fitted `rpart`
#'   objects).
#' @keywords internal
#' @noRd
port_search <- function(
  response,
  covariate_data,
  predictors,
  alpha,
  beta_value,
  gamma,
  n_total,
  control
) {
  remaining <- predictors
  rows <- list()
  trees <- list()

  for (k in seq_len(gamma)) {
    if (length(remaining) < k) {
      break
    }
    combos <- utils::combn(remaining, k, simplify = FALSE)
    flagged_predictors <- character(0)
    for (combo in combos) {
      tree <- fit_port_tree(response, covariate_data[combo], control)
      report <- read_port_tree(tree, n_total, alpha, beta_value)
      trees[[length(trees) + 1]] <- tree
      rows[[length(rows) + 1]] <- report$rows
      if (report$any_flagged) {
        flagged_predictors <- union(flagged_predictors, combo)
      }
    }
    remaining <- setdiff(remaining, flagged_predictors)
  }

  list(rows = rows, trees = trees)
}

#' Fit one PoRT regression tree
#'
#' @param response A 0/1 indicator response vector.
#' @param covariate_data A data frame of the covariates for this tree.
#' @param control The `rpart` control list.
#'
#' @return A fitted `rpart` object.
#' @keywords internal
#' @noRd
fit_port_tree <- function(response, covariate_data, control) {
  # The response column name must not collide with a covariate: rpart segfaults
  # on a formula whose response term also appears as a predictor.
  response_name <- make.unique(
    c(names(covariate_data), ".port_response")
  )[[ncol(covariate_data) + 1L]]
  tree_data <- covariate_data
  tree_data[[response_name]] <- response
  rpart::rpart(
    stats::reformulate(names(covariate_data), response = response_name),
    data = tree_data,
    method = "anova",
    control = control
  )
}

# ---- Reading the tree -----------------------------------------------------

#' Read the PoRT rule off one tree
#'
#' Walks the tree from the root and collects reported subgroups: a node that
#' meets the reading rule is reported as flagged and its descendants are
#' suppressed; a leaf that is not flagged is reported with a `FALSE` flag. The
#' whole-sample root is never reported as a subgroup unless the tree made no
#' split, in which case the single root node is returned.
#'
#' @param tree A fitted `rpart` object.
#' @param n_total The sample size the size gate is measured against.
#' @param alpha The reading-rule size gate.
#' @param beta The numeric prevalence threshold.
#'
#' @return A list with `rows` (a tibble of reported subgroups, without the
#'   `exposure_level` column) and `any_flagged` (a single logical).
#' @keywords internal
#' @noRd
read_port_tree <- function(tree, n_total, alpha, beta) {
  frame <- tree$frame
  node_ids <- as.integer(rownames(frame))
  is_leaf <- frame$var == "<leaf>"
  prevalence <- frame$yval
  size_n <- frame$n
  size_fraction <- size_n / n_total
  flagged <- (prevalence < beta | prevalence > 1 - beta) &
    (size_fraction >= alpha)

  collected <- integer(0)
  visit <- function(node) {
    idx <- match(node, node_ids)
    node_is_root <- node == 1L
    if (node_is_root && is_leaf[idx]) {
      collected <<- c(collected, node)
      return(invisible())
    }
    if (!node_is_root && (flagged[idx] || is_leaf[idx])) {
      collected <<- c(collected, node)
      return(invisible())
    }
    left <- 2L * node
    right <- 2L * node + 1L
    if (left %in% node_ids) {
      visit(left)
    }
    if (right %in% node_ids) {
      visit(right)
    }
  }
  visit(1L)

  rows <- port_tree_rows(
    tree,
    collected,
    node_ids,
    size_n,
    size_fraction,
    prevalence,
    flagged
  )
  list(rows = rows, any_flagged = any(rows$flagged))
}

#' Build the subgroup rows for a set of collected nodes
#'
#' @param tree A fitted `rpart` object.
#' @param nodes The collected node numbers.
#' @param node_ids The node numbers of every frame row.
#' @param size_n The per-node observation counts.
#' @param size_fraction The per-node sample-size fractions.
#' @param prevalence The per-node exposure prevalences.
#' @param flagged The per-node reading-rule outcomes.
#'
#' @return A tibble of subgroup rows.
#' @keywords internal
#' @noRd
port_tree_rows <- function(
  tree,
  nodes,
  node_ids,
  size_n,
  size_fraction,
  prevalence,
  flagged
) {
  if (length(nodes) == 0) {
    return(empty_port_tibble(sequential = FALSE))
  }
  paths <- rpart::path.rpart(tree, nodes = nodes, print.it = FALSE)
  described <- lapply(nodes, function(node) {
    idx <- match(node, node_ids)
    rule <- port_node_rule(paths[[as.character(node)]])
    tibble::tibble(
      subgroup = rule$subgroup,
      description = rule$description,
      n = as.integer(size_n[idx]),
      proportion = size_fraction[idx],
      prevalence = prevalence[idx],
      flagged = flagged[idx]
    )
  })
  vctrs::vec_rbind(!!!described)
}

#' Turn a node path into a subgroup label and a human-readable rule
#'
#' @param path The character vector returned for one node by
#'   [rpart::path.rpart()], beginning with `"root"`.
#'
#' @return A list with `subgroup` (the covariates in the rule) and `description`
#'   (the rule joined by ` & `).
#' @keywords internal
#' @noRd
port_node_rule <- function(path) {
  conditions <- path[path != "root"]
  if (length(conditions) == 0) {
    return(list(subgroup = "(overall)", description = "(overall)"))
  }
  parsed <- parse_port_conditions(conditions)
  variables <- unique(parsed$variable)
  simplified <- unlist(lapply(variables, function(variable) {
    simplify_variable_conditions(parsed[parsed$variable == variable, ])
  }))
  list(
    subgroup = paste(variables, collapse = ", "),
    description = paste(simplified, collapse = " & ")
  )
}

#' Split `rpart` path conditions into variable, operator, and value
#'
#' @param conditions The non-root elements of a [rpart::path.rpart()] path.
#'
#' @return A data frame with `variable`, `operator`, and `value` columns.
#' @keywords internal
#' @noRd
parse_port_conditions <- function(conditions) {
  matches <- regmatches(
    conditions,
    regexec("^(.*?)(>=|<=|>|<|=)(.*)$", conditions)
  )
  data.frame(
    variable = vapply(matches, function(m) trimws(m[[2]]), character(1)),
    operator = vapply(matches, function(m) m[[3]], character(1)),
    value = vapply(matches, function(m) trimws(m[[4]]), character(1)),
    stringsAsFactors = FALSE
  )
}

#' Reduce repeated conditions on one variable to the binding constraints
#'
#' For a numeric variable this keeps the largest lower bound and the smallest
#' upper bound, the only conditions that bind. For a factor variable it
#' intersects the level sets, since each split narrows the retained levels.
#'
#' @param cond A data frame of parsed conditions for a single variable.
#'
#' @return A character vector of at most two conditions.
#' @keywords internal
#' @noRd
simplify_variable_conditions <- function(cond) {
  variable <- cond$variable[[1]]
  if (all(cond$operator == "=")) {
    level_sets <- strsplit(cond$value, ",", fixed = TRUE)
    kept <- Reduce(intersect, level_sets)
    return(paste0(variable, "=", paste(kept, collapse = ",")))
  }
  bounds <- character(0)
  lower <- which(cond$operator %in% c(">=", ">"))
  upper <- which(cond$operator %in% c("<", "<="))
  if (length(lower) > 0) {
    pick <- lower[which.max(as.numeric(cond$value[lower]))]
    bounds <- c(bounds, paste0(variable, cond$operator[pick], cond$value[pick]))
  }
  if (length(upper) > 0) {
    pick <- upper[which.min(as.numeric(cond$value[upper]))]
    bounds <- c(bounds, paste0(variable, cond$operator[pick], cond$value[pick]))
  }
  bounds
}

# ---- Assembly -------------------------------------------------------------

#' An empty PoRT results tibble with the fixed columns and types
#'
#' @param sequential Whether to include the leading `time` column.
#'
#' @return A zero-row tibble.
#' @keywords internal
#' @noRd
empty_port_tibble <- function(sequential) {
  base <- tibble::tibble(
    subgroup = character(0),
    description = character(0),
    exposure_level = character(0),
    n = integer(0),
    proportion = double(0),
    prevalence = double(0),
    flagged = logical(0)
  )
  if (sequential) {
    return(vctrs::vec_cbind(tibble::tibble(time = integer(0)), base))
  }
  base
}

#' Assemble the per-subgroup rows into the fixed results tibble
#'
#' @param rows A list of subgroup tibbles.
#' @param sequential Whether the leading `time` column is present.
#'
#' @return A tibble with the fixed PoRT columns in order.
#' @keywords internal
#' @noRd
assemble_port_results <- function(rows, sequential) {
  base_cols <- c(
    "subgroup",
    "description",
    "exposure_level",
    "n",
    "proportion",
    "prevalence",
    "flagged"
  )
  cols <- if (sequential) c("time", base_cols) else base_cols
  populated <- rows[vapply(rows, nrow, integer(1)) > 0]
  if (length(populated) == 0) {
    return(empty_port_tibble(sequential))
  }
  combined <- vctrs::vec_rbind(!!!populated)
  combined[cols]
}

# ---- Methods --------------------------------------------------------------

method(print, port_result) <- function(x, ...) {
  flagged_count <- sum(x@results$flagged)
  cat_cli({
    cli::cli_h1("{S7::S7_class(x)@name}")
    cli::cli_text("Exposure: {.val {x@exposure}} ({x@exposure_type})")
    cli::cli_text("Observations: {x@n}")
    if ("time" %in% names(x@results)) {
      cli::cli_text("Time points: {length(unique(x@results$time))}")
    }
    cli::cli_text("Reading rule: alpha = {x@alpha}, gamma = {x@gamma}")
    beta_finite <- x@beta[is.finite(x@beta)]
    if (length(unique(beta_finite)) <= 1) {
      value <- if (length(beta_finite) == 0) {
        "not resolved"
      } else {
        signif(beta_finite[[1]], 3)
      }
      cli::cli_text("Prevalence threshold beta: {value}")
    } else {
      cli::cli_text(
        "Prevalence threshold beta: {length(beta_finite)} per-time values ({signif(min(beta_finite), 3)} to {signif(max(beta_finite), 3)})"
      )
    }
    if ("type" %in% names(x@results)) {
      is_censoring <- x@results$type == "censoring"
      exposure_res <- x@results[!is_censoring, , drop = FALSE]
      censoring_res <- x@results[is_censoring, , drop = FALSE]
      cli::cli_text(
        "Exposure subgroups: {nrow(exposure_res)} reported, {sum(exposure_res$flagged)} flagged"
      )
      cli::cli_text(
        "Censoring subgroups: {nrow(censoring_res)} reported, {sum(censoring_res$flagged)} flagged"
      )
    } else {
      cli::cli_text(
        "Subgroups: {nrow(x@results)} reported, {flagged_count} flagged"
      )
    }
  })
  invisible(x)
}
