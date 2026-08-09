# The pluggable conditional-density contract for the HDR diagnostic. A density
# estimator is an internal S7 object of class `hdr_density` carrying three
# function properties that check_hdr() calls in turn. S7 methods here are
# wired up at load time by S7::methods_register().

#' The conditional-density estimator class for HDR non-overlap
#'
#' `hdr_density` is the internal S7 class that supplies a pluggable conditional
#' density to [check_hdr()] and [check_hdr_seq()]. It is constructed through
#' [new_hdr_density()] or the default [hdr_density_normal()]; users do not build
#' it with `hdr_density()` directly.
#'
#' @usage NULL
#' @keywords internal
#' @noRd
hdr_density <- new_class(
  "hdr_density",
  properties = list(
    fit = class_function,
    density = class_function,
    hdr_threshold = NULL | class_function,
    label = class_character
  )
)

#' Construct a conditional-density estimator for the HDR diagnostic
#'
#' `new_hdr_density()` is the developer constructor for the pluggable
#' conditional-density contract used by [check_hdr()] and [check_hdr_seq()]. It
#' packages three functions into an estimator object that the HDR non-overlap
#' ratio of Bao and Schomaker (2025) is computed against. Supply an
#' `hdr_threshold` to give a closed-form density cutoff, or leave it `NULL` to
#' request the numeric-grid fallback.
#'
#' @details
#' The highest-density region (HDR) at a covariate profile \eqn{l} is the
#' smallest set of exposure values that captures probability mass `mass` of the
#' conditional density \eqn{f(a \mid l)}, namely
#' \eqn{A_\alpha(l) = \{ a : f(a \mid l) \ge f_\alpha(l) \}}. An estimator
#' therefore needs three pieces, which map to the three arguments:
#'
#' - `fit(formula, data)` fits the conditional-density model and returns any
#'   fitted state as a list. `check_hdr()` builds `formula` as
#'   `exposure ~ covariates` and passes the analysis `data`.
#' - `density(state, a, newdata)` returns \eqn{\hat{f}(a \mid l)} for the
#'   scalar or vector target `a`, recycled against the rows of `newdata`. It
#'   receives the state from `fit()`.
#' - `hdr_threshold(state, newdata, mass)` returns the density cutoff
#'   \eqn{f_\alpha(l_j)} whose HDR captures `mass`, one value per row of
#'   `newdata`. When it is `NULL`, [check_hdr()] instead evaluates `density`
#'   on a fine grid of exposure values per row and reads the cutoff off that
#'   grid, so the minimum contract is two functions.
#'
#' Membership of a target `a` in the HDR at row \eqn{j} is the test
#' \eqn{\hat{f}(a \mid l_j) \ge f_\alpha(l_j)}, and the non-overlap ratio
#' \eqn{\hat{\tau}(a)} is the fraction of rows whose HDR excludes `a`.
#'
#' `fit()` and `density()` receive the analysis data as a data frame holding the
#' exposure and the covariate columns as they were selected, so a factor or
#' character covariate arrives as a factor or character column. Encoding those
#' columns is the estimator's work. [hdr_density_normal()] leaves it to
#' [stats::lm()], which expands them into contrast terms from the formula; an
#' estimator that needs a numeric design matrix should build one with
#' [stats::model.matrix()].
#'
#' @param fit A function of `(formula, data)` returning fitted state as a list.
#' @param density A function of `(state, a, newdata)` returning the conditional
#'   density \eqn{\hat{f}(a \mid l)} for target `a` over the rows of `newdata`.
#' @param hdr_threshold An optional function of `(state, newdata, mass)`
#'   returning the per-row density cutoff. `NULL` (the default) requests the
#'   numeric-grid fallback. When it is `NULL`, [check_hdr()] evaluates `density`
#'   on a fixed grid of 1024 exposure values spanning the observed exposure
#'   range padded by its own width on each side; estimators whose conditional
#'   standard deviation is very small relative to that grid step should supply a
#'   closed-form `hdr_threshold` to keep the cutoff well resolved.
#' @param label A single string naming the estimator, stored on the object and
#'   reported by the [check_hdr()] result. Defaults to `"custom"`.
#'
#' @return An `hdr_density` estimator object.
#'
#' @references
#' Bao Y, Schomaker M (2025). Feasible Dose-Response Curves for Continuous
#' Treatments Under Positivity Violations.
#'
#' @examples
#' # A normal-equivalent estimator built through the public contract.
#' estimator <- new_hdr_density(
#'   fit = function(formula, data) {
#'     model <- lm(formula, data = data)
#'     list(model = model, sigma = sigma(model))
#'   },
#'   density = function(state, a, newdata) {
#'     mu <- predict(state$model, newdata = newdata)
#'     dnorm(a, mean = mu, sd = state$sigma)
#'   }
#' )
#' estimator
#'
#' @export
new_hdr_density <- function(
  fit,
  density,
  hdr_threshold = NULL,
  label = "custom"
) {
  if (!is.function(fit)) {
    abort(
      "{.arg fit} must be a function.",
      error_class = "positively_type_error"
    )
  }
  if (!is.function(density)) {
    abort(
      "{.arg density} must be a function.",
      error_class = "positively_type_error"
    )
  }
  if (!is.null(hdr_threshold) && !is.function(hdr_threshold)) {
    abort(
      "{.arg hdr_threshold} must be a function or {.code NULL}.",
      error_class = "positively_type_error"
    )
  }
  if (!is.character(label) || length(label) != 1 || is.na(label)) {
    abort(
      "{.arg label} must be a single string.",
      error_class = "positively_type_error"
    )
  }
  hdr_density(
    fit = fit,
    density = density,
    hdr_threshold = hdr_threshold,
    label = label
  )
}

#' The default normal conditional-density estimator for HDR non-overlap
#'
#' `hdr_density_normal()` builds the default conditional-density estimator used
#' by [check_hdr()] and [check_hdr_seq()]. It fits a linear model of the
#' exposure on the covariates, treats the conditional density as Gaussian with
#' the model's residual standard deviation, and supplies a closed-form HDR
#' threshold.
#'
#' @details
#' The estimator fits `lm(exposure ~ covariates)`, recovering fitted means
#' \eqn{\hat{\mu}(l)} and a homoskedastic residual standard deviation
#' \eqn{\hat{\sigma}}. The conditional density is
#' \eqn{\hat{f}(a \mid l) = \mathrm{dnorm}(a; \hat{\mu}(l), \hat{\sigma})}. Its
#' HDR at mass `mass` is a symmetric interval around \eqn{\hat{\mu}(l)}, so the
#' density cutoff is the closed form
#' \deqn{f_\alpha = \mathrm{dnorm}(z) / \hat{\sigma}, \qquad
#' z = \Phi^{-1}\!\left(\frac{1 + \mathrm{mass}}{2}\right).}
#' Membership of a target `a` in the HDR reduces to the interval test
#' \eqn{|a - \hat{\mu}(l)| \le z\,\hat{\sigma}}, and the non-overlap ratio at
#' `a` is the fraction of fitted means more than \eqn{z\,\hat{\sigma}} from `a`.
#'
#' Because the working model is a single Gaussian, this estimator detects
#' mean-shift support gaps, where a stratum's supported dose moves away from a
#' target, but not multimodal gaps: a hole between two modes of the true
#' conditional density is filled by the fitted normal and reported as
#' supported. Supply a flexible estimator through [new_hdr_density()] when
#' multimodal structure is expected.
#'
#' @return An `hdr_density` estimator object.
#'
#' @references
#' Bao Y, Schomaker M (2025). Feasible Dose-Response Curves for Continuous
#' Treatments Under Positivity Violations.
#'
#' @examples
#' estimator <- hdr_density_normal()
#' estimator
#'
#' @export
hdr_density_normal <- function() {
  fit <- function(formula, data) {
    model <- stats::lm(formula, data = data)
    # stats::sigma() returns the same number as summary(model)$sigma without
    # assembling the coefficient table, and so without the "essentially perfect
    # fit" warning summary.lm() raises when the exposure is close to a
    # deterministic function of the covariates. That warning belongs to
    # inference this estimator never reads, and it is unclassed, so letting it
    # through would put a bare base R warning in front of the user.
    list(model = model, sigma = stats::sigma(model))
  }
  density <- function(state, a, newdata) {
    mu <- stats::predict(state$model, newdata = newdata)
    stats::dnorm(a, mean = mu, sd = state$sigma)
  }
  hdr_threshold <- function(state, newdata, mass) {
    z <- stats::qnorm((1 + mass) / 2)
    rep(stats::dnorm(z) / state$sigma, nrow(newdata))
  }
  new_hdr_density(
    fit = fit,
    density = density,
    hdr_threshold = hdr_threshold,
    label = "normal"
  )
}

# ---- Threshold resolution -------------------------------------------------

#' Resolve the per-row HDR density cutoff for an estimator
#'
#' Returns the density cutoff \eqn{f_\alpha(l_j)} for every row of `newdata`.
#' When the estimator supplies a closed-form `hdr_threshold`, that is used
#' directly; otherwise the cutoff is computed numerically on a fine grid of
#' exposure values (the `hdr_threshold = NULL` fallback of the pluggable
#' contract).
#'
#' @param estimator An `hdr_density` estimator.
#' @param state The fitted state from `estimator@fit()`.
#' @param newdata The covariate rows to evaluate at.
#' @param mass The HDR probability mass.
#' @param exposure The observed exposure vector, used to place the numeric grid.
#' @param call The calling environment, used to build the error's call.
#'
#' @return A numeric vector of cutoffs, one per row of `newdata`.
#' @keywords internal
#' @noRd
hdr_thresholds <- function(
  estimator,
  state,
  newdata,
  mass,
  exposure,
  call = rlang::caller_env()
) {
  if (!is.null(estimator@hdr_threshold)) {
    estimator@hdr_threshold(state, newdata, mass)
  } else {
    numeric_hdr_thresholds(
      estimator,
      state,
      newdata,
      mass,
      exposure,
      call = call
    )
  }
}

#' Numeric-grid HDR density cutoffs
#'
#' Evaluates the estimator's conditional density on a shared fine grid of
#' exposure values, then for each row finds the highest-density cutoff whose
#' region captures probability mass `mass`. The grid spans the observed exposure
#' range padded by its own width on each side so that the conditional density's
#' tails are captured for rows whose mean sits near the edge of the range.
#'
#' @inheritParams hdr_thresholds
#'
#' @return A numeric vector of cutoffs, one per row of `newdata`.
#' @keywords internal
#' @noRd
numeric_hdr_thresholds <- function(
  estimator,
  state,
  newdata,
  mass,
  exposure,
  call = rlang::caller_env()
) {
  n_grid <- 1024L
  span <- diff(range(exposure))
  if (span == 0) {
    span <- 1
  }
  grid <- seq(
    min(exposure) - span,
    max(exposure) + span,
    length.out = n_grid
  )
  # grid_densities[j, k] = f_hat(grid[k] | l_j); vapply fills one column per
  # grid value, so each matrix row is a single covariate profile.
  grid_densities <- vapply(
    grid,
    function(a) estimator@density(state, a, newdata),
    numeric(nrow(newdata))
  )
  row_totals <- rowSums(grid_densities)
  bad <- !is.finite(row_totals) | row_totals <= 0
  n_bad <- sum(bad)
  if (n_bad > 0) {
    abort(
      c(
        "The numeric-grid HDR cutoff is undefined for {n_bad} observation{?s}.",
        x = "The estimator's conditional density is zero or non-finite at every grid value for {cli::qty(n_bad)}{?that observation/those observations}.",
        i = "Supply a closed-form `hdr_threshold` to `new_hdr_density()`, or use a density whose mass falls inside the padded exposure range."
      ),
      error_class = "positively_range_error",
      call = call
    )
  }
  apply(grid_densities, 1L, hdr_cutoff_from_grid, mass = mass)
}

#' The highest-density cutoff capturing a target mass on a density grid
#'
#' Sorts grid densities in decreasing order and returns the density value at
#' which the cumulative captured mass first reaches `mass`. Grid points with
#' density at or above the returned cutoff form the HDR.
#'
#' @param densities A numeric vector of density evaluations on a uniform grid.
#' @param mass The target probability mass.
#'
#' @return A single density cutoff.
#' @keywords internal
#' @noRd
hdr_cutoff_from_grid <- function(densities, mass) {
  ordered <- sort(densities, decreasing = TRUE)
  captured <- cumsum(ordered) / sum(ordered)
  ordered[which(captured >= mass)[1]]
}
