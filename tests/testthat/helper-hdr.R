# A structurally faithful normal-equivalent density estimator built through the
# public constructor. It mirrors the default estimator's geometry so the
# pluggable path can be checked against the closed-form path. With
# with_threshold = FALSE the cutoff is left NULL to request the numeric grid
# fallback.
make_normal_like <- function(with_threshold = FALSE) {
  fit <- function(formula, data) {
    model <- stats::lm(formula, data = data)
    list(model = model, sigma = stats::sigma(model))
  }
  density <- function(state, a, newdata) {
    mu <- stats::predict(state$model, newdata = newdata)
    stats::dnorm(a, mean = mu, sd = state$sigma)
  }
  if (!with_threshold) {
    return(new_hdr_density(fit = fit, density = density))
  }
  hdr_threshold <- function(state, newdata, mass) {
    z <- stats::qnorm((1 + mass) / 2)
    rep(stats::dnorm(z) / state$sigma, nrow(newdata))
  }
  new_hdr_density(fit = fit, density = density, hdr_threshold = hdr_threshold)
}
