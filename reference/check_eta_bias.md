# Diagnose positivity bias with the parametric bootstrap (ETA.Bias)

`check_eta_bias()` implements the ETA.Bias parametric-bootstrap
diagnostic of Petersen et al. (2012) for binary, categorical, and
continuous point exposures. It estimates the bias that a chosen causal
estimator would incur under the data's own fitted treatment and outcome
mechanisms, isolating the component of bias that stems from positivity
violations, weight truncation, and finite-sample sparsity.

## Usage

``` r
check_eta_bias(
  .data,
  .exposure,
  .outcome,
  .covariates,
  estimator = c("ipw", "gcomp", "aipw"),
  exposure_formula = NULL,
  outcome_formula = NULL,
  outcome_type = c("auto", "continuous", "binary"),
  truncation = NULL,
  truncation_grid = NULL,
  n_boot = 500,
  error_dist = c("normal", "empirical"),
  exposure_type = c("auto", "binary", "categorical", "continuous"),
  reference_level = NULL,
  msm_formula = NULL
)
```

## Arguments

- .data:

  A data frame.

- .exposure:

  The exposure column, selected with data-masking. A binary exposure has
  exactly two distinct values; a categorical exposure has as many levels
  as it has distinct values, and every non-reference level is contrasted
  against the reference; a continuous exposure has no levels, and is
  summarized by the working model `msm_formula` sets.

- .outcome:

  The outcome column, selected with data-masking.

- .covariates:

  The covariate columns, selected with tidyselect. Numeric, logical,
  factor, and character covariates are all supported. The fitted
  treatment and outcome models expand factor and character covariates
  into indicator terms, and the bootstrap resamples covariate rows, so
  factor levels are preserved throughout.

- estimator:

  The causal estimator to diagnose, either one of the built-in
  estimators or one of your own.

  The built-in estimators are named by a string: `"ipw"` (inverse
  probability weighting, matching
  [`propensity::ipw()`](https://r-causal.github.io/causalgenerics/reference/ipw.html)),
  `"gcomp"` (G-computation), or `"aipw"` (augmented, doubly robust). The
  augmented estimator's correction is a sum over the exposure levels, so
  it is an error for a continuous exposure.

  An estimator of your own is a function, or a length-one list whose
  name labels that function. A bare function is labeled `"custom"`. The
  label is what `@estimator`, `@params`, and
  [`generics::glance()`](https://generics.r-lib.org/reference/glance.html)
  report; the function itself is not kept on the result.

  The function is called once per bootstrap draw and receives one
  argument, the data frame of that draw: the resampled covariate rows,
  the exposure drawn from the fitted treatment mechanism, and the
  outcome drawn from the fitted outcome regression, under the column
  names they were selected by. A discrete exposure arrives as a factor
  over the levels the observed exposure holds, and a continuous one as
  the drawn numeric vector. The function must return one number per
  estimand term: named by the terms, in which case the values are read
  by name, or unnamed, in which case they are read in term order. Every
  draw's return is read on its own, so the two forms may be mixed from
  one draw to the next, and names that are anything other than the
  estimand terms are an error wherever they appear.

  The function has to aim at the estimand the run's terms define, which
  is the set of contrasts against `reference_level` for a discrete
  exposure and the coefficients of `msm_formula` for a continuous one.
  The truth is the same G-computation target whichever estimator is
  diagnosed, so an estimator aimed at anything else reports an estimand
  mismatch as part of its bias.

  `truncation` and `truncation_grid` bound a fitted nuisance model the
  function never receives, so supplying either alongside one is an
  error. A draw the function raises an error on is dropped exactly as a
  non-finite estimate is, and a run whose every draw was dropped is an
  error.

- exposure_formula:

  An optional model formula for the treatment mechanism. `NULL` (the
  default) fits a main-effects model of the exposure on the covariates,
  logistic at two levels, a multinomial logit past two, and a normal
  linear model for a continuous exposure.

- outcome_formula:

  An optional model formula for the outcome regression. `NULL` (the
  default) fits a main-effects model of the outcome on the exposure and
  covariates. Terms that recompute a data-dependent basis from the
  exposure column (for example
  [`poly()`](https://rdrr.io/r/stats/poly.html) of the exposure) are not
  supported, because the counterfactual predictions rebuild the design
  with the exposure set to each of its levels in turn.

- outcome_type:

  One of `"auto"` (detect from the outcome, the default),
  `"continuous"`, or `"binary"`.

- truncation:

  An optional length-two numeric vector `c(lower, upper)` bounding the
  fitted propensity for a single run. `NULL` (the default) applies no
  truncation. A pair of bounds describes the one-dimensional simplex of
  a two-level exposure, so it applies to binary exposures only and is an
  error past two levels or on a continuous exposure.

- truncation_grid:

  An optional numeric vector, each entry in `[0, 1/k)` for an exposure
  with `k` levels, which is `[0, 0.5)` for a binary exposure, and in
  `[0, 1)` for a continuous one. All levels share one set of bootstrap
  draws, and the grid overrides `truncation`. What an entry means
  depends on the exposure type: at two levels it is a lower bound and
  becomes the swept truncation level `c(lower, 1 - lower)`; past two
  levels it raises every fitted probability below it and renormalizes
  each row, leaving no upper bound to apply or report; and for a
  continuous exposure it is a quantile level `g` that caps the
  stabilized weight at its `1 - g` quantile, read once from the observed
  fit and held fixed across the bootstrap draws.

- n_boot:

  The number of bootstrap datasets. Defaults to `500`. Must be at least
  2, since the Monte Carlo standard error needs two draws.

- error_dist:

  The bootstrap error model for continuous outcomes, one of `"normal"`
  (add Gaussian noise matched to the residual standard deviation, the
  default) or `"empirical"` (resample residuals). Ignored for binary
  outcomes.

- exposure_type:

  One of `"auto"` (detect from the data, the default), `"binary"`,
  `"categorical"`, or `"continuous"`. A supplied type is authoritative
  and detection is not consulted, so it is rejected only when the
  exposure column cannot carry it: `"binary"` needs exactly two distinct
  values, `"continuous"` needs a numeric column, and `"categorical"`
  needs nothing beyond the two levels every discrete run needs. Within
  the discrete types the run is decided by the number of levels, so
  declaring `"categorical"` on a two-level column returns the binary run
  unchanged.

- reference_level:

  The exposure level every contrast is taken against. `NULL` (the
  default) takes the first level, which is the lower of the two values
  for a binary exposure and the first factor level otherwise. A
  continuous exposure has no levels to contrast, so supplying it there
  is an error.

- msm_formula:

  An optional one-sided formula for the marginal structural model a
  continuous exposure is summarized by, for example `~ a + I(a^2)`.
  `NULL` (the default) fits a model linear in the exposure, whose one
  coefficient is named after the exposure column. A discrete estimand is
  a set of contrasts and imposes no working model, so supplying it there
  is an error.

## Value

An `eta_bias_result` object, an S7 subclass of
[positivity_diagnostic](https://r-causal.github.io/positively/reference/positivity_diagnostic.md).
Its `@results` tibble has one row per estimand term and truncation
level, with the columns `term`, `truncation_lower`, `truncation_upper`,
`bias` (ETA.Bias), `mc_se` (the Monte Carlo standard error), and
`boot_mean` (the mean bootstrap estimate), in that order. `term` names
the estimand contrast as `<level> - <reference level>`, so an exposure
coded 0/1 has the single term `"1 - 0"`; for a continuous exposure it
names a working-model coefficient instead. A continuous run has no lower
cap to place on a density ratio, so `truncation_lower` is `NA_real_` and
`truncation_upper` holds the weight cap, which is `Inf` where no cap was
applied. The object also carries the properties `@estimator`, `@truth`
(a named numeric vector holding one truth per term), and
`@boot_estimates` (a list of bootstrap-estimate vectors, one per row of
`@results`).

[`generics::glance()`](https://generics.r-lib.org/reference/glance.html)
returns a one-row tibble with `n` (the sample size), `estimator`, and
`n_boot`. A run that dropped a bootstrap draw adds `n_boot_used`, the
number of draws retained, counted as the smallest number any reported
row was formed from; a run that dropped none omits the column. A
single-term estimand adds `truth`; an estimand of more than one term has
no single truth, so it reports `n_terms` in its place. A one-row result
then adds `bias` and `mc_se`, while a result of more rows replaces those
two with `n_levels`, `bias_min`, and `bias_max`. The row count is what
decides this, not the argument used, so a length-one `truncation_grid`
takes the one-row form.

## Details

The diagnostic (the paper's Section 4.1) treats the fitted
data-generating mechanism as a known truth and asks how a candidate
estimator behaves under it. The algorithm has three steps: fit the
treatment and outcome mechanisms, declare the fitted mechanism the
truth, then resample from it and refit. The exposure type changes what
each step does, and nothing about the order they come in.

First, it fits the treatment mechanism \\g_n\\ and the outcome
regression \\\bar{Q}\_n(A, W)\\, a main-effects linear or logistic model
of the outcome on the exposure and the covariates. Factor and character
covariates enter both models as indicator terms, exactly as
[`stats::glm()`](https://rdrr.io/r/stats/glm.html) expands them. The
treatment mechanism is a main-effects model of the exposure on the
covariates in the form its type calls for: the logistic model \\g_n(1
\mid W)\\ at two levels; the multinomial logit \\g_n(a \mid W)\\ past
two, which requires the nnet package; and for a continuous exposure the
normal linear model \\A \mid W \sim \mathcal{N}(m_n(W), \sigma_n^2)\\,
whose spread is the pooled maximum-likelihood residual standard
deviation \\\sigma_n = \sqrt{n^{-1} \sum_i (A_i - m_n(W_i))^2}\\.

The `truth` is the target the fitted mechanism implies, which is
available exactly because that fit is the truth of the bootstrap
data-generating mechanism. It depends on the outcome model alone, so it
is identical across estimators and constant across truncation levels. A
discrete exposure is summarized by the saturated factor model, whose
coefficients are the G-computation contrasts against the reference level
\\a_0\\, \\\psi_a = \frac{1}{n}\sum_i \[\bar{Q}\_n(a, W_i) -
\bar{Q}\_n(a_0, W_i)\]\\, one per non-reference level \\a\\. A saturated
model restricts nothing, so that summary imposes nothing; an exposure
coded 0/1 is the two-level case and its one contrast is \\\psi_1\\. A
continuous exposure has no levels to saturate. Its counterfactual mean
curve \\\frac{1}{n}\sum_i \bar{Q}\_n(a, W_i)\\ is evaluated over the
unique deciles of the observed exposure and projected by least squares
onto the working model `msm_formula` sets, and the truth is the
coefficients of that projection of the observed fit. Equal weights on a
decile grid weight the projection by the marginal density of the
exposure.

Second, it draws `n_boot` bootstrap datasets. Each resamples the
covariate rows with replacement, draws a bootstrap exposure \\A^\*\\
from the untruncated fitted treatment mechanism at the resampled
covariates, and draws a bootstrap outcome from \\\bar{Q}\_n(A^\*,
W^\*)\\: for continuous outcomes by adding normal or empirical-residual
error, and for binary outcomes by a Bernoulli draw. The exposure draw is
\\A^\* \sim \mathrm{Bernoulli}(g_n(1 \mid W^\*))\\ at two levels, a
categorical draw over the fitted probabilities \\g_n(a \mid W^\*)\\ past
two, and \\A^\* \sim \mathcal{N}(m_n(W^\*), \sigma_n^2)\\ for a
continuous exposure. Truncation is a property of the estimator, not of
the data-generating mechanism, so the bootstrap exposure is always drawn
from the untruncated mechanism.

Third, it refits both nuisance models on each bootstrap dataset and
applies the candidate estimator, truncating what that estimator divides
by. `ETA.Bias` is the mean of the bootstrap estimates minus `truth`, and
its Monte Carlo standard error is the standard deviation of the
bootstrap estimates divided by the square root of the number of retained
draws. A bootstrap draw whose estimate is non-finite carries no
information, and is dropped with a warning before the summaries are
formed.

What a `truncation_grid` entry \\g\\ truncates also follows the exposure
type. At two levels it bounds the fitted propensity into \\\[g, 1 -
g\]\\: the simplex is one-dimensional there, so a lower bound on one
level is an upper bound on the other and the rule is two-sided. Past two
levels it raises every fitted probability below \\g\\ up to \\g\\ and
renormalizes each row to sum to one, which is the matrix method of
[`propensity::ps_trunc()`](https://r-causal.github.io/propensity/reference/ps_trunc.html)
under `method = "ps"`. A lower pin on every level already bounds every
weight denominator, so no upper bound is applied and none is reported. A
continuous exposure divides by a density rather than by a probability,
and \\g\\ caps the stabilized weight \\f_n(A) / f_n(A \mid W)\\ at its
\\1 - g\\ quantile, read once from the observed fit and held fixed
across the bootstrap draws.

The three estimators behave differently under positivity violations.
G-computation ignores the treatment mechanism entirely, so its ETA.Bias
is Monte Carlo noise around zero and is flat across a truncation grid.
Inverse probability weighting divides by the fitted mechanism, so
extreme probabilities and heavy density ratios inflate its ETA.Bias, and
tightening the truncation bound trades bias for variance. The augmented
(doubly robust) estimator stays near zero in bias while carrying more
variance than G-computation, and it is defined for a discrete estimand
alone.

When `truncation_grid` is supplied, every grid point reuses one shared
set of bootstrap draws and nuisance refits, so the sweep isolates the
effect of the truncation bound alone: the inverse-probability-weighting
bias rises and its bootstrap spread shrinks as the bound tightens.

A continuous estimand is defined by the working model rather than
approximated by it, so two analyses of the same data under different
`msm_formula` values are aimed at different targets and report different
biases. That is a property of the estimand and not an inconsistency: the
bias reported is the bias of estimating the working model the analysis
actually uses.

ETA.Bias captures only the positivity, truncation, and sparsity
component of bias and excludes model-misspecification bias by
construction. It is a red flag for a fitted mechanism, not a bias
correction.

## References

Petersen ML, Porter KE, Gruber S, Wang Y, van der Laan MJ (2012).
Diagnosing and responding to violations in the positivity assumption.
*Statistical Methods in Medical Research*, 21(1):31–54.
[doi:10.1177/0962280210386207](https://doi.org/10.1177/0962280210386207)

## Examples

``` r
set.seed(1)
n <- 300
x1 <- rnorm(n)
x2 <- rnorm(n)
a <- rbinom(n, 1, plogis(x1 + x2))
y <- a + x1 + x2 + rnorm(n)
df <- data.frame(a = a, y = y, x1 = x1, x2 = x2)

# n_boot is small here to keep the example fast.
check_eta_bias(df, a, y, c(x1, x2), estimator = "ipw", n_boot = 25)
#> ℹ Treating `.exposure` as binary
#> 
#> ── ETA bias ────────────────────────────────────────────────────────────────────
#> Exposure: "a" (binary)
#> Observations: 300
#> Estimator: ipw
#> Bootstrap draws: 25
#> Truth: 1.113
#> ETA.Bias: 0.003 (MC SE 0.052)

# An estimator of your own reads one bootstrap dataset per draw. The exposure
# arrives as a factor of its own levels, so the fitted coefficient is named
# for the level it contrasts.
by_hand <- function(.data) {
  coef(lm(y ~ a + x1 + x2, data = .data))[["a1"]]
}
check_eta_bias(
  df,
  a,
  y,
  c(x1, x2),
  estimator = list(gformula = by_hand),
  n_boot = 25
)
#> ℹ Treating `.exposure` as binary
#> 
#> ── ETA bias ────────────────────────────────────────────────────────────────────
#> Exposure: "a" (binary)
#> Observations: 300
#> Estimator: gformula
#> Bootstrap draws: 25
#> Truth: 1.113
#> ETA.Bias: 0.014 (MC SE 0.027)
```
