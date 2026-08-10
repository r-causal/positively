# Diagnose positivity bias with the parametric bootstrap (ETA.Bias)

`check_eta_bias()` implements the ETA.Bias parametric-bootstrap
diagnostic of Petersen et al. (2012) for binary point exposures. It
estimates the bias that a chosen causal estimator would incur under the
data's own fitted treatment and outcome mechanisms, isolating the
component of bias that stems from positivity violations, weight
truncation, and finite-sample sparsity.

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
  exposure_type = c("auto", "binary")
)
```

## Arguments

- .data:

  A data frame.

- .exposure:

  The binary exposure column, selected with data-masking.
  `check_eta_bias()` aborts unless the exposure has exactly two distinct
  values, whether that type is detected or declared through
  `exposure_type`.

- .outcome:

  The outcome column, selected with data-masking.

- .covariates:

  The covariate columns, selected with tidyselect. Numeric, logical,
  factor, and character covariates are all supported. The fitted
  treatment and outcome models expand factor and character covariates
  into indicator terms, and the bootstrap resamples covariate rows, so
  factor levels are preserved throughout.

- estimator:

  The causal estimator to diagnose, one of `"ipw"` (inverse probability
  weighting, matching
  [`propensity::ipw()`](https://r-causal.github.io/propensity/reference/ipw.html)),
  `"gcomp"` (G-computation), or `"aipw"` (augmented, doubly robust).

- exposure_formula:

  An optional model formula for the treatment mechanism. `NULL` (the
  default) fits a main-effects logistic model of the exposure on the
  covariates.

- outcome_formula:

  An optional model formula for the outcome regression. `NULL` (the
  default) fits a main-effects model of the outcome on the exposure and
  covariates. Terms that recompute a data-dependent basis from the
  exposure column (for example
  [`poly()`](https://rdrr.io/r/stats/poly.html) of the exposure) are not
  supported, because the counterfactual predictions rebuild the design
  with the exposure set to 0 or 1.

- outcome_type:

  One of `"auto"` (detect from the outcome, the default),
  `"continuous"`, or `"binary"`.

- truncation:

  An optional length-two numeric vector `c(lower, upper)` bounding the
  fitted propensity for a single run. `NULL` (the default) applies no
  truncation.

- truncation_grid:

  An optional numeric vector of lower bounds, each in `[0, 0.5)`. Each
  becomes a swept truncation level `c(lower, 1 - lower)`, and all levels
  share one set of bootstrap draws. Overrides `truncation`.

- n_boot:

  The number of bootstrap datasets. Defaults to `500`. Must be at least
  2, since the Monte Carlo standard error needs two draws.

- error_dist:

  The bootstrap error model for continuous outcomes, one of `"normal"`
  (add Gaussian noise matched to the residual standard deviation, the
  default) or `"empirical"` (resample residuals). Ignored for binary
  outcomes.

- exposure_type:

  One of `"auto"` (detect from the data, the default) or `"binary"`. A
  supplied type is authoritative and detection is not consulted, so
  `"binary"` is rejected only when the exposure column does not have
  exactly two distinct values. Declaring it unlocks no call that
  `"auto"` would refuse, since detection reads any two-valued column as
  binary; it is accepted so that every diagnostic takes the same
  argument.

## Value

An `eta_bias_result` object, an S7 subclass of
[positivity_diagnostic](https://r-causal.github.io/positively/reference/positivity_diagnostic.md).
Its `@results` tibble has one row per truncation level with columns
`truncation_lower`, `truncation_upper`, `bias` (ETA.Bias), `mc_se` (the
Monte Carlo standard error), and `boot_mean` (the mean bootstrap
estimate). It also carries the properties `@estimator`, `@truth`, and
`@boot_estimates` (a list of bootstrap-estimate vectors, one per row of
`@results`).

[`generics::glance()`](https://generics.r-lib.org/reference/glance.html)
returns a one-row tibble with `n` (the sample size), `estimator`,
`n_boot`, and `truth`. A run at a single truncation level adds `bias`
and `mc_se`. A run over more than one level reports one bias per level,
so it replaces those two with `n_levels`, `bias_min`, and `bias_max`.
The level count is what decides this, not the argument used, so a
length-one `truncation_grid` takes the single-level form.

## Details

The diagnostic (the paper's Section 4.1) treats the fitted
data-generating mechanism as a known truth and asks how a candidate
estimator behaves under it. The algorithm has three steps.

First, it fits the treatment mechanism \\g_n(1 \mid W)\\, a main-effects
logistic model of the exposure on the covariates, and the outcome
regression \\\bar{Q}\_n(A, W)\\, a main-effects linear or logistic model
of the outcome on the exposure and covariates. Factor and character
covariates enter these models as indicator terms, exactly as
[`stats::glm()`](https://rdrr.io/r/stats/glm.html) expands them. The
target of inference is a single scalar, the G-computation estimate on
the observed data, \\\psi = \frac{1}{n}\sum_i \[\bar{Q}\_n(1, W_i) -
\bar{Q}\_n(0, W_i)\]\\. This `truth` depends on the outcome model alone.
It is identical across estimators and constant across truncation levels.

Second, it draws `n_boot` bootstrap datasets. Each resamples the
covariate rows with replacement, draws a bootstrap exposure \\A^\* \sim
\mathrm{Bernoulli}(g_n(1 \mid W^\*))\\ from the untruncated fitted
propensity, and draws a bootstrap outcome from \\\bar{Q}\_n(A^\*,
W^\*)\\: for continuous outcomes by adding normal or empirical-residual
error, and for binary outcomes by a Bernoulli draw. Truncation is a
property of the estimator, not of the data-generating mechanism, so the
bootstrap exposure is always drawn from the untruncated propensity.

Third, it refits both nuisance models on each bootstrap dataset and
applies the candidate estimator, bounding the fitted propensity into the
truncation interval where the estimator uses it. `ETA.Bias` is the mean
of the bootstrap estimates minus `truth`, and its Monte Carlo standard
error is the standard deviation of the bootstrap estimates divided by
the square root of the number of retained draws. A bootstrap draw whose
estimate is non-finite, which can happen when a resampled exposure lands
in a single arm or a refit propensity reaches exactly 0 or 1, is dropped
with a warning before the summaries are formed.

The three estimators behave differently under positivity violations.
G-computation ignores the propensity entirely, so its ETA.Bias is Monte
Carlo noise around zero and is flat across a truncation grid. Inverse
probability weighting divides by the fitted propensity, so extreme
scores inflate its ETA.Bias, and tightening the truncation bound trades
bias for variance. The augmented (doubly robust) estimator stays near
zero in bias while carrying more variance than G-computation.

When `truncation_grid` is supplied, every grid point reuses one shared
set of bootstrap draws and nuisance refits, so the sweep isolates the
effect of the truncation bound alone: the inverse-probability-weighting
bias rises and its bootstrap spread shrinks as the bound tightens.

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
#> 
#> ── ETA bias ────────────────────────────────────────────────────────────────────
#> Exposure: "a" (binary)
#> Observations: 300
#> Estimator: ipw
#> Bootstrap draws: 25
#> Truth: 1.113
#> ETA.Bias: 0.003 (MC SE 0.052)
```
