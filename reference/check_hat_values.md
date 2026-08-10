# Diagnose positivity for continuous exposures with hat values

`check_hat_values()` computes the hat-value (leverage) positivity
diagnostic of Moodie and Schulz (2025) for continuous exposures. It
measures how far a set of candidate exposure-covariate combinations sits
from the observed data by the leverage those combinations would exert on
a linear model, then compares the observed leverage profile against a
null in which the exposure is drawn independently of the covariates.

## Usage

``` r
check_hat_values(
  .data,
  .exposure,
  .covariates,
  probs = seq(0.05, 0.95, by = 0.05),
  threshold = 2,
  null_reps = 500,
  null_method = c("permutation", "bootstrap", "gaussian"),
  conf_level = 0.95,
  exposure_type = c("auto", "continuous")
)
```

## Arguments

- .data:

  A data frame.

- .exposure:

  The continuous exposure column, selected with data-masking. Under the
  default `exposure_type = "auto"`, `check_hat_values()` aborts when the
  exposure is detected as binary or categorical; declaring
  `exposure_type = "continuous"` passes any numeric column through the
  type gate.

- .covariates:

  The covariate columns, selected with tidyselect. Numeric, logical,
  factor, and character columns are accepted; a factor or character
  column is expanded internally to one fewer column than it has levels,
  as treatment-contrast indicators when the factor is unordered and as
  polynomial contrasts when it is ordered. The two encodings span the
  same space, so the hat values are the same either way.

- probs:

  A numeric vector of exposure percentiles, each in `[0, 1]`, at which
  candidate points are formed. Defaults to `seq(0.05, 0.95, by = 0.05)`.

- threshold:

  The leverage cutoff multiplier. A candidate is flagged when its hat
  value exceeds `threshold * p / n`. Defaults to `2`.

- null_reps:

  The number of null replicates used to calibrate \\\hat{\phi}\\.
  Defaults to `500`.

- null_method:

  The null-resampling scheme, one of `"permutation"` (permute observed
  exposures, the default), `"bootstrap"` (resample observed exposures
  with replacement), or `"gaussian"` (draw from a matched normal). A
  uniform null is deliberately not offered; see Details.

- conf_level:

  The null quantile, strictly between 0 and 1, against which the
  observed \\\hat{\phi}\\ is compared. Defaults to `0.95`.

- exposure_type:

  One of `"auto"` (detect from the data, the default) or `"continuous"`.
  A supplied type is authoritative and detection is not consulted, so
  `"continuous"` is rejected only when the exposure column is not
  numeric. Declaring it on a numeric column with few distinct values, a
  dose recorded at a handful of milligram levels for instance, is
  exactly the supported use: the unique-value heuristic reads such a
  column as categorical, so under `"auto"` the call aborts. Because
  numeric is the whole of the requirement, a two-valued numeric column
  is then accepted as well, which is the price of an authoritative
  declaration.

## Value

A `hat_values_result` object, an S7 subclass of
[positivity_diagnostic](https://r-causal.github.io/positively/reference/positivity_diagnostic.md).
Its `@results` tibble has one row per candidate point with columns `.id`
(the covariate row index), `prob` (the exposure percentile), `value`
(the candidate exposure value \\d_q\\), `hat_value` (the leverage), and
`high_leverage` (the logical flag). It also carries the scalar
properties `@phi_hat`, `@null_dist`, `@null_quantile`, `@exceeds_null`,
and `@p`, the number of model parameters \\p = q + 2\\ for `q` encoded
design columns, which sets the \\2p/n\\ leverage cutoff.

[`generics::glance()`](https://generics.r-lib.org/reference/glance.html)
returns a one-row tibble with `n` (the sample size), `phi_hat`,
`null_quantile`, `exceeds_null` (kept logical), the counts
`n_high_leverage` and `n_candidates` that `phi_hat` is the ratio of, and
`p`.

## Details

The method (the paper's Box 1) fits the linear design matrix \\M = \[1,
D, X\]\\, where `D` is the continuous exposure and `X` the encoded
covariates, so the number of parameters is \\p = q + 2\\ for `q` encoded
design columns. A numeric or logical covariate contributes one column; a
factor or character covariate contributes one fewer column than it has
levels, as treatment-contrast indicators when the factor is unordered
and as polynomial contrasts when it is ordered. The two encodings span
the same space, so `q` and the hat values are the same either way. For
each requested exposure percentile \\d_q\\ and each observed covariate
row \\x_i\\, it forms the candidate point \\x\_\* = (1, d_q, x_i)\\ and
computes its leverage (hat value) \\h\_\* = x\_\*^\top (M^\top M)^{-1}
x\_\*\\. A candidate is flagged as high leverage when \\h\_\* \>
\mathrm{threshold} \times p / n\\, with the default `threshold` of 2
reproducing the paper's \\2p/n\\ rule. The summary statistic
\\\hat{\phi}\\ is the proportion of candidates that are high leverage.

Because a high value of \\\hat{\phi}\\ can arise from the shape of the
exposure distribution alone, the diagnostic calibrates it against a null
in which the exposure is drawn independently of the covariates. Each of
`null_reps` replicates redraws the exposure, rebuilds the design matrix,
and recomputes \\\hat{\phi}\\ at the same candidate points. The observed
\\\hat{\phi}\\ is compared against the `conf_level` quantile of this
null distribution; when it exceeds that quantile, `exceeds_null` is
`TRUE`.

Three null-resampling schemes are available through `null_method`.
`"permutation"` (the default) permutes the observed exposures,
`"bootstrap"` resamples them with replacement, and `"gaussian"` draws
from a normal distribution matched to the observed exposure mean and
standard deviation. The permutation and bootstrap schemes are empirical:
they preserve the observed exposure distribution exactly. A uniform null
is deliberately not offered, because a uniform draw discards that
distribution: on a skewed exposure it spreads the draws more evenly than
the data are, which inflates the null quantile and costs power against a
real violation.

The per-candidate `high_leverage` flag is a conservative leverage
indicator that feeds \\\hat{\phi}\\, not a calibrated per-unit accept or
reject rule. The absolute value of \\\hat{\phi}\\ depends steeply on the
candidate percentile grid, so the null comparison, not the raw
magnitude, carries the diagnostic signal.

## References

Moodie EEM, Schulz J (2025). A Simple Diagnostic for the Positivity
Assumption for Continuous Exposures. *Statistics in Medicine*,
44:e70194. [doi:10.1002/sim.70194](https://doi.org/10.1002/sim.70194)

## Examples

``` r
set.seed(1)
n <- 200
x1 <- rnorm(n)
# The exposure depends on x1, so extreme doses are implausible for some
# covariate values.
dose <- rnorm(n, mean = x1)
df <- data.frame(dose = dose, x1 = x1)

# null_reps is small here to keep the example fast.
result <- check_hat_values(df, dose, x1, null_reps = 50)
result
#> 
#> ── Hat values ──────────────────────────────────────────────────────────────────
#> Exposure: "dose" (continuous)
#> Observations: 200
#> Null: permutation (50 reps), cutoff 2p/n
#> phi-hat: 0.19
#> Null 0.95 quantile: 0.052
#> Exceeds null: TRUE
#> High-leverage candidates: 723 of 3800
```
