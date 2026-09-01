# Diagnose positivity with effective data points

`check_edp()` computes the effective-data-points (EDP) positivity
diagnostic of Ring and Schomaker (2026) for binary, categorical, and
continuous exposures. For every observation moved to a candidate
intervention value, it sums a product kernel over the sample to measure
how much observed support surrounds that intervened-on point.

## Usage

``` r
check_edp(
  .data,
  .exposure,
  .covariates,
  .outcome_covariates = NULL,
  .treatment_covariates = NULL,
  values = NULL,
  variant = c("data", "estimator"),
  kernel = c("gaussian"),
  bw_exposure = NULL,
  bw_covariates = NULL,
  categorical_similarity = 0,
  exposure_type = c("auto", "binary", "categorical", "continuous")
)
```

## Arguments

- .data:

  A data frame.

- .exposure:

  The exposure column, selected with data-masking.

- .covariates:

  The covariate columns, selected with tidyselect.

- .outcome_covariates:

  The outcome-model covariate columns for the estimator variant,
  selected with tidyselect. Defaults to `.covariates`.

- .treatment_covariates:

  The treatment-model covariate columns for the estimator variant,
  selected with tidyselect. Defaults to `.covariates`.

- values:

  The interventions, one per element. `NULL` (the default) uses every
  level of a factor exposure and every observed value of a binary or
  character exposure, and the deciles of the observed exposure for
  continuous exposures.

  A vector holds one intervention value \\a^\*\\ per element, the value
  every observation is moved to. A list holds the same single values
  and, for a continuous exposure, functions of the observed exposure:
  `list(3, 5, lower = function(.x) .x - 5)`. A function is called once
  on the whole exposure vector and returns one intervened exposure
  \\d(a_i)\\ per observation, or a single value used for every
  observation. Functions need a continuous exposure; one supplied for a
  binary or categorical exposure is an error.

  Every intervention carries a label, which names it in the results and
  along a plot's axis. In a vector, every element is labeled by
  [`as.character()`](https://rdrr.io/r/base/character.html) of the value
  itself and any names the vector carries are ignored. In a list, a
  named element is labeled by its name, an unnamed single value by
  [`as.character()`](https://rdrr.io/r/base/character.html) of the value
  itself, and an unnamed function is an error, since it has nothing else
  to be called by. Two elements of a list cannot share a label.

- variant:

  One of `"data"` (report a single `edp` per point, the default) or
  `"estimator"` (report `edp_outcome`, `edp_treatment`, and
  `ideal_weight`).

- kernel:

  The continuous-dimension kernel. Only `"gaussian"` is currently
  available; categorical dimensions always use the match kernel.

- bw_exposure:

  The continuous-exposure half-distance. `NULL` (the default) uses
  `0.5 * sd(exposure)`. A half-distance of `0` counts exact matches.

- bw_covariates:

  The numeric-covariate half-distance, applied to every numeric
  covariate. `NULL` (the default) uses `1 * sd` per covariate.

- categorical_similarity:

  The kernel value in `[0, 1]` for non-matching categories, shared by
  categorical exposures and factor or character covariates. Defaults to
  `0`. A positive value masks categorical violations.

- exposure_type:

  One of `"auto"` (detect from the data, the default), `"binary"`,
  `"categorical"`, or `"continuous"`. A supplied type is authoritative
  and detection is not consulted, so it is rejected only when the
  exposure column cannot carry it: `"continuous"` needs a numeric column
  and `"binary"` needs exactly two distinct values, while
  `"categorical"` asks nothing of the column.

## Value

An `edp_result` object, an S7 subclass of
[positivity_diagnostic](https://r-causal.github.io/positively/reference/positivity_diagnostic.md).
Its `@results` tibble has one row per observation and intervention. The
data variant carries columns `.id` (the observation row index),
`intervention` (the intervention label), `value` (the intervened
exposure, the constant \\a^\*\\ for a single value and \\d(a_i)\\ for a
function), and `edp`. The estimator variant replaces `edp` with
`edp_outcome`, `edp_treatment`, and `ideal_weight`. It also carries the
`@variant` and `@bandwidths` properties.

[`generics::glance()`](https://generics.r-lib.org/reference/glance.html)
returns a one-row tibble with `n` (the sample size), `variant`,
`n_values` (the number of interventions), and the observed range of
every measure the variant reports. That is `edp_min` and `edp_max` for
the data variant, and `edp_outcome_min`, `edp_outcome_max`,
`edp_treatment_min`, `edp_treatment_max`, `ideal_weight_min`, and
`ideal_weight_max` for the estimator variant.

## Details

For each observation \\i\\ the diagnostic forms the intervened-on point
\\o^\*\_i = (l_i, a^\*)\\, holding the covariates at their observed
values and setting the exposure to an intervention value \\a^\*\\. The
effective data points at that point are \$\$\mathrm{EDP}\_i(a^\*) =
\sum\_{j=1}^{n} k(o_j, o^\*\_i),\$\$ a sum over the whole sample of a
product kernel that multiplies one factor per dimension. Each factor
lies in \\\[0, 1\]\\, so \\0 \le \mathrm{EDP} \le n\\. A continuous
dimension with half-distance \\h\\ contributes \\0.5^{(\delta / h)^2}\\,
which equals `1` at \\\delta = 0\\ and exactly `0.5` at \\\|\delta\| =
h\\; this is a reparameterized Gaussian kernel. A categorical dimension
contributes `1` on a match and `categorical_similarity` otherwise.
Numeric covariates and continuous exposures are continuous dimensions;
factor and character covariates, together with binary and categorical
exposures, are categorical dimensions.

An intervention need not be a single value every observation is moved
to. A function supplied through `values` is a rule \\d\\ applied to the
observed exposure, so observation \\i\\ is moved to \\d(a_i)\\ and the
intervened-on point is \\o^\*\_i = (l_i, d(a_i))\\. Lowering every
observed dose by five units is `function(.x) .x - 5`. EDP is then read
per natural observation under the shift rather than at one common
target; the algorithm is unchanged, and only the intervened-on points
come from the rule. Because those points differ across observations, the
exposure factor of the kernel compares every observed exposure against a
target that moves with the observation rather than against a single
\\a^\*\\. A run therefore holds one \\n \times n\\ matrix of
exposure-kernel values per function intervention, materialized before
any summing and held for the duration of the run.

The half-distances follow the paper's rule of thumb: `1 * sd` for each
numeric covariate and `0.5 * sd` for a continuous exposure. These are
marginal scales, so a mid-support gap in a strongly multimodal exposure
can be masked by an inflated standard deviation; supply `bw_exposure` or
`bw_covariates` directly when a robust scale is preferable. As the
half-distances grow without bound every kernel factor approaches `1` and
EDP approaches `n`, which renders the diagnostic uninformative; at a
half-distance of `0` a continuous dimension counts only exact matches.

EDP is a product kernel, so adding dimensions can only lower it. It is
meaningful relative across observations and strata for a fixed covariate
set, never against a universal threshold; compare a stratum against
another stratum or a candidate value against another candidate value
rather than against a fixed cutoff.

The estimator variant reframes EDP around the two nuisance models of a
causal estimator. `edp_outcome` conditions on the exposure and the
outcome-model covariates \\(A, L)\\; `edp_treatment` conditions on the
treatment-model covariates \\(L)\\ alone. Because the outcome measure
carries the extra exposure dimension, `edp_outcome <= edp_treatment`
whenever the two covariate sets are equal. The reported `ideal_weight`
is `edp_treatment / edp_outcome`, proportional to the inverse of the
estimated treatment density at the intervention, so it grows where
treatment support is thin. Where a target has no outcome-model support,
`edp_outcome` is zero and `ideal_weight` is infinite, the honest value
of one over an estimated density of zero.

## References

Ring C, Schomaker M (2026). A Diagnostic to Find and Help Combat
Stochastic Positivity Issues, with a Focus on Continuous Treatments.

## Examples

``` r
set.seed(1)
n <- 100
x1 <- rnorm(n)
dose <- rnorm(n, mean = x1)
df <- data.frame(dose = dose, x1 = x1)

result <- check_edp(df, dose, x1, values = c(0, 1), exposure_type = "continuous")
result
#> 
#> ── Effective data points ───────────────────────────────────────────────────────
#> Exposure: "dose" (continuous)
#> Observations: 100
#> Variant: data
#> Intervention values: 2
#> edp range: 0.793 to 26.809

# An intervention can also be a function of the observed exposure, read as a
# shift away from the dose each observation received
shifted <- check_edp(
  df,
  dose,
  x1,
  values = list(0, 1, lower = function(.x) .x - 1),
  exposure_type = "continuous"
)
tidy(shifted)
#> # A tibble: 300 × 4
#>      .id intervention value   edp
#>    <int> <chr>        <dbl> <dbl>
#>  1     1 0                0  17.9
#>  2     2 0                0  26.7
#>  3     3 0                0  14.6
#>  4     4 0                0  11.4
#>  5     5 0                0  26.8
#>  6     6 0                0  14.8
#>  7     7 0                0  26.1
#>  8     8 0                0  24.0
#>  9     9 0                0  25.5
#> 10    10 0                0  22.7
#> # ℹ 290 more rows
```
