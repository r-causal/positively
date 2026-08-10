# Summarize density ratios for positivity diagnosis

`check_density_ratios()` summarizes user-supplied density ratios, the
Radon-Nikodym derivatives that arise as weights in modified treatment
policy and longitudinal causal analyses. It is a bring-your-own-numbers
diagnostic: a numeric vector is treated as a point treatment, an n-by-T
matrix as a time-varying treatment whose cumulative products across time
points are summarized alongside the per-time ratios, and an lmtp fit is
read for its density-ratio component.

## Usage

``` r
check_density_ratios(ratios, ...)

# S3 method for class 'numeric'
check_density_ratios(
  ratios,
  probs = c(0.5, 0.9, 0.95, 0.99, 1),
  thresholds = c(10, 50),
  ...
)

# S3 method for class 'matrix'
check_density_ratios(
  ratios,
  probs = c(0.5, 0.9, 0.95, 0.99, 1),
  thresholds = c(10, 50),
  ...
)
```

## Arguments

- ratios:

  The density ratios. A numeric vector for a point treatment (a
  [propensity::psw](https://r-causal.github.io/propensity/reference/psw.html)
  object is accepted and read through
  [`vctrs::vec_data()`](https://vctrs.r-lib.org/reference/vec_data.html)),
  an n-by-T numeric matrix whose columns are time points for a
  time-varying treatment, or a fitted lmtp object whose density-ratio
  component is read out. Ratios must be non-negative, finite, and free
  of missing values.

- ...:

  Passed to methods.

- probs:

  A numeric vector of quantile probabilities in `[0, 1]`. Defaults to
  `c(0.5, 0.9, 0.95, 0.99, 1)`; the `1` quantile is reported as the
  maximum.

- thresholds:

  A numeric vector of positive exceedance thresholds. Defaults to
  `c(10, 50)`.

## Value

A `density_ratios_result` object, an S7 subclass of
[positivity_diagnostic](https://r-causal.github.io/positively/reference/positivity_diagnostic.md).
Its `@results` tibble is long, with columns `time` (the time point, `1`
for a point treatment), `statistic` (the summary name), and `value`. It
also carries the raw per-time ratios as the list property `@ratios`,
with one element per time point.

[`generics::glance()`](https://generics.r-lib.org/reference/glance.html)
returns a one-row tibble with `n` (the sample size), `n_times`, and the
headline statistics of the cumulative ratios at the final time point:
`mean`, `max`, `prop_zero`, and `ess_fraction`. For a point treatment
those are the statistics of the supplied ratios.

## Details

A density ratio is \\r = g^\*(a \mid h) / g(a \mid h)\\, the
intervention treatment law relative to the observed treatment law. It is
the weight that reweights the observed data to the interventional world,
and it is large exactly where the observed conditional treatment density
is small, that is, where positivity is threatened. Petersen and
colleagues (2012) recommend inspecting the distribution of these
weights, because near violations of positivity surface as a heavy upper
tail, while values of the observed density close to zero for some
treatment level signal a violation directly.

Because a density ratio is a Radon-Nikodym derivative, its mean is one
by construction whenever the intervention law is absolutely continuous
with respect to the observed law. The mean is therefore a specification
check, not a violation detector: well-behaved weights alone do not
guarantee positivity. The signal instead lives in the upper tail (the
quantiles, the maximum, and the proportion of ratios above the
thresholds) and in the Kish effective sample size, which collapses as a
few large weights come to dominate. A mean visibly below one is itself
informative: it marks a structural violation in which intervention
probability mass leaks outside the observed support, so that some ratios
are exactly zero.

For each time point the diagnostic reports the mean, the quantiles at
`probs` (named `quantile_50`, `quantile_90`, and so on, with the
`probs = 1` quantile reported as `max`), the proportion of ratios
exceeding each of `thresholds` (named `prop_gt_10`, `prop_gt_50`), the
proportion of exact-zero ratios (`prop_zero`), the Kish effective sample
size `ess`, and its fraction of the sample size `ess_fraction`. The Kish
effective sample size is \\(\sum_i r_i)^2 / \sum_i r_i^2\\; it equals
the sample size when the ratios are equal and is invariant to their
overall scale.

For a matrix input the same summaries are also applied to the cumulative
product of the ratios across time points, under the `cumulative_`
prefix. This is the sequential positivity signature: mild per-step
non-overlap compounds multiplicatively, so the cumulative effective
sample size can collapse even while each per-time fraction stays
healthy. For a single-column matrix the cumulative summaries coincide
with the per-time summaries. The default thresholds of 10 and 50 are
conventional weight-magnitude heuristics rather than values sanctioned
by either source paper.

## References

Petersen ML, Porter KE, Gruber S, Wang Y, van der Laan MJ (2012).
Diagnosing and responding to violations in the positivity assumption.
*Statistical Methods in Medical Research*, 21(1):31-54.
[doi:10.1177/0962280210386207](https://doi.org/10.1177/0962280210386207)

Ring C, Schomaker M (2026). A Diagnostic to Find and Help Combat
Stochastic Positivity Issues, with a Focus on Continuous Treatments.

## Examples

``` r
# Point treatment: a numeric vector of density ratios.
ratios <- c(0.5, 1, 1, 2, 8)
check_density_ratios(ratios)
#> 
#> ── Density ratios ──────────────────────────────────────────────────────────────
#> Density ratios: 5 observations across 1 time point
#> Mean ratio: 2.5
#> 50th percentile: 1
#> 90th percentile: 5.6
#> 95th percentile: 6.8
#> 99th percentile: 7.76
#> Maximum: 8
#> Proportion > 10: 0
#> Proportion > 50: 0
#> Proportion exactly zero: 0
#> Kish ESS fraction: 0.445

# Time-varying treatment: one column per time point. The cumulative products
# across columns are summarized alongside the per-time ratios.
m <- matrix(c(1, 1, 2, 1, 3, 1), nrow = 2)
check_density_ratios(m)
#> 
#> ── Density ratios ──────────────────────────────────────────────────────────────
#> Density ratios: 2 observations across 3 time points
#> Summaries shown for time point 3
#> Mean ratio: 2
#> 50th percentile: 2
#> 90th percentile: 2.8
#> 95th percentile: 2.9
#> 99th percentile: 2.98
#> Maximum: 3
#> Proportion > 10: 0
#> Proportion > 50: 0
#> Proportion exactly zero: 0
#> Kish ESS fraction: 0.8
#> Cumulative ESS fraction: 0.662
```
