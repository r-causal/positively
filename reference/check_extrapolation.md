# Diagnose extrapolation for binary exposures with Gower distances and the convex hull

`check_extrapolation()` implements the extrapolation diagnostics of King
and Zeng (2006) for binary point exposures. It measures, for every
observation, how well the opposite exposure group covers its position in
covariate space, using Gower distances that handle mixed numeric and
categorical covariates, and it optionally tests convex-hull membership
against the opposite group.

## Usage

``` r
check_extrapolation(
  .data,
  .exposure,
  .covariates,
  nearby = 1,
  hull = NULL,
  exposure_type = c("auto", "binary")
)
```

## Arguments

- .data:

  A data frame.

- .exposure:

  The binary exposure column, selected with data-masking.
  `check_extrapolation()` aborts unless the exposure has exactly two
  distinct values, whether that type is detected or declared through
  `exposure_type`.

- .covariates:

  The covariate columns, selected with tidyselect. Numeric and
  categorical covariates are both supported by the Gower distance.

- nearby:

  The reference radius as a multiple of the geometric variability. An
  observation counts an opposite-group unit as nearby when their Gower
  distance is at most `nearby` times the geometric variability. A single
  positive number, defaulting to `1`.

- hull:

  Whether to run the convex-hull membership test. `NULL` (the default)
  runs it automatically when there is at least one and at most ten
  numeric covariates and lpSolve is installed, and skips it otherwise.
  `TRUE` forces it to run, aborting when lpSolve is not installed and
  warning when there are no numeric covariates or the dimension is high.
  `FALSE` skips it.

- exposure_type:

  One of `"auto"` (detect from the data, the default) or `"binary"`. A
  supplied type is authoritative and detection is not consulted, so
  `"binary"` is rejected only when the exposure column does not have
  exactly two distinct values. Declaring it unlocks no call that
  `"auto"` would refuse, since detection reads any two-valued column as
  binary; it is accepted so that every diagnostic takes the same
  argument.

## Value

An `extrapolation_result` object, an S7 subclass of
[positivity_diagnostic](https://r-causal.github.io/positively/reference/positivity_diagnostic.md).
Its `@results` tibble has one row per observation with columns `.id`
(the row index), `exposure` (the exposure value), `frac_nearby`,
`gower_min`, `gower_mean`, `in_hull` (the hull membership flag, or `NA`
when the hull test did not run), and `low_support` (`TRUE` when the
nearest opposite-group unit lies farther away than one geometric
variability). It also carries the scalar properties
`@geometric_variability` and `@hull_run`.

[`summary()`](https://rdrr.io/r/base/summary.html) aggregates the
results by exposure group, returning one row per exposure level with the
columns `exposure`, `n`, `mean_gower_min`, `prop_supported` (the one
geometric variability fraction described in Details), and `prop_in_hull`
(the fraction inside the opposite group's hull, `NA` when the hull test
did not run).

[`generics::glance()`](https://generics.r-lib.org/reference/glance.html)
returns a one-row tibble with `n` (the sample size), `nearby`,
`geometric_variability`, `mean_frac_nearby`, `prop_supported`,
`hull_run` (kept logical), and `prop_in_hull`, which is `NA` when the
hull test did not run. The two readings, `prop_supported` and
`prop_in_hull`, answer different questions and are both reported.

## Details

Estimating a causal contrast requires comparing treated and untreated
units with similar covariates. Where one exposure group has no nearby
members of the other, the contrast rests on extrapolation rather than
data, the concern framed by Petersen and colleagues (2012) and by King
and Zeng (2006).

The Gower distance between two rows is the mean over covariates of a
per-covariate distance: for a numeric covariate, the absolute difference
divided by the pooled sample range; for a factor or character covariate,
one when the categories differ and zero when they agree. Distances
therefore lie in `[0, 1]` and mixed covariate types combine on a common
scale. Because the denominator is the sample range, a single extreme
covariate value rescales every distance on that axis.

The geometric variability is one half the mean of the full Gower
distance matrix, the Cuadras-Arenas quantity used by the WhatIf
software. It sets the reference radius: an observation is counted as
nearby-supported when an opposite-group unit lies within `nearby`
geometric variabilities. Three per-unit statistics are reported:
`gower_min` (distance to the nearest opposite-group unit), `gower_mean`
(mean distance to the opposite group), and `frac_nearby` (the fraction
of the opposite group within `nearby` geometric variabilities). A
planted support gap raises `gower_min` and lowers `frac_nearby` in both
directions. The per-group `prop_supported` reported by
[`summary()`](https://rdrr.io/r/base/summary.html), the `prop_supported`
statistic in
[`glance()`](https://generics.r-lib.org/reference/glance.html), and the
print method share one fixed definition: the fraction of units whose
nearest opposite-group unit lies within one geometric variability
(`gower_min <= geometric_variability`), independent of `nearby`. The
per-unit `low_support` column is the complement of that condition, so
`prop_supported` is the fraction of rows whose `low_support` is `FALSE`.

When every covariate is constant, all rows coincide, the geometric
variability is zero, and each observation is at distance zero from the
opposite group. Such an observation is fully supported: `gower_min` is
zero and `frac_nearby` is one. Constant numeric covariates likewise
contribute zero distance on their own axis without affecting the others.

The convex-hull test asks whether an observation lies inside the convex
hull of the opposite group, solved as a linear-programming feasibility
problem (weights that are non-negative, sum to one, and reproduce the
point). The hull test uses the numeric covariates only, since
categorical covariates carry no ordering to interpolate over. In-hull
membership collapses toward zero as the covariate dimension grows even
under perfect overlap, so the test is degenerate in high dimensions
(D'Amour and colleagues, 2021). It therefore runs automatically only
when there is at least one and at most ten numeric covariates and
lpSolve is installed. With no numeric covariates, more than ten of them,
or lpSolve absent, it is skipped and `in_hull` is `NA`. Setting
`hull = TRUE` forces the test to run and warns when the dimension is
high; `hull = FALSE` skips it.

## References

King G, Zeng L (2006). The Dangers of Extreme Counterfactuals.
*Political Analysis*, 14(2):131-159.
[doi:10.1093/pan/mpj004](https://doi.org/10.1093/pan/mpj004)

Petersen ML, Porter KE, Gruber S, Wang Y, van der Laan MJ (2012).
Diagnosing and responding to violations in the positivity assumption.
*Statistical Methods in Medical Research*, 21(1):31-54.
[doi:10.1177/0962280210386207](https://doi.org/10.1177/0962280210386207)

D'Amour A, Ding P, Feller A, Lei L, Sekhon J (2021). Overlap in
Observational Studies with High-Dimensional Covariates. *Journal of
Econometrics*, 221(2):644-654.
[doi:10.1016/j.jeconom.2019.10.014](https://doi.org/10.1016/j.jeconom.2019.10.014)

## Examples

``` r
set.seed(1)
n <- 100
x1 <- rnorm(n)
x2 <- rnorm(n)
exposure <- rbinom(n, 1, plogis(0.5 * x1 - 0.5 * x2))
df <- data.frame(exposure = exposure, x1 = x1, x2 = x2)

# hull = FALSE keeps the example free of the optional lpSolve dependency.
result <- check_extrapolation(df, exposure, c(x1, x2), hull = FALSE)
result
#> 
#> ── Extrapolation ───────────────────────────────────────────────────────────────
#> Exposure: "exposure" (binary)
#> Observations: 100
#> Geometric variability: 0.118
#> Nearby radius (1 x gv): 0.118
#> Mean fraction nearby: 0.156
#> Nearest opposite within one geometric variability: 97 of 100
#> Convex-hull test: not run
```
