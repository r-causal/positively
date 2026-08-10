# Diagnose sequential positivity for continuous exposures with HDR non-overlap

`check_hdr_seq()` applies the HDR non-overlap ratio of Bao and Schomaker
(2025) to a time-varying continuous exposure, one time point at a time.
At each time point it fits the conditional density of that exposure on
the history entering its conditioning set, following the
sequential-support logic of sPoRT (Chatton et al.), and reports the
non-overlap ratio at the target values.

## Usage

``` r
check_hdr_seq(
  .data,
  .exposures,
  .covariates,
  .baseline = NULL,
  mass = 0.95,
  values = NULL,
  density_estimator = hdr_density_normal(),
  lag = Inf,
  exposure_type = c("auto", "continuous")
)
```

## Arguments

- .data:

  A data frame in wide form, one row per subject.

- .exposures:

  An ordered tidyselect of exposure columns, one per time point.

- .covariates:

  A list of tidyselect expressions, one per time point, of the
  time-varying covariates. A length-one list is recycled across time
  points. Write it as a literal
  [`list()`](https://rdrr.io/r/base/list.html) call, for example
  `list(l0, l1, l2)`, rather than a pre-built list held in a variable.
  Numeric, logical, factor, and character columns are accepted, and each
  reaches the density estimator's formula as it was selected, as in
  [`check_hdr()`](https://r-causal.github.io/positively/reference/check_hdr.md).

- .baseline:

  A tidyselect of baseline covariates always included in every
  conditioning set, accepting the same column types as `.covariates`.
  Defaults to `NULL`.

- mass:

  The HDR probability mass, a single value strictly between 0 and 1.
  Defaults to `0.95`.

- values:

  A numeric vector of target exposure values. `NULL` (the default) uses
  a 100-point grid spanning the pooled exposure range.

- density_estimator:

  An `hdr_density` conditional-density estimator, built with
  [`hdr_density_normal()`](https://r-causal.github.io/positively/reference/hdr_density_normal.md)
  (the default) or
  [`new_hdr_density()`](https://r-causal.github.io/positively/reference/new_hdr_density.md).

- lag:

  The history window: the number of earlier time points whose covariates
  and exposures enter each conditioning set. Defaults to `Inf`, the full
  history.

- exposure_type:

  One of `"auto"` (detect from the data, the default) or `"continuous"`,
  applied to every exposure column in turn. A supplied type is
  authoritative and detection is not consulted, so `"continuous"` is
  rejected only when an exposure column is not numeric. Either way the
  error names every offending column at once, rather than the first one
  it meets. Declaring it on numeric columns with few distinct values,
  doses recorded at a handful of milligram levels for instance, is
  exactly the supported use: the unique-value heuristic reads such a
  column as categorical, so under `"auto"` the call aborts. Because
  numeric is the whole of what the type asks of a column, a two-valued
  numeric column is then accepted as well, which is the price of an
  authoritative declaration. A constant column is refused whichever way
  the type was decided, and the error names every constant column at
  once.

## Value

An `hdr_result` object, an S7 subclass of
[positivity_diagnostic](https://r-causal.github.io/positively/reference/positivity_diagnostic.md).
Its `@results` tibble has one row per time point and target value with
columns `time`, `value`, and `nonoverlap`.

[`generics::glance()`](https://generics.r-lib.org/reference/glance.html)
returns the columns
[`check_hdr()`](https://r-causal.github.io/positively/reference/check_hdr.md)
reports plus `n_times`, the number of time points. The non-overlap range
spans every wave rather than any one of them.

## Details

Data are supplied in wide form, one row per subject with one exposure
column per time point. For time point \\t\\ the diagnostic fits the
conditional density of the exposure at \\t\\ on its conditioning set,
then computes the point-in-time non-overlap ratio \\\hat{\tau}\_t(a)\\
exactly as
[`check_hdr()`](https://r-causal.github.io/positively/reference/check_hdr.md)
does. The conditioning set is the baseline covariates, the time-varying
covariates for each time point within the `lag` history window, and the
exposures at earlier time points within that window. With the default
`lag = Inf` the full history enters.

Because each time point is treated separately, a support gap that moves
a stratum's supported dose away from a target is isolated to the time
point where it occurs. The default normal estimator sees mean-shift gaps
only; a multimodal gap at a single time point is invisible to it, as for
[`check_hdr()`](https://r-causal.github.io/positively/reference/check_hdr.md).

A time point whose exposure is close to a deterministic function of its
conditioning set reads a ratio of one throughout that time point's
observed exposure range, for the reason
[`check_hdr()`](https://r-causal.github.io/positively/reference/check_hdr.md)
describes. Each time point is judged on a dense grid spanning its own
range, since the shared target grid spans the pooled range and can leave
a time point with a single target of its own. One warning names every
such time point.

## References

Bao Y, Schomaker M (2025). Feasible Dose-Response Curves for Continuous
Treatments Under Positivity Violations.

Chatton A, Schomaker M, Luque-Fernandez MA, Platt RW, Schnitzer ME
(2025). Is checking for sequential positivity violations getting you
down? Try sPoRT!

## Examples

``` r
set.seed(1)
n <- 200
df <- data.frame(
  l0 = rnorm(n),
  a1 = rnorm(n),
  l1 = rnorm(n),
  a2 = rnorm(n),
  l2 = rnorm(n),
  a3 = rnorm(n)
)

result <- check_hdr_seq(
  df,
  c(a1, a2, a3),
  list(l0, l1, l2),
  values = c(-2, 0, 2)
)
result
#> 
#> ── HDR non-overlap ─────────────────────────────────────────────────────────────
#> Exposure: "a1", "a2", and "a3" (continuous)
#> Observations: 200
#> HDR mass: 0.95
#> Density estimator: normal
#> Time points: 3
#> Non-overlap over 3 targets: 0 to 1
```
