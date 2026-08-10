# A container for a set of positivity diagnostics

`positivity_check` is the S7 container returned by
[`check_positivity()`](https://r-causal.github.io/positively/reference/check_positivity.md).
It holds one
[positivity_diagnostic](https://r-causal.github.io/positively/reference/positivity_diagnostic.md)
child per diagnostic, named for the diagnostic that produced it,
together with the exposure, covariates, exposure type, sample size, and
call the run resolved. Printing it gives a report: the exposure, its
type, the sample size, and the covariates stated once for the run, then
one section per diagnostic reading what that diagnostic found.
Diagnostics with a finding to report come first, and the rest follow in
the order they were requested.
[`sniff_violations()`](https://r-causal.github.io/positively/reference/sniff_violations.md)
is what those findings are.

## Usage

``` r
positivity_check(
  checks = list(),
  exposure = character(0),
  exposure_type = character(0),
  covariates = character(0),
  n = integer(0),
  call = quote({
 })
)
```

## Arguments

- checks:

  A named list of
  [positivity_diagnostic](https://r-causal.github.io/positively/reference/positivity_diagnostic.md)
  objects, in the order the diagnostics were requested. Each name is the
  diagnostic that produced the child.

- exposure:

  The exposure column name or names, time-ordered when the diagnostics
  are sequential.

- exposure_type:

  One of `"binary"`, `"categorical"`, or `"continuous"`.

- covariates:

  The covariate column names.

- n:

  The number of observations, an integer.

- call:

  The originating call.

## Value

A `positivity_check` object. Its
[`summary()`](https://rdrr.io/r/base/summary.html) method returns a
[tibble](https://tibble.tidyverse.org/reference/tibble.html) with one
row per statistic per diagnostic.

## Details

[`summary()`](https://rdrr.io/r/base/summary.html) gives the overview of
the run: one row per statistic per child, with the columns `diagnostic`,
`statistic`, `value`, and `threshold`. Each child contributes the
statistics its
[`generics::glance()`](https://generics.r-lib.org/reference/glance.html)
computed, so the overview reads the same way whichever diagnostics were
run.

[`sniff_violations()`](https://r-causal.github.io/positively/reference/sniff_violations.md)
reports what the run found, one row per finding across every child,
where [`summary()`](https://rdrr.io/r/base/summary.html) reports a
statistic per child whether or not anything was found.

Extract a child diagnostic with `x$port` or `x[["port"]]`, and list the
diagnostics the container holds with `names(x)`. Both extractors reject
a name the container does not hold rather than returning `NULL`. A child
carries its own
[`generics::tidy()`](https://generics.r-lib.org/reference/tidy.html) and
[`generics::glance()`](https://generics.r-lib.org/reference/glance.html)
methods, which are where a diagnostic's full results and its wide, typed
statistics are read.

[`autoplot()`](https://ggplot2.tidyverse.org/reference/autoplot.html)
draws each child's default view as one panel, and `autoplot(x, "port")`
draws a single child.

## Examples

``` r
set.seed(1)
n <- 150
x1 <- rnorm(n)
df <- data.frame(exposure = rbinom(n, 1, plogis(0.5 * x1)), x1 = x1)

res <- check_positivity(df, exposure, x1, diagnostics = "port")
#> ℹ Treating `.exposure` as binary

# The overview across every child.
summary(res)
#> # A tibble: 2 × 4
#>   diagnostic statistic     value threshold
#>   <chr>      <chr>         <dbl>     <dbl>
#> 1 port       n_subgroups      12     NA   
#> 2 port       n_low_support     0      0.05

# What the run found, one row per finding.
sniff_violations(res)
#> # A tibble: 0 × 7
#> # ℹ 7 variables: diagnostic <chr>, scope <chr>, label <chr>, n <int>,
#> #   statistic <chr>, value <dbl>, threshold <dbl>

# List the diagnostics and extract one child.
names(res)
#> [1] "port"
res$port
#> 
#> ── PoRT subgroups ──────────────────────────────────────────────────────────────
#> Exposure: "exposure" (binary)
#> Observations: 150
#> Reading rule: alpha = 0.05, gamma = 2
#> Prevalence threshold beta: 0.05
#> Subgroups: 12 reported, 0 with low support

# A child carries its own results and its own statistics.
tidy(res$port)
#> # A tibble: 12 × 7
#>    subgroup description   exposure_level     n proportion prevalence low_support
#>    <chr>    <chr>         <chr>          <int>      <dbl>      <dbl> <lgl>      
#>  1 x1       x1>=-0.6162 … 1                 16     0.107      0.0625 FALSE      
#>  2 x1       x1>=-1.265 &… 1                 13     0.0867     0.154  FALSE      
#>  3 x1       x1<-1.265     1                 13     0.0867     0.385  FALSE      
#>  4 x1       x1>=-0.6981 … 1                  9     0.06       0.889  FALSE      
#>  5 x1       x1>=0.8017 &… 1                 17     0.113      0.353  FALSE      
#>  6 x1       x1>=-0.1351 … 1                  6     0.04       0      FALSE      
#>  7 x1       x1>=-0.05811… 1                 19     0.127      0.474  FALSE      
#>  8 x1       x1>=-0.38 & … 1                 15     0.1        0.667  FALSE      
#>  9 x1       x1>=0.408 & … 1                 14     0.0933     0.5    FALSE      
#> 10 x1       x1>=0.6934 &… 1                  8     0.0533     0.75   FALSE      
#> 11 x1       x1>=0.2793 &… 1                 11     0.0733     0.909  FALSE      
#> 12 x1       x1>=1.489     1                  9     0.06       0.778  FALSE      
glance(res$port)
#> # A tibble: 1 × 6
#>       n n_subgroups n_low_support alpha  beta gamma
#>   <int>       <int>         <int> <dbl> <dbl> <dbl>
#> 1   150          12             0  0.05  0.05     2
```
