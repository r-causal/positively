
<!-- README.md is generated from README.Rmd. Please edit that file -->

# positively

<!-- badges: start -->

[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
<!-- badges: end -->

positively provides diagnostics for positivity violations and
extrapolation in causal inference. Causal effect estimates rest on the
positivity assumption: every unit must have a non-zero probability of
receiving each level of the exposure within every stratum of the
covariates. Where that probability is zero, or close to it, you face one
of two problems: 1) for an exposure-based model, like Inverse
Probability Weighting, you get a biased and unstable estimate or 2) for
an outcome-based model, like a regression model, an estimate comes from
the model extrapolating rather than from observed comparisons.

## Installation

positively is not yet on CRAN. You can install the development version
from [GitHub](https://github.com/r-causal/positively) with:

``` r
# install.packages("pak")
pak::pak("r-causal/positively")
```

## Example

`check_positivity()` detects the exposure type, or uses the one you
declare in `exposure_type`, and runs the diagnostics appropriate to that
type. The bundled `pos_violations` dataset has two positivity problems
planted in it: a structural subgroup that is never exposed, and a
practical near-violation where the fitted propensity approaches zero or
one in the tails of `x1`.

``` r
library(positively)

check <- check_positivity(pos_violations, exposure, c(x1, x2, region))
#> ℹ Treating `.exposure` as binary
check
#> 
#> ── Positivity check ────────────────────────────────────────────────────────────
#> Exposure: "exposure" (binary); 1000 observations; covariates x1, x2, and region
#> 
#> ── port ────────────────────────────────────────────────────────────────────────
#> 231 subgroups reported, 3 with low support
#> Rule: prevalence outside [0.05, 0.95] among subgroups of at least 5% of the
#> sample
#> 
#> ── extrapolation ───────────────────────────────────────────────────────────────
#> Geometric variability 0.144
#> 999 of 1000 have an opposite-exposure unit within one geometric variability;
#> 784 of 1000 fall in the opposite hull
#> 
#> ── edp ─────────────────────────────────────────────────────────────────────────
#> Data variant over 2 intervention values; edp 0.013 to 144.662
#> 
#> ℹ `sniff_violations()` for what was found, `$port` to extract a diagnostic,
#> `summary()` for every statistic.
```

`sniff_violations()` gathers what the run found into one table, across
every diagnostic, with the number behind each row in the `threshold`
column where there is one.

``` r
sniff_violations(check)
#> # A tibble: 5 × 7
#>   diagnostic    scope    label                      n statistic  value threshold
#>   <chr>         <chr>    <chr>                  <int> <chr>      <dbl>     <dbl>
#> 1 port          subgroup x1<-0.5027               305 prevalen… 0.0459     0.05 
#> 2 port          subgroup x1>=0.9626 & x1<1.154     53 prevalen… 0.962      0.05 
#> 3 port          subgroup x2>=1.063 & region=b      69 prevalen… 0          0.05 
#> 4 extrapolation overall  beyond one geometric …     1 prop_sup… 0.999      0.144
#> 5 extrapolation overall  outside the opposite-…   216 prop_in_… 0.784     NA
```

The positivity regression tree (`port`) reports each problem as a rule
over the covariates. Pull that diagnostic out with `$` and read its full
results with `tidy()`, one row per reported subgroup. The `low_support`
column marks those whose exposure prevalence falls below `beta` or above
`1 - beta` among subgroups at least `alpha` of the sample.

``` r
library(dplyr)

check$port |>
  tidy() |>
  filter(low_support)
#> # A tibble: 3 × 7
#>   subgroup   description  exposure_level     n proportion prevalence low_support
#>   <chr>      <chr>        <chr>          <int>      <dbl>      <dbl> <lgl>      
#> 1 x1         x1<-0.5027   1                305      0.305     0.0459 TRUE       
#> 2 x1         x1>=0.9626 … 1                 53      0.053     0.962  TRUE       
#> 3 x2, region x2>=1.063 &… 1                 69      0.069     0      TRUE
```

The three low-support subgroups match the planted ground truth. One rule
combining `region` and `x2` has an exposure prevalence of zero, the
structural violation, and two `x1` rules in the lower and upper tails
have prevalence near zero and one, the practical near-violation.

## Learn more

- Getting started: `vignette("positively")` walks through
  `check_positivity()` from end to end.
- The [package documentation](https://r-causal.github.io/positively/)
  holds the full reference and articles.
