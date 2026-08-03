
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
covariates. Where that probability is zero, or close to it, you face of
two problems: 1) for an exposure-based model, like Inverse Probability
Weighting, you get a biased and unstable estimate or 2) for an
outcome-based model, like a regression model, an estimate comes from the
model extrapolating rather than from observed comparisons.

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
#> 
#> ── edp ──
#> 
#> ── edp_result ──────────────────────────────────────────────────────────────────
#> Exposure: "exposure" (binary)
#> Observations: 1000
#> Variant: data
#> Intervention values: 2
#> edp range: 0.013 to 144.662
#> 
#> ── port ──
#> 
#> ── port_result ─────────────────────────────────────────────────────────────────
#> Exposure: "exposure" (binary)
#> Observations: 1000
#> Reading rule: alpha = 0.05, gamma = 2
#> Prevalence threshold beta: 0.05
#> Subgroups: 231 reported, 3 with low support
#> 
#> ── extrapolation ──
#> 
#> ── extrapolation_result ────────────────────────────────────────────────────────
#> Exposure: "exposure" (binary)
#> Observations: 1000
#> Geometric variability: 0.144
#> Nearby radius (1 x gv): 0.144
#> Mean fraction nearby: 0.284
#> Nearest opposite within one geometric variability: 999 of 1000
#> In opposite-group hull: 784 of 1000
```

The positivity regression tree (`port`) reports each problem as a rule
over the covariates. `tidy()` gives one row per reported subgroup, and
the `low_support` column marks those whose exposure prevalence falls
below `beta` or above `1 - beta` among subgroups at least `alpha` of the
sample.

``` r
library(dplyr)

check[["port"]] |>
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
