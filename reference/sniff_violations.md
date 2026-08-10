# Report what a set of diagnostics found

`sniff_violations()` answers what a run turned up, one row per finding,
across every diagnostic it holds. It is the row-level companion to
[`summary()`](https://rdrr.io/r/base/summary.html), which reports a
statistic per diagnostic whether or not anything was found.

## Usage

``` r
sniff_violations(x, scope = c("subgroup", "overall"), ...)
```

## Arguments

- x:

  A
  [positivity_diagnostic](https://r-causal.github.io/positively/reference/positivity_diagnostic.md)
  or a
  [positivity_check](https://r-causal.github.io/positively/reference/positivity_check.md).

- scope:

  The kinds of row to report, any of `"subgroup"`, `"overall"`, and
  `"unit"`. Defaults to the first two.

- ...:

  Passed to methods.

## Value

A [tibble](https://tibble.tidyverse.org/reference/tibble.html) with one
row per finding and the columns `diagnostic`, `scope`, `label`, `n`,
`statistic`, `value`, and `threshold`. `value` and `threshold` are
numeric, and `n` is an integer count that is `NA` where a row is not
about a set of rows.

## Details

A diagnostic reports a finding by declaring one. Each participating
class states which of its statistics it reports and under which label; a
class that declares nothing reports nothing, and no column is read as a
flag on its behalf. A diagnostic defined outside positively therefore
reports no rows until it defines a method, whatever its results hold.

`scope` names the kind of thing a row is about, and rows of different
kinds are not comparable:

- `"subgroup"`, a covariate region, which is something a reader can act
  on.

- `"overall"`, a reading taken across the whole run.

- `"unit"`, a single row of a distribution whose aggregate is the
  reading.

Unit rows are reachable but are never the default, because a diagnostic
that reports per unit produces far more rows than one that reports per
subgroup, and stacking the two buries the smaller set.

A reading contributes a row only when it fires, meaning at least one
unit falls on the wrong side of it. A run in which every unit is
supported reports nothing rather than reporting how well supported it
was.

`threshold` carries the one number behind a row, stated in the units of
the quantity it cuts, which are not always the units of `value`. Where
the statistic is itself a reading, the cut applies to that reading:
`phi_hat` sits beside the null quantile it was compared against, and a
subgroup's prevalence beside the prevalence it was cut at. Where the
statistic aggregates a per-row quantity, the cut applies to that
quantity rather than to the aggregate, so `prop_supported` sits beside a
Gower distance, which is what each unit was measured against, rather
than beside a proportion. The number is the run's own resolved value
rather than a package default, so a reader can see where a cut came
from. It is `NA` only where no number stands behind the reading at all,
as for the convex hull test, which reports containment rather than a
comparison.

## Examples

``` r
set.seed(1)
n <- 300
x1 <- rnorm(n)
df <- data.frame(exposure = rbinom(n, 1, plogis(3 * x1)), x1 = x1)

res <- check_positivity(
  df,
  exposure,
  x1,
  diagnostics = c("port", "extrapolation")
)
#> ℹ Treating `.exposure` as binary

sniff_violations(res)
#> # A tibble: 5 × 7
#>   diagnostic    scope    label                      n statistic  value threshold
#>   <chr>         <chr>    <chr>                  <int> <chr>      <dbl>     <dbl>
#> 1 port          subgroup x1<-0.828                 51 prevalen… 0         0.05  
#> 2 port          subgroup x1>=-0.695 & x1<-0.41…    33 prevalen… 0.0303    0.05  
#> 3 port          subgroup x1>=0.9318                53 prevalen… 0.981     0.05  
#> 4 extrapolation overall  beyond one geometric …    37 prop_sup… 0.877     0.0976
#> 5 extrapolation overall  outside the opposite-…    80 prop_in_… 0.733    NA     

# Unit rows are available but are never the default.
sniff_violations(res, scope = "unit")
#> # A tibble: 117 × 7
#>    diagnostic    scope label                         n statistic value threshold
#>    <chr>         <chr> <chr>                     <int> <chr>     <dbl>     <dbl>
#>  1 extrapolation unit  beyond one geometric var…    NA gower_min 0.252    0.0976
#>  2 extrapolation unit  beyond one geometric var…    NA gower_min 0.211    0.0976
#>  3 extrapolation unit  beyond one geometric var…    NA gower_min 0.117    0.0976
#>  4 extrapolation unit  beyond one geometric var…    NA gower_min 0.101    0.0976
#>  5 extrapolation unit  beyond one geometric var…    NA gower_min 0.139    0.0976
#>  6 extrapolation unit  beyond one geometric var…    NA gower_min 0.216    0.0976
#>  7 extrapolation unit  beyond one geometric var…    NA gower_min 0.178    0.0976
#>  8 extrapolation unit  beyond one geometric var…    NA gower_min 0.174    0.0976
#>  9 extrapolation unit  beyond one geometric var…    NA gower_min 0.127    0.0976
#> 10 extrapolation unit  beyond one geometric var…    NA gower_min 0.101    0.0976
#> # ℹ 107 more rows
```
