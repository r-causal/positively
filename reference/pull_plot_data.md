# Pull the data behind a plot

`pull_plot_data()` returns the tibble one view of a diagnostic result
draws. A reader can rebuild the figure with
[`ggplot2::ggplot()`](https://ggplot2.tidyverse.org/reference/ggplot.html)
directly, take the numbers somewhere else, or read what a view is
showing without rendering it.

## Usage

``` r
pull_plot_data(x, ...)
```

## Arguments

- x:

  A
  [positivity_diagnostic](https://r-causal.github.io/positively/reference/positivity_diagnostic.md)
  or a
  [positivity_check](https://r-causal.github.io/positively/reference/positivity_check.md).

- ...:

  Passed to methods. A class drawing more than one view takes `type`,
  the name of the view whose data to return, matched against that
  class's
  [`autoplot()`](https://ggplot2.tidyverse.org/reference/autoplot.html)
  menu. A
  [`check_port()`](https://r-causal.github.io/positively/reference/check_port.md)
  result takes `low_support_only` as
  [`autoplot()`](https://ggplot2.tidyverse.org/reference/autoplot.html)
  does, and a
  [positivity_check](https://r-causal.github.io/positively/reference/positivity_check.md)
  takes the name of the diagnostic to pull.

## Value

A [tibble](https://tibble.tidyverse.org/reference/tibble.html) holding
one view's data.

## Details

A view returns exactly one tibble, whatever the figure does with it.
Where a view splits its rows across layers, the whole frame comes back
and the split stays derivable from it.

`type` names the view. Each method mirrors its class's
[`autoplot()`](https://ggplot2.tidyverse.org/reference/autoplot.html)
menu, including the view that menu defaults to, so `pull_plot_data(x)`
returns the frame `autoplot(x)` draws. Every view is built from the
frame this returns, so the figure and the frame cannot disagree, and a
view a result cannot draw is refused here exactly as
[`autoplot()`](https://ggplot2.tidyverse.org/reference/autoplot.html)
refuses it.

For a
[`check_edp()`](https://r-causal.github.io/positively/reference/check_edp.md)
result:

- `"boxplot"` returns the box statistics themselves, one row per
  intervention, or one row per measure and intervention for the
  estimator variant. `ymin` and `ymax` are the fifth and the
  ninety-fifth percentile, `lower`, `middle`, and `upper` are the
  quartiles, `n` is the number of observations summarized, and
  `outliers` is a list-column of the values lying outside the whiskers,
  which the view draws as points. The estimator variant covers
  `edp_outcome` and `edp_treatment`; `ideal_weight` shares no scale with
  them and is left to the scatter view.

- `"histogram"`, `"ecdf"`, and `"density"` return the results at their
  own grain, one row per observation and intervention, with
  `intervention` as a factor whose levels follow the order the
  interventions were given in. The label is what separates the
  interventions rather than `value`, which a function intervention
  varies from observation to observation.

- `"scatter"` returns the whole estimator results tibble. The view draws
  the rows of finite `ideal_weight` in one layer and the infinite rows
  in another, and both are read off `ideal_weight`. It needs the
  estimator variant, and asking for it after a data-variant run is an
  error.

For a
[`check_density_ratios()`](https://r-causal.github.io/positively/reference/check_density_ratios.md)
result:

- `"distribution"` returns one row per ratio and time point, `time` and
  `ratio`. That is the grain the point-treatment histogram bins and the
  time-varying boxplot summarizes, and the raw ratios live nowhere else
  in the result.

- `"cumulative"` returns the cumulative-product summaries the series
  view draws, one row per time point and statistic, with `series` naming
  the summary in the words the key uses. It needs a time-varying input.

For a
[`check_eta_bias()`](https://r-causal.github.io/positively/reference/check_eta_bias.md)
result:

- `"bootstrap"` returns one row per bootstrap draw, keyed on the
  estimand `term`, on the truncation `level`, and on the `label` its
  facet strip reads. The view bins the `estimate`. Each draw also
  carries the `truth` its term was aimed at, which the view draws as a
  reference line, one per term rather than one per draw, since a term's
  truth is the same in every draw.

- `"sweep"` returns one row per term and truncation level, holding the
  `truncation` the sweep is drawn against, the `bias`, and the `lower`
  and `upper` ends of the band two Monte Carlo standard errors either
  side of it. It needs a sweep of more than one level.

For a
[`check_hat_values()`](https://r-causal.github.io/positively/reference/check_hat_values.md)
result:

- `"null"` returns the null replicates in `phi`. The three numbers the
  figure marks beside them, `null_quantile`, `phi_hat`, and
  `conf_level`, ride as repeated columns.

- `"profile"` returns one row per exposure percentile the candidates
  were built at, `prob` and the `fraction` of that percentile's
  candidates reading as high leverage. The results hold one row per
  candidate, so the aggregation is what this view adds.

For a
[`check_port()`](https://r-causal.github.io/positively/reference/check_port.md)
or
[`check_port_seq()`](https://r-causal.github.io/positively/reference/check_port_seq.md)
result there is one view, so the method takes `low_support_only` where a
menu would take `type`. It returns one row per bar: the reported
subgroups, with a leaf that more than one tree reported collapsed to a
single row, the display `label`, the `width` proportional to subgroup
size, and the pair of thresholds the row's prevalence was judged against
in `beta_lower` and `beta_upper`.

For a
[`check_hdr()`](https://r-causal.github.io/positively/reference/check_hdr.md)
or
[`check_hdr_seq()`](https://r-causal.github.io/positively/reference/check_hdr_seq.md)
result there is one view, and the frame is the results with `time` as a
factor where a sequential run carries one, because the wave separates
the curves rather than measuring anything.

For a
[`check_extrapolation()`](https://r-causal.github.io/positively/reference/check_extrapolation.md)
result:

- `"distribution"` returns the results with `exposure` as a factor,
  which is what panels the histograms.

- `"hull"` returns the same rows with `membership`, the reading of
  `in_hull` that the stacked bars and the key both use. It needs the
  hull test to have run.

A
[positivity_check](https://r-causal.github.io/positively/reference/positivity_check.md)
holds one frame per child rather than a frame of its own, so name the
diagnostic to pull: `pull_plot_data(check, "port")`. The frame is that
child's, unchanged, and a `type` passed alongside reaches it.

## Examples

``` r
set.seed(1)
n <- 100
x1 <- rnorm(n)
dose <- rnorm(n, mean = x1)
df <- data.frame(dose = dose, x1 = x1)
result <- check_edp(df, dose, x1, values = c(0, 1), exposure_type = "continuous")

pull_plot_data(result)
#> # A tibble: 2 × 8
#>   intervention     n  ymin lower middle upper  ymax outliers  
#>   <fct>        <int> <dbl> <dbl>  <dbl> <dbl> <dbl> <list>    
#> 1 0              100  5.59  16.6   22.7  25.5  26.7 <dbl [10]>
#> 2 1              100  4.43  14.4   20.9  24.4  25.0 <dbl [10]>
pull_plot_data(result, type = "ecdf")
#> # A tibble: 200 × 4
#>      .id intervention value   edp
#>    <int> <fct>        <dbl> <dbl>
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
#> # ℹ 190 more rows
```
