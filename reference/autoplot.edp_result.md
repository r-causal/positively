# Plot an effective-data-points diagnostic

Draws one view of a
[`check_edp()`](https://r-causal.github.io/positively/reference/check_edp.md)
result. `type = "boxplot"`, the default, boxes EDP by intervention. Its
whiskers reach the fifth and the ninety-fifth percentile, following Ring
and Schomaker (2026) rather than Tukey's fences, and the values outside
them are drawn as points. For the estimator variant it draws the
outcome-model and the treatment-model measure side by side in a facet
per intervention. `type = "histogram"` shows the distribution of EDP
faceted by intervention, `type = "ecdf"` shows its empirical cumulative
distribution colored by intervention, and `type = "density"` draws one
density ridgeline per intervention, which needs the ggridges package.
Each of those three draws one measure: `edp` for the data variant and
`edp_outcome` for the estimator variant, the measure carrying the
exposure dimension and so the one an intervention moves.
`type = "scatter"` plots `edp_outcome` against `edp_treatment` colored
by `ideal_weight`; it needs the estimator variant and aborts for the
data variant.
[`pull_plot_data()`](https://r-causal.github.io/positively/reference/pull_plot_data.md)
returns the tibble any of these views draws.

## Arguments

- object:

  An `edp_result` from
  [`check_edp()`](https://r-causal.github.io/positively/reference/check_edp.md).

- type:

  One of `"boxplot"`, `"histogram"`, `"ecdf"`, `"density"`, or
  `"scatter"`.

- ...:

  Not used.

## Value

A [ggplot2::ggplot](https://ggplot2.tidyverse.org/reference/ggplot.html)
object.

## Examples

``` r
set.seed(1)
n <- 100
x1 <- rnorm(n)
dose <- rnorm(n, mean = x1)
df <- data.frame(dose = dose, x1 = x1)
result <- check_edp(df, dose, x1, values = c(0, 1), exposure_type = "continuous")

autoplot(result)

autoplot(result, type = "histogram")

autoplot(result, type = "ecdf")

```
