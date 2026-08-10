# Plot an effective-data-points diagnostic

Draws one view of a
[`check_edp()`](https://r-causal.github.io/positively/reference/check_edp.md)
result. `type = "histogram"` shows the distribution of EDP faceted by
intervention value, and `type = "ecdf"` shows its empirical cumulative
distribution colored by intervention value. For the estimator variant,
`type = "scatter"` plots `edp_outcome` against `edp_treatment` colored
by `ideal_weight`; it aborts for the data variant.

## Arguments

- object:

  An `edp_result` from
  [`check_edp()`](https://r-causal.github.io/positively/reference/check_edp.md).

- type:

  One of `"histogram"`, `"ecdf"`, or `"scatter"`.

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

autoplot(result, type = "histogram")

autoplot(result, type = "ecdf")

```
