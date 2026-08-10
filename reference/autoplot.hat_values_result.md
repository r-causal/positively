# Plot a hat-value diagnostic

Draws one of two views of a
[`check_hat_values()`](https://r-causal.github.io/positively/reference/check_hat_values.md)
result. `type = "null"` shows the null distribution of \\\hat{\phi}\\ as
a histogram with the observed \\\hat{\phi}\\ and the `conf_level` null
quantile marked. `type = "profile"` shows the proportion of
high-leverage candidates at each exposure percentile, the leverage
profile across the exposure range.

## Arguments

- object:

  A `hat_values_result` from
  [`check_hat_values()`](https://r-causal.github.io/positively/reference/check_hat_values.md).

- type:

  One of `"null"` or `"profile"`.

- ...:

  Not used.

## Value

A [ggplot2::ggplot](https://ggplot2.tidyverse.org/reference/ggplot.html)
object.

## Examples

``` r
set.seed(1)
n <- 200
x1 <- rnorm(n)
dose <- rnorm(n, mean = x1)
df <- data.frame(dose = dose, x1 = x1)
# null_reps is small here to keep the example fast.
result <- check_hat_values(df, dose, x1, null_reps = 50)

autoplot(result, type = "null")

autoplot(result, type = "profile")

```
