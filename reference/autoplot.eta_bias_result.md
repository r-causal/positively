# Plot an ETA.Bias diagnostic

Draws one of two views of a
[`check_eta_bias()`](https://r-causal.github.io/positively/reference/check_eta_bias.md)
result. `type = "bootstrap"` shows the distribution of bootstrap
estimates as a histogram with the target `truth` marked, one facet per
truncation level. `type = "sweep"` shows ETA.Bias with a
two-Monte-Carlo-standard-error band against the truncation lower bound,
the bias-variance tradeoff of weight truncation.

## Arguments

- object:

  An `eta_bias_result` from
  [`check_eta_bias()`](https://r-causal.github.io/positively/reference/check_eta_bias.md).

- type:

  One of `"bootstrap"` or `"sweep"`.

- ...:

  Not used.

## Value

A [ggplot2::ggplot](https://ggplot2.tidyverse.org/reference/ggplot.html)
object.

## Examples

``` r
set.seed(1)
n <- 300
x1 <- rnorm(n)
x2 <- rnorm(n)
a <- rbinom(n, 1, plogis(2 * (x1 + x2)))
y <- a + x1 + x2 + rnorm(n)
df <- data.frame(a = a, y = y, x1 = x1, x2 = x2)
# n_boot is small here to keep the example fast.
result <- check_eta_bias(df, a, y, c(x1, x2), n_boot = 25)
autoplot(result, type = "bootstrap")


swept <- check_eta_bias(
  df,
  a,
  y,
  c(x1, x2),
  truncation_grid = c(0, 0.05, 0.1),
  n_boot = 25
)
autoplot(swept, type = "sweep")

```
