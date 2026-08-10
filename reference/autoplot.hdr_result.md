# Plot an HDR non-overlap diagnostic

Draws the non-overlap ratio \\\hat{\tau}(a)\\ against the target
exposure value `a`, with a reference line at zero. For a sequential
result the curves are colored by time point.

## Arguments

- object:

  An `hdr_result` from
  [`check_hdr()`](https://r-causal.github.io/positively/reference/check_hdr.md)
  or
  [`check_hdr_seq()`](https://r-causal.github.io/positively/reference/check_hdr_seq.md).

- ...:

  Not used.

## Value

A [ggplot2::ggplot](https://ggplot2.tidyverse.org/reference/ggplot.html)
object.

## Examples

``` r
set.seed(1)
n <- 300
l <- rnorm(n)
dose <- rnorm(n, mean = l)
df <- data.frame(dose = dose, l = l)
result <- check_hdr(df, dose, l)

autoplot(result)

```
