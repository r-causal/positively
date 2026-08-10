# Plot the diagnostics a positivity check holds

Draws each child's default view. With no `diagnostic`, the views are
composed into one panel; naming a diagnostic draws that child alone.

## Arguments

- object:

  A
  [positivity_check](https://r-causal.github.io/positively/reference/positivity_check.md)
  from
  [`check_positivity()`](https://r-causal.github.io/positively/reference/check_positivity.md).

- diagnostic:

  Optionally, the name of one diagnostic to draw. Defaults to `NULL`,
  which panels every child.

- ...:

  Passed to a child's own
  [`autoplot()`](https://ggplot2.tidyverse.org/reference/autoplot.html)
  method. The panel passes them to every child, so an argument only one
  diagnostic accepts belongs with a named `diagnostic` rather than with
  the panel.

## Value

A [ggplot2::ggplot](https://ggplot2.tidyverse.org/reference/ggplot.html)
object.

## Details

The panel needs the patchwork package, which composes the children into
a single object that is still a
[ggplot2::ggplot](https://ggplot2.tidyverse.org/reference/ggplot.html)
and so can be built on further. Without it, naming a diagnostic still
works, because that path draws one child's plot and composes nothing.

A diagnostic is named the way
[`names()`](https://rdrr.io/r/base/names.html) lists it and `$` extracts
it, which is also the name the printed report heads its sections with,
so a reader can go from a section straight to its plot.

## Examples

``` r
set.seed(1)
n <- 300
x1 <- rnorm(n)
df <- data.frame(exposure = rbinom(n, 1, plogis(3 * x1)), x1 = x1)

check <- check_positivity(df, exposure, x1, diagnostics = c("port", "edp"))
#> ℹ Treating `.exposure` as binary

autoplot(check)

autoplot(check, "port")
```
