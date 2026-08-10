# The abstract parent class for positivity diagnostics

`positivity_diagnostic` is the abstract S7 parent that every positively
diagnostic result inherits from. It cannot be instantiated directly. It
fixes the shared property set (a tidy results tibble, the exposure
column names and type, the sample size, the resolved method parameters,
and the originating call) and supplies the
[`generics::tidy()`](https://generics.r-lib.org/reference/tidy.html),
[`generics::glance()`](https://generics.r-lib.org/reference/glance.html),
[`summary()`](https://rdrr.io/r/base/summary.html), and
[`print()`](https://rdrr.io/r/base/print.html) behavior a diagnostic
inherits unless it overrides one. Package developers extend it when
adding a new diagnostic.

## Arguments

- results:

  A [tibble](https://tibble.tidyverse.org/reference/tibble.html) of tidy
  diagnostic output.

- exposure:

  The exposure column name or names, time-ordered when the diagnostic is
  sequential.

- exposure_type:

  One of `"binary"`, `"categorical"`, or `"continuous"`.

- n:

  The number of observations, an integer.

- params:

  A list of method parameters as supplied and resolved.

- call:

  The originating call.

## Value

An abstract class object. Construction of a subclass returns a
`positivity_diagnostic`.

## Details

[`generics::tidy()`](https://generics.r-lib.org/reference/tidy.html)
returns `@results` unchanged. The inherited
[`generics::glance()`](https://generics.r-lib.org/reference/glance.html)
returns a one-row
[tibble](https://tibble.tidyverse.org/reference/tibble.html) holding the
single column `n`, the sample size, because a subclass that adds no
statistics of its own has nothing else to report. Each diagnostic
overrides it to state its own statistics beside `n`, typed as they are
computed.

[`summary()`](https://rdrr.io/r/base/summary.html) reports those
statistics in long form, with the columns `statistic`, `value`, and
`threshold`. Not every
[`generics::glance()`](https://generics.r-lib.org/reference/glance.html)
column earns a row. `value` is numeric throughout, so the character and
logical statistics stay behind in the wide output, as do the sample
size, the resolved method parameters, and any column that another
statistic reports as its threshold, which would otherwise be stated
twice. A diagnostic may override
[`summary()`](https://rdrr.io/r/base/summary.html) to aggregate
something else entirely:
[`check_extrapolation()`](https://r-causal.github.io/positively/reference/check_extrapolation.md)
summarizes its per-unit results by exposure group instead.

`threshold` is the one number behind the row, stated in the units of the
quantity it cuts, which are not always the units of `value`. Where the
statistic is itself a reading, the cut applies to that reading, so
`phi_hat` sits beside the null quantile it was compared against. Where
the statistic counts how many rows crossed a cut, the cut applies to the
per-row quantity rather than to the count, so a count of low-support
subgroups sits beside a subgroup prevalence and a count of high-leverage
candidates beside a hat value.

`threshold` is `NA` wherever there is no one number to state. That
covers a statistic that is a raw magnitude, a statistic governed by
several parameters with no single cut, a run that resolved a different
cut at each wave rather than one throughout, and a statistic whose cut
is real but is not one of the pairings the package tracks. An `NA`
therefore reports that the summary has no cut to show, not that the
diagnostic compared nothing against anything. The pairings are fixed
inside positively; a diagnostic defined outside the package reports `NA`
throughout.
