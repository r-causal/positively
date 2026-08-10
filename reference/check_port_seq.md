# Locate sequential positivity violations with the sPoRT algorithm

`check_port_seq()` applies the PoRT reading rule one time point at a
time along a treatment regime, the sequential Positivity Regression Tree
(sPoRT) algorithm of Chatton et al. Under the stratified strategy it
fits a tree at each time point among the subjects still following the
regime, so that a violation masked by already-treated subjects is
exposed.

## Usage

``` r
check_port_seq(
  .data,
  .exposures,
  .covariates,
  .baseline = NULL,
  .censoring = NULL,
  strategy = c("stratified", "pooled"),
  lag = Inf,
  alpha = 0.05,
  beta = 0.05,
  gamma = 2,
  n_bins = 3,
  breaks = NULL,
  exposure_type = c("auto", "binary", "categorical", "continuous"),
  ...
)
```

## Arguments

- .data:

  A data frame in wide form, one row per subject.

- .exposures:

  An ordered tidyselect of exposure columns, one per time point.

- .covariates:

  A list of tidyselect expressions, one per time point, of the
  time-varying covariates. A length-one list is recycled across time
  points.

- .baseline:

  A tidyselect of baseline covariates always included in every
  conditioning set. Defaults to `NULL`.

- .censoring:

  An ordered tidyselect of censoring indicators, one per time point,
  each coded `1` when the subject is censored and `0` otherwise. When
  supplied, each indicator is analyzed as its own binary-exposure tree,
  in addition to the exposure trees. The selection must match
  `.exposures` in length and every indicator must be binary;
  monotonicity is not required of the data, since a subject censored at
  a time point is dropped from later censoring risk sets. An indicator
  that is `NA` at a time point is read as unknown censoring status and
  drops that subject from that time point's censoring risk set and every
  later one. Supplying `.censoring` also restricts each exposure risk
  set to the followers uncensored so far. Censoring is available under
  the stratified strategy, the only strategy currently implemented.
  Defaults to `NULL`, which runs the exposure trees alone.

- strategy:

  The sequential strategy. Only `"stratified"` (per-time trees among the
  regime followers, the default) is implemented; `"pooled"` is reserved
  for a later release.

- lag:

  The history window: the number of earlier time points whose covariates
  and exposures enter each conditioning set. Defaults to `Inf`, the full
  history.

- alpha, beta, gamma, n_bins, breaks:

  As in
  [`check_port()`](https://r-causal.github.io/positively/reference/check_port.md).
  `beta` is resolved per time point when it is `"gruber"`.

- exposure_type:

  One of `"auto"` (detect from the data, the default), `"binary"`,
  `"categorical"`, or `"continuous"`, resolved once from the exposure
  pooled across time points. A supplied type is authoritative and
  detection is not consulted, so it is rejected only when the pooled
  exposure cannot carry it: `"continuous"` needs a numeric column and
  `"binary"` needs exactly two distinct values, while `"categorical"`
  asks nothing of the column. The type also decides the risk sets, so
  declaring a non-binary type on a binary exposure changes the analysis
  rather than only the guardrail; see Details.

- ...:

  Passed to
  [`rpart::rpart.control()`](https://rdrr.io/pkg/rpart/man/rpart.control.html).

## Value

A `port_result` object, an S7 subclass of
[positivity_diagnostic](https://r-causal.github.io/positively/reference/positivity_diagnostic.md).
Its `@results` tibble carries the point columns plus a leading `time`
column, and a `type` column distinguishing exposure from censoring rows
when `.censoring` is supplied. Censoring rows report `exposure_level`
`"1"`, the censored level, so their `prevalence` is the subgroup's
censoring rate. `@beta` holds the per-time exposure thresholds ordered
by time point, and `@censoring_beta` holds the per-time censoring
thresholds when `.censoring` is supplied.

[`generics::glance()`](https://generics.r-lib.org/reference/glance.html)
returns the columns
[`check_port()`](https://r-causal.github.io/positively/reference/check_port.md)
reports plus `n_times`, the number of time points. When `.censoring` is
supplied it also reports `n_censoring_subgroups` and
`n_censoring_low_support` beside the exposure counts, which then cover
the exposure rows alone. The `beta` column is `NA_real_` whenever the
per-time thresholds differ, as `beta = "gruber"` makes them: `NA` there
means that the threshold varies by time point rather than that it is
unknown, and the full per-time vector stays available on `@beta`. It is
`NA_real_` in one further case, where the threshold genuinely is
unknown: no time point resolved a threshold, because every one of them
was skipped. Reading `@beta` tells the two apart, since it holds the
differing values in the first case and is entirely `NA` in the second.

## Details

Data are supplied in wide form, one row per subject with one exposure
column per time point. Under `strategy = "stratified"` the risk set at
time \\t\\ is the set of subjects still following the rule. For a
monotone binary treatment this is the subjects not yet initiated, that
is, those whose immediately preceding treatment is the untreated level;
every subject with an observed time-1 exposure is in the risk set at the
first time point. Restricting to these followers is what reveals a
violation that pooling over already-treated subjects would mask, because
the already-treated carry their treatment forward and inflate the
subgroup prevalence. A subject whose time-\\t\\ exposure is missing is
excluded from the time-\\t\\ risk set, so \\n_t\\ is always the size of
the sample the tree is fitted on.

At each time point the diagnostic runs
[`check_port()`](https://r-causal.github.io/positively/reference/check_port.md)
with the time-\\t\\ exposure as the response and a conditioning set of
the baseline covariates, the time-varying covariates within the `lag`
history window, and the earlier exposures within that window. The
reading rule and the combination search behave exactly as in the point
diagnostic.

The exposure type and the untreated level are both resolved once from
the non-missing pooled exposure across all time points, so the follower
set is defined against the regime's untreated level rather than
whichever levels happen to survive at a given wave, and an `NA`
exposure, for example one a censored subject carries after leaving the
study, distorts neither. The type is detected from the pooled exposure
unless you declare it through `exposure_type`, and a declaration is
authoritative.

Only a binary type carries the follower restriction; under any other
type every subject stays in every risk set. Declaring a non-binary type
on a binary exposure therefore changes the analysis rather than only the
guardrail: the already-treated are pooled back in, the diagnostic stops
exposing the violations that pooling over them masks, and the trees tend
to flag the trivial rule that an earlier exposure predicts the current
one. Declare a non-binary type on a binary exposure only when you intend
to drop the regime restriction.

When a risk set at some time point is empty or too small to support the
reading rule, for example a fully treated preceding wave or a single
follower under the Gruber bound, that time point is skipped with a
warning, contributes no rows, and records `NA` in `@beta`.

When `beta = "gruber"` the prevalence threshold is resolved per time
point from the risk-set size \\n_t\\, as \\5 / (\sqrt{n_t} \\ \ln
n_t)\\. Under a monotone rule the risk set shrinks over time, so the
per-time threshold grows: smaller follower sets are held to a looser
bound. The resolved thresholds are stored in `@beta`, ordered by time
point.

Supplying `.censoring` adds a second family of per-time trees for the
censoring process (Chatton et al. 2025). Each censoring indicator is
analyzed as a binary exposure with the same two-sided reading rule,
which flags a subgroup whose censoring prevalence is below `beta` or
above `1 - beta`. The upper tail flags a subgroup nearly certain to be
censored, just as the exposure trees flag a subgroup nearly certain to
be treated; the lower tail flags a subgroup nearly certain to remain
uncensored, the never-censored subgroups on which the counterfactual
under no censoring rests (Chatton et al. 2025, page 11). Because real
censoring is light in most subgroups, the lower tail usually accounts
for most of the low-support censoring rows. The censoring tree at time
\\t\\ runs on the subjects uncensored through the previous time point,
the never-censored risk set, which is the full uncensored cohort rather
than the treatment-regime followers, since the treatment rule plays no
role in the censoring process. Its conditioning set is the same baseline
and `lag`-window history the exposure trees use, with the current
treatment added as a covariate. Its prevalence threshold is resolved
from the size of the censoring risk set, which shrinks as earlier
censoring accrues. A censoring wave in which every subject in its risk
set shares the same censoring status offers no contrast to split on, so
it contributes no censoring rows while still recording its resolved
threshold in `@censoring_beta`. Censoring rows are marked in the results
by a `type` column with the value `"censoring"`, against `"exposure"`
for the exposure trees, and this column appears only when `.censoring`
is supplied.

Supplying `.censoring` also restricts the exposure risk sets. At time
\\t\\ the exposure risk set becomes the regime followers who are also
uncensored through time \\t - 1\\, since a subject censored earlier has
no observed exposure at time \\t\\ and cannot follow the regime (Chatton
et al. 2025, Figure 2). The exposure risk sets therefore shrink faster
than under the follower restriction alone, and the per-time thresholds
in `@beta` grow accordingly. A censored subject that carries an `NA`
exposure after it leaves the study falls outside the restricted exposure
risk sets, so its missing exposures never enter a tree.

## References

Danelian G, Foucher Y, Léger M, Le Borgne F, Chatton A (2023).
Identification of in-sample positivity violations using regression
trees: the PoRT algorithm.
[doi:10.1515/jci-2022-0032](https://doi.org/10.1515/jci-2022-0032)

Chatton A, Schomaker M, Luque-Fernandez MA, Platt RW, Schnitzer ME
(2025). Is checking for sequential positivity violations getting you
down? Try sPoRT!

## Examples

``` r
set.seed(1)
n <- 1000
c1 <- rnorm(n)
a1 <- rbinom(n, 1, 0.3)
# Monotone treatment: once treated, stays treated.
a2 <- a1
followers <- a1 == 0
# Subjects with c1 > 1 never initiate at time 2.
a2[followers] <- rbinom(sum(followers), 1, ifelse(c1[followers] > 1, 0, 0.4))
df <- data.frame(c1 = c1, a1 = a1, a2 = a2)

result <- check_port_seq(df, c(a1, a2), list(c1))
#> ℹ Treating `.exposures` as binary
result
#> 
#> ── sPoRT subgroups ─────────────────────────────────────────────────────────────
#> Exposure: "a1" and "a2" (binary)
#> Observations: 1000
#> Time points: 2
#> Reading rule: alpha = 0.05, gamma = 2
#> Prevalence threshold beta: 0.05
#> Subgroups: 144 reported, 1 with low support

# Add censoring indicators, coded 1 when censored. Censoring is near-certain
# at time 2 where c1 > 1, which the censoring trees flag.
cens1 <- rbinom(n, 1, 0.08)
cens2 <- rbinom(n, 1, ifelse(c1 > 1, 0.98, 0.1))
df$cens1 <- cens1
df$cens2 <- cens2

check_port_seq(df, c(a1, a2), list(c1), .censoring = c(cens1, cens2))
#> ℹ Treating `.exposures` as binary
#> 
#> ── sPoRT subgroups ─────────────────────────────────────────────────────────────
#> Exposure: "a1" and "a2" (binary)
#> Observations: 1000
#> Time points: 2
#> Reading rule: alpha = 0.05, gamma = 2
#> Prevalence threshold beta: 0.05
#> Exposure subgroups: 142 reported, 1 with low support
#> Censoring subgroups: 116 reported, 5 with low support
```
