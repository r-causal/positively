# Locate positivity-violating subgroups with the PoRT algorithm

`check_port()` runs the Positivity Regression Tree (PoRT) algorithm of
Danelian et al. (2023). It grows regression trees on the covariates with
the exposure as the response, then reads a fixed rule off each tree to
flag subgroups in which one exposure level is rare or absent.

## Usage

``` r
check_port(
  .data,
  .exposure,
  .covariates,
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

  A data frame.

- .exposure:

  The exposure column, selected with data-masking.

- .covariates:

  The covariate columns, selected with tidyselect. Numeric and factor
  covariates are both accepted; the trees split them directly.

- alpha:

  The minimum subgroup size for the reading rule, a fraction of the
  sample size strictly between 0 and 1. Defaults to `0.05`.

- beta:

  The prevalence threshold for the reading rule. A single number
  strictly between 0 and 0.5 (`0.05` by default), or the string
  `"gruber"` for the sample-size-adaptive bound \\5 / (\sqrt{n} \\ \ln
  n)\\. The upper limit is 0.5 because the rule flags a prevalence below
  `beta` or above `1 - beta`.

- gamma:

  The maximum number of covariates that may jointly define a subgroup, a
  whole number of at least one. Defaults to `2`.

- n_bins:

  For a continuous exposure, the number of quantile bins into which the
  exposure is categorized. A whole number of at least two, `3` by
  default; a single bin would place every observation at one exposure
  level and flag every subgroup. Ignored for binary and categorical
  exposures.

- breaks:

  For a continuous exposure, at least three distinct numeric cut points
  that override `n_bins`. `NULL` (the default) uses quantile bins.

- exposure_type:

  One of `"auto"` (detect from the data, the default), `"binary"`,
  `"categorical"`, or `"continuous"`. A supplied type is authoritative
  and detection is not consulted, so it is rejected only when the
  exposure column cannot carry it: `"continuous"` needs a numeric column
  and `"binary"` needs exactly two distinct values, while
  `"categorical"` asks nothing of the column. The type also selects how
  the exposure is read: only a continuous type is categorized into
  `n_bins` quantile bins.

- ...:

  Passed to
  [`rpart::rpart.control()`](https://rdrr.io/pkg/rpart/man/rpart.control.html),
  for example `cp`.

## Value

A `port_result` object, an S7 subclass of
[positivity_diagnostic](https://r-causal.github.io/positively/reference/positivity_diagnostic.md).
Its `@results` tibble has one row per reported subgroup with columns
`subgroup` (the covariates defining the rule), `description` (the
human-readable rule), `exposure_level` (the exposure level examined),
`n` (the subgroup size), `proportion` (its share of the sample),
`prevalence` (the exposure prevalence in the subgroup), and
`low_support` (the logical reading-rule outcome). It also carries the
`@trees`, `@alpha`, `@beta`, and `@gamma` properties.

[`generics::glance()`](https://generics.r-lib.org/reference/glance.html)
returns a one-row tibble with `n` (the sample size), `n_subgroups` (how
many subgroups were reported), `n_low_support` (how many of those the
reading rule marked), and the resolved rule parameters `alpha`, `beta`,
and `gamma`, which the counts cannot be read without.

## Details

PoRT searches for subgroups, defined by covariate splits, in which
positivity is threatened because an exposure level is (near) unobserved.
For a binary exposure the response is the indicator of the higher
exposure level; for a categorical exposure the algorithm runs once per
level with a one-versus-rest indicator; for a continuous exposure the
exposure is first categorized into `n_bins` quantile bins (or the bins
implied by `breaks`) and each bin is treated as a level. Covariates
enter the trees unchanged: `rpart` splits numeric covariates at
data-driven thresholds and groups factor levels, so no covariate
categorization is performed here.

The regression trees follow the control values quoted by Danelian et al.
(section 2.5): a minimum of 20 observations to attempt a split, at least
six observations in a leaf, a maximum depth of 30, and no
cost-complexity pruning (`cp = 0`), so that the largest number of
divisions is retained. The paper states a minimum leaf size of six;
`rpart`'s own default of `round(minsplit / 3)` is seven, so the value is
pinned to six to match the paper. This affects only the smallest leaves
and is immaterial at realistic sample sizes, where the `alpha` size gate
dominates.

The reading rule (Danelian et al., section 2.3) marks a node as low
support when its exposure prevalence is extreme, below `beta` or above
`1 - beta`, and its size is at least `alpha` of the sample. Once a node
is marked its descendants are suppressed, so the reading is reported at
the coarsest subgroup that isolates the violation. Because the rule is
deterministic given a tree, only the tree carries stochasticity.

The combination search controls how many covariates may jointly define a
subgroup (Danelian et al., section 2.4). Step one grows one tree per
single covariate; covariates that produce a low-support subgroup are
removed, and step `k` grows one tree per `k`-covariate combination of
the remaining covariates, up to `gamma`. Raising `gamma` uncovers joint
violations that no single covariate isolates, at the cost of more trees.

Setting `beta = "gruber"` resolves the threshold to \\5 / (\sqrt{n} \\
\ln n)\\, the sample-size-adaptive bound used by sPoRT. It is tighter
than `0.05` for samples larger than roughly 350 and looser below that,
so a small sample is held to a less stringent prevalence threshold.

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
n <- 400
x1 <- rnorm(n)
x2 <- factor(sample(c("a", "b"), n, replace = TRUE))
# No subject with x1 > 1 and x2 == "b" is ever treated.
ps <- ifelse(x1 > 1 & x2 == "b", 0, 0.5)
exposure <- rbinom(n, 1, ps)
df <- data.frame(exposure = exposure, x1 = x1, x2 = x2)

result <- check_port(df, exposure, c(x1, x2))
#> ℹ Treating `.exposure` as binary
result
#> 
#> ── PoRT subgroups ──────────────────────────────────────────────────────────────
#> Exposure: "exposure" (binary)
#> Observations: 400
#> Reading rule: alpha = 0.05, gamma = 2
#> Prevalence threshold beta: 0.05
#> Subgroups: 73 reported, 1 with low support
```
