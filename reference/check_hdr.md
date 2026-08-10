# Diagnose positivity for continuous exposures with HDR non-overlap

`check_hdr()` computes the highest-density-region (HDR) non-overlap
ratio of Bao and Schomaker (2025) for continuous exposures. For a common
target exposure value, it reports the fraction of covariate profiles
whose highest-density region excludes that value, that is, the share of
the population for which the target dose is not supported.

## Usage

``` r
check_hdr(
  .data,
  .exposure,
  .covariates,
  mass = 0.95,
  values = NULL,
  density_estimator = hdr_density_normal(),
  exposure_type = c("auto", "continuous")
)
```

## Arguments

- .data:

  A data frame.

- .exposure:

  The continuous exposure column, selected with data-masking. Under the
  default `exposure_type = "auto"`, `check_hdr()` aborts when the
  exposure is detected as binary or categorical; declaring
  `exposure_type = "continuous"` accepts any numeric column that takes
  more than one value.

- .covariates:

  The covariate columns, selected with tidyselect. Numeric, logical,
  factor, and character columns are accepted, and each reaches the
  density estimator's formula as it was selected. Encoding a non-numeric
  column is the estimator's work, which the default estimator leaves to
  [`stats::lm()`](https://rdrr.io/r/stats/lm.html).

- mass:

  The HDR probability mass, a single value strictly between 0 and 1.
  Defaults to `0.95`.

- values:

  A numeric vector of target exposure values at which to evaluate the
  non-overlap ratio. `NULL` (the default) uses a 100-point grid spanning
  the observed exposure range.

- density_estimator:

  An `hdr_density` conditional-density estimator, built with
  [`hdr_density_normal()`](https://r-causal.github.io/positively/reference/hdr_density_normal.md)
  (the default) or
  [`new_hdr_density()`](https://r-causal.github.io/positively/reference/new_hdr_density.md).

- exposure_type:

  One of `"auto"` (detect from the data, the default) or `"continuous"`.
  A supplied type is authoritative and detection is not consulted, so
  `"continuous"` is rejected only when the exposure column is not
  numeric. Declaring it on a numeric column with few distinct values, a
  dose recorded at a handful of milligram levels for instance, is
  exactly the supported use: the unique-value heuristic reads such a
  column as categorical, so under `"auto"` the call aborts. Because
  numeric is the whole of what the type asks of the column, a two-valued
  numeric column is then accepted as well, which is the price of an
  authoritative declaration. A constant column is refused whichever way
  the type was decided, since the ratio compares covariate profiles at a
  common target dose and one observed dose leaves nothing to compare.

## Value

An `hdr_result` object, an S7 subclass of
[positivity_diagnostic](https://r-causal.github.io/positively/reference/positivity_diagnostic.md).
Its `@results` tibble has one row per target value with columns `value`
(the target `a`) and `nonoverlap` (the ratio \\\hat{\tau}(a)\\). It also
carries the properties `@mass` and `@density_estimator`.

[`generics::glance()`](https://generics.r-lib.org/reference/glance.html)
returns a one-row tibble with `n` (the sample size), `mass`,
`density_estimator`, `n_values` (the number of target values), and
`nonoverlap_min` and `nonoverlap_max`, the range of the ratio across
those targets.

## Details

For a continuous exposure \\A\\ with covariates \\L\\, fix a probability
mass `mass` (the paper's \\\alpha\\). At covariate profile \\l\\ the HDR
is the smallest set of exposure values capturing that mass of the
conditional density, \$\$A\_\alpha(l) = \\ a : f(a \mid l) \ge
f\_\alpha(l) \\, \qquad P(A \in A\_\alpha(l) \mid L = l) =
\mathrm{mass}.\$\$ The non-overlap ratio at a target `a` is the fraction
of profiles whose HDR excludes it, \$\$\hat{\tau}(a) = \frac{1}{n}
\sum\_{j=1}^{n} \mathbf{1}\\ a \notin A\_\alpha(l_j) \\,\$\$ a value in
\\\[0, 1\]\\: zero when `a` is supported everywhere, one when it is
supported nowhere.

The conditional density is supplied by `density_estimator`. The default
[`hdr_density_normal()`](https://r-causal.github.io/positively/reference/hdr_density_normal.md)
fits `lm(exposure ~ covariates)`, treats the density as Gaussian with
the residual standard deviation \\\hat{\sigma}\\, and uses the
closed-form cutoff \\\mathrm{dnorm}(z) / \hat{\sigma}\\ with \\z =
\Phi^{-1}((1 + \mathrm{mass}) / 2)\\. Membership then reduces to the
interval test \\\|a - \hat{\mu}(l)\| \le z\\\hat{\sigma}\\, so
\\\hat{\tau}(a)\\ is the fraction of fitted means more than
\\z\\\hat{\sigma}\\ from `a`.

The non-overlap ratio is a population-level, common-target quantity. It
is not small merely because the density model fits well: when the
exposure is strongly covariate-driven, even the best central target has
a nonzero floor, because setting the whole population to one common dose
is infeasible. The default normal estimator detects mean-shift support
gaps, where a stratum's supported dose moves away from a target. It does
not detect multimodal gaps: a hole between two modes of the true
conditional density is filled by the fitted normal and reported as
supported. Supply a flexible estimator through
[`new_hdr_density()`](https://r-causal.github.io/positively/reference/new_hdr_density.md)
when multimodal structure is expected.

When the exposure is close to a deterministic function of the
covariates, the fitted conditional density is narrow enough that no
profile's region reaches another profile's dose, and the ratio is one
everywhere inside the observed exposure range. That reading is literally
correct, since a dose fixed by the covariates is supported at no other
profile, but it is set by the width of the fitted density rather than by
any comparison between profiles, so `check_hdr()` warns when it happens.
The warning is decided on a dense grid spanning the observed range, not
on `values`, so it reports a property of the fit and does not depend on
which doses were asked about: a ratio of one at a few probes aimed into
a genuine support gap is a finding rather than an artifact, and stays
unremarked.

## References

Bao Y, Schomaker M (2025). Feasible Dose-Response Curves for Continuous
Treatments Under Positivity Violations.

## Examples

``` r
set.seed(1)
n <- 300
l <- rnorm(n)
dose <- rnorm(n, mean = l)
df <- data.frame(dose = dose, l = l)

result <- check_hdr(df, dose, l, values = c(-2, 0, 2))
result
#> 
#> ── HDR non-overlap ─────────────────────────────────────────────────────────────
#> Exposure: "dose" (continuous)
#> Observations: 300
#> HDR mass: 0.95
#> Density estimator: normal
#> Non-overlap over 3 targets: 0.05 to 0.503
```
