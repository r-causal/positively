# Diagnose positivity across the applicable methods

`check_positivity()` resolves the exposure type, detecting it from the
data unless you declare one, runs the diagnostics that apply to that
type with their defaults, and returns a
[positivity_check](https://r-causal.github.io/positively/reference/positivity_check.md)
container. Printing it gives a report: the exposure, its type, the
sample size, and the covariates stated once for the run, then one
section per diagnostic reading what that diagnostic found.

## Usage

``` r
check_positivity(
  .data,
  .exposure,
  .covariates,
  diagnostics = NULL,
  exposure_type = c("auto", "binary", "categorical", "continuous"),
  args = list()
)
```

## Arguments

- .data:

  A data frame.

- .exposure:

  The exposure column, selected with data-masking.

- .covariates:

  The covariate columns, selected with tidyselect. Required, with no
  default, so that outcome columns are not swept in by accident. The
  exposure column must not be selected.

- diagnostics:

  The diagnostics to run. `NULL` (the default) selects the applicable
  set for the exposure type; otherwise a subset of `"edp"`, `"port"`,
  `"hat_values"`, `"hdr"`, and `"extrapolation"`, run in the order
  given.

- exposure_type:

  One of `"auto"` (detect from the data, the default), `"binary"`,
  `"categorical"`, or `"continuous"`. A supplied type is authoritative:
  detection is not consulted, and the resolved type is forwarded to
  every diagnostic that runs rather than re-derived from the data by
  each diagnostic in turn. It must be a type the exposure column can
  carry, and it must apply to every entry in `diagnostics`; both are
  checked before any diagnostic runs.

- args:

  A named list of per-diagnostic option lists, for example
  `list(port = list(alpha = 0.1))`. Each name must be a diagnostic being
  run.

## Value

A
[positivity_check](https://r-causal.github.io/positively/reference/positivity_check.md)
object holding one
[positivity_diagnostic](https://r-causal.github.io/positively/reference/positivity_diagnostic.md)
child per diagnostic, named for the diagnostic that produced it,
alongside the exposure, covariates, exposure type, and sample size the
call resolved.

## Details

`check_positivity()` resolves the exposure type a single time and hands
the resolved type to every diagnostic it composes, so none of them
derives a type of its own and the type is announced once, through an
informational message suppressed by `options(positively.quiet = TRUE)`.
Every other informational message a diagnostic emits still reaches you.

Detection is a default, never an override. Under
`exposure_type = "auto"` the type is inferred from the data; supply a
type instead and it is authoritative, with detection not consulted at
all. Either way the resolved type is checked against the exposure column
before any diagnostic runs: a continuous type needs a numeric column, a
binary type needs exactly two distinct values, and a categorical type
asks nothing of the column. A declared type therefore aborts only when
its math cannot run on the column as given, never because the
unique-value heuristic read the column differently than you did. A dose
recorded at eight distinct milligram levels is detected as categorical,
and `exposure_type = "continuous"` runs the continuous-exposure
diagnostics on it regardless.

The diagnostics that run must also apply to the resolved type. Under
`"auto"` that type is a guess, so an inapplicable-diagnostic error names
the `exposure_type` to declare in order to run the rejected diagnostics,
offering only a type the exposure column can carry. Under a declared
type the type is your premise, and the advice is to adjust `diagnostics`
instead.

The default diagnostic set depends on the exposure type:

- **binary**:
  [`check_edp()`](https://r-causal.github.io/positively/reference/check_edp.md),
  [`check_port()`](https://r-causal.github.io/positively/reference/check_port.md),
  [`check_extrapolation()`](https://r-causal.github.io/positively/reference/check_extrapolation.md)

- **categorical**:
  [`check_edp()`](https://r-causal.github.io/positively/reference/check_edp.md),
  [`check_port()`](https://r-causal.github.io/positively/reference/check_port.md)

- **continuous**:
  [`check_edp()`](https://r-causal.github.io/positively/reference/check_edp.md),
  [`check_port()`](https://r-causal.github.io/positively/reference/check_port.md)
  (with the exposure categorized),
  [`check_hat_values()`](https://r-causal.github.io/positively/reference/check_hat_values.md),
  [`check_hdr()`](https://r-causal.github.io/positively/reference/check_hdr.md)

Pass `diagnostics` to run an explicit subset instead; the requested
order is honored. Every name must be a diagnostic that
`check_positivity()` composes and that applies to the resolved exposure
type, otherwise the call aborts with the valid options listed.

[`check_eta_bias()`](https://r-causal.github.io/positively/reference/check_eta_bias.md)
and
[`check_density_ratios()`](https://r-causal.github.io/positively/reference/check_density_ratios.md)
are never composed here. The first needs an outcome and the second needs
user-supplied density ratios, so both are called directly rather than
through `check_positivity()`.

Per-method options are threaded through `args`, a named list of option
lists, for example `args = list(port = list(alpha = 0.1))`. Every name
in `args` must be a diagnostic that is actually being run.

## See also

[`check_edp()`](https://r-causal.github.io/positively/reference/check_edp.md),
[`check_port()`](https://r-causal.github.io/positively/reference/check_port.md),
[`check_hat_values()`](https://r-causal.github.io/positively/reference/check_hat_values.md),
[`check_hdr()`](https://r-causal.github.io/positively/reference/check_hdr.md),
and
[`check_extrapolation()`](https://r-causal.github.io/positively/reference/check_extrapolation.md)
for the individual diagnostics, and
[`check_eta_bias()`](https://r-causal.github.io/positively/reference/check_eta_bias.md)
and
[`check_density_ratios()`](https://r-causal.github.io/positively/reference/check_density_ratios.md)
for the two that `check_positivity()` does not call.

## Examples

``` r
set.seed(1)
n <- 150
x1 <- rnorm(n)
x2 <- rnorm(n)
ps <- 0.2 + 0.6 * plogis(0.5 * x1 - 0.5 * x2)
df <- data.frame(exposure = rbinom(n, 1, ps), x1 = x1, x2 = x2)

# The binary default set: edp, port, and extrapolation.
check_positivity(df, exposure, c(x1, x2))
#> ℹ Treating `.exposure` as binary
#> 
#> ── Positivity check ────────────────────────────────────────────────────────────
#> Exposure: "exposure" (binary); 150 observations; covariates x1 and x2
#> 
#> ── port ────────────────────────────────────────────────────────────────────────
#> 42 subgroups reported, 1 with low support
#> Rule: prevalence outside [0.05, 0.95] among subgroups of at least 5% of the
#> sample
#> 
#> ── extrapolation ───────────────────────────────────────────────────────────────
#> Geometric variability 0.107
#> 146 of 150 have an opposite-exposure unit within one geometric variability; 132
#> of 150 fall in the opposite hull
#> 
#> ── edp ─────────────────────────────────────────────────────────────────────────
#> Data variant over 2 intervention values; edp 1.558 to 36.088
#> 
#> ℹ `sniff_violations()` for what was found, `$port` to extract a diagnostic,
#> `summary()` for every statistic.

# Run a single diagnostic with a tuned option.
check_positivity(
  df,
  exposure,
  c(x1, x2),
  diagnostics = "port",
  args = list(port = list(alpha = 0.1))
)
#> ℹ Treating `.exposure` as binary
#> 
#> ── Positivity check ────────────────────────────────────────────────────────────
#> Exposure: "exposure" (binary); 150 observations; covariates x1 and x2
#> 
#> ── port ────────────────────────────────────────────────────────────────────────
#> 42 subgroups reported, 0 with low support
#> Rule: prevalence outside [0.05, 0.95] among subgroups of at least 10% of the
#> sample
#> 
#> ℹ `sniff_violations()` for what was found, `$port` to extract a diagnostic,
#> `summary()` for every statistic.
```
