# Binary-exposure data with planted positivity violations

A simulated point-treatment dataset that carries two known positivity
problems, so examples and tests can compare a diagnostic against ground
truth. It holds one structural violation (a subgroup that is never
exposed) and one practical near-violation (a region of covariate space
where the fitted propensity approaches zero or one while both exposure
levels are still observed overall). It serves
[`check_positivity()`](https://r-causal.github.io/positively/reference/check_positivity.md)
and the binary-exposure diagnostics
[`check_port()`](https://r-causal.github.io/positively/reference/check_port.md),
[`check_edp()`](https://r-causal.github.io/positively/reference/check_edp.md),
and
[`check_extrapolation()`](https://r-causal.github.io/positively/reference/check_extrapolation.md).

## Usage

``` r
pos_violations
```

## Format

A tibble with 1000 rows and 4 columns:

- exposure:

  Integer treatment indicator, 0 or 1.

- x1:

  Numeric covariate, standard normal. Its steep effect on the propensity
  drives the practical near-violation.

- x2:

  Numeric covariate, standard normal. It defines the structural subgroup
  together with `region`.

- region:

  Factor covariate with levels `"a"` and `"b"`.

## Source

Simulated by `data-raw/make-datasets.R` in the package source
repository, <https://github.com/r-causal/positively>.

## Details

The structural violation is exact: no subject with `region == "b"` and
`x2 > 1` is exposed, because the true propensity is set to zero in that
subgroup. Both exposure levels occur everywhere outside it, so the empty
cell is a property of the subgroup rather than of the sample size.

The practical near-violation comes from the strong dependence of
exposure on `x1`. Fitting the propensity model
`glm(exposure ~ x1 + x2 + region, family = binomial())` recovers fitted
values that fall below 0.01 in the lower tail of `x1` and rise above
0.99 in the upper tail, yet both exposure levels remain observed. This
is the signature of a practical violation: treatment is
near-deterministic in the tails without being structurally impossible.
