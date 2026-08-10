# Longitudinal continuous-exposure data with a time-2 support gap

A simulated wide-format longitudinal dataset with one row per subject
over three time points. It carries a known support gap in the time-2
exposure, so the sequential and continuous-exposure diagnostics can be
exercised against ground truth. It serves
[`check_hdr()`](https://r-causal.github.io/positively/reference/check_hdr.md),
[`check_hdr_seq()`](https://r-causal.github.io/positively/reference/check_hdr_seq.md),
and
[`check_port_seq()`](https://r-causal.github.io/positively/reference/check_port_seq.md).

## Usage

``` r
pos_violations_long
```

## Format

A tibble with 500 rows and 7 columns:

- id:

  Integer subject identifier, 1 to 500.

- l0:

  Numeric baseline covariate.

- a1:

  Numeric time-1 exposure.

- l1:

  Numeric time-1 covariate.

- a2:

  Numeric time-2 exposure, carrying the support gap.

- l2:

  Numeric time-2 covariate.

- a3:

  Numeric time-3 exposure.

## Source

Simulated by `data-raw/make-datasets.R` in the package source
repository, <https://github.com/r-causal/positively>.

## Details

The exposures `a1`, `a2`, and `a3` are continuous, each taking hundreds
of distinct values. The planted violation lives at time 2: `a2` is drawn
on `[0, 2]` or `[4, 6]`, so it never falls in the open interval
`(2, 4)`. The time-1 and time-3 exposures are unconstrained and place
values throughout that interval, which makes the gap specific to time 2
rather than a feature of the exposure distribution as a whole. The
covariates and exposures follow a simple time-ordered process in which
each variable depends on the one before it.
