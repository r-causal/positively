# check_hdr() rejects a Date covariate

    Code
      check_hdr(data, exposure, c(l, d))
    Condition
      Error in `check_hdr()`:
      ! `.covariates` must select numeric, logical, factor, or character columns.
      x "d" is another type.
      i Convert date or other ordered columns to numeric first.

# check_hdr_seq() rejects a Date covariate or baseline

    Code
      check_hdr_seq(data, c(a1, a2, a3), list(c(l1, d), l2, l3))
    Condition
      Error in `check_hdr_seq()`:
      ! `.covariates` must select numeric, logical, factor, or character columns.
      x "d" is another type.
      i Convert date or other ordered columns to numeric first.

---

    Code
      check_hdr_seq(data, c(a1, a2, a3), list(l1, l2, l3), .baseline = d)
    Condition
      Error in `check_hdr_seq()`:
      ! `.baseline` must select numeric, logical, factor, or character columns.
      x "d" is another type.
      i Convert date or other ordered columns to numeric first.

# the point print method is stable

    Code
      print(res)
    Output
      
      -- HDR non-overlap -------------------------------------------------------------
      Exposure: "exposure" (continuous)
      Observations: 300
      HDR mass: 0.95
      Density estimator: normal
      Non-overlap over 3 targets: 0.05 to 0.143

# the sequential print method is stable

    Code
      print(res)
    Output
      
      -- HDR non-overlap -------------------------------------------------------------
      Exposure: "a1", "a2", and "a3" (continuous)
      Observations: 300
      HDR mass: 0.95
      Density estimator: normal
      Time points: 3
      Non-overlap over 3 targets: 0 to 0.783

# check_hdr() argument validation messages are stable

    Code
      check_hdr(1:10, exposure, x1)
    Condition
      Error in `check_hdr()`:
      ! `.data` must be a data frame, not a <integer>.

---

    Code
      check_hdr(binary, exposure, c(x1, x2))
    Condition
      Error in `check_hdr()`:
      ! `check_hdr()` supports continuous exposures only.
      i `.exposure` was detected as "binary".
      i If `.exposure` is continuous, set `exposure_type = "continuous"`.

---

    Code
      check_hdr(data, exposure, tidyselect::starts_with("zzz"))
    Condition
      Error in `check_hdr()`:
      ! `.covariates` must select at least one column.

---

    Code
      check_hdr(data, c(exposure, l), l)
    Condition
      Error in `check_hdr()`:
      ! `.exposure` must select exactly one column, not 2.

---

    Code
      check_hdr(data, exposure, l, mass = 0)
    Condition
      Error in `check_hdr()`:
      ! `mass` must be between 0 and 1, exclusive.
      x Found 0.

---

    Code
      check_hdr(data, exposure, l, mass = 1.5)
    Condition
      Error in `check_hdr()`:
      ! `mass` must be between 0 and 1, exclusive.
      x Found 1.5.

---

    Code
      check_hdr(data, exposure, l, mass = c(0.9, 0.95))
    Condition
      Error in `check_hdr()`:
      ! `mass` must be a single number.

---

    Code
      check_hdr(data, exposure, l, values = "a")
    Condition
      Error in `check_hdr()`:
      ! `values` must be numeric, not a <character>.

---

    Code
      check_hdr(data, exposure, l, values = c(0, NA))
    Condition
      Error in `check_hdr()`:
      ! `values` must not contain missing values.

---

    Code
      check_hdr(data, exposure, l, values = numeric(0))
    Condition
      Error in `check_hdr()`:
      ! `values` must contain at least one value.

---

    Code
      check_hdr(data, exposure, l, values = c(0, Inf))
    Condition
      Error in `check_hdr()`:
      ! `values` must be finite.

---

    Code
      check_hdr(data, exposure, l, density_estimator = "not an estimator")
    Condition
      Error in `check_hdr()`:
      ! `density_estimator` must be an <hdr_density> object.
      i Build one with `hdr_density_normal()` or `new_hdr_density()`.
      x Received a <character>.

---

    Code
      check_hdr(na_covariate, exposure, l)
    Condition
      Error in `check_hdr()`:
      ! `.covariates` must not contain missing values.
      x Missing values in "l".

---

    Code
      check_hdr(one_row, exposure, l)
    Condition
      Error in `check_hdr()`:
      ! `.data` must have at least two observations, not 1.

---

    Code
      check_hdr(data, exposure, c(exposure, l))
    Condition
      Error in `check_hdr()`:
      ! `.covariates` must not include the exposure column.
      x "exposure" is also selected by `.exposure`.

---

    Code
      check_hdr(constant, exposure, l, exposure_type = "continuous")
    Condition
      Error in `check_hdr()`:
      ! `check_hdr()` needs `.exposure` to take more than one value.
      x `exposure` is constant.
      i The non-overlap ratio compares covariate profiles at a common target dose, and the target grid spans the observed exposure range. One observed dose leaves one target and a fitted density with no width, so the ratio would report that dose as supported at every profile.

# the flat-at-one non-overlap messages are stable

    Code
      res <- check_hdr(point, exposure, x, values = c(-1, 0, 1))
    Condition
      Warning in `check_hdr()`:
      The non-overlap ratio is 1 throughout the observed exposure range.
      i The exposure is close to a deterministic function of the columns it is conditioned on, so no covariate profile's highest-density region reaches another profile's dose.
      i The reading is correct, but it is set by the width of the fitted conditional density rather than by any comparison between profiles.

---

    Code
      res <- check_hdr_seq(sequential, c(a1, a2, a3), list(l0, l1, l2), values = c(-1,
        0, 1))
    Condition
      Warning in `check_hdr_seq()`:
      The non-overlap ratio is 1 throughout the observed exposure range of time points 2 and 3.
      i The exposure is close to a deterministic function of the columns it is conditioned on, so no covariate profile's highest-density region reaches another profile's dose.
      i The reading is correct, but it is set by the width of the fitted conditional density rather than by any comparison between profiles.

# check_hdr_seq() argument validation messages are stable

    Code
      check_hdr_seq(data, c(a1, a2, a3), list(l0, l1))
    Condition
      Error in `check_hdr_seq()`:
      ! `.covariates` must supply one selection per time point.
      i Found 2 selections for 3 exposures.

---

    Code
      check_hdr_seq(data, c(a1, a2, a3), list(l0, l1, l2), lag = -1)
    Condition
      Error in `check_hdr_seq()`:
      ! `lag` must be a single non-negative whole number or `Inf`.
      x Found -1.

---

    Code
      check_hdr_seq(data, c(a1, a2, a3), list(l0, l1, l2), lag = 1.5)
    Condition
      Error in `check_hdr_seq()`:
      ! `lag` must be a single non-negative whole number or `Inf`.
      x Found 1.5.

---

    Code
      check_hdr_seq(data, c(a1, a2, a3), list(l0, l1, l2), lag = "a")
    Condition
      Error in `check_hdr_seq()`:
      ! `lag` must be a single non-negative whole number or `Inf`.
      x Found "a".

---

    Code
      check_hdr_seq(binary, c(a1, a2, a3), list(l0, l1, l2))
    Condition
      Error in `check_hdr_seq()`:
      ! `check_hdr_seq()` supports continuous exposures only.
      i `a2` was detected as "binary".
      i If `a2` is continuous, set `exposure_type = "continuous"`.

---

    Code
      check_hdr_seq(binary_pair, c(a1, a2, a3), list(l0, l1, l2))
    Condition
      Error in `check_hdr_seq()`:
      ! `check_hdr_seq()` supports continuous exposures only.
      i `a2` and `a3` were detected as "binary".
      i If `a2` and `a3` are continuous, set `exposure_type = "continuous"`.

---

    Code
      check_hdr_seq(factors, c(a1, a2, a3), list(l0, l1, l2), exposure_type = "continuous")
    Condition
      Error in `check_hdr_seq()`:
      ! `check_hdr_seq()` needs a numeric `.exposures` for a continuous exposure.
      x `a2` and `a3` are <factor>.

---

    Code
      check_hdr_seq(mixed, c(a1, a2, a3), list(l0, l1, l2), exposure_type = "continuous")
    Condition
      Error in `check_hdr_seq()`:
      ! `check_hdr_seq()` needs a numeric `.exposures` for a continuous exposure.
      x `a1` is <factor>.
      x `a2` is <character>.
      x `a3` is <Date>.

---

    Code
      check_hdr_seq(ranked, c(a1, a2, a3), list(l0, l1, l2), exposure_type = "continuous")
    Condition
      Error in `check_hdr_seq()`:
      ! `check_hdr_seq()` needs a numeric `.exposures` for a continuous exposure.
      x `a1` is <ordered/factor>.
      x `a2` is <factor>.

---

    Code
      check_hdr_seq(data, c(a1, a2), list(c(l0, a1), l1), values = 0)
    Condition
      Error in `check_hdr_seq()`:
      ! `.covariates` and `.baseline` must not include an exposure at its own time point.
      x "a1" conditions its own model at time 1.

---

    Code
      check_hdr_seq(data, c(a1, a2, a3), list(l0, l1, l2), .baseline = a2, values = 0)
    Condition
      Error in `check_hdr_seq()`:
      ! `.covariates` and `.baseline` must not include an exposure at its own time point.
      x "a2" conditions its own model at time 2.

---

    Code
      check_hdr_seq(data, c(a1, a2, a3), list(tidyselect::starts_with("zzz")),
      values = 0)
    Condition
      Error in `check_hdr_seq()`:
      ! Every time point must have at least one conditioning column.
      x The conditioning set at time 1 is empty.
      i Supply time-varying covariates in `.covariates` or baseline covariates in `.baseline`.

---

    Code
      check_hdr_seq(constant, c(a1, a2, a3), list(l0, l1, l2), exposure_type = "continuous",
      values = 0)
    Condition
      Error in `check_hdr_seq()`:
      ! `check_hdr_seq()` needs `.exposures` to take more than one value.
      x `a1` and `a3` are constant.
      i The non-overlap ratio compares covariate profiles at a common target dose, and the target grid spans the observed exposure range. One observed dose leaves one target and a fitted density with no width, so the ratio would report that dose as supported at every profile.

# the HDR label and headline are stable

    Code
      diagnostic_label(res)
    Output
      [1] "HDR non-overlap"
    Code
      writeLines(diagnostic_headline(res))
    Output
      95% HDR, normal density; non-overlap 0.05 to 0.143 over 3 targets

