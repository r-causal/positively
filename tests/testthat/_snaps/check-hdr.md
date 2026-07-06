# the point print method is stable

    Code
      print(res)
    Output
      
      -- hdr_result ------------------------------------------------------------------
      Exposure: "exposure" (continuous)
      Observations: 300
      HDR mass: 0.95
      Density estimator: normal
      Non-overlap over 3 targets: 0.05 to 0.143

# the sequential print method is stable

    Code
      print(res)
    Output
      
      -- hdr_result ------------------------------------------------------------------
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
      x "a2" is not continuous.

