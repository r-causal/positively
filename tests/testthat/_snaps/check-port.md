# the degenerate Gruber bound message pins the bound and sample size

    Code
      check_port(data, exposure, x1, beta = "gruber")
    Condition
      Error in `check_port()`:
      ! The "gruber" prevalence threshold resolved to 0.580857101049332 at n = 12.
      i The reading rule needs a threshold below 0.5, but the sample-size-adaptive bound reaches 0.5 or more at small sample sizes.
      i Supply a numeric `beta` below 0.5 for a sample this small.

# check_port() rejects a dot that is not an rpart.control option

    Code
      check_port(data, exposure, c(x3, x1), alpa = 0.4)
    Condition
      Error in `check_port()`:
      ! `...` must name `rpart::rpart.control()` options only.
      x Unknown option: "alpa".
      i Valid options are "minsplit", "minbucket", "cp", "maxcompete", "maxsurrogate", "usesurrogate", "xval", "surrogatestyle", and "maxdepth".

# the point print method is stable

    Code
      print(res)
    Output
      
      -- port_result -----------------------------------------------------------------
      Exposure: "exposure" (binary)
      Observations: 1000
      Reading rule: alpha = 0.05, gamma = 2
      Prevalence threshold beta: 0.05
      Subgroups: 2 reported, 1 flagged

# the print method reports an unresolved beta

    Code
      print(res)
    Output
      
      -- port_result -----------------------------------------------------------------
      Exposure: "exposure" (binary)
      Observations: 1
      Reading rule: alpha = 0.05, gamma = 2
      Prevalence threshold beta: not resolved
      Subgroups: 0 reported, 0 flagged

# check_port() argument validation messages are stable

    Code
      check_port(1:10, exposure, x1)
    Condition
      Error in `check_port()`:
      ! `.data` must be a data frame, not a <integer>.

---

    Code
      check_port(data, exposure, tidyselect::starts_with("zzz"))
    Condition
      Error in `check_port()`:
      ! `.covariates` must select at least one column.

---

    Code
      check_port(data, c(exposure, x1), x2)
    Condition
      Error in `check_port()`:
      ! `.exposure` must select exactly one column, not 2.

---

    Code
      check_port(single_level, exposure, c(x1, x2))
    Condition
      Error in `check_port()`:
      ! `.exposure` must have at least two distinct values.
      x Found 1.

---

    Code
      check_port(data, exposure, c(x1, x2), alpha = 0)
    Condition
      Error in `check_port()`:
      ! `alpha` must be between 0 and 1, exclusive.
      x Found 0.

---

    Code
      check_port(data, exposure, c(x1, x2), alpha = 1.5)
    Condition
      Error in `check_port()`:
      ! `alpha` must be between 0 and 1, exclusive.
      x Found 1.5.

---

    Code
      check_port(data, exposure, c(x1, x2), alpha = c(0.05, 0.1))
    Condition
      Error in `check_port()`:
      ! `alpha` must be a single number.

---

    Code
      check_port(data, exposure, c(x1, x2), beta = 0)
    Condition
      Error in `check_port()`:
      ! `beta` must be between 0 and 1, exclusive.
      x Found 0.

---

    Code
      check_port(data, exposure, c(x1, x2), beta = 1.5)
    Condition
      Error in `check_port()`:
      ! `beta` must be between 0 and 1, exclusive.
      x Found 1.5.

---

    Code
      check_port(data, exposure, c(x1, x2), beta = 0.6)
    Condition
      Error in `check_port()`:
      ! `beta` must be below 0.5.
      i The reading rule flags prevalence below `beta` or above `1 - beta`, so a value of 0.5 or more flags every subgroup.
      x Found 0.6.

---

    Code
      check_port(data, exposure, c(x1, x2), beta = "banana")
    Condition
      Error in `check_port()`:
      ! `beta` must be a number in (0, 0.5) or "gruber".
      x Received "banana".

---

    Code
      check_port(data, exposure, c(x1, x2), gamma = 0)
    Condition
      Error in `check_port()`:
      ! `gamma` must be a single whole number of at least 1.
      x Found 0.

---

    Code
      check_port(data, exposure, c(x1, x2), gamma = 2.5)
    Condition
      Error in `check_port()`:
      ! `gamma` must be a single whole number of at least 1.
      x Found 2.5.

---

    Code
      check_port(data, exposure, c(x1, x2), n_bins = 0)
    Condition
      Error in `check_port()`:
      ! `n_bins` must be a single whole number of at least 2.
      x Found 0.

---

    Code
      check_port(data, exposure, c(x1, x2), breaks = "a")
    Condition
      Error in `check_port()`:
      ! `breaks` must be numeric or `NULL`, not a <character>.

# check_port() missing-exposure message is stable

    Code
      check_port(port_missing_exposure_data(), exposure, g)
    Condition
      Error in `check_port()`:
      ! `.exposure` must not contain missing values.

# check_port() degenerate binning messages are stable

    Code
      check_port(port_binning_data(), exposure, x1, breaks = 4, exposure_type = "continuous")
    Condition
      Error in `check_port()`:
      ! `breaks` must contain at least three distinct cut points.
      i The cut points bound the exposure bins, so three points define the minimum of two bins; a single bin would place every observation at one exposure level and flag every subgroup.
      x Found 1 distinct cut point.

---

    Code
      check_port(port_binning_data(), exposure, x1, breaks = c(0, 6), exposure_type = "continuous")
    Condition
      Error in `check_port()`:
      ! `breaks` must contain at least three distinct cut points.
      i The cut points bound the exposure bins, so three points define the minimum of two bins; a single bin would place every observation at one exposure level and flag every subgroup.
      x Found 2 distinct cut points.

---

    Code
      check_port(port_binning_data(), exposure, x1, breaks = c(0, NA, 6),
      exposure_type = "continuous")
    Condition
      Error in `check_port()`:
      ! `breaks` must not contain missing values.

---

    Code
      check_port(port_zero_inflated_data(), exposure, x1, n_bins = 3, exposure_type = "continuous")
    Condition
      Error:
      ! The exposure quantiles define fewer than two bins.
      i Ties in the exposure collapse the quantile cut points, so the requested bins reduce to one.
      i Supply explicit `breaks` that separate the tied values.

# check_port() missing-covariate message is stable

    Code
      check_port(data, exposure, c(x1, x2))
    Condition
      Error in `check_port()`:
      ! `.covariates` must not contain missing values.
      x Missing values in "x1".

