# check_extrapolation() rejects a Date covariate

    Code
      check_extrapolation(data, exposure, c(x1, d), hull = FALSE)
    Condition
      Error in `check_extrapolation()`:
      ! `.covariates` must select numeric, logical, factor, or character columns.
      x "d" is another type.
      i Convert date or other ordered columns to numeric first.

# the binary-only abort messages are stable

    Code
      check_extrapolation(continuous, exposure, x1)
    Condition
      Error in `check_extrapolation()`:
      ! `check_extrapolation()` supports binary exposures only.
      i `.exposure` was detected as "continuous".

---

    Code
      check_extrapolation(categorical, exposure, c(x1, x2))
    Condition
      Error in `check_extrapolation()`:
      ! `check_extrapolation()` supports binary exposures only.
      i `.exposure` was detected as "categorical".

# the argument validation messages are stable

    Code
      check_extrapolation(1:10, exposure, x1)
    Condition
      Error in `check_extrapolation()`:
      ! `.data` must be a data frame, not a <integer>.

---

    Code
      check_extrapolation(data, exposure, tidyselect::starts_with("zzz"))
    Condition
      Error in `check_extrapolation()`:
      ! `.covariates` must select at least one column.

---

    Code
      check_extrapolation(data, exposure, c(x1, x2), nearby = 0)
    Condition
      Error in `check_extrapolation()`:
      ! `nearby` must be positive.
      x Found 0.

---

    Code
      check_extrapolation(data, exposure, c(x1, x2), hull = "banana")
    Condition
      Error in `check_extrapolation()`:
      ! `hull` must be `NULL`, `TRUE`, or `FALSE`.

---

    Code
      check_extrapolation(data, c(exposure, x1), x2)
    Condition
      Error in `check_extrapolation()`:
      ! `.exposure` must select exactly one column, not 2.

# the missing-value and sample-size abort messages are stable

    Code
      check_extrapolation(na_exposure, exposure, c(x1, x2))
    Condition
      Error in `check_extrapolation()`:
      ! `.exposure` must not contain missing values.

---

    Code
      check_extrapolation(na_covariate, exposure, c(x1, x2))
    Condition
      Error in `check_extrapolation()`:
      ! `.covariates` must not contain missing values.
      x Missing values in "x1".

---

    Code
      check_extrapolation(too_small, exposure, c(x1, x2))
    Condition
      Error in `check_extrapolation()`:
      ! `.data` must have at least two observations, not 1.

---

    Code
      check_extrapolation(non_finite, exposure, c(x1, x2), hull = FALSE)
    Condition
      Error in `check_extrapolation()`:
      ! `.covariates` must not contain non-finite values.
      x Non-finite values in "x1".

# the hull warning and error messages are stable

    Code
      res <- check_extrapolation(mixed, exposure, g, hull = TRUE)
    Condition
      Warning in `check_extrapolation()`:
      The convex-hull test needs at least one numeric covariate; skipping it.

---

    Code
      check_extrapolation(gaussian, exposure, tidyselect::starts_with("x"), hull = TRUE)
    Condition
      Error in `check_extrapolation()`:
      ! The convex-hull test requires the lpSolve package.
      i Install lpSolve or set `hull = FALSE`.

# the hull-view abort without a hull run is stable

    Code
      ggplot2::autoplot(res, type = "hull")
    Condition
      Error in `autoplot.positively::extrapolation_result`:
      ! The convex-hull view needs the hull test to have run.
      i Rerun `check_extrapolation()` with `hull = TRUE`.

# the hull-view abort points to the numeric-covariate precondition when hull was skipped

    Code
      ggplot2::autoplot(res, type = "hull")
    Condition
      Error in `autoplot.positively::extrapolation_result`:
      ! The convex-hull view needs the hull test to have run.
      i The hull test needs at least one numeric covariate.

# the high-dimensional hull messages are stable

    Code
      check_extrapolation(data, exposure, tidyselect::starts_with("x"))
    Message
      i Treating `.exposure` as binary
      i Skipping the convex-hull test: 13 numeric covariates exceed the maximum of 10.
    Output
      
      -- extrapolation_result --------------------------------------------------------
      Exposure: "exposure" (binary)
      Observations: 60
      Geometric variability: 0.121
      Nearby radius (1 x gv): 0.121
      Mean fraction nearby: 0.002
      Nearest opposite within one geometric variability: 4 of 60
      Convex-hull test: not run

---

    Code
      check_extrapolation(data, exposure, tidyselect::starts_with("x"), hull = TRUE)
    Message
      i Treating `.exposure` as binary
    Condition
      Warning in `check_extrapolation()`:
      The convex-hull test is degenerate above 10 numeric covariates.
      i 13 numeric covariates were supplied; membership will be near zero.
    Output
      
      -- extrapolation_result --------------------------------------------------------
      Exposure: "exposure" (binary)
      Observations: 60
      Geometric variability: 0.121
      Nearby radius (1 x gv): 0.121
      Mean fraction nearby: 0.002
      Nearest opposite within one geometric variability: 4 of 60
      In opposite-group hull: 0 of 60

# the print method is stable

    Code
      print(res)
    Output
      
      -- extrapolation_result --------------------------------------------------------
      Exposure: "exposure" (binary)
      Observations: 200
      Geometric variability: 0.108
      Nearby radius (1 x gv): 0.108
      Mean fraction nearby: 0.076
      Nearest opposite within one geometric variability: 188 of 200
      In opposite-group hull: 129 of 200

