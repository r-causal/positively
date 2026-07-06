# the print method is stable

    Code
      print(res)
    Output
      
      -- density_ratios_result -------------------------------------------------------
      Density ratios: 200 observations across 1 time point
      Mean ratio: 0.981
      50th percentile: 0.577
      90th percentile: 2.059
      95th percentile: 3.004
      99th percentile: 5.513
      Maximum: 6.697
      Proportion > 10: 0
      Proportion > 50: 0
      Proportion exactly zero: 0
      Kish ESS fraction: 0.439

# the multi-time print method is stable

    Code
      print(res)
    Output
      
      -- density_ratios_result -------------------------------------------------------
      Density ratios: 100 observations across 3 time points
      Summaries shown for time point 3
      Mean ratio: 1.018
      50th percentile: 0.882
      90th percentile: 1.567
      95th percentile: 2.055
      99th percentile: 3.079
      Maximum: 3.319
      Proportion > 10: 0
      Proportion > 50: 0
      Proportion exactly zero: 0
      Kish ESS fraction: 0.783
      Cumulative ESS fraction: 0.604

# the percentile print labels use correct ordinals

    Code
      print(res)
    Output
      
      -- density_ratios_result -------------------------------------------------------
      Density ratios: 200 observations across 1 time point
      Mean ratio: 0.981
      1st percentile: 0.089
      2nd percentile: 0.114
      3rd percentile: 0.13
      21st percentile: 0.316
      99.5th percentile: 6.101
      Maximum: 6.697
      Proportion > 10: 0
      Proportion > 50: 0
      Proportion exactly zero: 0
      Kish ESS fraction: 0.439

# the teens percentiles take a plain 'th' ordinal

    Code
      print(res)
    Output
      
      -- density_ratios_result -------------------------------------------------------
      Density ratios: 200 observations across 1 time point
      Mean ratio: 0.981
      11th percentile: 0.198
      12th percentile: 0.212
      13th percentile: 0.219
      Maximum: 6.697
      Proportion > 10: 0
      Proportion > 50: 0
      Proportion exactly zero: 0
      Kish ESS fraction: 0.439

# print omits the quantile block when probs is the maximum alone

    Code
      print(res)
    Output
      
      -- density_ratios_result -------------------------------------------------------
      Density ratios: 200 observations across 1 time point
      Mean ratio: 0.981
      Maximum: 6.697
      Proportion > 10: 0
      Proportion > 50: 0
      Proportion exactly zero: 0
      Kish ESS fraction: 0.439

# the classed error messages are stable

    Code
      check_density_ratios(c(1, -0.5, 2))
    Condition
      Error in `check_density_ratios()`:
      ! `ratios` must be non-negative.
      x Density ratios are Radon-Nikodym derivatives and cannot be negative.

---

    Code
      check_density_ratios(c(1, NA, 2))
    Condition
      Error in `check_density_ratios()`:
      ! `ratios` must not contain missing values.

---

    Code
      check_density_ratios(numeric(0))
    Condition
      Error in `check_density_ratios()`:
      ! `ratios` must contain at least one value.

---

    Code
      check_density_ratios(c("a", "b"))
    Condition
      Error in `check_density_ratios()`:
      ! `ratios` must be a numeric vector, a numeric matrix, or an lmtp fit.
      x A <character> is not supported.

---

    Code
      check_density_ratios(tibble::tibble(a = 1:3))
    Condition
      Error in `check_density_ratios()`:
      ! `ratios` must be a numeric vector, a numeric matrix, or an lmtp fit.
      x A <tbl_df> is not supported.

---

    Code
      check_density_ratios(c(1, 1, 1), probs = c(0.5, 1.5))
    Condition
      Error in `check_density_ratios()`:
      ! `probs` must be between 0 and 1.
      x Found 1.5.

---

    Code
      check_density_ratios(c(1, 1, 1), probs = c(0.5, 0.5))
    Condition
      Error in `check_density_ratios()`:
      ! `probs` must not contain duplicate values.
      x Found 0.5 more than once.

---

    Code
      check_density_ratios(c(1, 1, 1), thresholds = c(-1, 50))
    Condition
      Error in `check_density_ratios()`:
      ! `thresholds` must be positive.
      x Found -1.

---

    Code
      check_density_ratios(c(1, 1, 1), thresholds = "10")
    Condition
      Error in `check_density_ratios()`:
      ! `thresholds` must be numeric, not a <character>.

---

    Code
      check_density_ratios(c(1, 1, 1), thresholds = c(10, 10))
    Condition
      Error in `check_density_ratios()`:
      ! `thresholds` must not contain duplicate values.
      x Found 10 more than once.

---

    Code
      check_density_ratios(malformed)
    Condition
      Error in `check_density_ratios()`:
      ! The lmtp fit has no density-ratio component.
      i Expected a density_ratios element in the fitted object.
      i Pass the density ratios directly as a numeric vector or matrix.

---

    Code
      ggplot2::autoplot(res_point, type = "cumulative")
    Condition
      Error in `autoplot.positively::density_ratios_result`:
      ! `type = "cumulative"` requires a time-varying (matrix) input.
      i This result summarizes a point treatment with a single time point.

---

    Code
      check_density_ratios(c(1, 1, 1), thresholds = numeric(0))
    Condition
      Error in `check_density_ratios()`:
      ! `thresholds` must contain at least one value.

---

    Code
      check_density_ratios(c(1, 1, 1), thresholds = c(10, NA))
    Condition
      Error in `check_density_ratios()`:
      ! `thresholds` must not contain missing values.

---

    Code
      check_density_ratios(matrix(c("a", "b", "c", "d"), nrow = 2))
    Condition
      Error in `check_density_ratios()`:
      ! `ratios` must be numeric, not a <matrix>.

