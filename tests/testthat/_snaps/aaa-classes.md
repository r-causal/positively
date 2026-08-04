# the invalid-index message is stable

    Code
      res[[1.9]]
    Condition
      Error in `res[[1.9]]`:
      ! `i` must be a diagnostic name or a whole-number position.

# the unknown-diagnostic message is stable

    Code
      res[["nonexistent"]]
    Condition
      Error in `res[["nonexistent"]]`:
      ! "nonexistent" is not a diagnostic in this container.
      i Available diagnostics are "edp" and "port".

---

    Code
      res$nonexistent
    Condition
      Error in `res$nonexistent`:
      ! "nonexistent" is not a diagnostic in this container.
      i Available diagnostics are "edp" and "port".

# printing a diagnostic is stable

    Code
      print(make_test_diagnostic())
    Output
      
      -- test_diagnostic -------------------------------------------------------------
      Exposure: "a" (binary)
      Observations: 3
      Results: 3 rows, 2 columns

# printing a container with two children is stable

    Code
      print(container)
    Output
      
      -- Positivity check ------------------------------------------------------------
      Exposure: "a" (binary); 3 observations; covariates x1 and x2
      
      -- edp -------------------------------------------------------------------------
      3 rows, 2 columns
      
      -- port ------------------------------------------------------------------------
      3 rows, 2 columns
      
      i `sniff_violations()` for what was found, `$edp` to extract a diagnostic,
      `summary()` for every statistic.

# the report is stable

    Code
      print(container)
    Output
      
      -- Positivity check ------------------------------------------------------------
      Exposure: "statin" (binary); 317 observations; covariates bmi and ldl
      
      -- port ------------------------------------------------------------------------
      45 subgroups reported, 4 with low support
      Rule: prevalence outside [0.05, 0.95] among subgroups of at least 5% of the
      sample
      
      -- extrapolation ---------------------------------------------------------------
      Geometric variability 0.105
      278 of 317 have an opposite-exposure unit within one geometric variability; 173
      of 317 fall in the opposite hull
      
      -- edp -------------------------------------------------------------------------
      Data variant over 2 intervention values; edp 0.238 to 78.554
      
      i `sniff_violations()` for what was found, `$port` to extract a diagnostic,
      `summary()` for every statistic.

