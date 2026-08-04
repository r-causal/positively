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
      
      -- edp --
      
      -- test_diagnostic -------------------------------------------------------------
      Exposure: "a" (binary)
      Observations: 3
      Results: 3 rows, 2 columns
      
      -- port --
      
      -- test_diagnostic -------------------------------------------------------------
      Exposure: "a" (binary)
      Observations: 3
      Results: 3 rows, 2 columns

