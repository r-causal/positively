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

