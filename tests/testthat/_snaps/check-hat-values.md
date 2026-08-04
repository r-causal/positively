# the continuous-only abort message is stable

    Code
      check_hat_values(binary, exposure, c(x1, x2))
    Condition
      Error in `check_hat_values()`:
      ! `check_hat_values()` supports continuous exposures only.
      i `.exposure` was detected as "binary".
      i If `.exposure` is continuous, set `exposure_type = "continuous"`.

---

    Code
      check_hat_values(categorical, exposure, x1)
    Condition
      Error in `check_hat_values()`:
      ! `check_hat_values()` supports continuous exposures only.
      i `.exposure` was detected as "categorical".
      i If `.exposure` is continuous, set `exposure_type = "continuous"`.

# the argument validation messages are stable

    Code
      check_hat_values(1:10, dose, x1)
    Condition
      Error in `check_hat_values()`:
      ! `.data` must be a data frame, not a <integer>.

---

    Code
      check_hat_values(data, dose, tidyselect::starts_with("zzz"))
    Condition
      Error in `check_hat_values()`:
      ! `.covariates` must select at least one column.

---

    Code
      check_hat_values(data, dose, x1, probs = c(0.5, 1.5))
    Condition
      Error in `check_hat_values()`:
      ! `probs` must be between 0 and 1.
      x Found 1.5.

---

    Code
      check_hat_values(data, dose, x1, conf_level = 1.5)
    Condition
      Error in `check_hat_values()`:
      ! `conf_level` must be between 0 and 1, exclusive.
      x Found 1.5.

---

    Code
      check_hat_values(data, dose, x1, null_method = "uniform")
    Condition
      Error in `check_hat_values()`:
      ! `null_method` must be one of "permutation", "bootstrap", or "gaussian", not "uniform".

---

    Code
      check_hat_values(data, c(dose, x1), x1)
    Condition
      Error in `check_hat_values()`:
      ! `.exposure` must select exactly one column, not 2.

---

    Code
      check_hat_values(data, dose, x1, conf_level = c(0.5, 0.9))
    Condition
      Error in `check_hat_values()`:
      ! `conf_level` must be a single number.

---

    Code
      check_hat_values(data, dose, x1, threshold = 0)
    Condition
      Error in `check_hat_values()`:
      ! `threshold` must be positive.
      x Found 0.

---

    Code
      check_hat_values(data, dose, x1, null_reps = 0)
    Condition
      Error in `check_hat_values()`:
      ! `null_reps` must be a single whole number of at least 1.
      x Found 0.

---

    Code
      check_hat_values(data, dose, x1, probs = numeric(0))
    Condition
      Error in `check_hat_values()`:
      ! `probs` must contain at least one value.

# the data-integrity error messages are stable

    Code
      check_hat_values(na_covariate, dose, x1, null_reps = 2)
    Condition
      Error in `check_hat_values()`:
      ! `.covariates` must not contain missing values.
      x Missing values in "x1".

---

    Code
      check_hat_values(factor_covariate, dose, group, null_reps = 2)
    Condition
      Error in `check_hat_values()`:
      ! `.covariates` must select numeric columns.
      x "group" is not numeric.
      i Encode factor or character covariates as numeric indicators first.

---

    Code
      check_hat_values(constant_covariate, dose, constant, null_reps = 2)
    Condition
      Error in `check_hat_values()`:
      ! The design matrix formed from `.exposure` and `.covariates` is not full rank.
      i Check for constant or collinear columns in `.covariates`.

# the print method is stable

    Code
      print(res)
    Output
      
      -- Hat values ------------------------------------------------------------------
      Exposure: "dose" (continuous)
      Observations: 150
      Null: permutation (50 reps), cutoff 2p/n
      phi-hat: 0.188
      Null 0.95 quantile: 0.059
      Exceeds null: TRUE
      High-leverage candidates: 537 of 2850

# the leverage label and headline are stable

    Code
      diagnostic_label(res)
    Output
      [1] "Hat values"
    Code
      writeLines(diagnostic_headline(res))
    Output
      phi-hat 0.188 exceeds the 0.95 permutation-null quantile of 0.059
      537 of 2850 unit-value pairs above the 2p/n cutoff (50 reps)

