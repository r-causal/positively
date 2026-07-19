# check_eta_bias() rejects an unsupported covariate type

    Code
      check_eta_bias(data, a, y, c(x1, bad), n_boot = 10)
    Condition
      Error in `check_eta_bias()`:
      ! `.covariates` must select numeric, logical, factor, or character columns.
      x "bad" is of an unsupported type.

# the sweep view aborts on a single-run result

    Code
      ggplot2::autoplot(single, type = "sweep")
    Condition
      Error in `autoplot_eta_bias_sweep()`:
      ! The sweep view needs a truncation sweep of more than one level.
      i Rerun `check_eta_bias()` with a `truncation_grid` of multiple lower bounds.

# the single-run print method is stable

    Code
      print(res)
    Output
      
      -- eta_bias_result -------------------------------------------------------------
      Exposure: "a" (binary)
      Observations: 300
      Estimator: ipw
      Bootstrap draws: 50
      Truth: 0.623
      ETA.Bias: -0.003 (MC SE 0.022)

# the truncation-sweep print method is stable

    Code
      print(res)
    Output
      
      -- eta_bias_result -------------------------------------------------------------
      Exposure: "a" (binary)
      Observations: 400
      Estimator: ipw
      Bootstrap draws: 50
      Truth: 1.064
      Truncation levels: 3
      ETA.Bias: 0.239 to 0.748

# the input and exposure validation messages are stable

    Code
      check_eta_bias(1:10, a, y, c(x1, x2), n_boot = 10)
    Condition
      Error in `check_eta_bias()`:
      ! `.data` must be a data frame, not a <integer>.

---

    Code
      check_eta_bias(continuous, a, y, c(x1, x2), n_boot = 10)
    Condition
      Error in `check_eta_bias()`:
      ! `check_eta_bias()` supports binary exposures only.
      i `.exposure` was detected as "continuous".

---

    Code
      check_eta_bias(categorical, a, y, c(x1, x2), n_boot = 10)
    Condition
      Error in `check_eta_bias()`:
      ! `check_eta_bias()` supports binary exposures only.
      i `.exposure` was detected as "categorical".

---

    Code
      check_eta_bias(good, a, y, tidyselect::starts_with("zzz"), n_boot = 10)
    Condition
      Error in `check_eta_bias()`:
      ! `.covariates` must select at least one column.

---

    Code
      check_eta_bias(good, a, not_a_column, c(x1, x2), n_boot = 10)
    Condition
      Error in `check_eta_bias()`:
      ! `.outcome` must select a column that exists in `.data`.
      Caused by error in `eval_select_column()`:
      ! Can't select columns that don't exist.
      x Column `not_a_column` doesn't exist.

---

    Code
      check_eta_bias(na_exposure, a, y, c(x1, x2), n_boot = 10)
    Condition
      Error in `check_eta_bias()`:
      ! `.exposure` must not contain missing values.

---

    Code
      check_eta_bias(good, a, y, c(x1, x2), outcome_type = "binary", n_boot = 10)
    Condition
      Error in `check_eta_bias()`:
      ! `.outcome` must contain only 0 and 1 when `outcome_type` is "binary".
      x Found other values, including -0.610598594863679, 1.13292680079015, -1.8162050739754, 1.50564523292254, and 0.370473934516183.

---

    Code
      check_eta_bias(good, a, y, c(a, x1, x2), n_boot = 10)
    Condition
      Error in `check_eta_bias()`:
      ! `.covariates` must not include the exposure or outcome column.
      x "a" is also selected by `.exposure` or `.outcome`.

# the truncation and n_boot validation messages are stable

    Code
      check_eta_bias(data, a, y, c(x1, x2), truncation = 0.05, n_boot = 10)
    Condition
      Error in `check_eta_bias()`:
      ! `truncation` must be a length-two numeric vector `c(lower, upper)`.

---

    Code
      check_eta_bias(data, a, y, c(x1, x2), truncation = c(-0.1, 0.95), n_boot = 10)
    Condition
      Error in `check_eta_bias()`:
      ! `truncation` must lie between 0 and 1.
      x Found -0.1 and 0.95.

---

    Code
      check_eta_bias(data, a, y, c(x1, x2), truncation = c(0.6, 0.4), n_boot = 10)
    Condition
      Error in `check_eta_bias()`:
      ! `truncation` must have its lower bound below its upper bound.
      x Found lower 0.6 and upper 0.4.

---

    Code
      check_eta_bias(data, a, y, c(x1, x2), truncation_grid = c(0.1, 0.6), n_boot = 10)
    Condition
      Error in `check_eta_bias()`:
      ! `truncation_grid` lower bounds must lie in `[0, 0.5)`.
      x Found 0.6.

---

    Code
      check_eta_bias(data, a, y, c(x1, x2), truncation_grid = c(0.1, NA), n_boot = 10)
    Condition
      Error in `check_eta_bias()`:
      ! `truncation_grid` must not contain missing values.

---

    Code
      check_eta_bias(data, a, y, c(x1, x2), n_boot = 2.5)
    Condition
      Error in `check_eta_bias()`:
      ! `n_boot` must be a single whole number of at least 1.
      x Found 2.5.

