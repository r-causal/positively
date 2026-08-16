# the n_boot floor message is stable

    Code
      check_eta_bias(data, a, y, x1, n_boot = 1)
    Condition
      Error in `check_eta_bias()`:
      ! `n_boot` must be a single whole number of at least 2.
      x Found 1.

# check_eta_bias() rejects an unsupported covariate type

    Code
      check_eta_bias(data, a, y, c(x1, bad), n_boot = 10)
    Condition
      Error in `check_eta_bias()`:
      ! `.covariates` must select numeric, logical, factor, or character columns.
      x "bad" is of an unsupported type.

# the degenerate bootstrap warning wording is stable

    Code
      res <- withr::with_seed(1, check_eta_bias(data, a, y, x1, estimator = "ipw",
        n_boot = 50))
    Condition
      Warning in `check_eta_bias()`:
      Dropped 2 of 50 bootstrap draws with a non-finite estimate.
      i Such draws arise when a bootstrap exposure leaves a level empty, when a fitted probability is exactly 0 or 1, or when a weight is degenerate.

# the continuous truncation grid messages are stable

    Code
      check_eta_bias(data, a, y, c(x1, x2), truncation_grid = c(0.1, 1.5), n_boot = 10)
    Condition
      Error in `check_eta_bias()`:
      ! `truncation_grid` quantile levels must lie in `[0, 1)`.
      i A level of 1 would cap every stabilized weight at the smallest one observed.
      x Found 1.5.

---

    Code
      check_eta_bias(data, a, y, c(x1, x2), truncation_grid = numeric(0), n_boot = 10)
    Condition
      Error in `check_eta_bias()`:
      ! `truncation_grid` must contain at least one quantile level.

# the sweep view aborts on a single-run result

    Code
      ggplot2::autoplot(single, type = "sweep")
    Condition
      Error in `autoplot.positively::eta_bias_result`:
      ! The sweep view needs a truncation sweep of more than one level.
      i Rerun `check_eta_bias()` with a `truncation_grid` of multiple lower bounds.

# the single-run print method is stable

    Code
      print(res)
    Output
      
      -- ETA bias --------------------------------------------------------------------
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
      
      -- ETA bias --------------------------------------------------------------------
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
      check_eta_bias(one_level, a, y, c(x1, x2), exposure_type = "categorical",
      n_boot = 10)
    Condition
      Error in `check_eta_bias()`:
      ! `.exposure` must have at least two levels.
      x It has 1 level.

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
      check_eta_bias(na_covariate, a, y, c(x1, x2), n_boot = 10)
    Condition
      Error in `check_eta_bias()`:
      ! `.covariates` must not contain missing values.
      x Missing values in "x1".

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
      i The ceiling is `1/k` for an exposure with 2 levels.
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
      ! `n_boot` must be a single whole number of at least 2.
      x Found 2.5.

# the ETA bias label and headline are stable

    Code
      diagnostic_label(res)
    Output
      [1] "ETA bias"
    Code
      writeLines(diagnostic_headline(res))
    Output
      ipw estimator against a truth of 0.623; ETA.Bias -0.003 (MC SE 0.022) from 50 bootstrap draws

# a three-level exposure without nnet aborts

    Code
      check_eta_bias(data, a, y, c(x1, x2), exposure_type = "categorical", n_boot = 10)
    Condition
      Error in `check_eta_bias()`:
      ! A categorical exposure of more than two levels requires the nnet package.
      i Install nnet to run `check_eta_bias()` on this exposure.

# reference_level rejects an unknown or non-scalar level

    Code
      check_eta_bias(data, a, y, c(x1, x2), exposure_type = "categorical",
      reference_level = "zzz", n_boot = 10)
    Condition
      Error in `check_eta_bias()`:
      ! `reference_level` must be one of the exposure's levels.
      x "zzz" is not among "a", "b", and "c".

---

    Code
      check_eta_bias(data, a, y, c(x1, x2), exposure_type = "categorical",
      reference_level = c("a", "b"), n_boot = 10)
    Condition
      Error in `check_eta_bias()`:
      ! `reference_level` must be a single value.

# the two-sided truncation bound is rejected past two levels

    Code
      check_eta_bias(data, a, y, c(x1, x2), truncation = c(0.05, 0.95), n_boot = 10)
    Condition
      Error in `check_eta_bias()`:
      ! `truncation` bounds a fitted probability from both sides, which describes a two-level exposure only.
      x `.exposure` has 3 levels.
      i Use `truncation_grid` to raise every level from below.

# a multi-term run prints one reading per term

    Code
      print(res)
    Output
      
      -- ETA bias --------------------------------------------------------------------
      Exposure: "a" (categorical)
      Observations: 300
      Estimator: ipw
      Bootstrap draws: 50
      Estimand terms: 2
      ETA.Bias (b - a): 0.093 (MC SE 0.066)
      ETA.Bias (c - a): 0.142 (MC SE 0.065)

# the continuous print method is stable

    Code
      print(res)
    Output
      
      -- ETA bias --------------------------------------------------------------------
      Exposure: "a" (continuous)
      Observations: 300
      Estimator: ipw
      Bootstrap draws: 30
      Truth: 0.945
      ETA.Bias: 0.213 (MC SE 0.034)

# the continuous guard messages are stable

    Code
      check_eta_bias(data, a, y, c(x1, x2), estimator = "aipw", n_boot = 10)
    Condition
      Error in `check_eta_bias()`:
      ! `estimator` "aipw" is not available for a continuous exposure.
      i The augmentation term integrates the outcome residual over the exposure grid, which this path does not compute.
      i Use "ipw" or "gcomp".

---

    Code
      check_eta_bias(data, a, y, c(x1, x2), truncation = c(0.05, 0.95), n_boot = 10)
    Condition
      Error in `check_eta_bias()`:
      ! `truncation` is a pair of probability bounds, and a continuous exposure weights by a density ratio rather than by a probability.
      i Use `truncation_grid` to cap the stabilized weight at a quantile of its own distribution.

---

    Code
      check_eta_bias(data, a, y, c(x1, x2), reference_level = 0, n_boot = 10)
    Condition
      Error in `check_eta_bias()`:
      ! `reference_level` names the level every contrast is taken against, which a continuous exposure does not have.
      i A continuous estimand is the coefficient set of the working model set by `msm_formula`.

# the function-estimator argument messages are stable

    Code
      check_eta_bias(data, a, y, c(x1, x2), estimator = list(estimate), n_boot = 10)
    Condition
      Error in `check_eta_bias()`:
      ! `estimator` must be the name of a built-in estimator, a function, or a length-one list whose name labels the function it holds.
      i The built-in estimators are "ipw", "gcomp", and "aipw".
      i A bare function is labeled "custom".
      x Found a <list>.

---

    Code
      check_eta_bias(data, a, y, c(x1, x2), estimator = estimate, truncation = c(0.05,
        0.95), n_boot = 10)
    Condition
      Error in `check_eta_bias()`:
      ! `truncation` bounds a fitted nuisance model, which an estimator of your own never receives.
      x `estimator` is the function labeled "custom".
      i Truncate inside the function, or diagnose a built-in estimator to sweep the bound.

---

    Code
      check_eta_bias(data, a, y, c(x1, x2), estimator = estimate, truncation_grid = c(
        0, 0.05), n_boot = 10)
    Condition
      Error in `check_eta_bias()`:
      ! `truncation_grid` bounds a fitted nuisance model, which an estimator of your own never receives.
      x `estimator` is the function labeled "custom".
      i Truncate inside the function, or diagnose a built-in estimator to sweep the bound.

# the estimator return messages are stable

    Code
      check_eta_bias(data, a, y, c(x1, x2), estimator = function(.data) c(0.25, 0.5),
      n_boot = 10)
    Condition
      Error in `check_eta_bias()`:
      ! The `estimator` function must return 1 numeric value, one per estimand term.
      x It returned 2 values of type <numeric>.
      i The estimand term is "1 - 0".

---

    Code
      check_eta_bias(data, a, y, c(x1, x2), estimator = function(.data) c(zzz = 0.25),
      n_boot = 10)
    Condition
      Error in `check_eta_bias()`:
      ! The `estimator` function's return must be named by the estimand term or not named at all.
      x It returned "zzz".
      i The estimand term is "1 - 0".
      i An unnamed return is read in term order.

---

    Code
      check_eta_bias(data, a, y, c(x1, x2), estimator = function(.data) stop(
        "no estimate here"), n_boot = 5)
    Condition
      Error in `check_eta_bias()`:
      ! Every bootstrap draw was dropped, leaving no estimate to compare against the truth.
      i A draw is dropped when its estimate is not finite, or when an estimator supplied as a function raises an error on it.
      i A non-finite estimate arises when a bootstrap exposure leaves a level empty, when a fitted probability is exactly 0 or 1, or when a weight is degenerate.

# the custom estimator print method is stable

    Code
      print(res)
    Output
      
      -- ETA bias --------------------------------------------------------------------
      Exposure: "a" (binary)
      Observations: 120
      Estimator: cbps
      Bootstrap draws: 10
      Truth: 0.581
      ETA.Bias: -0.331 (MC SE 0)

