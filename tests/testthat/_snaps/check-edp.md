# the scatter view-gate abort names the autoplot method

    Code
      ggplot2::autoplot(plain, type = "scatter")
    Condition
      Error in `autoplot.positively::edp_result`:
      ! The scatter view needs the estimator variant.
      i Rerun `check_edp()` with `variant = "estimator"`.

# the density view says what it needs when ggridges is absent

    Code
      ggplot2::autoplot(res, type = "density")
    Condition
      Error in `autoplot.positively::edp_result`:
      ! The density view requires the ggridges package.
      i Install ggridges, or draw another view: `autoplot(x, type = "histogram")`.

# the argument validation messages are stable

    Code
      check_edp(1:10, exposure, x1)
    Condition
      Error in `check_edp()`:
      ! `.data` must be a data frame, not a <integer>.

---

    Code
      check_edp(data, exposure, tidyselect::starts_with("zzz"), exposure_type = "continuous")
    Condition
      Error in `check_edp()`:
      ! `.covariates` must select at least one column.

---

    Code
      check_edp(data, c(exposure, x1), x1, exposure_type = "continuous")
    Condition
      Error in `check_edp()`:
      ! `.exposure` must select exactly one column, not 2.

---

    Code
      check_edp(data, exposure, x1, categorical_similarity = 1.5, exposure_type = "continuous")
    Condition
      Error in `check_edp()`:
      ! `categorical_similarity` must be between 0 and 1.
      x Found 1.5.

---

    Code
      check_edp(data, exposure, x1, categorical_similarity = "x", exposure_type = "continuous")
    Condition
      Error in `check_edp()`:
      ! `categorical_similarity` must be a single number.

---

    Code
      check_edp(data, exposure, x1, bw_exposure = -1, exposure_type = "continuous")
    Condition
      Error in `check_edp()`:
      ! `bw_exposure` must not be negative.
      x Found -1.

---

    Code
      check_edp(data, exposure, x1, bw_exposure = c(1, 2), exposure_type = "continuous")
    Condition
      Error in `check_edp()`:
      ! `bw_exposure` must be a single finite number.

---

    Code
      check_edp(data, exposure, x1, bw_covariates = -1, exposure_type = "continuous")
    Condition
      Error in `check_edp()`:
      ! `bw_covariates` must not be negative.
      x Found -1.

---

    Code
      check_edp(data, exposure, x1, variant = "bogus", exposure_type = "continuous")
    Condition
      Error in `check_edp()`:
      ! `variant` must be one of "data" or "estimator", not "bogus".

---

    Code
      check_edp(data, exposure, x1, kernel = "bogus", exposure_type = "continuous")
    Condition
      Error in `check_edp()`:
      ! `kernel` must be one of "gaussian", not "bogus".

# the exposure-type error messages are stable

    Code
      check_edp(character_exposure, exposure, x1, exposure_type = "continuous")
    Condition
      Error in `check_edp()`:
      ! `check_edp()` needs a numeric `.exposure` for a continuous exposure.
      x `.exposure` is <character>.

---

    Code
      check_edp(sim_edp_categorical(300), exposure, z2, exposure_type = "binary")
    Condition
      Error in `check_edp()`:
      ! `check_edp()` needs exactly two distinct values in `.exposure` for a binary exposure.
      x `.exposure` has 3 distinct values.

---

    Code
      check_edp(character_exposure, exposure, x1, exposure_type = "bogus")
    Condition
      Error in `check_edp()`:
      ! `exposure_type` must be one of "auto", "binary", "categorical", or "continuous", not "bogus".

# the data-integrity error messages are stable

    Code
      check_edp(one_row, exposure, x1, exposure_type = "continuous")
    Condition
      Error in `check_edp()`:
      ! `.data` must have at least two rows, not 1.

---

    Code
      check_edp(na_exposure, exposure, x1, exposure_type = "continuous")
    Condition
      Error in `check_edp()`:
      ! `.exposure` must not contain missing values.
      x Missing values in "exposure".

---

    Code
      check_edp(na_covariate, exposure, x1, exposure_type = "continuous")
    Condition
      Error in `check_edp()`:
      ! `.covariates` must not contain missing values.
      x Missing values in "x1".

# the value and estimator-selection error messages are stable

    Code
      check_edp(data, exposure, x1, values = "a", exposure_type = "continuous")
    Condition
      Error in `check_edp()`:
      ! `values` must be numeric for a continuous exposure, not a <character>.

---

    Code
      check_edp(data, exposure, x1, values = c(0, NA), exposure_type = "continuous")
    Condition
      Error in `check_edp()`:
      ! `values` must not contain missing values.

---

    Code
      check_edp(data, exposure, x1, values = numeric(0), exposure_type = "continuous")
    Condition
      Error in `check_edp()`:
      ! `values` must contain at least one value.

---

    Code
      check_edp(data, exposure, x1, .outcome_covariates = tidyselect::starts_with(
        "zzz"), variant = "estimator", values = 0, exposure_type = "continuous")
    Condition
      Error in `check_edp()`:
      ! `.outcome_covariates` must select at least one column.

---

    Code
      check_edp(data, exposure, x1, .treatment_covariates = tidyselect::starts_with(
        "zzz"), variant = "estimator", values = 0, exposure_type = "continuous")
    Condition
      Error in `check_edp()`:
      ! `.treatment_covariates` must select at least one column.

# the function-intervention error messages are stable

    Code
      check_edp(data, exposure, x1, values = list(shift), exposure_type = "continuous")
    Condition
      Error in `check_edp()`:
      ! Every function in `values` must be named.
      i The name labels the intervention in the results and on a plot.
      x Element 1 is an unnamed function.

---

    Code
      check_edp(data, exposure, x1, values = stats::setNames(list(shift, shift), c(
        "down", "down")), exposure_type = "continuous")
    Condition
      Error in `check_edp()`:
      ! `values` must label every intervention distinctly.
      x Duplicated label: "down".

---

    Code
      check_edp(data, exposure, x1, values = stats::setNames(list(5, 1), c("", "5")),
      exposure_type = "continuous")
    Condition
      Error in `check_edp()`:
      ! `values` must label every intervention distinctly.
      x Duplicated label: "5".
      i An element without a name is labeled by its own value.

---

    Code
      check_edp(data, exposure, x1, values = list(bad = function(.x) rep("a", length(
        .x))), exposure_type = "continuous")
    Condition
      Error in `check_edp()`:
      ! The `values` function "bad" must return a numeric vector, not a <character>.

---

    Code
      check_edp(data, exposure, x1, values = list(bad = function(.x) .x[1:2]),
      exposure_type = "continuous")
    Condition
      Error in `check_edp()`:
      ! The `values` function "bad" must return one value per observation or a single value.
      x It returned 2 values for 60 observations.

---

    Code
      check_edp(data, exposure, x1, values = list(bad = function(.x) replace(.x, 1,
        NA_real_)), exposure_type = "continuous")
    Condition
      Error in `check_edp()`:
      ! The `values` function "bad" must not return missing values.

---

    Code
      check_edp(data, exposure, x1, values = list(bad = function(.x) replace(.x, 1,
        Inf)), exposure_type = "continuous")
    Condition
      Error in `check_edp()`:
      ! The `values` function "bad" must not return non-finite values.

---

    Code
      check_edp(data, exposure, x1, values = list(bad = NULL), exposure_type = "continuous")
    Condition
      Error in `check_edp()`:
      ! `values` must contain single values or functions of the exposure.
      x Element 1 is a <NULL> of length 0.

---

    Code
      check_edp(data, exposure, x1, values = list(), exposure_type = "continuous")
    Condition
      Error in `check_edp()`:
      ! `values` must contain at least one value.

---

    Code
      check_edp(binary, exposure, c(x1, x2), values = list(down = shift))
    Condition
      Error in `check_edp()`:
      ! A function in `values` needs a continuous exposure.
      x The exposure is binary.
      i A function maps each observed exposure to its intervened value.

# the unused-bandwidth warnings are stable

    Code
      res <- check_edp(binary, exposure, c(x1, x2), bw_exposure = 0.1)
    Condition
      Warning in `check_edp()`:
      `bw_exposure` was supplied but the exposure is not continuous.
      i The half-distance is ignored for a binary exposure.

---

    Code
      res <- check_edp(categorical, exposure, s, bw_covariates = 0.5)
    Condition
      Warning in `check_edp()`:
      `bw_covariates` was supplied but no covariate is numeric.
      i The half-distance is ignored when every covariate is categorical.

# the data-variant unused estimator-covariate warning is stable

    Code
      res <- check_edp(data, exposure, x1, .outcome_covariates = c(x1, x2), values = 0,
      exposure_type = "continuous")
    Condition
      Warning in `check_edp()`:
      `.outcome_covariates` applies only to the estimator variant.
      i The selection is ignored when `variant` is "data".

# the print method is stable

    Code
      print(continuous)
    Output
      
      -- Effective data points -------------------------------------------------------
      Exposure: "exposure" (continuous)
      Observations: 150
      Variant: data
      Intervention values: 2
      edp range: 2.483 to 38.179

---

    Code
      print(binary)
    Output
      
      -- Effective data points -------------------------------------------------------
      Exposure: "exposure" (binary)
      Observations: 200
      Variant: data
      Intervention values: 2
      edp range: 2.038 to 48.847

---

    Code
      print(estimator)
    Output
      
      -- Effective data points -------------------------------------------------------
      Exposure: "exposure" (continuous)
      Observations: 150
      Variant: estimator
      Intervention values: 1
      edp_outcome range: 4.23 to 38.179
      edp_treatment range: 8.995 to 98.272
      ideal_weight range: 2.127 to 2.795

# the print method is stable for a function intervention

    Code
      print(res)
    Output
      
      -- Effective data points -------------------------------------------------------
      Exposure: "exposure" (continuous)
      Observations: 150
      Variant: data
      Intervention values: 2
      edp range: 0.039 to 38.715

# the EDP label and headline are stable

    Code
      diagnostic_label(res)
    Output
      [1] "Effective data points"
    Code
      writeLines(diagnostic_headline(res))
    Output
      Data variant over 7 intervention values; edp 0.001 to 38.179

