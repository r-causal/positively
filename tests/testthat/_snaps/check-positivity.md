# duplicate args names are rejected

    Code
      check_positivity(data, exposure, c(x1, x2), diagnostics = "port", args = list(
        port = list(alpha = 0.2), port = list(alpha = 0.4)))
    Condition
      Error in `check_positivity()`:
      ! `args` must not repeat a diagnostic name.
      x Duplicated: "port".

# the container print method is stable

    Code
      print(res)
    Output
      
      -- Positivity check ------------------------------------------------------------
      Exposure: "exposure" (binary); 200 observations; covariates x1 and x2
      
      -- port ------------------------------------------------------------------------
      32 subgroups reported, 1 with low support
      Rule: prevalence outside [0.05, 0.95] among subgroups of at least 5% of the
      sample
      
      -- extrapolation ---------------------------------------------------------------
      Geometric variability 0.108
      197 of 200 have an opposite-exposure unit within one geometric variability; 180
      of 200 fall in the opposite hull
      
      i `sniff_violations()` for what was found, `$port` to extract a diagnostic,
      `summary()` for every statistic.

# the classed errors are stable

    Code
      check_positivity(data, exposure, c(x1, x2), diagnostics = "bogus")
    Condition
      Error in `check_positivity()`:
      ! `diagnostics` names an unrecognised diagnostic: "bogus".
      i Valid diagnostics are "edp", "port", "hat_values", "hdr", and "extrapolation".

---

    Code
      check_positivity(data, exposure, c(x1, x2), diagnostics = "eta_bias")
    Condition
      Error in `check_positivity()`:
      ! `check_positivity()` does not run "eta_bias".
      i `check_eta_bias()` needs an outcome and `check_density_ratios()` needs user-supplied ratios.
      i Call it directly instead.

---

    Code
      check_positivity(data, exposure, c(x1, x2), diagnostics = "density_ratios")
    Condition
      Error in `check_positivity()`:
      ! `check_positivity()` does not run "density_ratios".
      i `check_eta_bias()` needs an outcome and `check_density_ratios()` needs user-supplied ratios.
      i Call it directly instead.

---

    Code
      check_positivity(data, exposure, c(x1, x2), diagnostics = "hat_values")
    Condition
      Error in `check_positivity()`:
      ! "hat_values" does not apply to a "binary" exposure.
      i Valid diagnostics for a "binary" exposure are "edp", "port", and "extrapolation".
      i Set `exposure_type = "continuous"` to run "hat_values".

---

    Code
      check_positivity(continuous, exposure, x1, diagnostics = "extrapolation")
    Condition
      Error in `check_positivity()`:
      ! "extrapolation" does not apply to a "continuous" exposure.
      i Valid diagnostics for a "continuous" exposure are "edp", "port", "hat_values", and "hdr".

---

    Code
      check_positivity(data, exposure, c(x1, x2), diagnostics = "port", args = list(
        bogus = list(alpha = 0.2)))
    Condition
      Error in `check_positivity()`:
      ! `args` names a diagnostic that is not being run: "bogus".
      i The diagnostics being run are "port".

---

    Code
      check_positivity(data, exposure, c(x1, x2), diagnostics = c("port", "port"))
    Condition
      Error in `check_positivity()`:
      ! `diagnostics` must not repeat a diagnostic.
      x Duplicated: "port".

---

    Code
      check_positivity(data, exposure, c(x1, x2), diagnostics = "port", args = list(
        port = "oops"))
    Condition
      Error in `check_positivity()`:
      ! Each element of `args` must be a list of options.
      x "port" is not a list.

---

    Code
      check_positivity(dgp_categorical(150), exposure, c(x1, x2), diagnostics = c(
        "port", "extrapolation"), exposure_type = "binary")
    Condition
      Error in `check_positivity()`:
      ! `check_positivity()` needs exactly two distinct values in `.exposure` for a binary exposure.
      x `.exposure` has 3 distinct values.

---

    Code
      check_positivity(data, exposure, c(x1, x2), diagnostics = "port", args = list(
        port = list(alpha = 1.5)))
    Condition
      Error in `check_positivity()`:
      ! `check_port()` failed while composing diagnostics.
      Caused by error in `check_port()`:
      ! `alpha` must be between 0 and 1, exclusive.
      x Found 1.5.

---

    Code
      check_positivity(data, c(exposure, x1), x2)
    Condition
      Error in `check_positivity()`:
      ! `.exposure` must select exactly one column, not 2.

---

    Code
      check_positivity(data, exposure, tidyselect::everything())
    Condition
      Error in `check_positivity()`:
      ! `.covariates` must not include the exposure column "exposure".
      i Exclude it from the selection, for example `c(everything(), -exposure)`.

---

    Code
      res[["nonexistent"]]
    Condition
      Error in `res[["nonexistent"]]`:
      ! "nonexistent" is not a diagnostic in this container.
      i Available diagnostics are "port".

---

    Code
      res[[5]]
    Condition
      Error in `res[[5]]`:
      ! Index 5 is out of bounds.
      i The container holds 1 diagnostic.

---

    Code
      res[[c("port", "extrapolation")]]
    Condition
      Error in `res[[c("port", "extrapolation")]]`:
      ! `i` must select a single diagnostic.
      x You supplied 2 values.

---

    Code
      check_positivity(data, exposure, c(x1, x2), diagnostics = c("hat_values", "hdr"))
    Condition
      Error in `check_positivity()`:
      ! "hat_values" and "hdr" do not apply to a "binary" exposure.
      i Valid diagnostics for a "binary" exposure are "edp", "port", and "extrapolation".
      i Set `exposure_type = "continuous"` to run "hat_values" and "hdr".

---

    Code
      check_positivity(dgp_categorical(150), exposure, c(x1, x2), diagnostics = c(
        "hat_values", "extrapolation"))
    Condition
      Error in `check_positivity()`:
      ! "hat_values" and "extrapolation" do not apply to a "categorical" exposure.
      i Valid diagnostics for a "categorical" exposure are "edp" and "port".

---

    Code
      check_positivity(data, exposure, c(x1, x2), diagnostics = "port", args = list(
        port = list(exposure_type = "categorical")))
    Condition
      Error in `check_positivity()`:
      ! `args` must not set `exposure_type`.
      x "port" sets it.
      i `check_positivity()` forwards its own `exposure_type` to every diagnostic, so declare the type there instead.

# the missing-exposure message is stable

    Code
      check_positivity(data, exposure, c(x1, x2))
    Condition
      Error in `check_positivity()`:
      ! `.exposure` must not contain missing values.

# autoplot() rejects a name the container does not hold

    Code
      autoplot(container, "hdr")
    Condition
      Error in `autoplot()`:
      ! "hdr" is not a diagnostic in this container.
      i Available diagnostics are "edp", "port", and "extrapolation".

# autoplot() says what it needs when patchwork is absent

    Code
      autoplot(container)
    Condition
      Error in `autoplot()`:
      ! Plotting a whole positivity check requires the patchwork package.
      i Install patchwork, or plot one diagnostic: `autoplot(check, "edp")`.

# an empty container says there is nothing to plot

    Code
      autoplot(container)
    Condition
      Error in `autoplot()`:
      ! There are no diagnostics to plot.

---

    Code
      plot(container)
    Condition
      Error in `autoplot()`:
      ! There are no diagnostics to plot.

