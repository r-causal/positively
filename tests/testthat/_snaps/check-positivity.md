# the container print method is stable

    Code
      print(res)
    Output
      
      -- Positivity check ------------------------------------------------------------
      
      -- port --
      
      -- port_result -----------------------------------------------------------------
      Exposure: "exposure" (binary)
      Observations: 200
      Reading rule: alpha = 0.05, gamma = 2
      Prevalence threshold beta: 0.05
      Subgroups: 32 reported, 1 flagged
      
      -- extrapolation --
      
      -- extrapolation_result --------------------------------------------------------
      Exposure: "exposure" (binary)
      Observations: 200
      Geometric variability: 0.108
      Nearby radius (1 x gv): 0.108
      Mean fraction nearby: 0.18
      Nearest opposite within one geometric variability: 197 of 200
      In opposite-group hull: 180 of 200

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
      check_positivity(sim_pos_categorical(150), exposure, c(x1, x2), diagnostics = c(
        "port", "extrapolation"), exposure_type = "binary")
    Condition
      Error in `check_positivity()`:
      ! `exposure_type` was set to "binary", but `.exposure` is detected as "categorical".
      x "extrapolation" cannot run on a detected "categorical" exposure.
      i Drop it from `diagnostics`, or call it directly.

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
      generics::tidy(res, diagnostic = "nonexistent")
    Condition
      Error in `generics::tidy()`:
      ! `diagnostic` must name a diagnostic in this container.
      i Available diagnostics are "port".

---

    Code
      res[["nonexistent"]]
    Condition
      Error in `res[["nonexistent"]]`:
      ! "nonexistent" is not a diagnostic in this container.
      i Available diagnostics are "port".

