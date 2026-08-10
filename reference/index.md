# Package index

## Check positivity

Diagnose positivity across every applicable method and work with the
resulting container object.

- [`check_positivity()`](https://r-causal.github.io/positively/reference/check_positivity.md)
  : Diagnose positivity across the applicable methods
- [`positivity_check()`](https://r-causal.github.io/positively/reference/positivity_check.md)
  : A container for a set of positivity diagnostics
- [`sniff_violations()`](https://r-causal.github.io/positively/reference/sniff_violations.md)
  : Report what a set of diagnostics found

## Diagnostics

Individual positivity and extrapolation diagnostics, including their
time-varying sequence variants.

- [`check_edp()`](https://r-causal.github.io/positively/reference/check_edp.md)
  : Diagnose positivity with effective data points
- [`check_port()`](https://r-causal.github.io/positively/reference/check_port.md)
  : Locate positivity-violating subgroups with the PoRT algorithm
- [`check_port_seq()`](https://r-causal.github.io/positively/reference/check_port_seq.md)
  : Locate sequential positivity violations with the sPoRT algorithm
- [`check_hat_values()`](https://r-causal.github.io/positively/reference/check_hat_values.md)
  : Diagnose positivity for continuous exposures with hat values
- [`check_extrapolation()`](https://r-causal.github.io/positively/reference/check_extrapolation.md)
  : Diagnose extrapolation for binary exposures with Gower distances and
  the convex hull
- [`check_eta_bias()`](https://r-causal.github.io/positively/reference/check_eta_bias.md)
  : Diagnose positivity bias with the parametric bootstrap (ETA.Bias)
- [`check_hdr()`](https://r-causal.github.io/positively/reference/check_hdr.md)
  : Diagnose positivity for continuous exposures with HDR non-overlap
- [`check_hdr_seq()`](https://r-causal.github.io/positively/reference/check_hdr_seq.md)
  : Diagnose sequential positivity for continuous exposures with HDR
  non-overlap
- [`check_density_ratios()`](https://r-causal.github.io/positively/reference/check_density_ratios.md)
  : Summarize density ratios for positivity diagnosis

## Plots

Autoplot methods that visualize each diagnostic result.

- [`autoplot.positivity_check`](https://r-causal.github.io/positively/reference/autoplot.positivity_check.md)
  [`plot.positivity_check`](https://r-causal.github.io/positively/reference/autoplot.positivity_check.md)
  : Plot the diagnostics a positivity check holds
- [`autoplot.edp_result`](https://r-causal.github.io/positively/reference/autoplot.edp_result.md)
  [`plot.edp_result`](https://r-causal.github.io/positively/reference/autoplot.edp_result.md)
  : Plot an effective-data-points diagnostic
- [`autoplot.port_result`](https://r-causal.github.io/positively/reference/autoplot.port_result.md)
  [`plot.port_result`](https://r-causal.github.io/positively/reference/autoplot.port_result.md)
  : Plot a PoRT subgroup diagnostic
- [`autoplot.hat_values_result`](https://r-causal.github.io/positively/reference/autoplot.hat_values_result.md)
  [`plot.hat_values_result`](https://r-causal.github.io/positively/reference/autoplot.hat_values_result.md)
  : Plot a hat-value diagnostic
- [`autoplot.extrapolation_result`](https://r-causal.github.io/positively/reference/autoplot.extrapolation_result.md)
  [`plot.extrapolation_result`](https://r-causal.github.io/positively/reference/autoplot.extrapolation_result.md)
  : Plot an extrapolation diagnostic
- [`autoplot.eta_bias_result`](https://r-causal.github.io/positively/reference/autoplot.eta_bias_result.md)
  [`plot.eta_bias_result`](https://r-causal.github.io/positively/reference/autoplot.eta_bias_result.md)
  : Plot an ETA.Bias diagnostic
- [`autoplot.hdr_result`](https://r-causal.github.io/positively/reference/autoplot.hdr_result.md)
  [`plot.hdr_result`](https://r-causal.github.io/positively/reference/autoplot.hdr_result.md)
  : Plot an HDR non-overlap diagnostic
- [`autoplot.density_ratios_result`](https://r-causal.github.io/positively/reference/autoplot.density_ratios_result.md)
  [`plot.density_ratios_result`](https://r-causal.github.io/positively/reference/autoplot.density_ratios_result.md)
  : Plot a density-ratio diagnostic

## Density estimators

Conditional-density estimators used by the highest density region
diagnostic.

- [`hdr_density_normal()`](https://r-causal.github.io/positively/reference/hdr_density_normal.md)
  : The default normal conditional-density estimator for HDR non-overlap
- [`new_hdr_density()`](https://r-causal.github.io/positively/reference/new_hdr_density.md)
  : Construct a conditional-density estimator for the HDR diagnostic

## Datasets

Example data with planted positivity violations.

- [`pos_violations`](https://r-causal.github.io/positively/reference/pos_violations.md)
  : Binary-exposure data with planted positivity violations
- [`pos_violations_long`](https://r-causal.github.io/positively/reference/pos_violations_long.md)
  : Longitudinal continuous-exposure data with a time-2 support gap

## Extending positively

Developer-facing classes for building new diagnostics.

- [`positivity_diagnostic`](https://r-causal.github.io/positively/reference/positivity_diagnostic.md)
  : The abstract parent class for positivity diagnostics
