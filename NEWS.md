# positively (development version)

* First development version of positively, the positivity and extrapolation
  member of the r-causal family. positively reports diagnostics for positivity
  violations and extrapolation in causal inference. Each function returns a
  tidy result object and an optional plot, and it leaves the interpretation to
  you: the package does not grade severity or recommend remediation.

* `check_positivity()` is the general entry point. It detects the exposure
  type, binary, categorical, or continuous, and runs the diagnostics
  appropriate to that type.

* Seven diagnostic families, each with its own `check_*()` function:
  effective data points (`check_edp()`), positivity regression trees
  (`check_port()`), hat values (`check_hat_values()`), highest density regions
  (`check_hdr()`), density ratios (`check_density_ratios()`), parametric
  bootstrap bias estimates (`check_eta_bias()`), and extrapolation
  (`check_extrapolation()`). The positivity regression tree and highest
  density region families add `check_port_seq()` and `check_hdr_seq()` for
  time-varying treatments.

* `tidy()`, `glance()`, and `autoplot()` methods for every result object, so
  each diagnostic yields a tidy data frame, a one-row summary, and a plot.

* Two example datasets with planted positivity problems: `pos_violations` for
  a single exposure and `pos_violations_long` for a longitudinal exposure.

* Five vignettes covering the general workflow, continuous exposures, density
  ratios, positivity regression trees, and estimator-focused diagnostics.

* `check_port_seq()` gains a `.censoring` argument, an ordered selection of
  per-time censoring indicators. Each indicator is analyzed as its own binary
  tree under a two-sided reading rule, and supplying it also restricts the
  exposure risk sets to the subjects uncensored so far.

* `check_eta_bias()` accepts factor and character covariates, which the fitted
  treatment and outcome models expand into indicator terms.

* The `positivity_check` container supports `x[["port"]]` to extract a child
  diagnostic and `names(x)` to list the diagnostics it holds.

* `check_hat_values()` records the number of fitted model parameters in a new
  `@p` property.

* The `check_density_ratios()` print method reports the proportion of density
  ratios exactly equal to zero.
