# positively (development version)

* First development version of positively, the positivity and extrapolation
  member of the r-causal family. positively reports diagnostics for positivity
  violations and extrapolation in causal inference. Each function returns a
  tidy result object and an optional plot, and it leaves the interpretation to
  you: the package does not grade severity or recommend remediation.

* `check_positivity()` is the general entry point. It resolves the exposure
  type once, binary, categorical, or continuous, runs the diagnostics
  appropriate to that type, and forwards the resolved type to every one of
  them, so the type is reported once rather than re-derived by each diagnostic
  in turn. Detection is a default and not an override: a supplied
  `exposure_type` is honored throughout, and the call aborts only when the
  exposure column cannot carry the type it was given. When the type was
  detected instead of declared, an error naming a diagnostic that does not
  apply to it also names the `exposure_type` to declare in order to run that
  diagnostic, provided the exposure column can carry that type.

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

* `exposure_type` is accepted and honored by every `check_*()` function that
  takes an exposure, including `check_port_seq()`, which previously had no such
  argument. A declared type is authoritative: it is checked against the
  exposure column rather than against the unique-value heuristic, so it is
  refused only when its math cannot run on the column as given. A continuous
  dose recorded at a handful of distinct levels is detected as categorical, and
  now runs under `exposure_type = "continuous"` in `check_hat_values()` and
  `check_hdr()` instead of being refused.

* `check_edp()` and `check_port()` validate a declared `exposure_type` against
  the exposure column. Previously `exposure_type = "binary"` on a three-level
  exposure was accepted, and it made `check_port()` search the top level alone,
  so the reported subgroups covered one level of three and the other two went
  unexamined. Such a declaration now aborts.

* `check_extrapolation()` and `check_eta_bias()` abort on a degenerate exposure
  that detection alone accepted. A constant factor exposure is detected as
  binary, and it made `check_extrapolation()` report infinite Gower distances
  behind a warning per observation. It made `check_eta_bias()` report an
  ETA.Bias of `NaN` under the default `estimator = "ipw"`, which dropped every
  bootstrap draw as non-finite, and a reported ETA.Bias of zero with no warning
  at all under `"gcomp"` and `"aipw"`.

* `check_port_seq()` gains `exposure_type`, resolved once from the exposure
  pooled across time points. Only a binary type carries the monotone-regime
  follower restriction, so declaring a non-binary type on a binary exposure
  keeps every subject in every risk set and stops the diagnostic from exposing
  the violations that pooling over the already-treated masks.
