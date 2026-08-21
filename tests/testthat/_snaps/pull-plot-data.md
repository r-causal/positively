# pull_plot_data() and autoplot() refuse the scatter view alike

    Code
      pull_plot_data(plain, type = "scatter")
    Condition
      Error in `pull_plot_data()`:
      ! The scatter view needs the estimator variant.
      i Rerun `check_edp()` with `variant = "estimator"`.

---

    Code
      ggplot2::autoplot(plain, type = "scatter")
    Condition
      Error in `autoplot.positively::edp_result`:
      ! The scatter view needs the estimator variant.
      i Rerun `check_edp()` with `variant = "estimator"`.

# pull_plot_data() and autoplot() refuse the cumulative view alike

    Code
      pull_plot_data(res, type = "cumulative")
    Condition
      Error in `pull_plot_data()`:
      ! `type = "cumulative"` requires a time-varying (matrix) input.
      i This result summarizes a point treatment with a single time point.

---

    Code
      ggplot2::autoplot(res, type = "cumulative")
    Condition
      Error in `autoplot.positively::density_ratios_result`:
      ! `type = "cumulative"` requires a time-varying (matrix) input.
      i This result summarizes a point treatment with a single time point.

# pull_plot_data() and autoplot() refuse the sweep view alike

    Code
      pull_plot_data(res, type = "sweep")
    Condition
      Error in `pull_plot_data()`:
      ! The sweep view needs a truncation sweep of more than one level.
      i Rerun `check_eta_bias()` with a `truncation_grid` of multiple lower bounds.

---

    Code
      ggplot2::autoplot(res, type = "sweep")
    Condition
      Error in `autoplot.positively::eta_bias_result`:
      ! The sweep view needs a truncation sweep of more than one level.
      i Rerun `check_eta_bias()` with a `truncation_grid` of multiple lower bounds.

# pull_plot_data() and autoplot() refuse a non-flag alike

    Code
      pull_plot_data(res, low_support_only = "yes")
    Condition
      Error in `pull_plot_data()`:
      ! `low_support_only` must be `TRUE` or `FALSE`.

---

    Code
      ggplot2::autoplot(res, low_support_only = "yes")
    Condition
      Error in `ggplot2::autoplot()`:
      ! `low_support_only` must be `TRUE` or `FALSE`.

# pull_plot_data() and autoplot() refuse the hull view alike

    Code
      pull_plot_data(res, type = "hull")
    Condition
      Error in `pull_plot_data()`:
      ! The convex-hull view needs the hull test to have run.
      i Rerun `check_extrapolation()` with `hull = TRUE`.

---

    Code
      ggplot2::autoplot(res, type = "hull")
    Condition
      Error in `autoplot.positively::extrapolation_result`:
      ! The convex-hull view needs the hull test to have run.
      i Rerun `check_extrapolation()` with `hull = TRUE`.

# a positivity check refuses to guess which diagnostic to pull

    Code
      pull_plot_data(check)
    Condition
      Error in `pull_plot_data()`:
      ! `diagnostic` must name a diagnostic in this container.
      i Available diagnostics are "port" and "edp".

---

    Code
      pull_plot_data(check, "bogus")
    Condition
      Error in `pull_plot_data()`:
      ! "bogus" is not a diagnostic in this container.
      i Available diagnostics are "port" and "edp".

