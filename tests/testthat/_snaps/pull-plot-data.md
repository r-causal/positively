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

