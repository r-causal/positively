# the numeric-grid undefined-cutoff message is stable

    Code
      check_hdr(df, exposure, l, values = c(0, 1), density_estimator = est)
    Condition
      Error in `numeric_hdr_thresholds()`:
      ! The numeric-grid HDR cutoff is undefined for 100 observations.
      x The estimator's conditional density is zero or non-finite at every grid value for those observations.
      i Supply a closed-form `hdr_threshold` to `new_hdr_density()`, or use a density whose mass falls inside the padded exposure range.

# new_hdr_density() validation messages are stable

    Code
      new_hdr_density(fit = "not a function", density = good_density)
    Condition
      Error in `new_hdr_density()`:
      ! `fit` must be a function.

---

    Code
      new_hdr_density(fit = good_fit, density = 1)
    Condition
      Error in `new_hdr_density()`:
      ! `density` must be a function.

---

    Code
      new_hdr_density(fit = good_fit, density = good_density, hdr_threshold = "not a function")
    Condition
      Error in `new_hdr_density()`:
      ! `hdr_threshold` must be a function or `NULL`.

---

    Code
      new_hdr_density(fit = good_fit, density = good_density, label = 1)
    Condition
      Error in `new_hdr_density()`:
      ! `label` must be a single string.

