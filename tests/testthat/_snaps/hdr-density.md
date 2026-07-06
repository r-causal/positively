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

