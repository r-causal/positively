# the point print method is stable

    Code
      print(res)
    Output
      
      -- port_result -----------------------------------------------------------------
      Exposure: "exposure" (binary)
      Observations: 1000
      Reading rule: alpha = 0.05, gamma = 2
      Prevalence threshold beta: 0.05
      Subgroups: 2 reported, 1 flagged

# the print method reports an unresolved beta

    Code
      print(res)
    Output
      
      -- port_result -----------------------------------------------------------------
      Exposure: "exposure" (binary)
      Observations: 1
      Reading rule: alpha = 0.05, gamma = 2
      Prevalence threshold beta: not resolved
      Subgroups: 0 reported, 0 flagged

# check_port() argument validation messages are stable

    Code
      check_port(1:10, exposure, x1)
    Condition
      Error in `check_port()`:
      ! `.data` must be a data frame, not a <integer>.

---

    Code
      check_port(data, exposure, tidyselect::starts_with("zzz"))
    Condition
      Error in `check_port()`:
      ! `.covariates` must select at least one column.

---

    Code
      check_port(data, c(exposure, x1), x2)
    Condition
      Error in `check_port()`:
      ! `.exposure` must select exactly one column, not 2.

---

    Code
      check_port(single_level, exposure, c(x1, x2))
    Condition
      Error in `check_port()`:
      ! `.exposure` must have at least two distinct values.
      x Found 1.

---

    Code
      check_port(data, exposure, c(x1, x2), alpha = 0)
    Condition
      Error in `check_port()`:
      ! `alpha` must be between 0 and 1, exclusive.
      x Found 0.

---

    Code
      check_port(data, exposure, c(x1, x2), alpha = 1.5)
    Condition
      Error in `check_port()`:
      ! `alpha` must be between 0 and 1, exclusive.
      x Found 1.5.

---

    Code
      check_port(data, exposure, c(x1, x2), alpha = c(0.05, 0.1))
    Condition
      Error in `check_port()`:
      ! `alpha` must be a single number.

---

    Code
      check_port(data, exposure, c(x1, x2), beta = 0)
    Condition
      Error in `check_port()`:
      ! `beta` must be between 0 and 1, exclusive.
      x Found 0.

---

    Code
      check_port(data, exposure, c(x1, x2), beta = 1.5)
    Condition
      Error in `check_port()`:
      ! `beta` must be between 0 and 1, exclusive.
      x Found 1.5.

---

    Code
      check_port(data, exposure, c(x1, x2), beta = 0.6)
    Condition
      Error in `check_port()`:
      ! `beta` must be below 0.5.
      i The reading rule flags prevalence below `beta` or above `1 - beta`, so a value of 0.5 or more flags every subgroup.
      x Found 0.6.

---

    Code
      check_port(data, exposure, c(x1, x2), beta = "banana")
    Condition
      Error in `check_port()`:
      ! `beta` must be a number in (0, 0.5) or "gruber".
      x Received "banana".

---

    Code
      check_port(data, exposure, c(x1, x2), gamma = 0)
    Condition
      Error in `check_port()`:
      ! `gamma` must be a single whole number of at least 1.
      x Found 0.

---

    Code
      check_port(data, exposure, c(x1, x2), gamma = 2.5)
    Condition
      Error in `check_port()`:
      ! `gamma` must be a single whole number of at least 1.
      x Found 2.5.

---

    Code
      check_port(data, exposure, c(x1, x2), n_bins = 0)
    Condition
      Error in `check_port()`:
      ! `n_bins` must be a single whole number of at least 2.
      x Found 0.

---

    Code
      check_port(data, exposure, c(x1, x2), breaks = "a")
    Condition
      Error in `check_port()`:
      ! `breaks` must be numeric or `NULL`, not a <character>.

