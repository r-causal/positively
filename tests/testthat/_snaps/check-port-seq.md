# the sequential print method is stable

    Code
      print(res)
    Output
      
      -- port_result -----------------------------------------------------------------
      Exposure: "a1", "a2", and "a3" (binary)
      Observations: 1000
      Time points: 3
      Reading rule: alpha = 0.05, gamma = 2
      Prevalence threshold beta: 3 per-time values (0.0229 to 0.0389)
      Subgroups: 262 reported, 1 with low support

# the sequential print is stable on the binary longitudinal fixture

    Code
      print(res)
    Output
      
      -- port_result -----------------------------------------------------------------
      Exposure: "a1", "a2", and "a3" (binary)
      Observations: 2000
      Time points: 3
      Reading rule: alpha = 0.05, gamma = 2
      Prevalence threshold beta: 0.05
      Subgroups: 45 reported, 2 with low support

# the sequential print reports censoring subgroups distinctly

    Code
      print(res)
    Output
      
      -- port_result -----------------------------------------------------------------
      Exposure: "a1", "a2", and "a3" (binary)
      Observations: 2000
      Time points: 3
      Reading rule: alpha = 0.05, gamma = 2
      Prevalence threshold beta: 0.05
      Exposure subgroups: 21 reported, 1 with low support
      Censoring subgroups: 32 reported, 1 with low support

# check_port_seq() censoring validation messages are stable

    Code
      check_port_seq(data, c(a1, a2, a3), list(l0, l1, l2), .censoring = c(c1, c2))
    Condition
      Error in `check_port_seq()`:
      ! `.censoring` must select one indicator per time point.
      i Found 2 indicators for 3 exposures.

---

    Code
      check_port_seq(data, c(a1, a2, a3), list(l0, l1, l2), .censoring = c(l0, l1, l2))
    Condition
      Error in `check_port_seq()`:
      ! `.censoring` must select binary 0/1 indicators.
      x "l0", "l1", and "l2" are not coded 0/1.
      i Code each censoring indicator as 1 when censored and 0 otherwise.

---

    Code
      check_port_seq(data, c(a1, a2, a3), list(l0, l1, l2), .censoring = c(c1, c2, c3),
      strategy = "pooled")
    Condition
      Error in `check_port_seq()`:
      ! `strategy` "pooled" is not yet implemented.
      i Use "stratified", the default.

# check_port_seq() argument validation messages are stable

    Code
      check_port_seq(1:10, c(a1, a2), list(c1))
    Condition
      Error in `check_port_seq()`:
      ! `.data` must be a data frame, not a <integer>.

---

    Code
      check_port_seq(data, tidyselect::starts_with("zzz"), list(c1))
    Condition
      Error in `check_port_seq()`:
      ! `.exposures` must select at least one column.

---

    Code
      check_port_seq(data, c(a1, a2, a3), list(c1, x))
    Condition
      Error in `check_port_seq()`:
      ! `.covariates` must supply one selection per time point.
      i Found 2 selections for 3 exposures.

---

    Code
      check_port_seq(data, c(a1, a2), list(tidyselect::starts_with("zzz")))
    Condition
      Error in `check_port_seq()`:
      ! Every time point must have at least one conditioning column.
      x The conditioning set at time 1 is empty.
      i Supply time-varying covariates in `.covariates` or baseline covariates in `.baseline`.

---

    Code
      check_port_seq(data, c(a1, a2, a3), list(c1), lag = -1)
    Condition
      Error in `check_port_seq()`:
      ! `lag` must be a single non-negative whole number or `Inf`.
      x Found -1.

---

    Code
      check_port_seq(data, c(a1, a2, a3), list(c1), lag = 1.5)
    Condition
      Error in `check_port_seq()`:
      ! `lag` must be a single non-negative whole number or `Inf`.
      x Found 1.5.

---

    Code
      check_port_seq(data, c(a1, a2, a3), list(c1), lag = "a")
    Condition
      Error in `check_port_seq()`:
      ! `lag` must be a single non-negative whole number or `Inf`.
      x Found "a".

---

    Code
      check_port_seq(data, c(a1, a2, a3), list(c1), alpha = 1.5)
    Condition
      Error in `check_port_seq()`:
      ! `alpha` must be between 0 and 1, exclusive.
      x Found 1.5.

---

    Code
      check_port_seq(data, c(a1, a2, a3), list(c1), beta = 1.5)
    Condition
      Error in `check_port_seq()`:
      ! `beta` must be between 0 and 1, exclusive.
      x Found 1.5.

---

    Code
      check_port_seq(data, c(a1, a2, a3), list(c1), beta = 0.6)
    Condition
      Error in `check_port_seq()`:
      ! `beta` must be below 0.5.
      i The reading rule flags prevalence below `beta` or above `1 - beta`, so a value of 0.5 or more flags every subgroup.
      x Found 0.6.

---

    Code
      check_port_seq(data, c(a1, a2, a3), list(c1), gamma = 0)
    Condition
      Error in `check_port_seq()`:
      ! `gamma` must be a single whole number of at least 1.
      x Found 0.

---

    Code
      check_port_seq(data, c(a1, a2, a3), list(c1), strategy = "pooled")
    Condition
      Error in `check_port_seq()`:
      ! `strategy` "pooled" is not yet implemented.
      i Use "stratified", the default.

# check_port_seq() missing-covariate message names the column and time

    Code
      check_port_seq(data, c(a1, a2, a3), list(c1))
    Condition
      Error in `purrr::map()`:
      i In index: 1.
      Caused by error in `port_seq_time()`:
      ! `.covariates` must not contain missing values in a risk set.
      x Missing values in "c1" at time point 1.

# check_port_seq() degenerate pooled-exposure messages are stable

    Code
      check_port_seq(constant, c(a1, a2), list(c1))
    Condition
      Error in `check_port_seq()`:
      ! `.exposures` must have at least two distinct values.
      x Found 1.

---

    Code
      check_port_seq(all_missing, c(a1, a2), list(c1))
    Condition
      Error in `check_port_seq()`:
      ! `.exposures` must have at least two distinct values.
      x Found 0.

