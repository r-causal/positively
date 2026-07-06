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
      Subgroups: 262 reported, 1 flagged

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
      Subgroups: 45 reported, 2 flagged

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

---

    Code
      check_port_seq(data, c(a1, a2, a3), list(c1), .censoring = c1)
    Condition
      Error in `check_port_seq()`:
      ! `.censoring` is not yet implemented.
      i Censoring-indicator trees are planned for a later release.

