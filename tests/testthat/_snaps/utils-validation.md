# validate_* failures produce stable messages

    Code
      validate_data_frame(1:10)
    Condition
      Error:
      ! `.data` must be a data frame, not a <integer>.

---

    Code
      validate_column_selection(integer(0))
    Condition
      Error:
      ! `.covariates` must select at least one column.

---

    Code
      validate_prob(1.5)
    Condition
      Error:
      ! `probs` must be between 0 and 1.
      x Found 1.5.

---

    Code
      validate_prob("a")
    Condition
      Error:
      ! `probs` must be numeric, not a <character>.

