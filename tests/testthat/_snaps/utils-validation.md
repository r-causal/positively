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

---

    Code
      validate_flag("yes", "flagged_only")
    Condition
      Error:
      ! `flagged_only` must be `TRUE` or `FALSE`.

# a selection error names the calling check in its message

    Code
      check_port(df, nope, x1)
    Condition
      Error in `check_port()`:
      ! `.exposure` must select a column that exists in `.data`.
      Caused by error in `eval_select_column()`:
      ! Can't select columns that don't exist.
      x Column `nope` doesn't exist.

