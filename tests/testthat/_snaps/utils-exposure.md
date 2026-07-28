# detect_exposure_type() announcement is stable

    Code
      exposure_type <- detect_exposure_type(c(0, 1, 0, 1))
    Message
      i Treating `.exposure` as binary

# resolve_exposure_type() errors point at the calling function

    Code
      check_fake(c(0, 1, 0, 1))
    Condition
      Error in `check_fake()`:
      ! `check_fake()` supports continuous exposures only.
      i `.exposure` was detected as "binary".
      i If `.exposure` is continuous, set `exposure_type = "continuous"`.

---

    Code
      check_fake(c(0, 1, 0, 1), exposure_type = "binary")
    Condition
      Error in `check_fake()`:
      ! `exposure_type` must be one of "auto" or "continuous", not "binary".

---

    Code
      check_fake(factor(c("a", "b", "c")), exposure_type = "continuous")
    Condition
      Error in `check_fake()`:
      ! `check_fake()` needs a numeric `.exposure` for a continuous exposure.
      x `.exposure` is <factor>.

