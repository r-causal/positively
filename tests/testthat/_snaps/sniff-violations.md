# an unrecognized scope is rejected

    Code
      sniff_violations(res, scope = "region")
    Condition
      Error in `sniff_violations()`:
      ! `scope` must be one of "subgroup", "overall", or "unit", not "region".

