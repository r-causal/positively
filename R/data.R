#' Binary-exposure data with planted positivity violations
#'
#' A simulated point-treatment dataset that carries two known positivity
#' problems, so examples and tests can compare a diagnostic against ground
#' truth. It holds one structural violation, a subgroup that is never exposed,
#' and one practical near-violation, a region of covariate space where the
#' fitted propensity approaches zero or one while both exposure levels are still
#' observed overall.
#'
#' @format A tibble with 1000 rows and 4 columns:
#' \describe{
#'   \item{exposure}{Integer treatment indicator, 0 or 1.}
#'   \item{x1}{Numeric covariate, standard normal. Its steep effect on the
#'     propensity drives the practical near-violation.}
#'   \item{x2}{Numeric covariate, standard normal. It defines the structural
#'     subgroup together with `region`.}
#'   \item{region}{Factor covariate with levels `"a"` and `"b"`.}
#' }
#'
#' @details
#' The structural violation is exact: no subject with `region == "b"` and
#' `x2 > 1` is exposed, because the true propensity is set to zero in that
#' subgroup. Both exposure levels occur everywhere outside it, so the empty cell
#' is a property of the subgroup rather than of the sample size.
#'
#' The practical near-violation comes from the strong dependence of exposure on
#' `x1`. Fitting the propensity model
#' `glm(exposure ~ x1 + x2 + region, family = binomial())` recovers fitted
#' values that fall below 0.01 in the lower tail of `x1` and rise above 0.99 in
#' the upper tail, yet both exposure levels remain observed. This is the
#' signature of a practical violation: treatment is near-deterministic in the
#' tails without being structurally impossible.
#'
#' @source Simulated by `data-raw/make-datasets.R`.
"pos_violations"

#' Longitudinal continuous-exposure data with a time-2 support gap
#'
#' A simulated wide-format longitudinal dataset with one row per subject over
#' three time points. It carries a known support gap in the time-2 exposure, so
#' the sequential and continuous-exposure diagnostics can be exercised against
#' ground truth. It serves [check_hdr()], [check_hdr_seq()], and
#' [check_port_seq()].
#'
#' @format A tibble with 500 rows and 7 columns:
#' \describe{
#'   \item{id}{Integer subject identifier, 1 to 500.}
#'   \item{l0}{Numeric baseline covariate.}
#'   \item{a1}{Numeric time-1 exposure.}
#'   \item{l1}{Numeric time-1 covariate.}
#'   \item{a2}{Numeric time-2 exposure, carrying the support gap.}
#'   \item{l2}{Numeric time-2 covariate.}
#'   \item{a3}{Numeric time-3 exposure.}
#' }
#'
#' @details
#' The exposures `a1`, `a2`, and `a3` are continuous, each taking hundreds of
#' distinct values. The planted violation lives at time 2: `a2` is drawn on
#' `[0, 2]` or `[4, 6]`, so it never falls in the open interval `(2, 4)`. The
#' time-1 and time-3 exposures are unconstrained and place values throughout
#' that interval, which makes the gap specific to time 2 rather than a feature
#' of the exposure distribution as a whole. The covariates and exposures follow
#' a simple time-ordered process in which each variable depends on the one
#' before it.
#'
#' @source Simulated by `data-raw/make-datasets.R`.
"pos_violations_long"
