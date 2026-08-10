# Construct a conditional-density estimator for the HDR diagnostic

`new_hdr_density()` is the developer constructor for the pluggable
conditional-density contract used by
[`check_hdr()`](https://r-causal.github.io/positively/reference/check_hdr.md)
and
[`check_hdr_seq()`](https://r-causal.github.io/positively/reference/check_hdr_seq.md).
It packages three functions into an estimator object that the HDR
non-overlap ratio of Bao and Schomaker (2025) is computed against.
Supply an `hdr_threshold` to give a closed-form density cutoff, or leave
it `NULL` to request the numeric-grid fallback.

## Usage

``` r
new_hdr_density(fit, density, hdr_threshold = NULL, label = "custom")
```

## Arguments

- fit:

  A function of `(formula, data)` returning fitted state as a list.

- density:

  A function of `(state, a, newdata)` returning the conditional density
  \\\hat{f}(a \mid l)\\ for target `a` over the rows of `newdata`.

- hdr_threshold:

  An optional function of `(state, newdata, mass)` returning the per-row
  density cutoff. `NULL` (the default) requests the numeric-grid
  fallback. When it is `NULL`,
  [`check_hdr()`](https://r-causal.github.io/positively/reference/check_hdr.md)
  evaluates `density` on a fixed grid of 1024 exposure values spanning
  the observed exposure range padded by its own width on each side;
  estimators whose conditional standard deviation is very small relative
  to that grid step should supply a closed-form `hdr_threshold` to keep
  the cutoff well resolved.

- label:

  A single string naming the estimator, stored on the object and
  reported by the
  [`check_hdr()`](https://r-causal.github.io/positively/reference/check_hdr.md)
  result. Defaults to `"custom"`.

## Value

An `hdr_density` estimator object.

## Details

The highest-density region (HDR) at a covariate profile \\l\\ is the
smallest set of exposure values that captures probability mass `mass` of
the conditional density \\f(a \mid l)\\, namely \\A\_\alpha(l) = \\ a :
f(a \mid l) \ge f\_\alpha(l) \\\\. An estimator therefore needs three
pieces, which map to the three arguments:

- `fit(formula, data)` fits the conditional-density model and returns
  any fitted state as a list.
  [`check_hdr()`](https://r-causal.github.io/positively/reference/check_hdr.md)
  builds `formula` as `exposure ~ covariates` and passes the analysis
  `data`.

- `density(state, a, newdata)` returns \\\hat{f}(a \mid l)\\ for the
  scalar or vector target `a`, recycled against the rows of `newdata`.
  It receives the state from `fit()`.

- `hdr_threshold(state, newdata, mass)` returns the density cutoff
  \\f\_\alpha(l_j)\\ whose HDR captures `mass`, one value per row of
  `newdata`. When it is `NULL`,
  [`check_hdr()`](https://r-causal.github.io/positively/reference/check_hdr.md)
  instead evaluates `density` on a fine grid of exposure values per row
  and reads the cutoff off that grid, so the minimum contract is two
  functions.

Membership of a target `a` in the HDR at row \\j\\ is the test
\\\hat{f}(a \mid l_j) \ge f\_\alpha(l_j)\\, and the non-overlap ratio
\\\hat{\tau}(a)\\ is the fraction of rows whose HDR excludes `a`.

`fit()` and [`density()`](https://rdrr.io/r/stats/density.html) receive
the analysis data as a data frame holding the exposure and the covariate
columns as they were selected, so a factor or character covariate
arrives as a factor or character column. Encoding those columns is the
estimator's work.
[`hdr_density_normal()`](https://r-causal.github.io/positively/reference/hdr_density_normal.md)
leaves it to [`stats::lm()`](https://rdrr.io/r/stats/lm.html), which
expands them into contrast terms from the formula; an estimator that
needs a numeric design matrix should build one with
[`stats::model.matrix()`](https://rdrr.io/r/stats/model.matrix.html).

## References

Bao Y, Schomaker M (2025). Feasible Dose-Response Curves for Continuous
Treatments Under Positivity Violations.

## Examples

``` r
# A normal-equivalent estimator built through the public contract.
estimator <- new_hdr_density(
  fit = function(formula, data) {
    model <- lm(formula, data = data)
    list(model = model, sigma = sigma(model))
  },
  density = function(state, a, newdata) {
    mu <- predict(state$model, newdata = newdata)
    dnorm(a, mean = mu, sd = state$sigma)
  }
)
estimator
#> <positively::hdr_density>
#>  @ fit          : function (formula, data)  
#>  @ density      : function (state, a, newdata)  
#>  @ hdr_threshold: NULL
#>  @ label        : chr "custom"
```
