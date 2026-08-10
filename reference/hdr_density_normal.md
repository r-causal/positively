# The default normal conditional-density estimator for HDR non-overlap

`hdr_density_normal()` builds the default conditional-density estimator
used by
[`check_hdr()`](https://r-causal.github.io/positively/reference/check_hdr.md)
and
[`check_hdr_seq()`](https://r-causal.github.io/positively/reference/check_hdr_seq.md).
It fits a linear model of the exposure on the covariates, treats the
conditional density as Gaussian with the model's residual standard
deviation, and supplies a closed-form HDR threshold.

## Usage

``` r
hdr_density_normal()
```

## Value

An `hdr_density` estimator object.

## Details

The estimator fits `lm(exposure ~ covariates)`, recovering fitted means
\\\hat{\mu}(l)\\ and a homoskedastic residual standard deviation
\\\hat{\sigma}\\. The conditional density is \\\hat{f}(a \mid l) =
\mathrm{dnorm}(a; \hat{\mu}(l), \hat{\sigma})\\. Its HDR at mass `mass`
is a symmetric interval around \\\hat{\mu}(l)\\, so the density cutoff
is the closed form \$\$f\_\alpha = \mathrm{dnorm}(z) / \hat{\sigma},
\qquad z = \Phi^{-1}\\\left(\frac{1 + \mathrm{mass}}{2}\right).\$\$
Membership of a target `a` in the HDR reduces to the interval test
\\\|a - \hat{\mu}(l)\| \le z\\\hat{\sigma}\\, and the non-overlap ratio
at `a` is the fraction of fitted means more than \\z\\\hat{\sigma}\\
from `a`.

Because the working model is a single Gaussian, this estimator detects
mean-shift support gaps, where a stratum's supported dose moves away
from a target, but not multimodal gaps: a hole between two modes of the
true conditional density is filled by the fitted normal and reported as
supported. Supply a flexible estimator through
[`new_hdr_density()`](https://r-causal.github.io/positively/reference/new_hdr_density.md)
when multimodal structure is expected.

## References

Bao Y, Schomaker M (2025). Feasible Dose-Response Curves for Continuous
Treatments Under Positivity Violations.

## Examples

``` r
estimator <- hdr_density_normal()
estimator
#> <positively::hdr_density>
#>  @ fit          : function (formula, data)  
#>  @ density      : function (state, a, newdata)  
#>  @ hdr_threshold: function (state, newdata, mass)  
#>  @ label        : chr "normal"
```
