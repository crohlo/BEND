
<!-- README.md is generated from README.Rmd. Please edit that file -->

# BEND

<!-- badges: start -->
<!-- badges: end -->

The goal of BEND is to provide a set of models to estimate nonlinear
longitudinal data using Bayesian estimation methods. These models
include the:

1.  Bayesian Piecewise Random Effects Model (`Bayes_PREM()`) which
    estimates a piecewise random effects (mixture) model for a given
    number of latent classes and a latent number of possible
    changepoints in each class, and can incorporate class and outcome
    predictive covariates (see Lamm, 2022 and Lock et al., 2018 for more
    details).

2.  Bayesian Crossed Random Effects Model (`Bayes_CREM()`) which
    estimates a linear, quadratic, exponential, or piecewise crossed
    random effects models where individuals are changing groups over
    time (e.g., students and schools; see Rohloff et al., 2024 for more
    details).

3.  Bayesian Bivariate Piecewise Random Effects Model (`Bayes_BPREM()`)
    which estimates a bivariate piecewise random effects model to
    jointly model two related outcomes (e.g., reading and math
    achievement; see Peralta et al., 2022 for more details).

This package requires Just Another Gibbs Sampler (JAGS) to be installed
on your computer (<https://mcmc-jags.sourceforge.io/>), and depends on
the packages `rjags` and `label.switching`.

### References

Lamm, R. (2022). Incorporation of covariates in Bayesian piecewise
growth mixture models. <https://hdl.handle.net/11299/252533>

Lock, E. F., Kohli, N., & Bose, M. (2018). Detecting multiple random
changepoints in Bayesian piecewise growth mixture models. Psychometrika,
83(3), 733–750. <https://doi.org/10.1007/s11336-017-9594-5>

Peralta, Y., Kohli, N., Lock, E. F., & Davison, M. L. (2022). Bayesian
modeling of associations in bivariate piecewise linear mixed-effects
models. Psychological Methods, 27(1), 44–64.
<https://doi.org/10.1037/met0000358>

Rohloff, C. T., Kohli, N., & Lock, E. F. (2024). Identifiability and
estimability of Bayesian linear and nonlinear crossed random effects
models. British Journal of Mathematical and Statistical Psychology.
<https://doi.org/10.1111/bmsp.12334>

## Installation

You can install the development version of BEND from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
library(devtools)
install_github("crohlo/BEND")
```
