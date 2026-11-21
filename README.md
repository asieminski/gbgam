# gbgam

Gradient boosting of generalised additive models in R: component-wise
boosting with penalised regression spline base learners built from 'mgcv'
smooth terms, including scale-invariant tensor product interactions. The
effective degrees of freedom are tracked during fitting and used for
adaptive smoothness relaxation and AIC-based stopping. Uncertainty
intervals come from a parametric bootstrap or Bayesian backfitting.

## Installation

```r
# install.packages("remotes")
remotes::install_github("asieminski/gbgam")
```

## Example

```r
library(gbgam)
dat <- gamSim(1, n = 400)
fit <- gbgam(y ~ t2(x0, k = 20) + t2(x1, k = 20) + t2(x2, k = 20) +
               t2(x3, k = 20), data = dat)
summary(fit)
plot(fit)
```

## Licence

GPL (>= 3)
