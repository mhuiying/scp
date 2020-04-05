
<!-- README.md is generated from README.Rmd. Please edit that file -->

# Spatial Conformal Prediction

<!-- badges: start -->

<!-- badges: end -->

The goal of scp is to provide spatial prediction intervals using Global
Spatial Conformal Prediction (GSCP) and Local Spatial Conformal
Prediction (LSCP).

## Installation

You can install the released version of scp from
[CRAN](https://CRAN.R-project.org) with:

``` r
<!-- install.packages("scp") -->
install_github("mhuiying/scp",auth_token = "53232d302e23533a1219e3fcdbde5d22c105dfc0")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(scp)

N = 21; n = N^2
S = seq(0,1,length=N)
s = expand.grid(S,S)
d = as.matrix(dist(s))

theta        = c(0,3,0.1,0.7)
names(theta) = c("Nugget","PartialSill","Range","Smoothness")
C = mat_cov(d,theta)
X = t(chol(C))%*%rnorm(n)
Y = X^3 + rnorm(n)

# Estimate spatial covariance parameters
bins     = seq(0.01,0.2,0.01)
thetaHat = get_theta(s,Y,dists=bins)

# spatial prediction
s0  = c(0.5, 0.5)
idx = which(s[,1]==s0[1] & s[,2]==s0[2])
PI  = conformal_pred(s0,s[-idx,],Y[-idx],thetaHat,eta=0.1)
cat(paste("True value: ", Y[idx], "\n"))
#> True value:  0.163383695736793
cat(paste("Prediction Interval: [ ", PI[1], ",", PI[2], "]"))
#> Prediction Interval: [  -10.4945921329898 , 29.9807915632017 ]
```

A visualization of the spatial process:

<img src="man/figures/README-visual-1.png" width="100%" />

User can also customize the nonconformity measure, (running time too
long, so I commented it out for
now)

``` r
# PI2 = scp(s0,s[-idx,],Y[-idx],pred_fun=krige_pred, dfun="abs_residual",precision=0.1)
# cat(paste("Prediction Interval: [ ", PI2[1], ",", PI2[2], "]"))
```

or the predictive function.

``` r
pred_fun = function(s0,s,Y,alpha) return(mean(Y))
PI3 = scp(s0,s[-idx,],Y[-idx],pred_fun=pred_fun, dfun="abs_residual",precision=0.1)
cat(paste("Prediction Interval: [ ", PI3[1], ",", PI3[2], "]"))
#> Prediction Interval: [  -33.8332433015048 , 30.0667566984952 ]
```
