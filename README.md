
<!-- README.md is generated from README.Rmd. Please edit that file -->

# R Package for Spatial Conformal Prediction

<!-- badges: start -->

<!-- [![Travis build status](https://travis-ci.com/mhuiying/scp.svg?branch=master)](https://travis-ci.com/mhuiying/scp) -->

<!-- [![AppVeyor build status](https://ci.appveyor.com/api/projects/status/github/mhuiying/scp?branch=master&svg=true)](https://ci.appveyor.com/project/mhuiying/scp) -->

<!-- [![Codecov test coverage](https://codecov.io/gh/mhuiying/scp/branch/master/graph/badge.svg)](https://codecov.io/gh/mhuiying/scp?branch=master) -->

<!-- badges: end -->

The goal of “scp” is to provide valid model-free spatial prediction
intervals.

## Installation

The current development version can be installed from source using
devtools.

``` r
devtools::install_github("mhuiying/scp", build_vignettes = TRUE)
```

## Example

``` r
library(scp)

# an example sample data
data('sample_data')
s  = sample_data$s
Y  = sample_data$Y

# locations to predict
s0  = c(0.5,0.5)
s0s = rbind(c(0.4, 0.4), c(0.5,0.5), c(0.6, 0.6))

# default prediction interval
scp(s0=s0,s=s,Y=Y)
scp(s0=s0s,s=s,Y=Y)

# user define eta=0.1, where LSCP is considered
scp(s0=s0,s=s,Y=Y,eta=0.1)

# user define non-conformity measure
scp(s0=s0,s=s,Y=Y,dfun="std_residual2")

# user define prediction function
fun = function(s0,s,Y) return(mean(Y))
scp(s0=s0,s=s,Y=Y,pred_fun=fun)
```

Want more example, please check our `vignettes`.

``` r
browseVignettes('scp')
```

## References

Mao, Huiying, Ryan Martin, and Brian Reich. **Valid model-free spatial
prediction**, 2020. [\[arxiv\]](https://arxiv.org/pdf/2006.15640.pdf)
