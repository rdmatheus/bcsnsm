
# bcsnsm

<!-- badges: start -->
[![R-CMD-check](https://github.com/rdmatheus/bcsnsm/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/rdmatheus/bcsnsm/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

The `bcsnsm` package is dedicated to the estimation and analysis of goodness-of-fit
within the class of the multivariate BCS-NSM distributions and regression models. The
package includes implementations of some Box-Cox symmetric (BCS) and log-symmetric distributions for
marginal modeling. It provides the Gaussian, $t$, slash, and hyperbolic copulas
with some options of association matrices for dependence modeling. 

## Installation

The development version of the `bcnsm` package is hosted on
[GitHub](https://github.com/rdmatheus/bcnsm). The installation of the
package can be done using the `devtools` package. If it has not been
previously installed, it can be added by executing the following command
in the `R` environment: `install.packages("devtools")`.

After installing the `devtools` package, if one is using *Windows*, it
may be necessary to install the current version of the `RTools` program,
available for download [at
here](https://cran.r-project.org/bin/windows/Rtools). Finally, the
development version of the `bcnsm` can be installed by executing the
following command:

``` r
devtools::install_github("rdmatheus/bcsnsm")
```
