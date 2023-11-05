afttest: Model Diagnostics for AFT Models
================

[![CRAN_Status_Badge](https://www.r-pkg.org/badges/version/afttest)](https://CRAN.R-project.org/package=afttest)
[![Total_Downloads](https://cranlogs.r-pkg.org/badges/grand-total/afttest)](https://CRAN.R-project.org/package=afttest)
[![Downloads from the RStudio CRAN
mirror](https://cranlogs.r-pkg.org/badges/afttest)](https://CRAN.R-project.org/package=afttest)
[![BuildStatus](https://github.com/WoojungBae/afttest/workflows/R-CMD-check/badge.svg)](https://github.com/WoojungBae/afttest/actions)
[![arXiv](https://img.shields.io/badge/arXiv-10.48550/arXiv.2305.11445-<COLOR>.svg)](https://doi.org/10.48550/arXiv.2305.11445)

The R package **afttest** is intended to be a user-friendly *supplementary*
package to the base package **aftgee**.

## Features

In addition to the R interface, **afttest** provides a C++ header-only library 
integrated with **Rcpp**, which allows the construction of spline basis 
functions directly in C++ with the help of **Rcpp** and **RcppArmadillo**. 
Thus, it can also be treated as one of the **Rcpp\** packages. 

## Installation of CRAN Version

You can install the released version from
[CRAN](https://CRAN.R-project.org/package=afttest).

``` r
install.packages("afttest")
```

## Development

The latest version of the package is under development at
[GitHub](https://github.com/WoojungBae/afttest). If it is able to pass
the automated package checks, one may install it by

``` r
if (!require(devtools)) install.packages("devtools")
devtools::install_github("WoojungBae/afttest")
```

## Performance
The implementation of the main functions has been written in C++ with the 
help of the **Rcpp** and **RcppArmadillo** packages. The computational 
performance has thus been boosted.

## License
[GNU General Public License](https://www.gnu.org/licenses/) (≥ 3)
