
<!-- README.md is generated from README.Rmd. Please edit that file -->

# weakARMA

<!-- badges: start -->
<!-- badges: end -->

The goal of weakARMA is to allows the study of nonlinear time series
models through weak ARMA representations.

## Installation (Gitlab)

### Current released

You can install the released version of weakARMA from
[PLMlab](https://plmlab.math.cnrs.fr) with:

``` r
install.packages("remotes")
remotes::install_gitlab("jrolland/weakARMA", host="https://plmlab.math.cnrs.fr")
```

### Development version

You can install the currently developped version of weakARMA from
[PLMlab](https://plmlab.math.cnrs.fr) with:

``` r
install.packages("remotes")
remotes::install_git("https://plmlab.math.cnrs.fr/jrolland/weakARMA.git", ref="develop")
```

## Installation (CRAN)

Once accepted, you’ll be able to install the released version of
weakARMA from [CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("weakARMA")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(weakARMA)
## basic example code
```

<!--
What is special about using `README.Rmd` instead of just `README.md`? You can include R chunks like so:


```r
summary(cars)
#>      speed           dist       
#>  Min.   : 4.0   Min.   :  2.00  
#>  1st Qu.:12.0   1st Qu.: 26.00  
#>  Median :15.0   Median : 36.00  
#>  Mean   :15.4   Mean   : 42.98  
#>  3rd Qu.:19.0   3rd Qu.: 56.00  
#>  Max.   :25.0   Max.   :120.00
```

You'll still need to render `README.Rmd` regularly, to keep `README.md` up-to-date. `devtools::build_readme()` is handy for this. You could also use GitHub Actions to re-render `README.Rmd` every time you push. An example workflow can be found here: <https://github.com/r-lib/actions/tree/master/examples>.

You can also embed plots, for example:

<img src="man/figures/README-pressure-1.png" width="100%" />

In that case, don't forget to commit and push the resulting figure files, so they display on GitHub and CRAN.
-->
