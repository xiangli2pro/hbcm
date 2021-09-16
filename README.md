
<!-- README.md is generated from README.Rmd. Please edit that file -->

# hbcm package

<!-- badges: start -->
<!-- badges: end -->

hbcm package contains all the necessary functions to conduct cluster
analysis on network data by using the proposed heterogeneous block
covariance model (HBCM).

## Installation

To the date (09/12/2021) the package is still under development. You can
install the development version from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("xiangli2pro/hbcm")

# load package
library('hbcm')
```

## Examples

1.  Create a matrix data `x` of dimension `500x500`, with columns
    belonging to three non-overlapping groups (groups labeled as 1, 2,
    and 3). `x` values are determined by three things: parameter vector
    `alpha` (follows multivariate normal distribution), heterogeneous
    parameters vector `hlambda` and `hsigma`.

``` r
# check function arguments and return values documentation
?hbcm::data_gen
```

``` r
# x dimension
n <- 500
p <- 500

# cluster number
centers <- 3
# cluster labels follow a multinulli distribution with probability ppi
ppi <- rep(1/centers, centers)
# simulate a vector of labels from the multinulli distribution
labels <- sample(c(1:centers), size = p, replace = TRUE, prob = ppi) 

# specify the (mu, sigma) of the MVN distribution of alpha
mu <- rep(0, centers)

off_diag <- 0.5
sigma <- diag(rep(1, centers))
for (i in 1:centers) {
  for (j in 1:centers) {
    if (i!=j){
      sigma[i,j] = off_diag
    } 
  }
}

# set up the generating function of hlambda and hsigma
hparam_func <- list(
  lambda_func = function(p) stats::rnorm(p, 0, 1),
  sigma_func = function(p) stats::rchisq(p, 2) + 1
)

# set up the number of simulation data
size <- 1

# generate data
data_list <- hbcm::data_gen(n, p, centers, mu, sigma, labels, size, hparam_func)
x <- data_list$x_list[[1]]
```

2.  Use heterogeneous block covariance model (HBCM) to cluster the
    columns of data `x`. Need to provide a starting label guess and the
    number of clusters.

``` r
# check function arguments and return values documentation
?hbcm::heterogbcm()
```

``` r
# use spectral clustering to make a label guess
system.time(
start_labels <- kernlab::specc(abs(cor(x)), centers = 3)@.Data
)
#>    user  system elapsed 
#>   4.157   0.082   4.275

# use hbcm to perform clustering
system.time(
  res <- hbcm::heterogbcm(x, centers = 3, 
                  tol = 1e-3, iter = 100, iter_init = 3, 
                  labels = start_labels, 
                  verbose = FALSE)
)
#>    user  system elapsed 
#>   2.702   0.059   2.892
```

3.  Use metric [Rand-Index](https://en.wikipedia.org/wiki/Rand_index)
    and adjusted Rand-Index to compare the estimated label assignment
    with the true label assignment. The higher the value, the better the
    performance.

``` r
# check function arguments and return values documentation
?hbcm::matchLabel()
```

``` r
# hbcm::heterogbcm() returns the optimal posterior distribution of the latent label variables

# need to first convert distribution to label assignment
hbcm_labels <- apply(res$qc, 2, which.max)

# evaluate clustering performance
library('dplyr')
specc_eval <- hbcm::matchLabel(labels, start_labels) %>% 
  unlist() %>% round(3)
hbcm_eval <- hbcm::matchLabel(labels, hbcm_labels) %>% 
  unlist() %>% round(3)

print(specc_eval)
#>    Rand adjRand 
#>   0.526   0.075
print(hbcm_eval)
#>    Rand adjRand 
#>   0.788   0.527
```
