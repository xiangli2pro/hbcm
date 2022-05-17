
<!-- README.md is generated from README.Rmd. Please edit that file -->

# hbcm package

[![Project Status: Active - The project has reached a stable, usable
state and is being actively
developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![R build
status](https://github.com/rossellhayes/ipa/workflows/R-CMD-check/badge.svg)](https://github.com/rossellhayes/ipa/actions)

<!-- badges: start -->
<!-- badges: end -->

Community detection is a clustering method based on objects’ pairwise
relationships such that objects classified in the same group are more
densely connected than objects from different groups and correlations
within the same cluster are homogeneous. Most of the model-based
community detection methods such as the stochastic block model and its
variants are designed for networks with binary (yes/no) connectivity
values, which ignores the practical scenarios where the pairwise
relationships are continuous, reflecting different degrees of
connectivity. The heterogeneous block covariance model (HBCM) proposes a
novel clustering structure applicable on signed and continuous
connections, such as the covariances or correlations between objects.
Furthermore, it takes into account the heterogeneous property of each
object within a community. A novel variational EM algorithm is employed
to estimate the optimal group membership. HBCM has provable consistent
estimation of clustering memberships and its practical performance is
demonstrated by numerical simulations. The HBCM is illustrated on yeast
gene expression data and gene groups with correlated expression levels
responding to the same transcription factors are detected.

## Installation

To the date (2022-05-17) the package is still under development. You can
install the development version from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("xiangli2pro/hbcm")

# load package
library('hbcm')
```

## Examples

1.  Create a matrix data `x` of dimension `NxP=500x500`, with columns
    belonging to three (`K=3`) non-overlapping groups (groups labeled as
    1, 2, and 3). `x` values are determined by three things: parameter
    vector `alpha` of size `NxK` (follows multivariate normal
    distribution), heterogeneous parameters vector `hlambda` and
    `hsigma` of sizes `Px1` respectively.

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

# specify the (mu, omega) of the MVN distribution of alpha
mu <- rep(0, centers)

off_diag <- 0.5
omega <- diag(rep(1, centers))
for (i in 1:centers) {
  for (j in 1:centers) {
    if (i!=j){
      omega[i,j] = off_diag
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
data_list <- hbcm::data_gen(n, p, centers, mu, omega, labels, size, hparam_func)
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
#>   4.179   0.070   4.272

# use hbcm to perform clustering
system.time(
  res <- hbcm::heterogbcm(x, centers = centers, 
                  tol = 1e-3, iter = 100, iter_init = 3, 
                  labels = start_labels, 
                  verbose = FALSE)
)
#>    user  system elapsed 
#>   3.724   0.094   3.825
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

# result shows that hbcm model is better than spectral-clustering model in terms of rand index.
print(specc_eval)
#>    Rand adjRand 
#>   0.540   0.074
print(hbcm_eval)
#>    Rand adjRand 
#>   0.599   0.150
```
