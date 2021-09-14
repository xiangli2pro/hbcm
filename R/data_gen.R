#' Generate simulation data and parameters
#'
#' @description 
#' `data_gen( )` generates simulation data and parameters for implementation HBCM model
#' 
#' @return A list of lists, each sub-list contains a set of data X, parameter Alpha, Lambda 
#' and Sigma.
#' 
#' @param n An integer specifying the number of observations (rows) of data X. 
#' @param p An integer specifying the number of variables (columns) of data X. 
#' @param centers An integer specifying the number of clusters.
#' @param mu A vector of size `p` specifying the mean of the multivariate normal distribution of Alpha.
#' @param sigma A matrix of size `p x p` specifying the variance matrix of the 
#' multivariate normal distribution of Alpha.
#' @param labels A vector of size `p` specifying the cluster labels of the variables.
#' @param size An integer specifying the number of simulated data sets.
#' @param hpara_func A list of size 2 specifying the function of generating random numbers for
#' parameter Lambda and Sigma.
#' 
#' @importFrom stats rnorm rchisq
#' @importFrom MASS mvrnorm
#' 
#' @export   
#' @examples
#' n <- 100
#' p <- 2
#' centers <- 2
#' mu <- rep(0, p)
#' sigma <- matrix(c(1,0,0,1), nrow=p)
#' labels <- c(1,2)
#' size <- 5
#' hpara_func <- list(
#'   lambda_func = function(p) stats::rnorm(p, 0, 1),
#'   sigma_func = function(p) stats::rchisq(p, 2) + 1
#' )
#' data_list <- data_gen(n, p, centers, mu, sigma, labels, size, hpara_func)
data_gen <- function(n, p, centers, mu, sigma, labels, size, hpara_func) {
  sample_data <- lapply(1:size, function(k) sample_gen(n, p, mu, sigma, labels, hpara_func))

  Alpha_list <- lapply(sample_data, function(x) x[["Alpha"]])
  X_list <- lapply(sample_data, function(x) x[["X"]])
  HLambda_list <- lapply(sample_data, function(x) x[["HLambda"]])
  HSigma_list <- lapply(sample_data, function(x) x[["HSigma"]])

  list(Alpha_list = Alpha_list, X_list = X_list, HLambda_list = HLambda_list, HSigma_list = HSigma_list)
}

# Helper function
sample_gen <- function(n, p, mu, sigma, labels, hpara_func) {
  HLambda <- (hpara_func$lambda_func)(p)
  HSigma <- (hpara_func$sigma_func)(p)
  
  Alpha <- MASS::mvrnorm(n, mu, sigma)
  X <- matrix(rep(0, n * p), nrow = n)
  for (i in 1:n) {
    for (j in 1:p) {
      X[i, j] <- HLambda[j] * (Alpha[i, labels[j]]) + HSigma[j] * rnorm(1, 0, 1)
    }
  }
  list(Alpha = Alpha, X = X, HLambda = HLambda, HSigma = HSigma)
}




