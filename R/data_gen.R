#' Generate simulation data and parameters
#'
#' @description
#' `data_gen( )` generates simulation data X and parameters (alpha, Lambda, Sigma) for experimenting
#' with the HBCM model.
#'
#' @return A list of lists.
#' \item{X_list}{A list of data X.}
#' \item{alpha_list}{A list of alpha.}
#' \item{hlambda_list}{A list of Lambda.}
#' \item{hsigma_list}{A list of Sigma.}
#'
#' @param n An integer specifying the number of observations (rows) of data X.
#' @param p An integer specifying the number of variables (columns) of data X.
#' @param centers An integer specifying the number of clusters.
#' @param mu A vector of size `centers` specifying the mean vector of the 
#' multivariate normal distribution of alpha.
#' @param sigma A matrix of size `centers x centers` specifying the covariance matrix of the
#' multivariate normal distribution of alpha.
#' @param labels A vector of size `p` specifying the cluster labels of the variables.
#' @param size An integer specifying the number of simulation data sets.
#' @param hparam_func A list of size two specifying the function of generating random numbers of
#' parameter Lambda and Sigma.
#'
#' @importFrom stats rnorm rchisq
#' @importFrom MASS mvrnorm
#'
#' @export
#' @examples
#' n <- 500
#' p <- 500
#' centers <- 3
#' 
#' mu <- rep(0, centers)
#' 
#' off_diag <- 0.5
#' sigma <- diag(rep(1, centers))
#' for (i in 1:centers) {
#'  for (j in 1:centers) {
#'    if (i!=j){
#'      sigma[i,j] = off_diag
#'    } 
#'  }
#' }
#' 
#' ppi <- rep(1/centers, centers)
#' labels <- sample(c(1:centers), size = p, replace = TRUE, prob = ppi) 
#' 
#' size <- 5
#' hparam_func <- list(
#'   lambda_func = function(p) stats::rnorm(p, 0, 1),
#'   sigma_func = function(p) stats::rchisq(p, 2) + 1
#' )
#' data_list <- data_gen(n, p, centers, mu, sigma, labels, size, hparam_func)
data_gen <- function(n, p, centers, mu, sigma, labels, size, hparam_func) {
  sample_data <- lapply(1:size, function(k) sample_gen(n, p, mu, sigma, labels, hparam_func))

  X_list <- lapply(sample_data, function(x) x[["X"]])
  alpha_list <- lapply(sample_data, function(x) x[["alpha"]])
  hlambda_list <- lapply(sample_data, function(x) x[["hlambda"]])
  hsigma_list <- lapply(sample_data, function(x) x[["hsigma"]])

  list(
    X_list = X_list, alpha_list = alpha_list,
    hlambda_list = hlambda_list, hsigma_list = hsigma_list
  )
}

# Helper function
sample_gen <- function(n, p, mu, sigma, labels, hparam_func) {
  hlambda <- (hparam_func$lambda_func)(p)
  hsigma <- (hparam_func$sigma_func)(p)

  alpha <- MASS::mvrnorm(n, mu, sigma)
  X <- matrix(rep(0, n * p), nrow = n)
  for (i in 1:n) {
    for (j in 1:p) {
      X[i, j] <- hlambda[j] * (alpha[i, labels[j]]) + hsigma[j] * rnorm(1, 0, 1)
    }
  }

  list(
    X = X, alpha = alpha,
    hlambda = hlambda, hsigma = hsigma
  )
}
