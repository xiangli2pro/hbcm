#' Generate simulation data and parameters
#'
#' @description
#' `data_gen( )` generates matrix data (x) and parameters (alpha, Lambda, Sigma) for experimenting
#' with the HBCM model.
#'
#' @return A list of lists.
#' \item{x_list}{a list of data x.}
#' \item{alpha_list}{a list of alpha.}
#' \item{hlambda_list}{a list of Lambda.}
#' \item{hsigma_list}{a list of Sigma.}
#'
#' @param n an integer specifying the number of observations (rows) of data x.
#' @param p an integer specifying the number of variables (columns) of data x.
#' @param centers an integer specifying the number of clusters.
#' @param mu a vector of size `centers` specifying the mean vector of the 
#' multivariate normal distribution of alpha.
#' @param omega a matrix of size `centers x centers` specifying the covariance matrix of the
#' multivariate normal distribution of alpha.
#' @param labels a vector of size `p` specifying the cluster labels of the variables.
#' @param size an integer specifying the number of simulation data sets.
#' @param hparam_func a list of size two specifying the function of generating random numbers of
#' parameter Lambda and Sigma.
#'
#' @importFrom stats rnorm rchisq
#' @importFrom MASS mvrnorm
#'
#' @rdname simulation_data
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
#' omega <- diag(rep(1, centers))
#' for (i in 1:centers) {
#'  for (j in 1:centers) {
#'    if (i!=j){
#'      omega[i,j] = off_diag
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
#' data_list <- data_gen(n, p, centers, mu, omega, labels, size, hparam_func)
data_gen <- function(n, p, centers, mu, omega, labels, size, hparam_func) {
  sample_data <- lapply(1:size, function(k) sample_gen(n, p, mu, omega, labels, hparam_func))

  x_list <- lapply(sample_data, function(data) data[["x"]])
  alpha_list <- lapply(sample_data, function(data) data[["alpha"]])
  hlambda_list <- lapply(sample_data, function(data) data[["hlambda"]])
  hsigma_list <- lapply(sample_data, function(data) data[["hsigma"]])

  list(
    x_list = x_list, alpha_list = alpha_list,
    hlambda_list = hlambda_list, hsigma_list = hsigma_list
  )
}

# Helper function
sample_gen <- function(n, p, mu, omega, labels, hparam_func) {
  hlambda <- (hparam_func$lambda_func)(p)
  hsigma <- (hparam_func$sigma_func)(p)

  alpha <- MASS::mvrnorm(n, mu, omega)
  x <- matrix(rep(0, n * p), nrow = n)
  for (i in 1:n) {
    for (j in 1:p) {
      x[i, j] <- hlambda[j] * (alpha[i, labels[j]]) + hsigma[j] * rnorm(1, 0, 1)
    }
  }

  list(
    x = x, alpha = alpha,
    hlambda = hlambda, hsigma = hsigma
  )
}


