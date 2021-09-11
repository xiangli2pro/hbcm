#' Title
#'
#' @param n 
#' @param p 
#' @param mu 
#' @param sigma 
#' @param labels 
#' @param hpara_func 
#' @return
#' @export
sample_gen <- function(n, p, mu, sigma, labels, hpara_func) {
  PLambda <- (hpara_func$lambda_func)(p)
  PSigma <- (hpara_func$sigma_func)(p)

  Alpha <- MASS::mvrnorm(n, mu, sigma)
  X <- matrix(rep(0, n * p), nrow = n)
  for (i in 1:n) {
    for (j in 1:p) {
      X[i, j] <- PLambda[j] * (Alpha[i, labels[j]]) + PSigma[j] * rnorm(1, 0, 1)
    }
  }
  list(Alpha = Alpha, X = X, PLambda = PLambda, PSigma = PSigma)
}

#' @export
data_gen <- function(n, p, centers, mu, sigma, labels, size, hpara_func) {
  sample_data <- lapply(1:size, function(k) sample_gen(n, p, mu, sigma, labels, hpara_func))

  Alpha_list <- lapply(sample_data, function(x) x[["Alpha"]])
  X_list <- lapply(sample_data, function(x) x[["X"]])
  PLambda_list <- lapply(sample_data, function(x) x[["PLambda"]])
  PSigma_list <- lapply(sample_data, function(x) x[["PSigma"]])

  list(Alpha_list = Alpha_list, X_list = X_list, PLambda_list = PLambda_list, PSigma_list = PSigma_list)
}







