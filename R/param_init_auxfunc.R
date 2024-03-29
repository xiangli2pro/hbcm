#' @rdname param_init
#' @description `init_omega( )` gives the initial estimation of group-correlation matrix omega.
#' @export
#' 
#' @param x a numeric matrix data.
#' @param centers an integer specifying the number of clusters.
#' @param labels a vector specifying the cluster labels of the columns of x.
#' @param hlambda heterogeneous parameter vector Lambda.
#' @param hsigma heterogeneous parameter vector Sigma.
init_omega <- function(x, centers, labels, hlambda, hsigma) {
  
  x <- as.matrix(x)
  n <- nrow(x)
  p <- ncol(x)
  
  covx <- t(x) %*% x / n
  hlambda_mat <- hlambda %*% t(hlambda)
  S <- covx / hlambda_mat
  
  omega <- matrix(0, centers, centers)
  for (k in 1:centers) {
    for (l in k:centers)
    {
      omega[k, l] <- mean(S[labels == k, labels == l])
    }
  }
  omega <- omega + t(omega) - diag(diag(omega))
  
  ## find the best estimate of omega
  
  # ## option 1
  # omega_v1 <- matrix(0, centers, centers)
  # S <- (covx - diag(hsigma^2, nrow = p, ncol = p)) / hlambda_mat
  # 
  # for (k in 1:centers) {
  #   for (l in k:centers)
  #   {
  #     # similar to linear regression
  #     # S(i,j) = lambda_i*lambda_j*W_ij
  #     # W_ij = (Xt*X)^(-1)*X*Y
  #     # W_ij = (hlambda*halmbda)^(-1)*halmbda*S
  #     omega_v1[k, l] <- sum(S[labels == k, labels == l] * hlambda_mat[labels == k, labels == l]) / sum(hlambda_mat[labels == k, labels == l]^2)
  #   }
  # }
  # omega_v1 <- omega_v1 + t(omega_v1) - diag(diag(omega_v1))
  # 
  # ## option 2
  # omega_v2 <- matrix(0, centers, centers)
  # S <- covx / hlambda_mat
  # 
  # for (k in 1:centers) {
  #   for (l in k:centers)
  #   {
  #     omega_v2[k, l] <- mean(S[labels == k, labels == l])
  #   }
  # }
  # omega_v2 <- omega_v2 + t(omega_v2) - diag(diag(omega_v2))
  # 
  # # return omega
  # if(det(omega_v1) > 0){
  #   omega <- omega_v1
  # } else {
  #   omega <- omega_v2
  # }
  
  omega
}


#' @rdname param_init
obj_init_qalpha <- function(x, hlambda, hsigma) {
  
  x <- as.matrix(x)
  n <- nrow(x)
  p <- ncol(x)

  # D_j_sum <- sum(hlambda^2 / hsigma^2)
  # Dbi_j_sum <- x %*% (hlambda / hsigma^2)

  alpha_cov <- 1 / (1 + sum(hlambda^2 / hsigma^2))
  alpha_mu <- alpha_cov * (x %*% (hlambda / hsigma^2))
  alpha_mu_inter <- alpha_cov + alpha_mu^2

  # return
  list(
    alpha_mu = alpha_mu,
    alpha_mu_inter = alpha_mu_inter,
    alpha_cov = alpha_cov
  )
}

#' @rdname param_init
obj_qalpha_logL <- function(x, qalpha, hlambda, hsigma) {
  
  x <- as.matrix(x)
  n <- nrow(x)
  p <- ncol(x)

  # given qalpha, entropy of alpha is constant
  entropy_alpha <- -n / 2 * log(2 * pi) - n / 2 * log(qalpha$alpha_cov) - n / 2

  logL <- -(n / 2) * log(2 * pi) - n / 2 * log(1) +
    -1 / 2 * 1 * sum(qalpha$alpha_mu_inter) +
    -n / 2 * sum(log(2 * pi * (hsigma^2))) +
    -1 / 2 * sum((x^2) %*% diag(1 / (hsigma^2), nrow = p, ncol = p)) +
    -1 / 2 * sum(qalpha$alpha_mu_inter) * sum(hlambda^2 / hsigma^2) +
    sum(diag(c(qalpha$alpha_mu)) %*% x %*% diag(hlambda / hsigma^2, nrow = p, ncol = p)) +
    -entropy_alpha

  # return -logL
  -logL
}

#' @rdname param_init
obj_init_hlambda <- function(x, qalpha) {
  
  x <- as.matrix(x)
  n <- nrow(x)
  p <- ncol(x)

  # return
  colSums(diag(c(qalpha$alpha_mu)) %*% x) / sum(qalpha$alpha_mu_inter)
}

#' @rdname param_init
obj_init_hsigma <- function(x, qalpha, hlambda) {
  
  x <- as.matrix(x)
  n <- nrow(x)
  p <- ncol(x)

  # return
  sqrt(1 / n * (colSums(x^2) +
    hlambda^2 * sum(qalpha$alpha_mu_inter) +
    -2 * colSums(diag(c(qalpha$alpha_mu)) %*% x %*% diag(hlambda, nrow = p, ncol = p))))
}


