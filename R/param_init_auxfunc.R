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
  S <- covx - diag(hsigma^2, nrow = p, ncol = p)
  
  for (i in 1:p) {
    for (j in 1:p) {
      S[i, j] <- S[i, j] / (hlambda[i] * hlambda[j])
    }
  }
  
  omega <- matrix(0, centers, centers)
  hlambda_mat <- hlambda %*% t(hlambda)
  for (k in 1:(centers - 1)) {
    for (l in (k + 1):centers)
    {
      omega[k, l] <- sum(S[labels == k, labels == l] * hlambda_mat[labels == k, labels == l]) / sum(hlambda_mat[labels == k, labels == l]^2)
    }
  }
  
  omega <- omega + t(omega)
  diag(omega) <- 1
  
  # return
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


