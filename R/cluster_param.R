#' @rdname cluster_mod
obj_qc <- function(x, centers, ppi, omega, qalpha, hlambda, hsigma) {
  n <- nrow(x)
  p <- ncol(x)

  if (centers == 1) {
    qc <- rep(1, p)
  } else {
    n <- nrow(x)
    p <- ncol(x)

    qc0 <- rcpp_qc(
      n, p, centers, ppi,
      hsigma, hlambda,
      qalpha$alpha_mu, x, qalpha$alpha_mu_inter
    )
    qc0 <- exp(apply(qc0, 2, function(qc_col) qc_col - max(qc_col)))

    # normalize to be multi-nulli distribution
    qc <- t(t(qc0) / colSums(qc0))
  }

  # return
  qc
}

#' @rdname cluster_mod
obj_qalpha <- function(x, centers, omega, qc, hlambda, hsigma) {
  
  x <- as.matrix(x)
  n <- nrow(x)
  p <- ncol(x)

  if (centers == 1) {
    dsum_j <- qc %*% (hlambda^2 / hsigma^2)
    dbsumi_j <- qc %*% t(x %*% diag((hlambda / hsigma^2)))

    alpha_cov <- 1 / (1 + dsum_j)
    alpha_mu <- t(alpha_cov %*% dbsumi_j)
    alpha_mu_inter <- lapply(as.list(alpha_mu), function(mu) mu^2 + alpha_cov)
  } else {
    dsum_j <- diag(c(qc %*% (hlambda^2 / hsigma^2)))
    dbsumi_j <- qc %*% t(x %*% diag((hlambda / hsigma^2)))

    alpha_cov <- Matrix::solve(Matrix::solve(omega) + dsum_j) # covariance is the same for all i
    alpha_mu <- t(alpha_cov %*% dbsumi_j)
    alpha_mu_inter <- rcpp_qalpha_mu_inter(n, p, centers, alpha_mu, alpha_cov)
  }

  # return
  list(
    alpha_mu = alpha_mu,
    alpha_mu_inter = alpha_mu_inter,
    alpha_cov = alpha_cov
  )
}

#' @rdname cluster_mod
obj_ppi <- function(centers, qc) {
  if (centers == 1) {
    return(1)
  } else {
    return(colSums(t(qc)) / sum(qc))
  }
}

#' @rdname cluster_mod
obj_omega <- function(centers, qalpha) {
  n <- nrow(qalpha$alpha_mu)
  
  if (centers == 1) {
    return(1)
  } else {
    return(1 / n * base::Reduce("+", qalpha$alpha_mu_inter))
  }
}

#' @rdname cluster_mod
obj_hlambda <- function(x, centers, qc, qalpha) {
  n <- nrow(x)
  p <- ncol(x)

  if (centers == 1) {
    alpha_mu_inter <- unlist(qalpha$alpha_mu_inter)

    hlambda <- colSums(diag(c(qalpha$alpha_mu)) %*% x) / sum(alpha_mu_inter)
  } else {
    alpha_mu_inter <- t(sapply(qalpha$alpha_mu_inter, function(mu_inter) diag(mu_inter)))

    hlambda <- colSums(x * (qalpha$alpha_mu %*% qc)) / colSums(alpha_mu_inter %*% qc)
  }

  # return
  hlambda
}

#' @rdname cluster_mod
obj_hsigma <- function(x, centers, qc, qalpha, hlambda) {
  n <- nrow(x)
  p <- ncol(x)

  if (centers == 1) {
    alpha_mu_inter <- unlist(qalpha$alpha_mu_inter)
    hsigma <- sqrt(1 / n * (colSums(x^2) + hlambda^2 * rep(sum(alpha_mu_inter), p) +
      -2 * hlambda * colSums(diag(c(qalpha$alpha_mu)) %*% x)))
  } else {
    alpha_mu_inter <- t(sapply(qalpha$alpha_mu_inter, function(mu_inter) diag(mu_inter)))
    hsigma <- sqrt(1 / n * (colSums(x^2) + hlambda^2 * colSums(alpha_mu_inter %*% qc) +
      -2 * hlambda * colSums(x * (qalpha$alpha_mu %*% qc))))
  }

  # return
  hsigma
}
