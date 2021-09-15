
#' @rdname optim_vem
obj_logL <- function(X, centers, ppi, sigma, qc, qalpha, hlambda, hsigma) {
  n <- nrow(X)
  p <- ncol(X)

  if (any(hsigma == 0) | any(ppi == 0)) {
    if (any(hsigma == 0)) {
      stop("Error: hsigma has 0 value!")
    } else if (any(ppi == 0)) {
      stop("Error: ppi has 0 value!")
    }
  }

  if (centers == 1) {
    return(obj_qalpha_logL(X, qalpha, hlambda, hsigma))
  } else {
    qc_mat <- qc * log(qc)
    qc_mat[is.na(qc_mat)] <- 0
    entropy_c <- sum(qc_mat)

    entropy_alpha <- -centers / 2 * n * log(2 * pi) - n / 2 * log(det(qalpha$alpha_cov)) - centers / 2 * n

    if (det(sigma) <= 0) {
      det_sigma <- 1
    } else {
      det_sigma <- det(sigma)
    }

    # use vapply to specify the output type
    logL <- sum(t(qc) %*% log(ppi)) +
      -(n / 2) * centers * log(2 * pi) - n / 2 * log(det_sigma) +
      -1 / 2 * sum(sapply(qalpha$alpha_mu_inter, function(x) sum(diag(x %*% solve(sigma))))) +
      -n / 2 * sum(log(2 * pi * hsigma^2)) +
      -1 / 2 * sum((X^2) %*% (1 / (hsigma^2))) +
      -1 / 2 * sum((t(colSums(t(sapply(qalpha$alpha_mu_inter, function(x) diag(x))))) %*% qc) * (hlambda^2 / hsigma^2)) +
      sum((hlambda / hsigma^2) * colSums(X * (qalpha$alpha_mu %*% qc))) +
      -entropy_c - entropy_alpha
  }

  # return
  -logL
}