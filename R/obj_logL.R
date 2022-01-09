#' log-likelihood function
#' @rdname cluster_mod
obj_logL <- function(x, centers, ppi, omega, qc, qalpha, hlambda, hsigma) {
  n <- nrow(x)
  p <- ncol(x)

  if (any(hsigma == 0) | any(ppi == 0)) {
    if (any(hsigma == 0)) {
      stop("Error: hsigma has 0 value!")
    } else if (any(ppi == 0)) {
      stop("Error: ppi has 0 value!")
    }
  }

  if (centers == 1) {
    return(obj_qalpha_logL(x, qalpha, hlambda, hsigma))
  } else {
    qc_mat <- qc * log(qc)
    qc_mat[is.na(qc_mat)] <- 0
    entropy_c <- sum(qc_mat)

    entropy_alpha <- -centers / 2 * n * log(2 * pi) - n / 2 * log(det(qalpha$alpha_cov)) - centers / 2 * n

    if (det(omega) <= 0) {
      det_omega <- 1
    } else {
      det_omega <- det(omega)
    }

    # use vapply to specify the output type
    logL <- sum(t(qc) %*% log(ppi)) +
      -(n / 2) * centers * log(2 * pi) - n / 2 * log(det_omega) +
      -1 / 2 * sum(sapply(qalpha$alpha_mu_inter, function(mu_inter) sum(diag(mu_inter %*% solve(omega))))) +
      -n / 2 * sum(log(2 * pi * hsigma^2)) +
      -1 / 2 * sum((x^2) %*% (1 / (hsigma^2))) +
      -1 / 2 * sum((t(colSums(t(sapply(qalpha$alpha_mu_inter, function(mu_inter) diag(mu_inter))))) %*% qc) * (hlambda^2 / hsigma^2)) +
      sum((hlambda / hsigma^2) * colSums(x * (qalpha$alpha_mu %*% qc))) +
      -entropy_c - entropy_alpha
  }

  # return
  -logL
}
