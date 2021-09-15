#' Estimate the optimal posterior distribution of the data column labels.
#'
#' @description
#' `hbcm( )` gives the optimal posterior distribution of the labels, which can be used
#' to derive the optimal label assignment of the data columns.
#'
#' @return A list of values.
#' \item{sigma}{estimated optimal group-correlation matrix.}
#' \item{hlambda}{estimated optimal heterogeneous parameter Lambda.}
#' \item{hsigma}{estimated optimal heterogeneous parameter Sigma.}
#' \item{obj_logL_val}{-logL after each iteration.}
#' \item{qc}{estimated optimal posterior distribution of the column labels.}
#'
#' @rdname optim_vem
#' 
#' @param X matrix data.
#' @param centers An integer specifying the number of clusters.
#' @param labels A vector specifying the cluster labels of the columns of X.
#' @param tol numerical tolerance of the iteration updates.
#' @param iter number of iterations.
#' @param verbose if TRUE, print parameters estimation on each iteration.
#' @param hlambda heterogeneous parameter Lambda.
#' @param hsigma heterogeneous parameter Sigma.
#' @param qalpha distribution of parameter alpha.
#' @param ppi probability of multi-nulli distribution.
#' @param sigma group-correlation matrix.
#' @param qc posterior distribution of labels.
#'
#' @export
hbcm <- function(X, centers, tol, iter, labels, verbose = FALSE) {
  n <- nrow(X)
  p <- ncol(X)

  # initial values of hlambda and hsigma
  init_hparameters <- init_hparam(X, centers, labels, tol, iter = 3, verbose)
  hlambda <- init_hparameters$hlambda
  hsigma <- init_hparameters$hsigma

  if (centers == 1) {
    sigma <- 1
    ppi <- 1
    qc0 <- rep(1, p)
  } else {

    # initial estimate of group-correlation matrix sigma
    sigma <- init_sigma(X, centers, labels, hlambda, hsigma)
    # initial estimate of the probablity of the multi-nulli distribution
    ppi <- table(labels) / p
    # initial distribution of c based on ppi
    qc0 <- sapply(labels, function(x) x == c(1:centers)) * 1
  }

  # initial distribution of alpha based on sigma, qc0, hlambda, hsigma
  qalpha <- obj_qalpha(X, centers, sigma, qc0, hlambda, hsigma)

  # initial distribution of c based on sigma, qalpha
  qc <- obj_qc(X, centers, ppi, sigma, qalpha, hlambda, hsigma)

  # initial -logL
  obj_logL_val <- vector()
  obj_logL_val[1] <- obj_logL(
    X, centers, ppi, sigma,
    qc, qalpha,
    hlambda, hsigma
  )

  min_val <- obj_logL_val[1]
  iiter <- 1

  while (iiter <= iter) {

    # update ppi
    ppi_new <- obj_ppi(centers, qc)

    min_val <- verbose_print(
      verbose, "ppi", min_val,
      X, centers, ppi_new, sigma,
      qc, qalpha,
      hlambda, hsigma
    )


    # update sigma
    sigma_new <- obj_sigma(centers, qalpha)

    min_val <- verbose_print(
      verbose, "sigma", min_val,
      X, centers, ppi_new, sigma_new,
      qc, qalpha,
      hlambda, hsigma
    )


    # update hlambda
    hlambda_new <- obj_hlambda(X, centers, qc, qalpha)

    min_val <- verbose_print(
      verbose, "hlambda", min_val,
      X, centers, ppi_new, sigma_new,
      qc, qalpha,
      hlambda_new, hsigma
    )


    # update hsigma
    hsigma_new <- obj_hsigma(X, centers, qc, qalpha, hlambda_new)

    min_val <- verbose_print(
      verbose, "hsigma", min_val,
      X, centers, ppi_new, sigma_new,
      qc, qalpha,
      hlambda_new, hsigma_new
    )

    # update qalpha
    qalpha_new <- obj_qalpha(X, centers, sigma_new, qc, hlambda_new, hsigma_new)

    min_val <- verbose_print(
      verbose, "alpha", min_val,
      X, centers, ppi_new, sigma_new,
      qc, qalpha_new,
      hlambda_new, hsigma_new
    )


    # update qc
    qc_new <- obj_qc(
      X, centers, ppi_new, sigma_new, qalpha_new,
      hlambda_new, hsigma_new
    )

    min_val <- verbose_print(
      verbose, "ppi", min_val,
      X, centers, ppi_new, sigma_new,
      qc_new, qalpha_new,
      hlambda_new, hsigma_new
    )

    obj_logL_val[iiter + 1] <- obj_logL(
      X, centers, ppi_new, sigma_new,
      qc_new, qalpha_new,
      hlambda_new, hsigma_new
    )

    if (abs(abs(obj_logL_val[iiter + 1] - obj_logL_val[iiter]) / obj_logL_val[iiter]) < tol) break

    iiter <- iiter + 1

    sigma <- sigma_new
    hlambda <- hlambda_new
    hsigma <- hsigma_new
    qc <- qc_new
    qalpha <- qalpha_new
  }

  list(
    sigma = sigma,
    hlambda = hlambda, hsigma = hsigma,
    obj_logL_val = obj_logL_val,
    qc = qc
  )
}

# Helper function -------------------------------------

## print iteration info
verbose_print <- function(verbose, param_name, min_val,
                          X, centers, ppi, sigma, qc, qalpha, hlambda, hsigma) {
  if (verbose == TRUE) {
    obj_logL_val_new <- obj_logL(
      X, centers, ppi, sigma,
      qc, qalpha,
      hlambda, hsigma
    )

    if (obj_logL_val_new <= min_val) {
      min_val <- obj_logL_val_new
      cat(paste0("Update ", param_name, " : ", obj_logL_val_new, " --"), "\n")
    } else {
      cat(paste0("Update ", param_name, " : ", obj_logL_val_new, " ++"), "\n")
    }
  }

  # return
  min_val
}
