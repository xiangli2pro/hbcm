#' Estimate the optimal posterior distribution of the data column labels.
#'
#' @description
#' `heterogbcm( )` gives the optimal posterior distribution of the labels, which can be used
#' to derive the optimal label assignment of the data columns.
#'
#' @return A list of values.
#' \item{sigma}{estimated optimal group-correlation matrix.}
#' \item{hlambda}{estimated optimal heterogeneous parameter Lambda.}
#' \item{hsigma}{estimated optimal heterogeneous parameter Sigma.}
#' \item{obj_logL_val}{vector of -logL from each iteration.}
#' \item{qc}{estimated optimal posterior distribution of the column labels.}
#'
#' @rdname cluster_mod
#'
#' @param x a numeric matrix data.
#' @param centers an integer specifying the number of clusters.
#' @param labels a vector specifying the cluster labels of the columns of x.
#' @param tol numerical tolerance of the iteration updates.
#' @param iter number of iterations.
#' @param iter_init number of iterations of parameters initial estimation, default is 3.
#' @param verbose if TRUE, print iteration information.
#' @param hlambda heterogeneous parameter vector Lambda.
#' @param hsigma heterogeneous parameter vector Sigma.
#' @param qalpha posterior distribution of parameter vector alpha.
#' @param ppi probability of multi-nulli distribution.
#' @param sigma group-correlation matrix.
#' @param qc posterior distribution of labels.
#'
#' @export
heterogbcm <- function(x, centers, tol, iter, iter_init = 3, labels, verbose = FALSE) {
  n <- nrow(x)
  p <- ncol(x)

  # initial values of hlambda and hsigma
  init_hparameters <- init_hparam(x, centers, labels, tol, iter_init, verbose)
  hlambda <- init_hparameters$hlambda
  hsigma <- init_hparameters$hsigma

  # if centers == 1, sigma, ppi and qc0 are fixed
  if (centers == 1) {
    
    sigma <- 1
    ppi <- 1
    qc0 <- rep(1, p)
    
  } else {

    # initial estimate of group-correlation matrix sigma
    sigma <- init_sigma(x, centers, labels, hlambda, hsigma)
    # initial estimate of the probablity of the multi-nulli distribution
    ppi <- table(labels) / p
    # initial distribution of c based on ppi
    qc0 <- sapply(labels, function(grp) grp == c(1:centers)) * 1
  }

  # initial distribution of alpha based on sigma, qc0, hlambda, hsigma
  qalpha <- obj_qalpha(x, centers, sigma, qc0, hlambda, hsigma)

  # initial distribution of c based on sigma, qalpha
  qc <- obj_qc(x, centers, ppi, sigma, qalpha, hlambda, hsigma)

  # initial -logL
  obj_logL_val <- vector()
  obj_logL_val[1] <- obj_logL(
    x, centers, ppi, sigma,
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
      x, centers, ppi_new, sigma,
      qc, qalpha,
      hlambda, hsigma
    )


    # update sigma
    sigma_new <- obj_sigma(centers, qalpha)

    min_val <- verbose_print(
      verbose, "sigma", min_val,
      x, centers, ppi_new, sigma_new,
      qc, qalpha,
      hlambda, hsigma
    )


    # update hlambda
    hlambda_new <- obj_hlambda(x, centers, qc, qalpha)

    min_val <- verbose_print(
      verbose, "hlambda", min_val,
      x, centers, ppi_new, sigma_new,
      qc, qalpha,
      hlambda_new, hsigma
    )


    # update hsigma
    hsigma_new <- obj_hsigma(x, centers, qc, qalpha, hlambda_new)

    min_val <- verbose_print(
      verbose, "hsigma", min_val,
      x, centers, ppi_new, sigma_new,
      qc, qalpha,
      hlambda_new, hsigma_new
    )

    # update qalpha
    qalpha_new <- obj_qalpha(x, centers, sigma_new, qc, hlambda_new, hsigma_new)

    min_val <- verbose_print(
      verbose, "qalpha", min_val,
      x, centers, ppi_new, sigma_new,
      qc, qalpha_new,
      hlambda_new, hsigma_new
    )


    # update qc
    qc_new <- obj_qc(
      x, centers, ppi_new, sigma_new, qalpha_new,
      hlambda_new, hsigma_new
    )

    min_val <- verbose_print(
      verbose, "qc", min_val,
      x, centers, ppi_new, sigma_new,
      qc_new, qalpha_new,
      hlambda_new, hsigma_new
    )

    obj_logL_val[iiter + 1] <- obj_logL(
      x, centers, ppi_new, sigma_new,
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
                          x, centers, ppi, sigma, qc, qalpha, hlambda, hsigma) {
  if (verbose == TRUE) {
    obj_logL_val_new <- obj_logL(
      x, centers, ppi, sigma,
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
