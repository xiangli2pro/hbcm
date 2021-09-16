#' Initial estimation of heterogeneous parameters Lambda and Sigma
#'
#' @description
#' `init_hparam( )` gives the initial estimation of parameters Lambda and Sigma.
#'
#' @return A list of values.
#' \item{hlambda}{optimal initial estimation of Lambda.}
#' \item{hsigma}{optimal initial estimation of Sigma.}
#' 
#' @param X matrix data.
#' @param centers An integer specifying the number of clusters.
#' @param labels A vector specifying the cluster labels of the columns of X.
#' @param tol numerical tolerance of the iteration updates.
#' @param iter number of iterations.
#' @param verbose if TRUE, print iteration information.
#' @param hlambda heterogeneous parameter vector Lambda.
#' @param hsigma heterogeneous parameter vector Sigma.
#' @param qalpha distribution of parameter vector alpha.
#'
#' @importFrom RSpectra eigs_sym
#'
#' @export
#' @rdname param_init
init_hparam <- function(X, centers, labels,
                         tol, iter, verbose = FALSE) {

  # if all the columns belong to the same cluster
  if (centers == 1) {
    result <- init_hparam0(X, tol, iter, verbose = verbose)
    hsigma <- result$hsigma
    hlambda <- result$hlambda
  } else {
    n <- nrow(X)
    p <- ncol(X)

    # eig_max is the eigen-pair with largest eigenvalue from covariance matrix of X
    # take eigvector/sqrt(eigenvalue) as the start estimation of hlambda
    covX <- (t(X) %*% X) / n
    eig_max <- RSpectra::eigs_sym(covX - diag(diag(covX)), 1, which = "LM")
    hlambda_sign <- eig_max$vectors[, 1] * sqrt(eig_max$values[1])

    hlambda <- rep(0, p)
    hsigma <- rep(0, p)

    # update parts of hlamba that correspond to different clusters
    for (k in 1:centers)
    {
      result <- init_hparam0(X[, labels == k], tol, iter, verbose = verbose)

      hlambda_temp <- result$hlambda
      if (hparam_sign(hlambda_temp, hlambda_sign[labels == k]) <
        hparam_sign(-hlambda_temp, hlambda_sign[labels == k])) {
        hlambda_temp <- -hlambda_temp
      }

      hlambda[labels == k] <- hlambda_temp
      hsigma[labels == k] <- result$hsigma
    }
  }

  list(hsigma = hsigma, hlambda = hlambda)
}


#' @rdname param_init
init_hparam0 <- function(X, tol, iter, verbose = FALSE) {
  n <- nrow(X)
  p <- ncol(X)
  covX <- t(X) %*% X / n

  # eig_max is the eigen-pair with largest eigenvalue from covariance matrix of X
  # take eigvector/sqrt(eigenvalue) as the start estimation of hlambda
  eig_max <- RSpectra::eigs_sym(covX - diag(diag(covX)), 1, which = "LM")
  hlambda <- eig_max$vectors[, 1] * sqrt(eig_max$values[1])
  
  hsigma <- diag(covX) - diag(hlambda %*% t(hlambda))
  hsigma[hsigma < 0] <- 0
  hsigma <- sqrt(hsigma)

  # minimum of logL after each iteraction sub-step (qalpha, hsigma, hlambda)
  min_val <- Inf
  # logL after each iteration
  obj_val <- vector()
  obj_val[1] <- Inf

  iiter <- 1
  while (iiter <= iter) {

    # update alpha
    qalpha_new <- obj_init_qalpha(X, hlambda, hsigma) 
    
    min_val <- verbose_print_alpha(verbose, "alpha", min_val,
                                   X, qalpha_new, hlambda, hsigma)

    # update hlambda
    hlambda_new <- obj_init_hlambda(X, qalpha_new) 

    min_val <- verbose_print_alpha(verbose, "hlambda", min_val,
                                   X, qalpha_new, hlambda_new, hsigma)

    # update hsigma
    hsigma_new <- obj_init_hsigma(X, qalpha_new, hlambda_new) 

    min_val <- verbose_print_alpha(verbose, "hsigma", min_val,
                                   X, qalpha_new, hlambda_new, hsigma_new)

    obj_val[iiter + 1] <- obj_qalpha_logL(X, qalpha_new, hlambda_new, hsigma_new)
    if (abs(obj_val[iiter + 1] - obj_val[iiter]) < tol) break

    # Continue update
    iiter <- iiter + 1

    hlambda <- hlambda_new
    hsigma <- hsigma_new
  }

  return(list(hlambda = hlambda, hsigma = hsigma))
}


# Helper functions ---------------------------------------------------------------

## calculate the number of pairs that hlambda and hsigma have the same sign
hparam_sign <- function(x, y) {
  sum((x > 0) & (y > 0)) + sum((x < 0) & (y < 0))
}

## print iteration info
verbose_print_alpha <- function(verbose, param_name, min_val,
                                X, qalpha, hlambda, hsigma){
  
  if (verbose == TRUE) {
    
    obj_val_new <- obj_qalpha_logL(X, qalpha, hlambda, hsigma) 
    
    if (obj_val_new <= min_val) {
      
      min_val <- obj_val_new
      cat(paste0("Update ", param_name, " : ", obj_val_new, " --"), "\n")
    } else {
      cat(paste0("Update ", param_name, " : ", obj_val_new, " ++"), "\n")
    }
  }
  
  # return 
  min_val
}






