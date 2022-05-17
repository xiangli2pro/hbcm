#' Initial estimation of heterogeneous parameters Lambda and Sigma
#'
#' @description
#' `init_hparam( )` gives the initial estimation of parameters Lambda and Sigma.
#'
#' @return A list of values.
#' \item{hlambda}{optimal initial estimation of Lambda.}
#' \item{hsigma}{optimal initial estimation of Sigma.}
#' 
#' @param x a numeric matrix data.
#' @param centers an integer specifying the number of clusters.
#' @param labels a vector specifying the cluster labels of the columns of x.
#' @param tol numerical tolerance of the iteration updates.
#' @param iter number of iterations.
#' @param verbose if TRUE, print iteration information.
#' @param hlambda heterogeneous parameter vector Lambda.
#' @param hsigma heterogeneous parameter vector Sigma.
#' @param qalpha distribution of parameter vector alpha.
#' @param vec1 a vector of hlambda estimation.
#' @param vec2 a vector of hlambda estimation.
#' @param param_name name of the parameter being updated.
#' @param min_val minimum objective-value being calculated.
#'
#' @importFrom RSpectra eigs_sym
#'
#' @export
#' @rdname param_init
init_hparam <- function(x, centers, labels,
                         tol, iter, verbose = FALSE) {

  # if all the columns belong to the same cluster
  if (centers == 1) {
    result <- init_hparam0(x, tol, iter, verbose = verbose)
    hsigma <- result$hsigma
    hlambda <- result$hlambda
  } else {
    n <- nrow(x)
    p <- ncol(x)

    # eig_max is the eigen-pair with largest eigenvalue from covariance matrix of x
    # take eigvector/sqrt(eigenvalue) as the start estimation of hlambda
    covx <- (t(x) %*% x) / n
    eig_max <- RSpectra::eigs_sym(covx - diag(diag(covx)), 1, which = "LM")
    hlambda_sign <- eig_max$vectors[, 1] * sqrt(eig_max$values[1])

    hlambda <- rep(0, p)
    hsigma <- rep(0, p)

    # update parts of hlambda that correspond to different clusters
    for (k in 1:centers)
    {
      result <- init_hparam0(x[, labels == k], tol, iter, verbose = verbose)

      hlambda_temp <- result$hlambda
      if (hparam_sign(hlambda_temp, hlambda_sign[labels == k]) <
        hparam_sign(-hlambda_temp, hlambda_sign[labels == k])) {
        hlambda_temp <- -hlambda_temp
      }

      hlambda[labels == k] <- hlambda_temp
      hsigma[labels == k] <- result$hsigma
    }
    
    ## update 05/17
    hsigma[hsigma == 0] <- 1
  }

  list(hsigma = hsigma, hlambda = hlambda)
}


#' @rdname param_init
init_hparam0 <- function(x, tol, iter, verbose = FALSE) {
  
  x <- as.matrix(x)
  
  n <- nrow(x)
  p <- ncol(x)
  covx <- t(x) %*% x / n
  
  # if (p <= 2){
  #   # ??????????????
  #   # Each cluster needs to have at least 3 objects
  #   # update 2022/02/22 
  #   # if there is only 1 column in the cluster, take the constant 
  #   hlambda = ifelse(p == 1, sd(x), apply(x, 2, sd))
  #   hsigma = rep(0.1, p)
  #   
  # } else {
  #   
  # }
  
  # eig_max is the eigen-pair with largest eigenvalue from covariance matrix of x
  # take eigvector/sqrt(eigenvalue) as the start estimation of hlambda
  eig_max <- RSpectra::eigs_sym(covx - diag(diag(covx)), 1, which = "LM")
  # update 2022/02/22
  # If the eigenvalue less than 0, take it to be 1
  hlambda <- ifelse(eig_max$values[1] >= 0, 
                    eig_max$vectors[, 1] * sqrt(eig_max$values[1]),
                    eig_max$vectors[, 1])
  
  hsigma <- diag(covx) - diag(hlambda %*% t(hlambda))
  # updated 2022/02/22, hsigma is in the denominator, cannot be 0
  # hsigma[hsigma < 0] <- 0.1 
  hsigma <- abs(hsigma)
  hsigma <- sqrt(hsigma)

  # minimum of logL after each iteraction sub-step (qalpha, hsigma, hlambda)
  min_val <- Inf
  # logL after each iteration
  obj_val <- vector()
  obj_val[1] <- Inf

  iiter <- 1
  while (iiter <= iter) {

    # update alpha
    qalpha_new <- obj_init_qalpha(x, hlambda, hsigma) 
    
    min_val <- verbose_print_alpha(verbose, "alpha", min_val,
                                   x, qalpha_new, hlambda, hsigma)

    # update hlambda
    hlambda_new <- obj_init_hlambda(x, qalpha_new) 

    min_val <- verbose_print_alpha(verbose, "hlambda", min_val,
                                   x, qalpha_new, hlambda_new, hsigma)

    # update hsigma
    hsigma_new <- obj_init_hsigma(x, qalpha_new, hlambda_new) 

    min_val <- verbose_print_alpha(verbose, "hsigma", min_val,
                                   x, qalpha_new, hlambda_new, hsigma_new)

    obj_val[iiter + 1] <- obj_qalpha_logL(x, qalpha_new, hlambda_new, hsigma_new)
    if (abs(obj_val[iiter + 1] - obj_val[iiter]) < tol) break

    # Continue update
    iiter <- iiter + 1

    hlambda <- hlambda_new
    hsigma <- hsigma_new
  }

  return(list(hlambda = hlambda, hsigma = hsigma))
}


# Helper functions ---------------------------------------------------------------

#' @rdname param_init
#' @description calculate the number of pairs that hlambda and hsigma have the same sign
hparam_sign <- function(vec1, vec2) {
  sum((vec1 > 0) & (vec2 > 0)) + sum((vec1 < 0) & (vec2 < 0))
}

#' @rdname param_init
#' @description print iteration info
verbose_print_alpha <- function(verbose, param_name, min_val,
                                x, qalpha, hlambda, hsigma){
  
  obj_val_new <- obj_qalpha_logL(x, qalpha, hlambda, hsigma) 
  
  if (obj_val_new <= min_val) {
    
    min_val <- obj_val_new
    
    if (verbose == TRUE) cat(paste0("Update Initial ", param_name, " : ", obj_val_new, " --"), "\n")
  } else {
    
    if (verbose == TRUE) cat(paste0("Update Initial ", param_name, " : ", obj_val_new, " ++"), "\n")
  }
  
  # return 
  min_val
}






