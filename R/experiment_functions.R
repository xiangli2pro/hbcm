#' Estimate the optimal posterior distribution of the data column labels.
#'
#' @description
#' `heterogbcm_hparam( )` gives the optimal posterior distribution of the labels, use random initial values for hlambda and hsigma.
#'
#' @return A list of values.
#' \item{omega}{estimated optimal group-correlation matrix.}
#' \item{hlambda}{estimated optimal heterogeneous parameter Lambda.}
#' \item{hsigma}{estimated optimal heterogeneous parameter Sigma.}
#' \item{obj_logL_val}{vector of -logL from each iteration.}
#' \item{qc}{estimated optimal posterior distribution of the column labels.}
#' \item{cluster}{a vector of integers (from 1:k) indicating the cluster to which each column is allocated.}
#' 
#' @rdname experiment_functions
#'
#' @param x a numeric matrix data.
#' @param centers an integer specifying the number of clusters.
#' @param labels a vector specifying the cluster labels of the columns of x.
#' @param tol numerical tolerance of the iteration updates.
#' @param iter number of iterations.
#' @param verbose if TRUE, print iteration information.
#' @param init_hlambda initial values for parameter vector Lambda.
#' @param init_hsigma initial values for parameter vector Sigma.
#' @param iter_init iteration times of initial parameter estimation
#'
#' @export
heterogbcm_hparam <- function(x, centers, tol, iter, labels, verbose = FALSE,
                       init_hlambda, init_hsigma) {
  
  # make the data to be matrix
  x <- as.matrix(x)
  
  n <- nrow(x)
  p <- ncol(x)
  
  # initial values of hlambda and hsigma
  # init_hparameters <- init_hparam(x, centers, labels, tol, iter_init, verbose)
  hlambda <- init_hlambda
  hsigma <- init_hsigma
  
  # if centers == 1, omega, ppi and qc0 are fixed
  if (centers == 1) {
    
    omega <- 1
    ppi <- 1
    qc0 <- rep(1, p)
    
  } else {
    
    # initial estimate of group-correlation matrix omega
    omega <- init_omega(x, centers, labels, hlambda, hsigma)
    # initial estimate of the probablity of the multi-nulli distribution
    ppi <- table(labels) / p
    # initial distribution of c based on ppi
    qc0 <- sapply(labels, function(grp) grp == c(1:centers)) * 1
  }
  
  # initial distribution of alpha based on omega, qc0, hlambda, hsigma
  qalpha <- obj_qalpha(x, centers, omega, qc0, hlambda, hsigma)
  
  # initial distribution of c based on omega, qalpha
  qc <- obj_qc(x, centers, ppi, omega, qalpha, hlambda, hsigma)
  
  # initial -logL
  obj_logL_val <- vector()
  obj_logL_val[1] <- obj_logL(
    x, centers, ppi, omega,
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
      x, centers, ppi_new, omega,
      qc, qalpha,
      hlambda, hsigma
    )
    
    
    # update omega
    omega_new <- obj_omega(centers, qalpha)
    
    min_val <- verbose_print(
      verbose, "omega", min_val,
      x, centers, ppi_new, omega_new,
      qc, qalpha,
      hlambda, hsigma
    )
    
    
    # update hlambda
    hlambda_new <- obj_hlambda(x, centers, qc, qalpha)
    
    min_val <- verbose_print(
      verbose, "hlambda", min_val,
      x, centers, ppi_new, omega_new,
      qc, qalpha,
      hlambda_new, hsigma
    )
    
    
    # update hsigma
    hsigma_new <- obj_hsigma(x, centers, qc, qalpha, hlambda_new)
    
    min_val <- verbose_print(
      verbose, "hsigma", min_val,
      x, centers, ppi_new, omega_new,
      qc, qalpha,
      hlambda_new, hsigma_new
    )
    
    # update qalpha
    qalpha_new <- obj_qalpha(x, centers, omega_new, qc, hlambda_new, hsigma_new)
    
    min_val <- verbose_print(
      verbose, "qalpha", min_val,
      x, centers, ppi_new, omega_new,
      qc, qalpha_new,
      hlambda_new, hsigma_new
    )
    
    
    # update qc
    qc_new <- obj_qc(
      x, centers, ppi_new, omega_new, qalpha_new,
      hlambda_new, hsigma_new
    )
    
    min_val <- verbose_print(
      verbose, "qc", min_val,
      x, centers, ppi_new, omega_new,
      qc_new, qalpha_new,
      hlambda_new, hsigma_new
    )
    
    obj_logL_val[iiter + 1] <- obj_logL(
      x, centers, ppi_new, omega_new,
      qc_new, qalpha_new,
      hlambda_new, hsigma_new
    )
    
    if (abs(abs(obj_logL_val[iiter + 1] - obj_logL_val[iiter]) / obj_logL_val[iiter]) < tol) break
    
    iiter <- iiter + 1
    
    omega <- omega_new
    hlambda <- hlambda_new
    hsigma <- hsigma_new
    qc <- qc_new
    qalpha <- qalpha_new
  }
  
  cluster <- apply(qc, 2, which.max)
  
  list(
    omega = omega,
    hlambda = hlambda, hsigma = hsigma,
    obj_logL_val = -obj_logL_val, # minimize -> maximize
    qc = qc,
    cluster = cluster
  )
}

#' @rdname experiment_functions
#' @export
#' @description
#' `heterogbcm_constant_hlambda_hsigma( )` experiment function to keep hlambda and hsigma constant across simulation
heterogbcm_constant_hlambda_hsigma <- function(x, centers, tol, iter, iter_init = 3, labels, verbose = FALSE) {
  
  # make the data to be matrix
  x <- as.matrix(x)
  
  n <- nrow(x)
  p <- ncol(x)
  
  # initial values of hlambda and hsigma
  hlambda_const <- colSums(t(x) %*% x) / (n * p)
  hsigma_const <- rep(1, p)
  hlambda <- hlambda_const
  hsigma <- hsigma_const
  # official version
  # init_hparameters <- init_hparam(x, centers, labels, tol, iter_init, verbose)
  # hlambda <- init_hparameters$hlambda
  # hsigma <- init_hparameters$hsigma
  
  # if centers == 1, omega, ppi and qc0 are fixed
  if (centers == 1) {
    
    omega <- 1
    ppi <- 1
    qc0 <- rep(1, p)
    
  } else {
    
    # initial estimate of group-correlation matrix omega
    omega <- init_omega(x, centers, labels, hlambda, hsigma)
    # initial estimate of the probablity of the multi-nulli distribution
    ppi <- table(labels) / p
    # initial distribution of c based on ppi
    qc0 <- sapply(labels, function(grp) grp == c(1:centers)) * 1
  }
  
  # initial distribution of alpha based on omega, qc0, hlambda, hsigma
  qalpha <- obj_qalpha(x, centers, omega, qc0, hlambda, hsigma)
  
  # initial distribution of c based on omega, qalpha
  qc <- obj_qc(x, centers, ppi, omega, qalpha, hlambda, hsigma)
  
  # initial -logL
  obj_logL_val <- vector()
  obj_logL_val[1] <- obj_logL(
    x, centers, ppi, omega,
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
      x, centers, ppi_new, omega,
      qc, qalpha,
      hlambda, hsigma
    )
    
    
    # update omega
    omega_new <- obj_omega(centers, qalpha)
    
    min_val <- verbose_print(
      verbose, "omega", min_val,
      x, centers, ppi_new, omega_new,
      qc, qalpha,
      hlambda, hsigma
    )
    
    
    # update hlambda
    # official version
    # hlambda_new <- obj_hlambda(x, centers, qc, qalpha)
    hlambda_new <- hlambda_const
      
    min_val <- verbose_print(
      verbose, "hlambda", min_val,
      x, centers, ppi_new, omega_new,
      qc, qalpha,
      hlambda_new, hsigma
    )
    
    
    # update hsigma
    # official version
    # hsigma_new <- obj_hsigma(x, centers, qc, qalpha, hlambda_new)
    hsigma_new <- hsigma_const
    
    min_val <- verbose_print(
      verbose, "hsigma", min_val,
      x, centers, ppi_new, omega_new,
      qc, qalpha,
      hlambda_new, hsigma_new
    )
    
    # update qalpha
    qalpha_new <- obj_qalpha(x, centers, omega_new, qc, hlambda_new, hsigma_new)
    
    min_val <- verbose_print(
      verbose, "qalpha", min_val,
      x, centers, ppi_new, omega_new,
      qc, qalpha_new,
      hlambda_new, hsigma_new
    )
    
    
    # update qc
    qc_new <- obj_qc(
      x, centers, ppi_new, omega_new, qalpha_new,
      hlambda_new, hsigma_new
    )
    
    min_val <- verbose_print(
      verbose, "qc", min_val,
      x, centers, ppi_new, omega_new,
      qc_new, qalpha_new,
      hlambda_new, hsigma_new
    )
    
    obj_logL_val[iiter + 1] <- obj_logL(
      x, centers, ppi_new, omega_new,
      qc_new, qalpha_new,
      hlambda_new, hsigma_new
    )
    
    if (abs(abs(obj_logL_val[iiter + 1] - obj_logL_val[iiter]) / obj_logL_val[iiter]) < tol) break
    
    iiter <- iiter + 1
    
    omega <- omega_new
    hlambda <- hlambda_new
    hsigma <- hsigma_new
    qc <- qc_new
    qalpha <- qalpha_new
  }
  
  cluster <- apply(qc, 2, which.max)
  
  list(
    omega = omega,
    hlambda = hlambda, hsigma = hsigma,
    obj_logL_val = -obj_logL_val, # minimize -> maximize
    qalpha = qalpha,
    qc = qc,
    cluster = cluster
  )
}


#' @rdname experiment_functions
#' @export
#' @description
#' `heterogbcm_logL( )` experiment function to compare the initial logL and last logL.
heterogbcm_logL <- function(x, centers, tol, iter, iter_init = 3, labels, verbose = FALSE) {
  
  # make the data to be matrix
  x <- as.matrix(x)
  
  n <- nrow(x)
  p <- ncol(x)
  
  # initial values of hlambda and hsigma
  init_hparameters <- init_hparam(x, centers, labels, tol, iter_init, verbose)
  hlambda <- init_hparameters$hlambda
  hsigma <- init_hparameters$hsigma
  
  # if centers == 1, omega, ppi and qc0 are fixed
  if (centers == 1) {
    
    omega <- 1
    ppi <- 1
    qc0 <- rep(1, p)
    
  } else {
    
    # initial estimate of group-correlation matrix omega
    omega <- init_omega(x, centers, labels, hlambda, hsigma)
    # initial estimate of the probablity of the multi-nulli distribution
    ppi <- table(labels) / p
    # initial distribution of c based on ppi
    qc0 <- sapply(labels, function(grp) grp == c(1:centers)) * 1
  }
  
  # initial distribution of alpha based on omega, qc0, hlambda, hsigma
  qalpha <- obj_qalpha(x, centers, omega, qc0, hlambda, hsigma)
  
  # initial distribution of c based on omega, qalpha
  # update 0301
  # qc <- obj_qc(x, centers, ppi, omega, qalpha, hlambda, hsigma)
  qc <- qc0
    
  # initial -logL
  obj_logL_val <- vector()
  obj_logL_val[1] <- obj_logL(
    x, centers, ppi, omega,
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
      x, centers, ppi_new, omega,
      qc, qalpha,
      hlambda, hsigma
    )
    
    
    # update omega
    omega_new <- obj_omega(centers, qalpha)
    
    min_val <- verbose_print(
      verbose, "omega", min_val,
      x, centers, ppi_new, omega_new,
      qc, qalpha,
      hlambda, hsigma
    )
    
    
    # update hlambda
    hlambda_new <- obj_hlambda(x, centers, qc, qalpha)
    
    min_val <- verbose_print(
      verbose, "hlambda", min_val,
      x, centers, ppi_new, omega_new,
      qc, qalpha,
      hlambda_new, hsigma
    )
    
    
    # update hsigma
    hsigma_new <- obj_hsigma(x, centers, qc, qalpha, hlambda_new)
    
    min_val <- verbose_print(
      verbose, "hsigma", min_val,
      x, centers, ppi_new, omega_new,
      qc, qalpha,
      hlambda_new, hsigma_new
    )
    
    # update qalpha
    qalpha_new <- obj_qalpha(x, centers, omega_new, qc, hlambda_new, hsigma_new)
    
    min_val <- verbose_print(
      verbose, "qalpha", min_val,
      x, centers, ppi_new, omega_new,
      qc, qalpha_new,
      hlambda_new, hsigma_new
    )
    
    
    # update qc
    qc_new <- obj_qc(
      x, centers, ppi_new, omega_new, qalpha_new,
      hlambda_new, hsigma_new
    )
    
    min_val <- verbose_print(
      verbose, "qc", min_val,
      x, centers, ppi_new, omega_new,
      qc_new, qalpha_new,
      hlambda_new, hsigma_new
    )
    
    obj_logL_val[iiter + 1] <- obj_logL(
      x, centers, ppi_new, omega_new,
      qc_new, qalpha_new,
      hlambda_new, hsigma_new
    )
    
    if (abs(abs(obj_logL_val[iiter + 1] - obj_logL_val[iiter]) / obj_logL_val[iiter]) < tol) break
    
    iiter <- iiter + 1
    
    omega <- omega_new
    hlambda <- hlambda_new
    hsigma <- hsigma_new
    qc <- qc_new
    qalpha <- qalpha_new
  }
  
  cluster <- apply(qc, 2, which.max)
  
  list(
    omega = omega,
    hlambda = hlambda, hsigma = hsigma,
    obj_logL_val = -obj_logL_val, # minimize -> maximize
    qc = qc,
    cluster = cluster
  )
}


#' @rdname experiment_functions
#' @export
#' @description
#' `heterogbcm_iterStep( )` experiment function with different update algorithm
#' first update q2 and parameters Theta, then update q1 and parameters Theta.
heterogbcm_iterStep <- function(x, centers, tol, iter, iter_init = 3, 
                                labels, verbose = FALSE) {
  
  # make the data to be matrix
  x <- as.matrix(x)
  
  n <- nrow(x)
  p <- ncol(x)
  
  # initial values of hlambda and hsigma
  init_hparameters <- init_hparam(x, centers, labels, tol, iter_init, verbose)
  hlambda <- init_hparameters$hlambda
  hsigma <- init_hparameters$hsigma
  
  # if centers == 1, omega, ppi and qc0 are fixed
  if (centers == 1) {
    
    omega <- 1
    ppi <- 1
    qc0 <- rep(1, p)
    
  } else {
    
    # initial estimate of group-correlation matrix omega
    omega <- init_omega(x, centers, labels, hlambda, hsigma)
    # initial estimate of the probablity of the multi-nulli distribution
    ppi <- table(labels) / p
    # initial distribution of c based on ppi
    qc0 <- sapply(labels, function(grp) grp == c(1:centers)) * 1
  }
  
  # initial distribution of alpha based on omega, qc0, hlambda, hsigma
  qalpha <- obj_qalpha(x, centers, omega, qc0, hlambda, hsigma)
  
  # initial distribution of c based on omega, qalpha
  qc <- obj_qc(x, centers, ppi, omega, qalpha, hlambda, hsigma)
  
  # initial -logL
  obj_logL_val <- vector()
  obj_logL_val[1] <- obj_logL(
    x, centers, ppi, omega,
    qc, qalpha,
    hlambda, hsigma
  )
  
  min_val <- obj_logL_val[1]
  iiter <- 1
  
  while (iiter <= iter) {
    
    # ----------------------- update qalpha and theta
    # update qalpha
    qalpha_new <- obj_qalpha(x, centers, omega, qc, hlambda, hsigma)
    
    min_val <- verbose_print(
      verbose, "qalpha", min_val,
      x, centers, ppi, omega,
      qc, qalpha_new,
      hlambda, hsigma
    )
    
    # update ppi (first iteration)
    ppi_new_i1 <- obj_ppi(centers, qc)
    
    min_val <- verbose_print(
      verbose, "ppi", min_val,
      x, centers, ppi_new_i1, omega,
      qc, qalpha_new,
      hlambda, hsigma
    )
    
    # update omega (first iteration)
    omega_new_i1 <- obj_omega(centers, qalpha_new)
    
    min_val <- verbose_print(
      verbose, "omega", min_val,
      x, centers, ppi_new_i1, omega_new_i1,
      qc, qalpha_new,
      hlambda, hsigma
    )
    
    # update hlambda (first iteration)
    hlambda_new_i1 <- obj_hlambda(x, centers, qc, qalpha_new)
    
    min_val <- verbose_print(
      verbose, "hlambda", min_val,
      x, centers, ppi_new_i1, omega_new_i1,
      qc, qalpha_new,
      hlambda_new_i1, hsigma
    )
    
    # update hsigma (first iteration)
    hsigma_new_i1 <- obj_hsigma(x, centers, qc, qalpha_new, hlambda_new_i1)
    
    min_val <- verbose_print(
      verbose, "hsigma", min_val,
      x, centers, ppi_new_i1, omega_new_i1,
      qc, qalpha_new,
      hlambda_new_i1, hsigma_new_i1
    )
    
    # ----------------------- update qc and theta
    
    # update qc
    qc_new <- obj_qc(
      x, centers, ppi_new_i1, omega_new_i1, qalpha_new,
      hlambda_new_i1, hsigma_new_i1
    )
    
    min_val <- verbose_print(
      verbose, "qc", min_val,
      x, centers, ppi_new_i1, omega_new_i1,
      qc_new, qalpha_new,
      hlambda_new_i1, hsigma_new_i1
    )
    
    # update ppi (second iteration)
    ppi_new_i2 <- obj_ppi(centers, qc_new)
    
    min_val <- verbose_print(
      verbose, "ppi", min_val,
      x, centers, ppi_new_i2, omega_new_i1,
      qc_new, qalpha_new,
      hlambda_new_i1, hsigma_new_i1
    )
    
    # omega (no need to update omega, only depend on qalpha)
    omega_new_i2 <- omega_new_i1
    
    # update hlambda (first iteration)
    hlambda_new_i2 <- obj_hlambda(x, centers, qc_new, qalpha_new)
    
    min_val <- verbose_print(
      verbose, "hlambda", min_val,
      x, centers, ppi_new_i2, omega_new_i2,
      qc_new, qalpha_new,
      hlambda_new_i2, hsigma_new_i1
    )
    
    # update hsigma (first iteration)
    hsigma_new_i2 <- obj_hsigma(x, centers, qc_new, qalpha_new, hlambda_new_i2)
    
    min_val <- verbose_print(
      verbose, "hsigma", min_val,
      x, centers, ppi_new_i2, omega_new_i2,
      qc_new, qalpha_new,
      hlambda_new_i2, hsigma_new_i2
    )
    
    
    # ----------------------- check if convergence
    obj_logL_val[iiter + 1] <- obj_logL(
      x, centers, ppi_new_i2, omega_new_i2,
      qc_new, qalpha_new,
      hlambda_new_i2, hsigma_new_i2
    )
    
    if (abs(abs(obj_logL_val[iiter + 1] - obj_logL_val[iiter]) / obj_logL_val[iiter]) < tol) break
    
    iiter <- iiter + 1
    
    omega <- omega_new_i2
    hlambda <- hlambda_new_i2
    hsigma <- hsigma_new_i2
    qc <- qc_new
    qalpha <- qalpha_new
  }
  
  cluster <- apply(qc, 2, which.max)
  
  list(
    omega = omega,
    hlambda = hlambda, hsigma = hsigma,
    obj_logL_val = -obj_logL_val, # minimize -> maximize
    qc = qc,
    cluster = cluster
  )
}

#' @rdname experiment_functions
#' @export
#' @description
#' `heterogbcm_noInitLabel( )` experiment function with no initial guess of the labels.
heterogbcm_noInitLabel <- function(x, centers, tol, iter, iter_init = 3, 
                                   labels=NA, verbose = FALSE) {
  
  # make the data to be matrix
  x <- as.matrix(x)
  
  n <- nrow(x)
  p <- ncol(x)
  
  # assign random cluster probability
  qc0 <- t(MCMCpack::rdirichlet(p, rep(1, centers)))
  labels <- apply(qc0, 2, which.max)
  
  # initial values of hlambda and hsigma
  init_hparameters <- init_hparam(x, centers, labels, tol, iter_init, verbose)
  hlambda <- init_hparameters$hlambda
  hsigma <- init_hparameters$hsigma
  
  # if centers == 1, omega, ppi and qc0 are fixed
  if (centers == 1) {
    
    omega <- 1
    ppi <- 1
    qc0 <- rep(1, p)
    
  } else {
    
    # initial estimate of group-correlation matrix omega
    omega <- init_omega(x, centers, labels, hlambda, hsigma)
    # initial estimate of the probablity of the multi-nulli distribution
    ppi <- table(labels) / p
    # initial distribution of c based on ppi
    # qc0 <- sapply(labels, function(grp) grp == c(1:centers)) * 1
  }
  
  # initial distribution of alpha based on omega, qc0, hlambda, hsigma
  qalpha <- obj_qalpha(x, centers, omega, qc0, hlambda, hsigma)
  
  # initial distribution of c based on omega, qalpha
  qc <- qc0
  
  # initial -logL
  obj_logL_val <- vector()
  obj_logL_val[1] <- obj_logL(
    x, centers, ppi, omega,
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
      x, centers, ppi_new, omega,
      qc, qalpha,
      hlambda, hsigma
    )
    
    
    # update omega
    omega_new <- obj_omega(centers, qalpha)
    
    min_val <- verbose_print(
      verbose, "omega", min_val,
      x, centers, ppi_new, omega_new,
      qc, qalpha,
      hlambda, hsigma
    )
    
    
    # update hlambda
    hlambda_new <- obj_hlambda(x, centers, qc, qalpha)
    
    min_val <- verbose_print(
      verbose, "hlambda", min_val,
      x, centers, ppi_new, omega_new,
      qc, qalpha,
      hlambda_new, hsigma
    )
    
    
    # update hsigma
    hsigma_new <- obj_hsigma(x, centers, qc, qalpha, hlambda_new)
    
    min_val <- verbose_print(
      verbose, "hsigma", min_val,
      x, centers, ppi_new, omega_new,
      qc, qalpha,
      hlambda_new, hsigma_new
    )
    
    # update qalpha
    qalpha_new <- obj_qalpha(x, centers, omega_new, qc, hlambda_new, hsigma_new)
    
    min_val <- verbose_print(
      verbose, "qalpha", min_val,
      x, centers, ppi_new, omega_new,
      qc, qalpha_new,
      hlambda_new, hsigma_new
    )
    
    # update qc
    qc_new <- obj_qc(
      x, centers, ppi_new, omega_new, qalpha_new,
      hlambda_new, hsigma_new
    )
    
    min_val <- verbose_print(
      verbose, "qc", min_val,
      x, centers, ppi_new, omega_new,
      qc_new, qalpha_new,
      hlambda_new, hsigma_new
    )
    
    obj_logL_val[iiter + 1] <- obj_logL(
      x, centers, ppi_new, omega_new,
      qc_new, qalpha_new,
      hlambda_new, hsigma_new
    )
    
    if (abs(abs(obj_logL_val[iiter + 1] - obj_logL_val[iiter]) / obj_logL_val[iiter]) < tol) break
    
    iiter <- iiter + 1
    
    omega <- omega_new
    hlambda <- hlambda_new
    hsigma <- hsigma_new
    qc <- qc_new
    qalpha <- qalpha_new
  }
  
  cluster <- apply(qc, 2, which.max)
  
  list(
    omega = omega,
    hlambda = hlambda, hsigma = hsigma,
    obj_logL_val = -obj_logL_val, # minimize -> maximize
    qc = qc,
    cluster = cluster
  )
}



#' @rdname experiment_functions
#' @export
#' @description
#' `heterogbcm_qcDiscrete( )` experiment function to update qc with discrete values.
heterogbcm_qcDiscrete <- function(x, centers, tol, iter, iter_init = 3, labels, verbose = FALSE) {
  
  # make the data to be matrix
  x <- as.matrix(x)
  
  n <- nrow(x)
  p <- ncol(x)
  
  # initial values of hlambda and hsigma
  init_hparameters <- init_hparam(x, centers, labels, tol, iter_init, verbose)
  hlambda <- init_hparameters$hlambda
  hsigma <- init_hparameters$hsigma
  
  # if centers == 1, omega, ppi and qc0 are fixed
  if (centers == 1) {
    
    omega <- 1
    ppi <- 1
    qc0 <- rep(1, p)
    
  } else {
    
    # initial estimate of group-correlation matrix omega
    omega <- init_omega(x, centers, labels, hlambda, hsigma)
    # initial estimate of the probablity of the multi-nulli distribution
    ppi <- table(labels) / p
    # initial distribution of c based on ppi
    qc0 <- sapply(labels, function(grp) grp == c(1:centers)) * 1
  }
  
  # initial distribution of alpha based on omega, qc0, hlambda, hsigma
  qalpha <- obj_qalpha(x, centers, omega, qc0, hlambda, hsigma)
  
  # initial distribution of c based on omega, qalpha
  qc <- obj_qc(x, centers, ppi, omega, qalpha, hlambda, hsigma)
  
  # initial -logL
  obj_logL_val <- vector()
  obj_logL_val[1] <- obj_logL(
    x, centers, ppi, omega,
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
      x, centers, ppi_new, omega,
      qc, qalpha,
      hlambda, hsigma
    )
    
    
    # update omega
    omega_new <- obj_omega(centers, qalpha)
    
    min_val <- verbose_print(
      verbose, "omega", min_val,
      x, centers, ppi_new, omega_new,
      qc, qalpha,
      hlambda, hsigma
    )
    
    
    # update hlambda
    hlambda_new <- obj_hlambda(x, centers, qc, qalpha)
    
    min_val <- verbose_print(
      verbose, "hlambda", min_val,
      x, centers, ppi_new, omega_new,
      qc, qalpha,
      hlambda_new, hsigma
    )
    
    
    # update hsigma
    hsigma_new <- obj_hsigma(x, centers, qc, qalpha, hlambda_new)
    
    min_val <- verbose_print(
      verbose, "hsigma", min_val,
      x, centers, ppi_new, omega_new,
      qc, qalpha,
      hlambda_new, hsigma_new
    )
    
    # update qalpha
    qalpha_new <- obj_qalpha(x, centers, omega_new, qc, hlambda_new, hsigma_new)
    
    min_val <- verbose_print(
      verbose, "qalpha", min_val,
      x, centers, ppi_new, omega_new,
      qc, qalpha_new,
      hlambda_new, hsigma_new
    )
    
    
    # update qc
    qc_new0 <- obj_qc(
      x, centers, ppi_new, omega_new, qalpha_new,
      hlambda_new, hsigma_new
    )
    
    # discrete qc
    qc_newLabel <- apply(qc_new0, 2, function(x) which.max(x)[1])
    qc_new <- sapply(qc_newLabel, function(grp) grp == c(1:centers)) * 1
    
    min_val <- verbose_print(
      verbose, "qc", min_val,
      x, centers, ppi_new, omega_new,
      qc_new, qalpha_new,
      hlambda_new, hsigma_new
    )
    
    obj_logL_val[iiter + 1] <- obj_logL(
      x, centers, ppi_new, omega_new,
      qc_new, qalpha_new,
      hlambda_new, hsigma_new
    )
    
    if (abs(abs(obj_logL_val[iiter + 1] - obj_logL_val[iiter]) / obj_logL_val[iiter]) < tol) break
    
    iiter <- iiter + 1
    
    omega <- omega_new
    hlambda <- hlambda_new
    hsigma <- hsigma_new
    qc <- qc_new
    qalpha <- qalpha_new
  }
  
  cluster <- apply(qc, 2, which.max)
  
  list(
    omega = omega,
    hlambda = hlambda, hsigma = hsigma,
    obj_logL_val = -obj_logL_val, # minimize -> maximize
    qc = qc,
    cluster = cluster
  )
}





