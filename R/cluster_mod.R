#' Estimate the optimal posterior distribution of the data column labels.
#'
#' @description
#' `heterogbcm( )` gives the optimal posterior distribution of the labels, which can be used
#' to derive the optimal label assignment of the data columns.
#'
#' @return A list of values.
#' \item{omega}{estimated optimal group-correlation matrix.}
#' \item{hlambda}{estimated optimal heterogeneous parameter Lambda.}
#' \item{hsigma}{estimated optimal heterogeneous parameter Sigma.}
#' \item{obj_logL_val}{vector of -logL from each iteration.}
#' \item{qalpha}{estimated optimal posterior distribution of the alpha.}
#' \item{qc}{estimated optimal posterior distribution of the column labels.}
#' \item{cluster}{a vector of integers (from 1:k) indicating the cluster to which each column is allocated.}
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
#' @param omega group-correlation matrix.
#' @param qc posterior distribution of labels.
#' @export
heterogbcm <- function(x, centers, tol, iter, iter_init = 3, labels, verbose = FALSE) {
  
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

# Helper function -------------------------------------

## print iteration info
verbose_print <- function(verbose, param_name, min_val,
                          x, centers, ppi, omega, qc, qalpha, hlambda, hsigma) {
  
  obj_logL_val_new <- obj_logL(
    x, centers, ppi, omega,
    qc, qalpha,
    hlambda, hsigma
  )
  
  if (obj_logL_val_new <= min_val) {
    min_val <- obj_logL_val_new
    if (verbose == TRUE) cat(paste0("Update ", param_name, " : ", obj_logL_val_new, " --"), "\n")
  } else {
    if (verbose == TRUE) cat(paste0("Update ", param_name, " : ", obj_logL_val_new, " ++"), "\n")
  }

  # return
  min_val
}

## alternative functions
#' @rdname cluster_mod
#' @export
#' @description
#' `get_random_qc( )` assign initial prob between 0 and 1 to cluster matrix
get_random_qc <- function(grp, centers){
  
  # get discrete 1 or 0 for each cluster
  qc_vec <- (grp == c(1:centers))*1
  # assign maximum probability for the cluster equal to 1
  idx_max <- qc_vec == 1
  max_prob <- runif(1, min=0.5, max = 1)
  # max_prob <- runif(1, min=1/centers, max = 2/centers)
  # assign probabilities to other cluster with sum equal to 1
  qc_vec[idx_max] <- max_prob
  qc_vec[!idx_max] <- get_prob(centers-1, 1 - max_prob)
  
  qc_vec
}

#' @rdname cluster_mod
#' @export
#' @description
#' `get_prob( )` get num_prob with sum equal to total_prob
get_prob <- function(num_prob, total_prob) {
  # Generate k-1 random values between 0 and 1
  random_values <- runif(num_prob - 1)
  # Append 0 and 1 to the list of random values and sort
  sorted_values <- sort(c(0, random_values, 1))
  # Compute the differences between successive elements
  prob_vec <- diff(sorted_values) * total_prob
  return(prob_vec)
}

#' @rdname cluster_mod
#' @export
#' @description
#' `heterogbcm_prob_qc( )` set initial distribution of label to be probability rather than discrete value.
heterogbcm_prob_qc <- function(x, centers, tol, iter, iter_init = 3, labels, verbose = FALSE) {
  
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
    qc0 <- sapply(labels, get_random_qc, centers=centers)
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
    qalpha = qalpha,
    qc = qc,
    cluster = cluster
  )
}


#' @rdname cluster_mod
#' @export
#' @description
#' `heterogbcm_converge_qc( )` set initial distribution of label to be probability rather than discrete value.
heterogbcm_converge_qc <- function(x, centers, tol, iter, iter_init = 3, labels, verbose = FALSE) {
  
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
    qc0 <- sapply(labels, get_random_qc, centers=centers)
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
  
  qc_mat <- list()
  qc_mat[[1]] <- qc
  
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
    
    qc_mat[[iiter+1]] <- qc_new
    # use max. abs. diff. between two probability matrix as criteria
    if (max(abs(qc_mat[[iiter+1]] - qc_mat[[iiter]])) < tol) break
    
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


