#' Calculate the rand-index of two label assignments
#'
#' @description 
#' `matchLabel( )` calculate the rand-index and adjusted rand-index of estimated label assignment
#' compared to true label assignment, which can be used to evaluate the performance of the estimated
#' label. The metric takes value between 0 and 1, and higher value indicates better performance.
#' 
#' @param reference true label assignment.
#' @param label estimated label assignment.
#'
#' @return
#' \item{Rand}{rand-index}
#' \item{adjRand}{adjusted rand-index}
#' @export
#' @rdname cluster_eval
matchLabel <- function(reference, label) {
  max_ref <- max(reference)
  max_lab <- max(label)
  
  confmat <- matrix(0, max_ref, max_lab)
  for (i in 1:max_ref) {
    for (j in 1:max_lab)
    {
      confmat[i, j] <- sum((reference == i) & (label == j))
    }
  }
  
  # return
  calRand(confmat)
}

#' @rdname cluster_eval
#' @description 
#' `calRand( )` calculates the rand-index from a confusion matrix.
#' @param confmat a 2-dimensional confusion matrix.
#' @export
calRand <- function(confmat) {
  sum_total <- sum(confmat)
  sum_rows <- apply(confmat, 1, sum)
  sum_cols <- apply(confmat, 2, sum)
  
  pair_confmat <- choose(confmat, 2)
  pair_rows <- choose(sum_rows, 2)
  pair_cols <- choose(sum_cols, 2)
  adjRand <- (sum(pair_confmat) - (sum(pair_rows) * sum(pair_cols)) / choose(sum_total, 2)) / 
    (1 / 2 * (sum(pair_rows) + sum(pair_cols)) - 
       (sum(pair_rows) * sum(pair_cols)) / choose(sum_total, 2))
  
  a <- sum(pair_confmat)
  b <- sum(pair_rows) - sum(pair_confmat)
  c <- sum(pair_cols) - sum(pair_confmat)
  d <- choose(sum_total, 2) - a - b - c
  Rand <- (a + d) / choose(sum_total, 2)
  
  # return
  list(Rand = Rand, adjRand = adjRand)
}

#' @rdname cluster_eval
#' @description 
#' `rSpecc( )` a customized spectral clustering model.
#' @param x numeric matrix of data.
#' @param centers the number of clusters.
#' @param iter.max the maximum number of iterations allowed.
#' @param nstart how many random sets in the kmeans step should be chosen? default is 10.
#' 
#' @return 
#' \item{.Data}{A vector of integers indicating the cluster to which each point is allocated.}
#' \item{size}{The number of points in each cluster.}
#' \item{totss}{The total sum of squares.}
#' \item{withinss}{Vector of within-cluster sum of squares, one component per cluster.}
#' \item{tot.withinss}{Total within-cluster sum of squares, i.e. `sum(withinss)`.}
#' \item{betweenss}{The between-cluster sum of squares, i.e. `totss-tot.withinss`.}
#' @export
rSpecc <- function(x, centers, iter.max = 100, nstart = 10) {
  
  d <- 1 / sqrt(rowSums(x))
  l <- d * x %*% diag(d)
  
  eig_x <- RSpectra::eigs_sym(l, centers, which = "LM")$vectors
  eig_y <- eig_x / sqrt(rowSums(eig_x^2))
  
  res <- stats::kmeans(eig_y, centers, iter.max, nstart)
  
  # return
  list(.Data = res$cluster, 
       size = res$size,
       totss = res$totss,
       withinss = res$withinss,
       tot.withinss = res$tot.withinss,
       betweenss = res$betweenss
  )
}




