#' Visulize the clustering results
#'
#' @description 
#' `colMat_reorder( )` reorder the columns of affinity matrix by specified groups. 
#'  
#' @param affMatrix a square affinity matrix.
#' @param centers number of clusters.
#' @param labels label assignment of the columns.
#'
#' @return
#' \item{matReorder}{affinity matrix reordered by groups.}
#' \item{groupLens}{a vector of group lengths.}
#' 
#' @importFrom dplyr %>%
#' @importFrom reshape2 melt
#' @importFrom stats cov
#' @importFrom graphics title
#' @importFrom ggplot2 aes element_blank
#' 
#' 
#' @export
#' @rdname cluster_visual
colMat_reorder <- function(affMatrix, centers, labels){
  
  # group-level matrix (average)
  groupMatrix <- matrix(NA, nrow = centers, ncol = centers)
  for (i in c(1:centers)){
    
    for(j in c(1:centers))
      
      groupMatrix[i, j] <- mean(affMatrix[labels == i, labels == j])
  }
  orderDiag <- order(abs(diag(groupMatrix)))
  
  # column index of each group elements
  colOrder <- list()
  for (i in c(1:centers)){
    groupIdx <- which(labels == i)
    colOrder[[i]] <- groupIdx
  }
  
  # number of elements of each group
  groupLens <- sapply(colOrder[orderDiag], length)
  
  # order columns by ascending group-level value
  colOrder2 <- colOrder[orderDiag] %>% unlist()
  # reorder the columns
  matReorder <- affMatrix[colOrder2, colOrder2]
  
  return(list(matReorder = matReorder, 
              groupLens = groupLens))
}


#' @description 
#' `colMat_heatMap( )` plots the heatmap of affMatrix ordered by groups.
#' @param margin a numeric value adjusting the position of group lines.
#' @param midpoint a numeric value of the midpoint of heatmap range.
#' @param limit a numeric range of the values.
#' @param legendName a character name of the legend.
#' @param title a character name of the heatmap.
#' @param size size of the groupline.
#' 
#' @export
#' @rdname cluster_visual
colMat_heatMap <- function(affMatrix, centers, labels,
                           margin = 0.5,
                           midpoint = 0, limit = c(-1,1), size = 0.2,
                           legendName = "", title = ""){
  
  mat <- colMat_reorder(affMatrix, centers, labels)
  matReorder <- mat$matReorder
  groupLens <- mat$groupLens
  
  groupLensCum <- cumsum(groupLens)
  axis_idxStart <- c(margin, groupLensCum[-length(groupLensCum)] + margin)
  axis_idxEnd <- c(groupLensCum + margin)
  matLines <- data.frame(x = c(axis_idxStart, axis_idxEnd),
                         xend = c(axis_idxEnd, axis_idxEnd),
                         y = c(axis_idxEnd, axis_idxStart),
                         yend = c(axis_idxEnd, axis_idxEnd))
  
  matReorder %>%
    reshape2::melt(na.rm = TRUE) %>%
    ggplot2::ggplot(ggplot2::aes(x=Var1, y=Var2, fill=value)) +
    ggplot2::geom_tile(color = "white")+
    ggplot2::scale_fill_gradient2(low = "blue", high = "red",
                         midpoint = midpoint, limit = limit, space = "Lab",
                         name=legendName) +
    ggplot2::theme(
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      axis.text.x = element_blank(),
      axis.text.y = element_blank(),
      panel.grid.major = element_blank(),
      panel.border = element_blank(),
      panel.background = element_blank(),
      axis.ticks = element_blank())+
    ggplot2::ggtitle(title) +
    ggplot2::coord_fixed()+
    ggplot2::geom_segment(data=matLines, 
                 aes(x, y,xend=xend, yend=yend), 
                 size=size, inherit.aes=F)
  
}

## global define the variable .
utils::globalVariables(".")


#' @description 
#' `crossValid( )` 2-fold crossValidation selects the number of clusters.
#' @param x an input matrix.
#' @param centers number of clusters.
#' @param pt partition times.
#' 
#' @export
#' @rdname cluster_visual
crossValid_func_adjR <- function(x, centers, pt){
  
  n <- nrow(x)
  p <- ncol(x)
  
  adjR <- sapply(c(1:pt), function(i){
    
    idx <- sample(c(1:n), size = floor(n/2), replace = FALSE)
    x_test <- x[idx,]
    x_valid <- x[-idx,]
    
    spec_test <- rSpecc(abs(cov(x_test)), centers = centers)$.Data
    hbcm_test <- hbcm::heterogbcm(scale(x_test, center = TRUE, scale = FALSE),
                                  centers = centers,
                                  tol = 1e-3, iter = 100, iter_init = 3,
                                  labels = spec_test,
                                  verbose = FALSE)
    
    spec_valid <- rSpecc(abs(cov(x_valid)), centers = centers)$.Data
    hbcm_valid <- hbcm::heterogbcm(scale(x_valid, center = TRUE, scale = FALSE),
                                   centers = centers,
                                   tol = 1e-3, iter = 100, iter_init = 3,
                                   labels = spec_valid,
                                   verbose = FALSE)
    
    matchLabel(hbcm_test$cluster, hbcm_valid$cluster)$adjRand
  })
  
  mean(adjR)
}


