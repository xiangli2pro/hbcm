#' Transcript Factor (TF) regulators of yeast genes.
#'
#' The raw data is downloaded from the supplementary file of paper https://www.tandfonline.com/doi/suppl/10.1080/15384101.2019.1570655?scroll=top.
#' OR, see github https://github.com/xiangli2pro/hbcm/tree/main/inst/extdata/ to access the raw data.
#' Check github https://github.com/xiangli2pro/hbcm/tree/main/data-raw/yeast_gene_data.R to see how data is processed.
#'
#' @format A data frame with three variables: \code{Probe}, \code{Gene}, \code{TF_regulator}.
"data_gene_group"


#' Yeast gene expression from the experiment.
#'
#' The raw data is downloaded the paper https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE75694
#' OR, see github https://github.com/xiangli2pro/hbcm/tree/main/inst/extdata/ to access the raw data.
#' Check github https://github.com/xiangli2pro/hbcm/tree/main/data-raw/yeast_gene_data.R to see how data is processed.
#'
#' @format A data frame with 241 variables (genes).
"data_gene_sample"