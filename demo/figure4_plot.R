#################################################################################################
## Manuscript: Community Detection with Heterogeneous Block Covariance Model
## Plot code for manuscript Figure 4 (heatmap of gene expression cluster)
#################################################################################################


# ------------ Install packages
library(devtools)
devtools::install_github("xiangli2pro/hbcm")
library(hbcm)
library(tidyverse)
library(ggpubr)
library(scales)
library(ggforce)
library(patchwork)

set.seed(202205)


# ------------ Load gene expression data and true cluster info
hbcm::data_gene_sample?
hbcm::data_gene_group?


# ------------ Gene cluster estimated by HBCM and spectral cluster

# assume K=5
centers <- 5
# spectral cluster estimation
spec_labels <- rSpecc(abs(cov(data_gene_sample)), centers = centers)$.Data
# HBCM cluster estimation
hbcm_res <- hbcm::heterogbcm(scale(data_gene_sample, center = TRUE, scale = FALSE),
  centers = centers,
  tol = 1e-3, iter = 100, iter_init = 3,
  labels = spec_labels,
  verbose = FALSE
)
# adjRand between HBCM and spectral cluster
matchLabel(spec_labels, hbcm_res$cluster)


# ------------ Compare the estimation cluster with recommended cluster

# get recommended cluster label
matchProbe <- match(names(data_gene_sample), data_gene_group$Probe)
matchTF <- data_gene_group$TF_regulator[matchProbe]
uniTFs <- unique(matchTF)

# convert to numerical label
trueLabel <- rep(0, ncol(data_gene_sample))
for (i in c(1:length(uniTFs))) {
  trueLabel[which(matchTF == uniTFs[i])] <- i
}
table(trueLabel)

# adjRand of recommended label VS HBCM
matchLabel(trueLabel, hbcm_res$cluster)$adjRand

# adjRand of recommended label VS spectral clustering
matchLabel(trueLabel, spec_labels)$adjRand


# ------------ Make heatmap plot

colheatMap_spectral <- hbcm::colMat_heatMap(
  affMatrix = cor(data_gene_sample), centers, labels = hbcm_res$cluster,
  margin = 0.5, midpoint = 0, limit = c(-1,1), size = 0.2,
  legendName = "Correlation", title = "HeatMap of correlation by HBCM groups")

ggsave("plots/heatmap_HBCM.png", dpi = 300)
