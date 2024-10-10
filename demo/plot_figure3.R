#################################################################################################
## Manuscript: Community Detection with Heterogeneous Block Covariance Model
## Plot code for manuscript Figure 3 (cross-validation)
#################################################################################################


#************************** Install packages
library(devtools)
devtools::install_github("xiangli2pro/hbcm")
# load_all()
# document()
library(hbcm)
library(tidyverse)
library(ggpubr)
library(scales)
library(ggforce)
library(patchwork)

set.seed(2024)
output_dir = '/Users/xiangli/Downloads/hbcm_jcgs_revise/simulation/final/' 

#************************** Simulate data

# make simulation data

# data size and cluster number
centers <- 5 # number of Classes
n <- 1500 # number of Observations
p <- 500 # number of Genes

# mean vector
mu <- rep(0, centers)

# cluster level covariance matrix
off_diag <- 0.2
omega <- diag(rep(1, centers))
for (i in 1:centers) {
  for (j in 1:centers) {
    if (i != j) {
      omega[i, j] <- off_diag
    }
  }
}

# equally distributed class
ppi <- rep(1 / centers, centers)
true_labels <- sample(c(1:centers), size = p, replace = TRUE, prob = ppi)

# individual hlambda and hsigma
hparam_func <- list(
  lambda_func = function(p) rep(1, p),
  sigma_func = function(p) rep(6, p)
)
# sample data
data_list <- hbcm::data_gen(n, p, centers, mu, omega, true_labels, size=1, hparam_func)
X <- data_list$x_list[[1]]


#************************** Use cross-validation to identify cluster number K in simulation data

# parallel process
registerDoParallel(detectCores())

# candidate K
kVec <- c(2:9)
# get adjRand from cross-validation 
res <- foreach(
  K = kVec, .errorhandling = "pass",
  .packages = c("MASS", "Matrix", "matrixcalc", "kernlab", "RSpectra")
) %dopar%
  hbcm::crossValid_func_adjR(X, K, pt = 20)

# remove NaN
cv_res <- sapply(res, function(x) {
  if (length(x) == 1) {
    x
  } else {
    NA
  }
})

# prepare data for plot
cv_data <- cbind(kVec, cv_res) %>%
  `colnames<-`(c("K", "adjRand")) %>%
  as.data.frame()

# make plot
pSim <- ggplot(cv_data) +
  geom_line(aes(x = K, y = adjRand)) +
  xlab("K") +
  ylab("Adjusted Rand Index") +
  theme_bw() +
  theme(panel.grid.minor = element_blank()) +
  scale_x_continuous(breaks = c(2:9)) +
  geom_vline(xintercept = 5, linetype = "dashed", color = "#f37735")
pSim
# save data
load(paste0(output_dir, "data/figure3_sim.rda"))
save(true_labels, data_list, cv_data, pSim, file = paste0(output_dir, "data/figure3_sim.rda"))
ggsave(paste0(output_dir, "plot/figure3_sim.png"), dpi = 600)
ggsave(paste0(output_dir, "plot/figure3_sim.eps"), dpi = 600)

#************************** Use cross-validation to identify cluster number K in gene expression data

# load data from the package, see more in help page
hbcm::data_gene_sample?
data("data_gene_sample")
data_gene_sample = data_gene_sample
sum(data_gene_sample>0)
sum(data_gene_sample<0)

data_centered <- scale(data_gene_sample, center = TRUE, scale = FALSE)
c(data_centered)

library(plot.matrix)
par(mar=c(5.1, 4.1, 4.1, 4.1)) # adapt margins
plot(spec_grp_cor_abs, digits=1, text.cell=list(cex=0.8))
plot(spec_grp_cor, digits=1, text.cell=list(cex=0.8))
plot(spec_grp_cov_abs, digits=1, text.cell=list(cex=0.8))
plot(spec_grp_cov, digits=1, text.cell=list(cex=0.8))
plot(hbcm_grp_cor_abs, digits=1, text.cell=list(cex=0.8))
plot(hbcm_grp_cov_abs, digits=1, text.cell=list(cex=0.8))
plot(hbcm_grp_cor, digits=1, text.cell=list(cex=0.8))
plot(hbcm_grp_cov, digits=1, text.cell=list(cex=0.8))
plot(spec_grp_cor-hbcm_grp_cor, digits=1, text.cell=list(cex=0.8))
plot(abs(spec_grp_cor_abs-hbcm_grp_cor_abs)/spec_grp_cor_abs, digits=1, text.cell=list(cex=0.8))


# parallel computing
registerDoParallel(detectCores())

# cross-validation to calculate adjRand
kVec <- c(2:20)
res_adjR <- foreach(
  K = kVec, .errorhandling = "pass",
  .packages = c("MASS", "Matrix", "matrixcalc", "kernlab", "RSpectra")
) %dopar%
  hbcm::crossValid_func_adjR(data_gene_sample, K, pt = 20)

spec_labels_gene <- rSpecc(abs(cov(data_gene_sample)), centers =5)$.Data

# -----------------HBCM cluster estimation with Spectral initial
hbcm_spec_gene <-  hbcm::heterogbcm_converge_qc(
  scale(data_gene_sample, center = TRUE, scale = FALSE), # standardize the data
  centers = 5, # set number of clusters
  tol = 1e-2, # iteration stop criteria
  iter = 200, # max iteration steps
  iter_init = 3, # iteration steps for estimation of hsigma and hlambda
  labels = spec_labels_gene, # initial label estimation
  verbose = FALSE # whether or not to print out the iteration message
)
length(hbcm_spec_gene$hlambda[hbcm_spec_gene$hlambda<0])
length(hbcm_spec_gene$hlambda[hbcm_spec_gene$hlambda>=0])

res_adjR_clean <- unlist(res_adjR[1:13])

# make plot
pReal <- cbind(kVec[1:13], res_adjR_clean) %>%
  `colnames<-`(c("K", "adjR_mean")) %>%
  as.data.frame() %>%
  ggplot() +
  geom_line(aes(x = K, y = adjR_mean)) +
  xlab("K") +
  ylab("Adjusted Rand Index") +
  theme_bw() +
  theme(panel.grid.minor = element_blank()) +
  scale_x_continuous(breaks = kVec) +
  geom_vline(xintercept = 5, linetype = "dashed", color = "#f37735")
ggtitle("Real data\n2-fold crossValidation of adjR")
pReal

# save data
load(paste0(output_dir, "data/figure3_real.rda"))
save(res_adjR, file = paste0(output_dir, "data/figure3_real.rda"))
ggsave(paste0(output_dir, "plot/figure3_real.png"), dpi = 600)
ggsave(paste0(output_dir, "plot/figure3_real.eps"), dpi = 600)

# ------------ Combine plots to a single plot

ggarrange(pSim, pReal,
          ncol = 2, nrow = 1, align = "hv",
          common.legend = TRUE
)
ggsave(paste0(output_dir, "plot/figure3_sim_real.png"), dpi = 600, height = 4, width = 8)
ggsave(paste0(output_dir, "plot/figure3_sim_real.eps"), dpi = 600, height = 4, width = 8)
