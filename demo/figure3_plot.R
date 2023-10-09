#################################################################################################
## Manuscript: Community Detection with Heterogeneous Block Covariance Model
## Plot code for manuscript Figure 3 (cross-validation)
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


# ------------ Simulate data

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
labels <- sample(c(1:centers), size = p, replace = TRUE, prob = ppi)

# individual hlambda and hsigma
hlambda <- rep(1, p)
hsigma <- rep(6, p)

# sample data
data_list <- sample_gen(n, p, mu, omega, labels, hlambda, hsigma)
X <- data_list$x


# ------------ Use cross-validation to identify cluster number K in simulation data

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
  `colnames<-`(c("K", "logLikelihood")) %>%
  as.data.frame()

# make plot
pSim <- ggplot(cv_data) +
  geom_line(aes(x = K, y = logLikelihood)) +
  xlab("K") +
  ylab("Adjusted Rand Index") +
  theme_bw() +
  theme(panel.grid.minor = element_blank()) +
  scale_x_continuous(breaks = c(2:9)) +
  geom_vline(xintercept = 5, linetype = "dashed", color = "#f37735")

# save data
save(cv_data, file = "data/K5_SimData_0523.rda")
ggsave("plots/K5_Sim.png", dpi = 300)


# ------------ Use cross-validation to identify cluster number K in gene expression data

# load data from the package, see more in help page
hbcm::data_gene_sample?
  
# parallel computing
registerDoParallel(detectCores())

# cross-validation to calculate adjRand
hbcm::data_gene_sample?
kVec <- c(2:20)
res_adjR <- foreach(
  K = kVec, .errorhandling = "pass",
  .packages = c("MASS", "Matrix", "matrixcalc", "kernlab", "RSpectra")
) %dopar%
  crossValid_func_adjR(X = gene_cdkData, K, pt = 20)

# make plot
pReal <- cbind(kVec, res_adjR %>% unlist()) %>%
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

# save plots
save(res_adjR, file = "data/K5_RealData.rda")
ggsave("plots/K5_Real.png", dpi = 300)


# ------------ Combine plots to a single plot

ggarrange(pSim, pReal,
  ncol = 2, nrow = 1, align = "hv",
  common.legend = TRUE
)

ggsave("./plots/KbyCV.eps", dpi = 300, width = 8, height = 4)
# ggsave("./plots/KbyCV.png", dpi = 300, width = 8, height = 4)
