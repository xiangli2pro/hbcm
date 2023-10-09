#################################################################################################
## Manuscript: Community Detection with Heterogeneous Block Covariance Model
## Simulation code for manuscript Figure 1 (left)
#################################################################################################


# ------------ Install packages
library(devtools)
devtools::install_github("xiangli2pro/hbcm")
packages <- c("hbcm", "parallel", "foreach", "doParallel", "kernlab", "matrixcalc")
lapply(packages, require, char = TRUE)


# ------------ Simulate data
# set seed to replicate the result
set.seed(2022)


###### -----------------------------------------------------------------------------------######
# to replicate the manuscript Figure 1 (left),
# run the script for hsigma_value values in {2,3,4,5,6,7,8,9,10}
hsigma_value <- 2
hparam_func <- list(
  lambda_func = function(p) rep(1, p),
  sigma_func = function(p) rep(hsigma_value, p)
)
###### -----------------------------------------------------------------------------------######


# data size and cluster number
centers <- 3 # number of Classes
n <- 1000 # number of Observations
p <- 1000 # number of Genes

# mean vector
mu <- rep(0, centers)

# class-level covariance matrix
off_diag <- 0.5
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

# set up the number of simulation data
size <- 100

# generate data
data_list <- hbcm::data_gen(n, p, centers, mu, omega, labels, size, hparam_func)
# save data
save(data_list, labels,
  file = paste0("hlambda", 1, "_hsigma", hsigma_value, ".rda")
)


# ------------ Fit cluster model

# get simulated data
X_list <- data_list$x_list
# spectral clustering class estimation
spec_tm <- system.time(
  spec_labels <- lapply(X_list, function(x) rSpecc(abs(cor(x)), centers = centers)$.Data),
  gcFirst = FALSE
)


# HBCM class estimation
# run in parallel
registerDoParallel(detectCores())

hbcm_tm <- system.time(
  hbcm_res <- foreach(
    m = c(1:size), .errorhandling = "pass",
    .packages = c("MASS", "Matrix", "matrixcalc", "kernlab", "RSpectra")
  ) %dopar%
    hbcm::heterogbcm(scale(X_list[[m]], center = TRUE, scale = FALSE),
      centers = centers,
      tol = 1e-3, iter = 100, iter_init = 3,
      labels = spec_labels[[m]],
      verbose = FALSE
    ),
  gcFirst = FALSE
)

parallel::stopCluster(cl)

# save results
save(spec_tm, hbcm_tm,
  spec_labels, hbcm_res,
  labels,
  file = paste0("res_", "hlambda", 1, "_hsigma", hsigma_value, ".rda")
)
