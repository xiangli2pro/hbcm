#################################################################################################
## Manuscript: Community Detection with Heterogeneous Block Covariance Model
## Simulation code for manuscript Figure 1 (middle)
#################################################################################################


# ------------ Install packages
library(devtools)
devtools::install_github("xiangli2pro/hbcm")
packages <- c("hbcm", "parallel", "foreach", "doParallel", "kernlab", "matrixcalc")
lapply(packages, require, char = TRUE)


# ------------ Simulate data
# set seed to replicate the result
set.seed(2022)

# mean vector
mu <- rep(0, centers)


###### -----------------------------------------------------------------------------------######
# to replicate the Figure 1 (middle),
# run the script for off_diag values in {0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9}
off_diag <- 0.1
omega <- diag(rep(1, centers))
for (i in 1:centers) {
  for (j in 1:centers) {
    if (i != j) {
      omega[i, j] <- off_diag
    }
  }
}
###### -----------------------------------------------------------------------------------######


# specify data size and cluster number
centers <- 3 # number of Classes
n <- 1000 # number of Observations
p <- 1000 # number of Genes

# equally distributed class
ppi <- rep(1 / centers, centers)
labels <- sample(c(1:centers), size = p, replace = TRUE, prob = ppi)

# take hlambda=1 and hsigma=6
hparam_func <- list(
  lambda_func = function(p) rep(1, p),
  sigma_func = function(p) rep(6, p)
)

# set up the number of simulation data
size <- 100

# generate data
data_list <- hbcm::data_gen(n, p, centers, mu, omega, labels, size, hparam_func)
# save data
save(data_list, labels,
  file = paste0("omega_", off_diag, ".rda")
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
  file = paste0("res_", "omega_", off_diag, ".rda")
)
