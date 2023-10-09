#################################################################################################
## Manuscript: Community Detection with Heterogeneous Block Covariance Model
## Simulation code for manuscript Table 1
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
# to replicate the Table 1 result,
# run the script for all the combinations of centers={3,5,7,10} X n={1000, 3000} X p={500, 1000, 3000, 5000}
centers <- 3 # number of Classes
n <- 3000 # number of Observations
p <- 500 # number of Genes
###### -----------------------------------------------------------------------------------######


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

# take random values for hlambda and hsigma
hparam_func <- list(
  lambda_func = function(p) rnorm(p, 0, 1),
  sigma_func = function(p) rchisq(p, 2) + 1
)

# set up the number of simulation data
size <- 100

# generate data
data_list <- hbcm::data_gen(n, p, centers, mu, omega, labels, size, hparam_func)
# save data
save(data_list, labels,
  file = paste0("n", n, "_p", p, "_K", centers, ".rda")
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

# save data
save(spec_tm, hbcm_tm,
  spec_labels, hbcm_res,
  labels,
  file = paste0("res_", "n", n, "_p", p, "_K", centers, ".rda")
)
