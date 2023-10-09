#################################################################################################
## Manuscript: Community Detection with Heterogeneous Block Covariance Model
## Simulation code for manuscript Figure 1 (right) t-distribution with correction
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
# to replicate the Figure 1 (right) t-distribution without correction
# run the script for degree of freedom (df) values in {3, 4, 5, 6, 7}
df <- 3
sample_gen <- function(n, p, mu, omega, labels, hlambda, hsigma) {
  alpha <- MASS::mvrnorm(n, mu, omega)
  x <- matrix(rep(0, n * p), nrow = n)
  for (i in 1:n) {
    for (j in 1:p) {
      x[i, j] <- hlambda[j] * (alpha[i, labels[j]]) + hsigma[j] * rt(1, df) / (sqrt(df / (df - 2))) # t-dist correction
    }
  }

  list(
    x = x, alpha = alpha,
    hlambda = hlambda, hsigma = hsigma
  )
}
###### -----------------------------------------------------------------------------------######


# cluster number and data size
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

# take hlambda=1 and hsigma=6
hlambda <- rep(1, p)
hsigma <- rep(6, p)

# set up the number of simulation data
size <- 100

# generate data
data_list <- lapply(c(1:size), function(i) sample_gen(n, p, mu, omega, labels, hlambda, hsigma))
# save data
save(data_list, labels,
  file = paste0("t_corrected_", df, ".rda")
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
  file = paste0("res_", "t_corrected_", df, ".rda")
)
