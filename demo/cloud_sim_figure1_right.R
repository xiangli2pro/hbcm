#************************** set random seed
set.seed(2024)


#************************** generate data
# data size and cluster number
centers <- 3 # number of Classes
n <- 1000 # number of Rows
p <- 1000 # number of Columns

# mean vector of normal distribution
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

# set probability of clusters
ppi <- rep(1/centers, centers)
# set true labels
true_labels <- sample(c(1:centers), size = p, replace = TRUE, prob = ppi)

# set function to generate hlambda and hsigma
hlambda <- rep(1, p)
hsigma <- rep(6, p)

# generate data
sample_gen <- function(n, p, mu, omega, labels, hlambda, hsigma) {
  alpha <- MASS::mvrnorm(n, mu, omega)
  x <- matrix(rep(0, n * p), nrow = n)
  for (i in 1:n) {
    for (j in 1:p) {
      x[i, j] <- hlambda[j] * (alpha[i, labels[j]]) + hsigma[j] * rt(1, df_value) # update degree
    }
  }
  
  list(
    x_list = x, alpha_list = alpha,
    hlambda = hlambda, hsigma = hsigma
  )
}

data_list <- lapply(
  c(1:size), 
  function(i) sample_gen(n, p, mu, omega, true_labels, hlambda, hsigma)
)

save(data_list, true_labels,
     file = paste0(output_dir, "data_", "n",n,"_p",p,"_k",centers, '_tdist_df', df_value,".rda"))


#************************** model fit
## Spectral Clustering
X_list <- lapply(data_list, function(x) x$x_list)

start <- Sys.time()
spec_systime <- system.time(
  spec_labels <- lapply(X_list, 
                        function(x) rSpecc(abs(cor(x)), centers=centers)$.Data),
  gcFirst = FALSE
)
end <- Sys.time()
run_time_spec <- difftime(end, start, units="mins")

## HBCM
# cl <- parallel::makeCluster(2)
# doParallel::registerDoParallel(cl)
registerDoParallel(detectCores()-2)

start <- Sys.time()
hbcm_systime <- system.time(
  hbcm_res <- foreach(m=c(1:size),.errorhandling = 'pass',
                      .packages = c("MASS","Matrix","matrixcalc","kernlab", "RSpectra")) %dopar%
    hbcm::heterogbcm_converge_qc(scale(X_list[[m]], center = TRUE, scale = FALSE),
                                 centers = centers,
                                 tol = 1e-3, 
                                 iter = 200, 
                                 iter_init = 3,
                                 labels = spec_labels[[m]],
                                 verbose = FALSE),
  gcFirst = FALSE
)
end <- Sys.time()
run_time_hbcm <- difftime(end, start, units="mins")
# parallel::stopCluster(cl)

#************************** summary
## Generate Summary
get_ARIsummary <- function(
    spec_labels, hbcm_res, 
    run_time_spec, run_time_hbcm,
    true_labels, len, sim_name
){
  
  # get result where HBCM is not NULL
  nonna_indx <- sapply(hbcm_res, function(x) length(x)==len)
  nonna_num <- sum(nonna_indx)
  # glue('{nonna_num} number of simulations return non-Null value')
  
  # summary of spectral
  summary_spec <- c(mean(sapply(spec_labels[nonna_indx], function(x) matchLabel(true_labels, x)$adjRand)),
                    sd(sapply(spec_labels[nonna_indx], function(x) matchLabel(true_labels, x)$adjRand), na.rm=TRUE))
  
  # summary of HBCM
  summary_hbcm <- c(mean(sapply(hbcm_res[nonna_indx], function(x) matchLabel(true_labels, x$cluster)$adjRand)),
                    sd(sapply(hbcm_res[nonna_indx], function(x) matchLabel(true_labels, x$cluster)$adjRand), na.rm=TRUE))
  
  # combine results
  cbind(sim_name, 
        paste0(round(summary_spec[1],2), " (", round(summary_spec[2],2), ")"), 
        run_time_spec,
        paste0(round(summary_hbcm[1],2), " (", round(summary_hbcm[2],2), ")"),
        run_time_hbcm) %>% 
    as.data.frame() %>%
    `colnames<-`(c('Simulation', "Spectral", 'Spectral Time', "HBCM", "HBCM Time"))
}
sim_summary <- get_ARIsummary(spec_labels, hbcm_res, 
                              run_time_spec, run_time_hbcm,
                              true_labels, 
                              len=7, sim_name=paste0('sim_k',centers,'_n',n,'_p',p, '_tdist_df', df_value))

save(spec_systime, hbcm_systime,
     run_time_spec, run_time_hbcm,
     spec_labels, hbcm_res, true_labels, 
     sim_summary,
     file = paste0(output_dir,"res_", "n",n,"_p",p,"_k", centers, '_tdist_df', df_value, ".rda"))


