#************************** set random seed
set.seed(2024)


#************************** generate data
# set alpha mean
mu <- rep(0, centers)

# set omega
off_diag <- 0.5
omega <- diag(rep(1, centers))
for (i in 1:centers) {
  for (j in 1:centers) {
    if (i!=j){
      omega[i,j] = off_diag
    }
  }
}

# set probability of clusters
ppi <- rep(1/centers, centers)
# set true labels
true_labels <- sample(c(1:centers), size = p, replace = TRUE, prob = ppi)

# set function to generate hlambda and hsigma
hparam_func <- list(
  lambda_func = function(p) rnorm(p, 0, 1),
  sigma_func = function(p) rchisq(p, 2) + 1
)

# generate data
data_list <- hbcm::data_gen(n, p, centers, mu, omega, true_labels, size, hparam_func)

save(data_list, true_labels,
     file = paste0(output_dir, "data_", "n",n,"_p",p,"_k",centers,".rda"))


#************************** model fit
## Spectral Clustering
X_list <- data_list$x_list

# start <- Sys.time()
# spec_systime <- system.time(
#   spec_labels <- lapply(X_list, 
#                         function(x) rSpecc(abs(cor(x)), centers=centers)$.Data),
#   gcFirst = FALSE
# )
# end <- Sys.time()
# run_time_spec <- difftime(end, start, units="mins")

## HBCM
# cl <- parallel::makeCluster(2)
# doParallel::registerDoParallel(cl)
registerDoParallel(detectCores()-2)

start <- Sys.time()
wgcna_systime <- system.time(
 wgcna_res <- foreach(m=c(1:size),.errorhandling = 'pass',
                      .packages = c("WGCNA")) %dopar%
    blockwiseModules(
      X_list[[m]], power = power, 
      TOMType = "unsigned",
      minModuleSize = minModuleSize, 
      mergeCutHeight = mergeCutHeight, 
      # reassignThreshold = 0, 
      numericLabels = TRUE, 
      # pamRespectsDendro = FALSE,
      saveTOMs = FALSE,
      # saveTOMFileBase = "femaleMouseTOM",
      verbose = 0
    ),
  gcFirst = FALSE
)
end <- Sys.time()
run_time_wgcna <- difftime(end, start, units="mins")
# parallel::stopCluster(cl)

#************************** summary
## Generate Summary
get_ARIsummary <- function(
    wgcna_res, run_time_wgcna,
    true_labels, sim_name
    ){
  
  # summary of WGCNA
  summary_wgcna <- c(mean(sapply(wgcna_res, function(x) matchLabel(true_labels, x$colors+1)$adjRand)),
                     sd(sapply(wgcna_res, function(x) matchLabel(true_labels, x$colors+1)$adjRand), na.rm=TRUE))
  
  # combine results
  cbind(sim_name, 
        paste0(round(summary_wgcna[1],2), " (", round(summary_wgcna[2],2), ")"), 
        run_time_wgcna) %>% 
    as.data.frame() %>%
    `colnames<-`(c('Simulation', "WGCNA", "WGCNA Time"))
}
sim_summary <- get_ARIsummary(wgcna_res, run_time_wgcna,
                              true_labels, 
                              sim_name=paste0('sim_k',centers,'_n',n,'_p',p))

save(wgcna_systime, run_time_wgcna, wgcna_res, true_labels, sim_summary,
     file = paste0(output_dir,"res_", "n",n,"_p",p,"_k", centers,"_wgcna.rda"))


