#************************** setup environment
library(devtools)
load_all()
packages <- c("hbcm", "parallel", "foreach", "doParallel", "kernlab", "matrixcalc")
lapply(packages, require, char=TRUE)
output_dir = '/Users/xiangli/Downloads/hbcm_jcgs_revise/simulation/final/data/' 


#************************** set random seed
set.seed(2024)


#************************** generate data
# replicate times
size = 100

# n=500, p=300
centers <- 3 # number of Classes
n <- 500  # number of Observations
p <- 300  # number of Genes

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
     file = paste0(output_dir, "data_", "n",n,"_p",p,"_k",centers,"_wgcna.rda"))


#************************** model fit
## Spectral Clustering
X_list <- data_list$x_list

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

#************************** WGCNA
# # download packages
# install.packages("WGCNA")
# install.packages("BiocManager")
# BiocManager::install("GO.db")
# BiocManager::install("impute")
# BiocManager::install("preprocessCore")
# BiocManager::install("DESeq2")
library(WGCNA)

wgcna_power6_size10_cut0.6 = list()
start <- Sys.time()
for (m in c(1:size)){
  wgcna_power6_size10_cut0.6[[m]] <- blockwiseModules(
    X_list[[m]], power = 6,
    TOMType = "unsigned",
    minModuleSize = 10,
    mergeCutHeight = 0.60,
    # reassignThreshold = 0, 
    numericLabels = TRUE, 
    # pamRespectsDendro = FALSE,
    saveTOMs = FALSE,
    # saveTOMFileBase = "femaleMouseTOM",
    verbose = 0
  )
}
end <- Sys.time()
run_time_wgcna_power6_size10_cut0.6 <- difftime(end, start, units="mins")
run_time_wgcna_power6_size10_cut0.6


#************************** summary
## Generate Summary
get_ARIsummary <- function(
    spec_labels, hbcm_res, wgcna_res,
    run_time_spec, run_time_hbcm, run_time_wgcna,
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
  
  # summary of WGCNA
  summary_wgcna <- c(mean(sapply(wgcna_res, function(x) matchLabel(true_labels, x$colors+1)$adjRand)),
                             sd(sapply(wgcna_res, function(x) matchLabel(true_labels, x$colors+1)$adjRand), na.rm=TRUE))
  
  # combine results
  cbind(sim_name, 
        paste0(round(summary_spec[1],2), " (", round(summary_spec[2],2), ")"), 
        run_time_spec,
        paste0(round(summary_hbcm[1],2), " (", round(summary_hbcm[2],2), ")"),
        run_time_hbcm,
        paste0(round(summary_wgcna[1],2), " (", round(summary_wgcna[2],2), ")"),
        run_time_wgcna) %>% 
    as.data.frame() %>%
    `colnames<-`(c('Simulation', "Spectral", 'Spectral Time', "HBCM", "HBCM Time", "WGCNA", "WGCNA Time"))
}
sim_summary <- get_ARIsummary(spec_labels, hbcm_res, wgcna_power6_size10_cut0.6,
                              run_time_spec, run_time_hbcm, run_time_wgcna_power6_size10_cut0.6,
                              true_labels, 
                              len=7, sim_name=paste0('sim_k',centers,'_n',n,'_p',p))

save(spec_systime, hbcm_systime, 
     run_time_spec, run_time_hbcm, run_time_wgcna_power6_size10_cut0.6,
     spec_labels, hbcm_res, true_labels, wgcna_power6_size10_cut0.6,
     sim_summary,
     file = paste0(output_dir,"res_", "n",n,"_p",p,"_k", centers,"_wgcna.rda"))

# analysis
wgcna_cluster_num = unlist(lapply(wgcna_power6_size10_cut0.6, function(x) max(x$colors+1)))
table(wgcna_cluster_num)
# wgcna_cluster_num
# 2  3  4 
# 29 68  3 
wgcna_cluster_num_max = unlist(lapply(wgcna_power6_size10_cut0.6, function(x) max(table(x$colors+1))))
quantile(wgcna_cluster_num_max)
# 0%  25%  50%  75% 100% 
# 226  240  245  252  281 



  