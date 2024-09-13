# ------------ SET UP
# !!! if you want to save R packages in a different directory other than default, put it in the first argument

# .libPaths(c("/CCAS/home/xiangli/Rpackages", .libPaths()))
# .libPaths()

# !!! make sure the required packages are already installed in the cloud server R environment

# install hbcm from github
# devtools::install_github("xiangli2pro/hbcm")
# library(devtools)
# load_all()
# document()
packages <- c("hbcm", "parallel", "foreach", "doParallel", "kernlab", "matrixcalc")
lapply(packages, require, char=TRUE)

# !!! update the input and output path in cloud server
input_dir = '/Users/xiangli/Dropbox/GWU/hbcm/demo/cloud_sim_figure2.R'  
#!!!must keep the / by the end
output_dir = '/Users/xiangli/Downloads/hbcm_jcgs_revise/simulation/' 
# specify how many replications for each simulation scenario
size = 100

# ------------ RUN SIMULATION
# execute simulation for figure1 (left)
for (hsigma_value in seq(10, 15, 1)){
  source(input_dir)
}


