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
input_dir_f1_left = '/Users/xiangli/Dropbox/GWU/hbcm/demo/cloud_sim_figure1_left.R' 
input_dir_f1_middle = '/Users/xiangli/Dropbox/GWU/hbcm/demo/cloud_sim_figure1_middle.R' 
input_dir_f1_right = '/Users/xiangli/Dropbox/GWU/hbcm/demo/cloud_sim_figure1_right.R' 
input_dir_f1_right_corrected = '/Users/xiangli/Dropbox/GWU/hbcm/demo/cloud_sim_figure1_right_corrected.R' 

#!!!must keep the / by the end
output_dir = '/Users/xiangli/Downloads/hbcm_jcgs_revise/simulation/' 
# specify how many replications for each simulation scenario
size = 100

# ------------ RUN SIMULATION
# execute simulation for figure1 (left)
for (hsigma_value in seq(2, 10, 1)){
  source(input_dir_f1_left)
}

# execute simulation for figure1 (middle)
for (off_diag_value in seq(0.1, 0.9, 0.1)){
  source(input_dir_f1_middle)
}

# execute simulation for figure1 (right)
for (df_value in seq(3, 10, 1)){
  source(input_dir_f1_right)
}

# execute simulation for figure1 (right, standardized t-dist)
for (df_value in seq(3, 10, 1)){
  source(input_dir_f1_right_corrected)
}

