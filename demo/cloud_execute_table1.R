# ------------ SET UP
# !!! if you want to save R packages in a different directory other than default, put it in the first argument

# .libPaths(c("/CCAS/home/xiangli/Rpackages", .libPaths()))
# .libPaths()

# !!! make sure the required packages are already installed in the cloud server R environment

# install hbcm from github
# devtools::install_github("xiangli2pro/hbcm")
packages <- c("hbcm", "parallel", "foreach", "doParallel", "kernlab", "matrixcalc")
lapply(packages, require, char=TRUE)

# !!! update the input and output path in cloud server
input_dir = '/Users/xiangli/cloud_sim_table1.R'
output_dir = '/Users/xiangli/table1_sim_results/' #!!!must keep the / by the end
# specify how many replications for each simulation scenario
size = 100

# ------------ RUN SIMULATION
# execute simulation for table 1
# n=500
centers <- 3 # number of Classes
n <- 500  # number of Observations
p <- 300  # number of Genes
source(input_dir)

centers <- 5   # number of Classes
n <- 500  # number of Observations
p <- 500  # number of Genes
source(input_dir)

centers <- 7   # number of Classes
n <- 500  # number of Observations
p <- 1000  # number of Genes
source(input_dir)


# n = 1000
centers <- 3 # number of Classes
n <- 1000  # number of Observations
p <- 500  # number of Genes
source(input_dir)

centers <- 5   # number of Classes
n <- 1000  # number of Observations
p <- 1000  # number of Genes
source(input_dir)

centers <- 7   # number of Classes
n <- 1000  # number of Observations
p <- 1500  # number of Genes
source(input_dir)




