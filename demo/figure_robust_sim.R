# install helper package
packages_helper <- c(
  'devtools', "parallel",
  "foreach", "doParallel",
  "kernlab", "matrixcalc",
  "tidyverse", "ggpubr",
  "scales", "ggforce",
  "patchwork"
)
install.packages(packages_helper, repos = "http://cran.us.r-project.org")

# install hbcm package
devtools::install_github("xiangli2pro/hbcm")

# load all required packages
lapply(c('hbcm', packages_helper), require, char = TRUE)

# set random seed to reproduce the result
set.seed(2024)
output_dir = '/Users/xiangli/Downloads/hbcm_jcgs_revise/simulation/final/plot/' 

df <- data.frame(
  Value = c(0:5),
  Spectral = c(0.67, 0.68, 0.67, 0.63, 0.52, 0.30),
  HBCM = c(0.83, 0.83, 0.82, 0.82, 0.81, 0.71)
)

# make manuscript Figure 1 (left)
plot_hsigma <- df %>%
  as.data.frame() %>%
  `colnames<-`(c("Value", "Spectral", "HBCM")) %>%
  ggplot() +
  geom_line(aes(x = Value, y = Spectral, color = "Spectral")) +
  geom_line(aes(x = Value, y = HBCM, color = "HBCM")) +
  scale_color_manual(name = "Model", values = c("HBCM" = "#00aedb", "Spectral" = "#f37735")) +
  theme_bw() +
  theme(
    panel.grid.major = element_blank()
    # , legend.position = ""
  ) +
  xlab(expression(paste("Mis-specification of Within Cluster Heterogeneous Variance"))) +
  ylab("Adjusted Rand Index") 

plot_hsigma
ggsave(paste0(output_dir, "figure_robust.eps"), dpi = 600, width = 6, height = 4)
ggsave(paste0(output_dir, "figure_robust.png"), dpi = 600, width = 6, height = 4)
