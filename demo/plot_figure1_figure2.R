#################################################################################################
## Manuscript: Community Detection with Heterogeneous Block Covariance Model
## Plot code for manuscript Figure 1 and Figure 2
#################################################################################################


#************************** Install packages
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

#**************************  Read simulation output
df_f1_left = data.frame(
  Value = c(2:10),
  Spectral = c(1,1,1,0.99,0.91,0.68,0.33,0.06,0.01),
  HBCM = c(1,1,1,0.99,0.95,0.86,0.72,0.48,0.17)
)

df_f1_middle = data.frame(
  Value = seq(0.1, 0.9, 0.1),
  Spectral = c(0.98, 0.98, 0.97, 0.94, 0.91, 0.82, 0.67, 0.37, 0.02),
  HBCM = c(1, 0.99, 0.99, 0.98, 0.95, 0.89, 0.77, 0.47, 0.03)
)

df_f1_right = data.frame(
  Value = seq(3, 10, 1),
  Spectral = c(0, 0.16, 0.43, 0.56, 0.64, 0.70, 0.74, 0.76),
  HBCM = c(0.13, 0.62, 0.76, 0.81, 0.85, 0.87, 0.88, 0.89)
)

df_f1_right_corrected = data.frame(
  Value = seq(3, 10, 1),
  Spectral = c(0.91, 0.90, 0.90, 0.90, 0.90, 0.91, 0.91, 0.91),
  HBCM = c(0.95, 0.95, 0.95, 0.95, 0.95, 0.95, 0.95, 0.95)
)

df_f2 = data.frame(
  Value = c(10:15),
  Spectral = c(0.89, 0.85, 0.81, 0.51, 0.43, 0.32),
  HBCM = c(0.90, 0.87, 0.83, 0.80, 0.74, 0.45),
  Sepctral_Mis = c(0.01, 0, 0.01, 0.18, 0.27, 0.47),
  HBCM_Mis = c(0, 0, 0, 0, 0, 0)
)


#**************************  Make plots of manuscript Figure 1

# make manuscript Figure 1 (left)
plot_hsigma <- df_f1_left %>%
  as.data.frame() %>%
  `colnames<-`(c("Value", "Spectral", "HBCM")) %>%
  ggplot() +
  geom_line(aes(x = Value, y = Spectral, color = "Spectral")) +
  geom_line(aes(x = Value, y = HBCM, color = "HBCM")) +
  scale_color_manual(name = "Model", values = c("HBCM" = "#00aedb", "Spectral" = "#f37735")) +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    legend.position = ""
  ) +
  xlab(expression(paste("Heterogeneous ", sigma))) +
  ylab("Adjusted Rand Index") +
  annotate("text",
           x = 2.5, y = c(0.6, 0.55, 0.50),
           label = c(
             expression(paste(n, "=", 1000, " ", p, "=", 1000)),
             expression(paste(lambda, "=", 1)),
             expression(paste(K, "=", 3))
           ),
           size = 3, hjust = 0
  )


# make manuscript Figure 1 (middle)
plot_omega <- df_f1_middle %>%
  as.data.frame() %>%
  `colnames<-`(c("Value", "Spectral", "HBCM")) %>%
  ggplot() +
  geom_line(aes(x = Value, y = Spectral, color = "Spectral")) +
  geom_line(aes(x = Value, y = HBCM, color = "HBCM")) +
  scale_color_manual(name = "Model", values = c("HBCM" = "#00aedb", "Spectral" = "#f37735")) +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    legend.position = ""
  ) +
  xlab(expression(paste(Omega, " off-diagonal"))) +
  ylab("") +
  annotate("text",
           x = 0.2, y = c(0.6, 0.55, 0.50),
           label = c(
             expression(paste(n, "=", 1000, " ", p, "=", 1000)),
             expression(paste(lambda, "=", 1, " ", sigma, "=", 6, )),
             expression(paste(K, "=", 3))
           ),
           size = 3, hjust = 0
  ) +
  scale_x_continuous(breaks = seq(0.1, 0.9, by = 0.2))


# make manuscript Figure 1 (right)
# make plot
plot_t <- cbind(df_f1_right, df_f1_right_corrected[, -1]) %>%
  as.data.frame() %>%
  `colnames<-`(c("Value", "Spectral", "HBCM", "SpectralC", "HBCMC")) %>%
  ggplot() +
  geom_line(aes(x = Value, y = Spectral, color = "Spectral", linetype = "t")) +
  geom_line(aes(x = Value, y = HBCM, color = "HBCM", linetype = "t")) +
  geom_line(aes(x = Value, y = SpectralC, color = "Spectral", linetype = "t-corrected")) +
  geom_line(aes(x = Value, y = HBCMC, color = "HBCM", linetype = "t-corrected")) +
  scale_color_manual(
    name = "Model",
    values = c("HBCM" = "#00aedb", "Spectral" = "#f37735")
  ) +
  scale_linetype_manual(
    name = "Distribution",
    values = c("t" = 1, "t-corrected" = 2)
  ) +
  theme_bw() +
  theme(panel.grid.major = element_blank()) +
  xlab(expression(paste("df of t-distribution"))) +
  ylab("") +
  ylim(c(0, 1)) +
  annotate("text",
           x = 4.5, y = c(0.6, 0.55, 0.50) - 0.4,
           label = c(
             expression(paste(n, "=", 1000, " ", p, "=", 1000)),
             expression(paste(lambda, "=", 1, " ", sigma, "=", 6, )),
             expression(paste(K, "=", 3))
           ),
           size = 3, hjust = 0
  )


# combine plots to a single plot

ggarrange(plot_hsigma, plot_omega, plot_t,
          ncol = 3, nrow = 1, align = "hv",
          common.legend = TRUE
)

ggsave(paste0(output_dir, "figure1.eps"), dpi = 600, width = 8, height = 4)
ggsave(paste0(output_dir, "figure1.png"), dpi = 600, width = 8, height = 4)


#**************************  Make plots of manuscript Figure 2

plot_f2 <-  df_f2%>%
  as.data.frame() %>%
  `colnames<-`(c("Value", "Spectral", "HBCM", "Spectral_mis", "HBCM_mis")) %>%
  ggplot() +
  geom_line(aes(x = Value, y = Spectral, color = "Spectral", linetype = "True")) +
  geom_line(aes(x = Value, y = HBCM, color = "HBCM", linetype = "True")) +
  geom_line(aes(x = Value, y = Spectral_mis, color = "Spectral", linetype = "Mislead")) +
  geom_line(aes(x = Value, y = HBCM_mis, color = "HBCM", linetype = "Mislead")) +
  scale_color_manual(
    name = "Model",
    values = c("HBCM" = "#00aedb", "Spectral" = "#f37735")
  ) +
  scale_linetype_manual(
    name = "Reference Label",
    values = c("True" = 1, "Mislead" = 2)
  ) +
  theme_bw() +
  theme(panel.grid.major = element_blank()) +
  xlab(expression(paste("Heterogeneous ", sigma))) +
  ylab("Adjusted Rand Index") +
  ylim(c(0, 1)) +
  annotate("text",
           x = 10.5, y = c(0.6, 0.55, 0.50) - 0.1,
           label = c(
             expression(paste(n, "=", 1000, " ", p, "=", 1000)),
             expression(paste(lambda, "= c(1,5,25)")),
             expression(paste(K, "=", 3))
           ),
           size = 3, hjust = 0
  )

# save plots
ggsave(paste0(output_dir, "figure2.eps"), dpi = 600)
ggsave(paste0(output_dir, "figure2.png"), dpi = 600)
