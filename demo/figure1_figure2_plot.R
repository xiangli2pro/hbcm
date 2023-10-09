#################################################################################################
## Manuscript: Community Detection with Heterogeneous Block Covariance Model
## Plot code for manuscript Figure 1 and Figure 2
#################################################################################################


# ------------ Install packages
library(devtools)
devtools::install_github("xiangli2pro/hbcm")
library(hbcm)
library(tidyverse)
library(ggpubr)
library(scales)
library(ggforce)
library(patchwork)


# ------------ Define functions

# function to read data and make a table of AdjRand of spectral clustering and HBCM
tb_func <- function(dataset, values, len) {
  
  # use lapply to process in parallel
  ls <- lapply(
    # list of file name
    paste0(dataset, values, ".rda"),
    # function to apply on each file
    function(f) {
      # load data
      load(f)
      # calculate mean adjRand
      mean <- c(
        # mean of spectral cluster
        mean(sapply(spec_labels, function(x) matchLabel(labels, x)$adjRand)),
        # mean of HBCM
        mean(sapply(hbcm_res, function(x) {
          # if condition because HBCM returns NaN for some simulated data
          if (length(x) == len) {
            matchLabel(labels, x$cluster)$adjRand
          } else {
            NA
          }
        }), na.rm = TRUE)
      )
      # calculate standard deviation of adjRand
      sd <- c(
        # sd of spectral cluster
        sd(sapply(spec_labels, function(x) matchLabel(labels, x)$adjRand)),
        # sd of HBCM
        sd(sapply(hbcm_res, function(x) {
          if (length(x) == len) {
            matchLabel(labels, x$cluster)$adjRand
          } else {
            NA
          }
        }), na.rm = TRUE)
      )
      # return mean and sd for each file
      return(list(
        mean = mean,
        sd = sd
      ))
    }
  )

  # combine all the results
  cbind(
    values,
    sapply(ls, function(x) x$mean) %>% round(2) %>% t()
  ) %>%
    as.data.frame() %>%
    `colnames<-`(c("Value", "Spectral", "HBCM"))
}


# function to read data and make a table of AdjRand of spectral clustering and HBCM
# use a customized labels
tb_func_mislead <- function(dataset, values, len) {
  
  # for manuscript Figure 2 plot
  # assume the mislead label is
  misLead_labels <- c(rep(1, 330), rep(2, 330), rep(3, 340))
  
  # calculate the adjRand using the mislead label
  ls <- lapply(
    paste0(dataset, values, ".rda"),
    function(f) {
      load(f)
      mean <- c(
        mean(sapply(spec_labels, function(x) matchLabel(misLead_labels, x)$adjRand)),
        mean(sapply(hbcm_res, function(x) {
          if (length(x) == len) {
            matchLabel(misLead_labels, x$cluster)$adjRand
          } else {
            NA
          }
        }), na.rm = TRUE)
      )
      sd <- c(
        sd(sapply(spec_labels, function(x) matchLabel(misLead_labels, x)$adjRand)),
        sd(sapply(hbcm_res, function(x) {
          if (length(x) == len) {
            matchLabel(misLead_labels, x$cluster)$adjRand
          } else {
            NA
          }
        }), na.rm = TRUE)
      )
      return(list(
        mean = mean,
        sd = sd
      ))
    }
  )

  cbind(
    values,
    sapply(ls, function(x) x$mean) %>% round(2) %>% t()
  ) %>%
    as.data.frame() %>%
    `colnames<-`(c("Value", "Spectral", "HBCM"))
}


# ------------ Make plots of manuscript Figure 1

# make manuscript Figure 1 (left)

# read data
res_prefix = "./simulation_data/res_hlambda1_hsigma"
res_hsigma <- tb_func(res_prefix, c(2:10), len = 7)
# make plot
plot_hsigma <- res_hsigma %>%
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

# read data
res_prefix = "./simulation_data/res_omega_"
res_omega <- tb_func(res_prefix, seq(0.1, 0.9, by = 0.1), len = 7)
# make plot
plot_omega <- res_omega %>%
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

# read data
res_prefix = "./simulation_data/res_t_"
res_prefix_tc = "./simulation_data/res_t_correct_"
res_t <- tb_func(res_prefix, c(3:7), len = 7)
res_tcorrect <- tb_func(res_prefix_tc, c(3:7), len = 7)
# make plot
plot_t <- cbind(res_t, res_tcorrect[, -1]) %>%
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

ggsave("plots/simParameter.eps", dpi = 600, width = 8, height = 4)


# ------------ Make plots of manuscript Figure 2

# read data
res_prefix = "./simulation_data2/res_hsigma"
res_hlambda_margin1 <- tb_func(res_prefix, c(10:15), len = 7)
res_hlambda_margin2 <- tb_func_mislead(res_prefix, c(10:15), len = 7)
# make plot
plot_t <- cbind(c(10:15), res_hlambda_margin1[, 2:3], res_hlambda_margin2[2:3]) %>%
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
ggsave("lambdaMargin.eps", dpi = 300)
