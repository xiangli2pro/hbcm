## Code for manuscript *Community Detection with Heterogeneous Block Covariance Model*


The folder provides comprehensive code scripts to replicate the numerical results and graphs presented in manuscript *Community Detection with Heterogeneous Block Covariance Model*. The readme doc gives instructions to help reader navigate the files.


### Data

The yeast gene expression data sets used in manuscript Figure 3 and Figure 4 are embedded in the `hbcm` package, and can be called by their names (`data_gene_sample`, `data_gene_group`) in R. The data sets can also be downloaded from the [package github webpage](https://github.com/xiangli2pro/hbcm/tree/main/data). Those data sets are processed following the [script](https://github.com/xiangli2pro/hbcm/blob/main/data-raw/yeast_gene_data.R) here. Raw data pre-process is available in [package github webpage](https://github.com/xiangli2pro/hbcm/tree/main/inst/extdata).


Simulation data can be replicated by running the R scripts provided in this directory.

### Script

To replicate manuscript Table 1 result:
- Run `table1_simulation.R` for different combinations of K, N, P.

To replicate manuscript Figure 1 and Figure 2 results:
- Run `figure1_left_simulation.R` for different values of hsigma;
- Run `figure1_middle_simulation.R` for different values of omega;
- Run `figure1_right_tdist_simulation.R` for different degree of freedom (df) for t-distribution;
- Run `figure1_right_tdist_corrected_simulation.R` for different degree of freedom (df) for t-distribution (divided by their standard deviation);
- After getting the results from above scripts, run `figure1_figure2_plot.R` to make plots of Figure 1 and Figure 2.


To replicate manuscript Figure 3 result:
- Run `figure3_plot.R` to generate the cross-validation plots of simulation data and gene expression data.

To replicate manuscript Figure 4 result:
- Run `figure4_plot.R` to generate the heatmap plot of gene expression data clustered by HBCM.






