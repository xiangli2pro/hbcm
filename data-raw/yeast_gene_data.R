# code to prepare `yeast_gene_data` dataset goes here


# ------------ Data source
## GSE75694_matrix_378_samples.txt: supplementary file downloaded from the bottom of webpage https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE75694
## Table S1.xlsx: supplementary table downloaded from https://www.tandfonline.com/doi/suppl/10.1080/15384101.2019.1570655?scroll=top


# ------------ Data process 

# use load_all() to activate the interactive call of system.file, otherwise it returns empty string
library(devtools)
library(readxl)
library(tidyverse)
load_all()

# add data path
path2sample = system.file("extdata", "GSE75694_matrix_378_samples.txt", package = "hbcm")
path2genegroup = system.file("extdata", "Table_S1.xlsx", package = "hbcm")

# gene expression data 5900x378 from the paper
gene_sample <- read.delim(path2sample)
# pairs of gene, transcript factor (TF) form the paper (241 genes)
gene_group <- read_xlsx(path2genegroup)
# get unique genes
gene_unique <- intersect(gene_group$Probe, gene_sample$ID_REF)

# keep unique genes in gene sample data
gene_sample_uni <- gene_sample %>% 
  filter(gene_sample$ID_REF %in% gene_unique) 
# transpose matrix
gene_sample_uni_t <- gene_sample_uni %>% 
  t() %>% 
  `colnames<-`(gene_sample_uni$ID_REF) %>%
  as.data.frame() %>% 
  slice(-1)

# keep unique genes in gene group data
gene_group_uni <- gene_group %>% 
  filter(gene_group$Probe %in% gene_unique) %>% 
  rename("TF_regulator" = "TF regulator")


# ------------ Save data
data_gene_sample <- gene_sample_uni_t
data_gene_group <- gene_group_uni
usethis::use_data(data_gene_sample, overwrite = TRUE)
usethis::use_data(data_gene_group, overwrite = TRUE)

