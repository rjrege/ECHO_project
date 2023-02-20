library(tidyverse)
Raw_Annotation <- read_rds(file = "Wilberding_BPD_noannot_chr_hg38_hgmd.rds")
Genotypes <- Raw_Annotation[["gt"]][["gt_GT"]]
Allele_Frequency <- Raw_Annotation[["fix"]][["AF"]]
