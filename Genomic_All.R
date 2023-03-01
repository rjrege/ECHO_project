library(tidyverse)
library(FREGAT)

WASHU_all_weeks <- read_csv(file = "WASHU_all_Weeks.csv")
Raw_Annotation <- read_rds(file = "Wilberding_BPD_noannot_chr_hg38_cardio.rds")
bed <- read.table("gene_list_50kbpad.bed",header = FALSE, sep="\t",stringsAsFactors=FALSE, quote="")

Genotypes <- Raw_Annotation[["gt"]]

table(Genotypes$gt_GT)
#dots might represent insertions/deletions? exclude for now
dotted <- Genotypes %>%
  filter(gt_GT == '0/.')

length(unique(dotted$POS))
View(dotted)
head(dotted)
#test stuff, ignore
data(example.data) 
View(genodata)
