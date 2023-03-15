library(tidyverse)
library(FREGAT)

WASHU_all_weeks <- read_csv(file = "WASHU_all_Weeks.csv")
Raw_Annotation <- read_rds(file = "Wilberding_BPD_noannot_chr_hg38_cardio.rds")
bed <- read.table("gene_list_50kbpad.bed",header = FALSE, sep="\t",stringsAsFactors=FALSE, quote="")

Genotypes <- Raw_Annotation[["gt"]]
Fix <- Raw_Annotation[["fix"]] %>%
  filter(POS == 237456575)

%>%
  select(ChromKey, CHROM) %>%
  unique()
Genotypes <- left_join(Genotypes, Fix, by = "ChromKey") 

Genotypes <- Genotypes %>%
  select(CHROM, -ChromKey, everything())
table(Genotypes$gt_GT)
#dots might represent insertions/deletions? exclude for now
keep_this <- Genotypes %>%
  #filter(gt_GT == '0/.' | gt_GT == './1' | gt_GT == '1/.') %>%
  filter(POS == 237456575) 

%>%
  filter(CHROM == 'chr1')
  filter(gt_GT == '0/0' | gt_GT == '0/1' | gt_GT == '1/1') %>%
  mutate(geno = getGene(gt_GT)) %>%
  mutate(id = substr(Indiv, 1, 6)) %>%
  select(id, geno, POS, CHROM) %>%
  group_by(POS, CHROM) %>%
  summarize(num = n()) %>%
  arrange(num)

table(Genotypes$gt_GT)

#test stuff, ignore
data(example.data) 
View(genodata)
View(phenodata)
View(kin)

write_csv(keep_this, file = "temp_genotypes_table.csv")


getGene <- function(x){
  
  tranTable <- rbind( c("0/0", 0), c("0/1", 1), c("1/1", 2))
  geno <- gsub(":.+", "", x)
  for (i in 1:nrow(tranTable)){
    
    ind <- which(geno == tranTable[i,1])
    geno[ind] <- as.numeric(tranTable[i,2])
  }
  return(geno)
}
