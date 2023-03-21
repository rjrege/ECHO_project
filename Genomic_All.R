library(tidyverse)
library(hablar)
library(FREGAT)

WASHU_all_weeks <- read_csv(file = "WASHU_all_Weeks.csv")
Raw_Annotation <- read_rds(file = "Wilberding_BPD_noannot_chr_hg38_cardio.rds")
bed <- read.table("gene_list_50kbpad.bed",header = FALSE, sep="\t",stringsAsFactors=FALSE, quote="")

Genotypes <- Raw_Annotation[["gt"]]
Fix <- Raw_Annotation[["fix"]] %>%
  #select(POS, AF) %>%
  convert(num(AF)) %>%
  mutate(homo_alt = AF*AF) %>%
  mutate(homo_ref = (1-AF)*(1-AF)) %>%
  mutate(hetero = 1-homo_alt-homo_ref) %>%
  #select(-AF) %>%
  arrange(POS)

Genotypes <- left_join(Genotypes, Fix, by = "ChromKey") 

Genotypes <- Genotypes %>%
  select(CHROM, -ChromKey, everything())
table(Genotypes$gt_GT)
#dots might represent insertions/deletions? exclude for now
keep_this <- Genotypes %>%
  filter(is.na(gt_GT) | gt_GT != '0/.') %>% 
  filter(is.na(gt_GT) | gt_GT != '1/.') %>%
  filter(is.na(gt_GT) | gt_GT != './1') %>%
  mutate(id = substr(Indiv, 1, 6)) %>%
  select(id, everything()) %>%
  select(-Indiv) %>%
  #filter(POS == 4697315) %>%
  #filter(!is.na(gt_GT)) %>%
  select(-gt_GT_alleles) %>%
  distinct() 
# list of POS that have an NA greater than 20%
removed_POS_list <- keep_this %>%
  group_by(POS) %>%
  summarize(prop = sum(is.na(gt_GT)) / n()) %>%
  arrange(prop) %>%
  filter(prop <= 0.2) %>%
  select(POS)
#takes out NA POS greater than 20%
updated_keep_this <- keep_this %>%
  filter(POS %in% removed_POS_list$POS) %>%
  left_join(Fix, by = "POS") %>%
  mutate(new_geno = 12345)

# imputes the NA values
for (i in 1:length(updated_keep_this$new_geno)) {
  updated_keep_this$new_geno[i] <- sample(c(0, 1, 2), size = 1, replace = TRUE, prob = c(updated_keep_this$homo_ref[i], updated_keep_this$hetero[i], updated_keep_this$homo_alt[i]))
  
  if(is.na(updated_keep_this$gt_GT[i])) {
    updated_keep_this$gt_GT[i] <- updated_keep_this$new_geno[i]
  }
}

# Why are there still POS with insane amounts of entries? Also what are the duplicates in POS?
final_keep_this <- updated_keep_this %>%
  mutate(geno = getGene(gt_GT)) %>%
  select(-gt_GT, -new_geno) %>%
  filter(POS == 55463977)


%>%
  group_by(POS) %>%
  summarize(num = n()) %>%
  arrange(num)

colMeans(is.na(keep_this))
  
  #filter(POS == 237456575) 

#%>%
 # filter(CHROM == 'chr1')
 # filter(gt_GT == '0/0' | gt_GT == '0/1' | gt_GT == '1/1') %>%
  #mutate(geno = getGene(gt_GT)) %>%
  #mutate(id = substr(Indiv, 1, 6)) %>%
  #select(geno, everything()) %>%
  #select(id, geno, POS, CHROM) %>%
  group_by(gt_GT) %>%
  summarize(num = n())
  #arrange(num)

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
