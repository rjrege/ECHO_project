library(tidyverse)
library(hablar)
library(famSKATRC)

Raw_Annotation <- read_rds(file = "Wilberding_BPD_noannot_chr_hg38_cardio.rds")

Genotypes <- Raw_Annotation[["gt"]]
Fix <- Raw_Annotation[["fix"]] %>%
  select(POS, AF, ChromKey) %>%
  group_by(POS, ChromKey) %>%
  mutate(dupe = max(row_number()) > 1) %>%
  ungroup() %>%
  filter(dupe == FALSE) %>%
  distinct() %>%
  select(-dupe) %>%
  convert(num(AF)) %>%
  mutate(homo_alt = AF*AF) %>%
  mutate(homo_ref = (1-AF)*(1-AF)) %>%
  mutate(hetero = 1-homo_alt-homo_ref) %>%
  #select(-AF) %>%
  arrange(POS)

Genotypes_new <- left_join(Fix, Genotypes, by = c("ChromKey", "POS")) %>%
  filter(is.na(gt_GT) | gt_GT != '0/.') %>% 
  filter(is.na(gt_GT) | gt_GT != '1/.') %>%
  filter(is.na(gt_GT) | gt_GT != './1')


keep_this <- Genotypes_new %>%
  mutate(id = substr(Indiv, 1, 6)) %>%
  select(id, everything()) %>%
  select(-Indiv)
  #filter(!is.na(gt_GT)) %>%
  #select(-gt_GT_alleles) 

# list of POS that have an NA less than 20%
removed_POS_list <- keep_this %>%
  group_by(POS, ChromKey) %>%
  summarize(prop = sum(is.na(gt_GT)) / n()) %>%
  arrange(prop) %>%
  filter(prop <= 0.2) %>%
  select(POS, ChromKey)
#takes out NA POS greater than 20%
updated_keep_this <- keep_this %>%
  filter(POS %in% removed_POS_list$POS) %>%
  mutate(new_geno = 12345)

# imputes the NA values
for (i in 1:length(updated_keep_this$gt_GT)) {
  
  if(is.na(updated_keep_this$gt_GT[i])) {
    updated_keep_this$new_geno[i] <- sample(c(0, 1, 2), size = 1, replace = TRUE, prob = c(updated_keep_this$homo_ref[i], updated_keep_this$hetero[i], updated_keep_this$homo_alt[i]))
    updated_keep_this$gt_GT[i] <- updated_keep_this$new_geno[i]
  }
}

# Why are there still POS with insane amounts of entries? Also what are the duplicates in POS?
variants_last_210059 <- updated_keep_this %>%
  mutate(geno = getGene(gt_GT)) %>%
  select(id, POS, ChromKey, geno) %>%
  filter(id == 210059) %>%
  distinct(POS, .keep_all = TRUE)

variants_last_not_210059 <- updated_keep_this %>%
  mutate(geno = getGene(gt_GT)) %>%
  select(id, POS, ChromKey, geno) %>%
  filter(id != 210059)

final_keep_this <- rbind(variants_last_not_210059, variants_last_210059) %>%
  select(-ChromKey) %>%
  pivot_wider(names_from = c(POS), values_from = geno) %>%
  arrange(id)

geno_data <- final_keep_this %>%
  select(-id)

final_geno_data <- apply(geno_data, 2, as.numeric)

pheno_data <- read_csv("pheno_data.csv") %>%
  filter(PROP %in% final_keep_this$id) %>%
  mutate_all( ~ ifelse(is.na(.x), median(.x, na.rm = T), .)) %>%
  rename(id = PROP) %>%
  as.data.frame()

#write csv for kinships
id_nums <- pheno_data %>%
  select(id)
write_csv(id_nums, "id_nums.csv")

kinship_data <- read_csv("kinships.csv") %>%
  remove_rownames %>% 
  column_to_rownames(var="...1") %>%
  as.matrix()

wuweights_r <- function(maf) ifelse(maf>0, dbeta(maf,1,25), 0)
wuweights_c <- function(maf) ifelse(maf>0, dbeta(maf,0.5,0.5), 0)

famSKAT_TAPSE <- famSKAT_RC(PHENO = pheno_data[,"TAPSE_36"], genotypes = as.matrix(final_geno_data), binomialimpute=TRUE, id = id_nums$id, fullkins = kinship_data, maf = 0.05, sqrtweights_c=wuweights_c, sqrtweights_r=wuweights_r, phi = 1)

#SKAT
TAPSE_36_data <- as.matrix(pheno_data%>%select(TAPSE_36))
covariates <- as.matrix(pheno_data%>%select(GA:Sex))
SKAT_TAPSE_36_data <- list(TAPSE_36_data, covariates, final_geno_data, kinship_data)
names(SKAT_TAPSE_36_data) <- c('pheno', 'covars', 'geno', 'kin')
SKAT_TAPSE_36_test <- SKAT_NULL_emmaX(pheno ~ covars, K=kin, data=SKAT_TAPSE_36_data)

#test stuff, ignore

data("SKAT.fam.example")
K = SKAT.fam.example$K
Z = SKAT.fam.example$Z
sample_data_traits <- SKAT.fam.example$X
obj<-SKAT_NULL_emmaX(y ~ X, K=K, data=SKAT.fam.example)
SKAT(Z, obj)$p.value

library(kinship2)
sample.ped.geno <- process_data()
KIN = kinship(sample.ped.geno$IID, sample.ped.geno$FA, sample.ped.geno$MO)
IID = sample.ped.geno$IID
wuweights_r <- function(maf) ifelse(maf>0, dbeta(maf,1,25), 0)
wuweights_c <- function(maf) ifelse(maf>0, dbeta(maf,0.5,0.5), 0)

sample.ped.geno[,"Phenotype"]
P_VALUES <- famSKAT_RC(PHENO=sample.ped.geno[,"Phenotype"],genotypes=as.matrix(
  sample.ped.geno[,7:ncol(sample.ped.geno)]), binomialimpute=TRUE,
  id=IID,fullkins=KIN,maf=0.05, sqrtweights_c=wuweights_c,
  sqrtweights_r=wuweights_r, phi = c(0,0.2,0.5,0.9))
print(P_VALUES)
test_geno <- as.matrix(
  sample.ped.geno[,7:ncol(sample.ped.geno)])
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
