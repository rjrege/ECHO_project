library(tidyverse)
library(hablar)
library(FREGAT)

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
  select(-id) %>%
  as.snp.data()

pheno_data <- read_csv("pheno_data.csv") %>%
  filter(PROP %in% final_keep_this$id) %>%
  mutate_all( ~ ifelse(is.na(.x), median(.x, na.rm = T), .)) %>%
  rename(id = PROP)

#write csv for kinships
id_nums <- pheno_data %>%
  select(PROP)
write_csv(id_nums, "id_nums.csv")

kinship_data <- read_csv("kinships.csv") %>%
  remove_rownames %>% 
  column_to_rownames(var="...1") %>%
  as.matrix()

pheno_data_TAPSE_36 <- pheno_data %>%
  select(-PAAT_36, -PAAT_RVET_36) %>%
  as.data.frame()

out_TAPSE_36 <- FFBSKAT(formula = TAPSE_36 ~ GA + BW + Race + Ethnicity + Sex, phenodata = pheno_data_TAPSE_36, genodata = final_keep_this, kin = kinship_data)
out <- FFBSKAT(trait ~ age + sex, phenodata, genodata, kin)
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
