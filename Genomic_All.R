library(tidyverse)
library(hablar)
library(famSKATRC)
library(SKAT)
library(purrr)
Raw_Annotation <- read_rds(file = "Wilberding_BPD_noannot_chr_hg38_cardio.rds")
bed <- read.table(file = "gene_list_50kbpad.bed", as.is = TRUE, sep = '\t')

Genotypes <- Raw_Annotation[["gt"]]
Fix <- Raw_Annotation[["fix"]] %>%
  select(POS, AF, ChromKey, CHROM) %>%
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

updated_Fix <- addGene(Fix, bed)

Genotypes_new <- left_join(updated_Fix, Genotypes, by = c("ChromKey", "POS")) %>%
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
  select(id, POS, ChromKey, geno, GENE) %>%
  filter(id == 210059) %>%
  distinct(POS, .keep_all = TRUE)

variants_last_not_210059 <- updated_keep_this %>%
  mutate(geno = getGene(gt_GT)) %>%
  select(id, POS, ChromKey, geno, GENE) %>%
  filter(id != 210059)

final_keep_this <- rbind(variants_last_not_210059, variants_last_210059) %>%
  select(-ChromKey) 

snp_info <- final_keep_this %>%
  select(-id, -geno) %>%
  distinct()

#%>%
  pivot_wider(names_from = c(POS), values_from = geno) %>%
  arrange(id)

geno_data <- final_keep_this %>%
  select(-id)

final_geno_data <- apply(geno_data, 2, as.numeric)

pheno_data <- read_csv("pheno_data.csv") %>%
  filter(PROP %in% final_keep_this$id) %>%
  mutate_all( ~ ifelse(is.na(.x), median(.x, na.rm = T), .)) %>%
  rename(id = PROP) %>%
  as.data.frame() %>%
  mutate(Race = as.factor(Race), Sex = as.factor(Sex), Ethnicity = as.factor(Ethnicity))

#write csv for kinships
# id_nums <- pheno_data %>%
#   select(id)
# write_csv(id_nums, "id_nums.csv")
# 
# kinship_data <- read_csv("kinships.csv") %>%
#   remove_rownames %>% 
#   column_to_rownames(var="...1") %>%
#   as.matrix()
# 
SKAT_results <- runSKAT(final_geno_data, snp_info)
 SKAT_results <- SKAT_results %>% as.data.frame()
 names(SKAT_results) <- c("Gene", "TAPSE36 p-value", "TAPSE36 n marker test", "PAAT36 p-value", "PAAT36 n marker test", "PAATRVET36 p-value", "PAATRVET n marker test")

geno_data_num <- apply(geno_data, 2, as.numeric)
af <- calc_af(geno_data_num, 0.05)

combined_pheno_geno_data <- cbind(pheno_data, af$geno_data_filt) %>% as.data.frame(stringsAsFactors=FALSE)

pheno_names <- c('TAPSE_36', 'PAAT_36', 'PAAT_RVET_36')
out_table <- c()
for (pheno in pheno_names) {
  combined_pheno_geno_data <- combined_pheno_geno_data %>%
    mutate_(phenotype = pheno)
  
  if(pheno=='TAPSE_36') {
    covariate = 'BW'
  }  else {
    covariate = 'Sex'
  }
  combined_pheno_geno_data <- combined_pheno_geno_data %>%
    mutate_(cov = covariate)
  
  lm_1 <- combined_pheno_geno_data %>% 
    select(-id, -TAPSE_36, -PAAT_36, -PAAT_RVET_36, -BW, -Sex) %>%
    map(~lm(combined_pheno_geno_data$phenotype ~ .x + cov, data = combined_pheno_geno_data)) %>% 
    map(summary)
  
  for (var in names(lm_1)) {
    out_table = rbind(out_table, c(pheno,var,lm_1[[var]]$coefficients[2,4]))
  
  }
}

out_table <- out_table %>% as.data.frame() 
names(out_table) <- c("pheno", "var", "pval")
out_table <- out_table %>% mutate(pval = as.numeric(pval)) %>% filter(var != 'GA') %>% filter(var != 'Race') %>% filter(var != 'Ethnicity') %>% filter(var != 'phenotype') %>% filter(var != 'cov') %>% rename(POS = var) %>% mutate(POS = as.integer(POS))

pos_chrom <- updated_keep_this %>% select(POS, CHROM, GENE) %>% unique() %>% filter(POS %in% out_table$POS)
out_table <- left_join(out_table, pos_chrom, by = 'POS')

out_table <- out_table %>% 
  arrange(GENE, POS) %>%
  mutate(log_pval = -log10(pval), POS_2 = row_number())


p <- ggplot(data = out_table, aes(x = POS_2, y = log_pval)) +
  geom_point(aes(color = GENE))

p = p + facet_grid(.~pheno)

print(p)
ggplot(data = out_table %>% filter(pheno == 'PAAT_36'), mapping = aes(x = POS_2, y = log_pval)) +
  geom_point(aes(color = GENE))

ggplot(data = out_table %>% filter(pheno == 'PAAT_RVET_36'), mapping = aes(x = POS_2, y = log_pval)) +
  geom_point(aes(color = GENE))

#Covariate Analysis - REMOVE ETHNICITY

# BW should be kept for TAPSE
l_TAPSE = lm(TAPSE_36 ~ GA + BW + Race + Sex, data = pheno_data)
l_TAPSE_2 = lm(TAPSE_36 ~ BW, data = pheno_data)

# Sex should be kept for PAAT
l_PAAT = lm(PAAT_36 ~ GA + BW + Race + Sex, data = pheno_data)
l_PAAT_2 = lm(PAAT_36 ~ Sex, data = pheno_data)

# Sex should be kept for PAAT RVET
l_PAAT_RVET = lm(PAAT_RVET_36 ~ GA + BW + Race + Sex, data = pheno_data)
l_PAAT_RVET_2 = lm(PAAT_RVET_36 ~ Sex, data = pheno_data)

# #SKAT
# TAPSE_36_data <- as.matrix(pheno_data%>%select(TAPSE_36))
# covariates <- as.matrix(pheno_data%>%select(BW, Sex))
# adjusted_final_geno_data <- apply(final_geno_data, 2, as.integer)
# adjusted_kinship_data <- kinship_data / 2
# SKAT_TAPSE_36_data <- list(TAPSE_36_data, covariates, adjusted_final_geno_data, adjusted_kinship_data)
# names(SKAT_TAPSE_36_data) <- c('y', 'X', 'Z', 'K')
# SKAT_TAPSE_36_test <- SKAT_NULL_emmaX(y ~ X, K=adjusted_kinship_data, data=SKAT_TAPSE_36_data)
# SKAT(Z, SKAT_TAPSE_36_test)$p.value
# 
# PAAT_36_data <- as.matrix(pheno_data%>%select(PAAT_36))
# SKAT_PAAT_36_data <- list(PAAT_36_data, covariates, adjusted_final_geno_data, adjusted_kinship_data)
# names(SKAT_PAAT_36_data) <- c('y', 'X', 'Z', 'K')
# SKAT_PAAT_36_test <- SKAT_NULL_emmaX(y ~ X, K=adjusted_kinship_data, data=SKAT_PAAT_36_data)
# SKAT(Z, SKAT_PAAT_36_test)$p.value
# 
# PAAT_RVET_36_data <- as.matrix(pheno_data%>%select(PAAT_RVET_36))
# SKAT_PAAT_RVET_36_data <- list(PAAT_RVET_36_data, covariates, adjusted_final_geno_data, adjusted_kinship_data)
# names(SKAT_PAAT_RVET_36_data) <- c('y', 'X', 'Z', 'K')
# SKAT_PAAT_RVET_36_test <- SKAT_NULL_emmaX(y ~ X, K=adjusted_kinship_data, data=SKAT_PAAT_RVET_36_data)
# SKAT(Z, SKAT_PAAT_RVET_36_test)$p.value
# b <- SKAT(Z, SKAT_PAAT_RVET_36_test)
# #test stuff, ignore
# 
# data("SKAT.fam.example")
# K = SKAT.fam.example$K
# Z = SKAT.fam.example$Z
# sample_data_traits <- SKAT.fam.example$X
# obj<-SKAT_NULL_emmaX(y ~ X, K=K, data=SKAT.fam.example)
# SKAT(Z, obj)$p.value
# 
# library(kinship2)
# sample.ped.geno <- process_data()
# KIN = kinship(sample.ped.geno$IID, sample.ped.geno$FA, sample.ped.geno$MO)
# IID = sample.ped.geno$IID
# wuweights_r <- function(maf) ifelse(maf>0, dbeta(maf,1,25), 0)
# wuweights_c <- function(maf) ifelse(maf>0, dbeta(maf,0.5,0.5), 0)
# 
# sample.ped.geno[,"Phenotype"]
# P_VALUES <- famSKAT_RC(PHENO=sample.ped.geno[,"Phenotype"],genotypes=as.matrix(
#   sample.ped.geno[,7:ncol(sample.ped.geno)]), binomialimpute=TRUE,
#   id=IID,fullkins=KIN,maf=0.05, sqrtweights_c=wuweights_c,
#   sqrtweights_r=wuweights_r, phi = c(0,0.2,0.5,0.9))
# print(P_VALUES)
# test_geno <- as.matrix(
#   sample.ped.geno[,7:ncol(sample.ped.geno)])
# write_csv(keep_this, file = "temp_genotypes_table.csv")
# 
# wuweights_r <- function(maf) ifelse(maf>0, dbeta(maf,1,25), 0)
# wuweights_c <- function(maf) ifelse(maf>0, dbeta(maf,0.5,0.5), 0)
# 
# famSKAT_TAPSE <- famSKAT_RC(PHENO = pheno_data[,"TAPSE_36"], genotypes = as.matrix(final_geno_data), binomialimpute=TRUE, id = id_nums$id, fullkins = kinship_data, maf = 0.05, sqrtweights_c=wuweights_c, sqrtweights_r=wuweights_r, phi = 1)

getGene <- function(x){
  
  tranTable <- rbind( c("0/0", 0), c("0/1", 1), c("1/1", 2))
  geno <- gsub(":.+", "", x)
  for (i in 1:nrow(tranTable)){
    
    ind <- which(geno == tranTable[i,1])
    geno[ind] <- as.numeric(tranTable[i,2])
  }
  return(geno)
}

addGene <- function(data, bed){
  data$GENE = ""
  names(bed) = c("CHROM", "start", "end", "gene")
  for (i in 1:nrow(bed)) {
    chr = bed$CHROM[i]
    start = bed$start[i]
    end = bed$end[i]
    gene = bed$gene[i]
    data <- data %>%
      mutate(GENE = ifelse(CHROM==chr & POS>=start & POS<= end, gene, GENE))
  }
  return(data)
}

runSKAT <- function(geno, snp_info){
  genes = unique(snp_info$GENE)
  results_table <- c("", "", "", "", "", "", "")
  names(results_table) = c("Gene", "TAPSE36 p-value", "TAPSE36 n marker test", "PAAT36 p-value", "PAAT36 n marker test", "PAATRVET36 p-value", "PAATRVET n marker test")
  for (gene in genes) {
    snps = as.character(snp_info$POS[snp_info$GENE==gene])
    geno_gene <- as.data.frame(geno) %>%
      select(snps)
    final_geno_gene <- apply(geno_gene, 2, as.integer)
    covariates_temp <- as.data.frame(covariates)
    covariates_TAPSE <- matrix(covariates_temp$BW, nrow(covariates_temp), 1)
    covariates_PAAT <- matrix(covariates_temp$Sex, nrow(covariates_temp), 1)
    SKAT_TAPSE_36_data <- list(TAPSE_36_data, covariates_TAPSE, final_geno_gene, adjusted_kinship_data)
    names(SKAT_TAPSE_36_data) <- c('y', 'X', 'Z', 'K')
    SKAT_TAPSE_36_test <- SKAT_NULL_emmaX(y ~ X, K=adjusted_kinship_data, data=SKAT_TAPSE_36_data)
    TAPSE_36 <- SKAT(SKAT_TAPSE_36_data$Z, SKAT_TAPSE_36_test)
    TAPSE36p <- TAPSE_36$p.value
    TAPSE36n <- TAPSE_36$param$n.marker.test
    
    SKAT_PAAT_36_data <- list(PAAT_36_data, covariates_PAAT, final_geno_gene, adjusted_kinship_data)
    names(SKAT_PAAT_36_data) <- c('y', 'X', 'Z', 'K')
    SKAT_PAAT_36_test <- SKAT_NULL_emmaX(y ~ X, K=adjusted_kinship_data, data=SKAT_PAAT_36_data)
    PAAT_36 <- SKAT(SKAT_PAAT_36_data$Z, SKAT_PAAT_36_test)
    PAAT36p <- PAAT_36$p.value
    PAAT36n <- PAAT_36$param$n.marker.test
    
    SKAT_PAAT_RVET_36_data <- list(PAAT_RVET_36_data, covariates_PAAT, final_geno_gene, adjusted_kinship_data)
    names(SKAT_PAAT_RVET_36_data) <- c('y', 'X', 'Z', 'K')
    SKAT_PAAT_RVET_36_test <- SKAT_NULL_emmaX(y ~ X, K=adjusted_kinship_data, data=SKAT_PAAT_RVET_36_data)
    PAAT_RVET_36 <- SKAT(SKAT_PAAT_RVET_36_data$Z, SKAT_PAAT_RVET_36_test)
    PAAT_RVET36p <- PAAT_RVET_36$p.value
    PAAT_RVET36n <- PAAT_RVET_36$param$n.marker.test
    
    results_table <- rbind(results_table, c(gene, TAPSE36p, TAPSE36n, PAAT36p, PAAT36n, PAAT_RVET36p, PAAT_RVET36n, ncol(final_geno_gene)))
  }
  return(results_table)
}

calc_af <- function(geno_data, af_threshold) {
  af = colSums(geno_data_num, na.rm = T) / (2 * colSums(is.na(geno_data_num) == FALSE))
  af <- data.frame(snp = colnames(geno_data_num), af = af)
  #ind = max(af, )
  af_filt = af %>% rowwise() %>% filter(min(af, 1-af) >= af_threshold)
  geno_data_filt = geno_data[,colnames(geno_data)%in%af_filt$snp]
  #geno_data <- geno_data %>% as.matrix()
  #snp_names = colnames(geno_data)[af > af_threshold]
  return(list(af = af, geno_data_filt = geno_data_filt))
}
