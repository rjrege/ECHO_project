library(tidyverse)
WASHU_RawSet <- read_csv(file = "WASHU_all.csv")
WASHU_ALL_Clean <- read_csv("WASHU_all.csv", 
col_types = cols(`Shenan (0=No BPD, 1= BPD, 4=dead, transfer, withdrew)` = col_factor(levels = c("0", 
"1", "4")), `2004 (0=No BPD, 1=BPD)` = col_factor(levels = c("0", 
"1")), MODIFEID = col_factor(levels = c("0", 
"1")), `NIH (0=No, 1=mild, 2= moderate, 3=severe` = col_factor(levels = c("0", 
"1", "2", "3")), PH_32 = col_factor(levels = c("0", 
"1")), PH_36 = col_factor(levels = c("0", 
"1")), PH_One = col_factor(levels = c("0", 
"1")), `WEIGHT AT 32 WEEKS` = col_double(), 
Weight...21 = col_double(), Weight...30 = col_double()))

WASHU_tidied <- write_csv(file = "WASHU_tidied.csv")
WASHU_tidied <- write_csv(WASHU_ALL_CLEAN, file = "WASHU_tidied.csv")
WASHU_tidied <- write_csv(WASHU_ALL_Clean, file = "WASHU_tidied.csv")
