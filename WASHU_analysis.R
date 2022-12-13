
library(tidyverse)

WASHU_all_weeks <- read_csv(file = "WASHU_all_Weeks.csv")

## PAAT RVET 32

paat_rvet_32 <- WASHU_all_weeks %>%
  select(`Shenan (0=No BPD, 1= BPD, 4=dead, transfer, withdrew)`, PAAT_RVET_32) %>%
  filter(!is.na(PAAT_RVET_32)) %>%
  filter(!(`Shenan (0=No BPD, 1= BPD, 4=dead, transfer, withdrew)` == 4))


paat_rvet_32 %>%
  group_by(`Shenan (0=No BPD, 1= BPD, 4=dead, transfer, withdrew)`) %>%
  summarize(mean = mean(PAAT_RVET_32))

paat_rvet_32 %>%
  table()

paat_rvet_32_aov <- aov(PAAT_RVET_32 ~ `Shenan (0=No BPD, 1= BPD, 4=dead, transfer, withdrew)`, data = paat_rvet_32)

summary(paat_rvet_32_aov)  

## PAAT RVET 36

paat_rvet_36 <- WASHU_all_weeks %>%
  select(`Shenan (0=No BPD, 1= BPD, 4=dead, transfer, withdrew)`, PAAT_RVET_36) %>%
  filter(!is.na(PAAT_RVET_36)) %>%
  filter(!(`Shenan (0=No BPD, 1= BPD, 4=dead, transfer, withdrew)` == 4))


paat_rvet_36 %>%
  group_by(`Shenan (0=No BPD, 1= BPD, 4=dead, transfer, withdrew)`) %>%
  summarize(mean = mean(PAAT_RVET_36))

paat_rvet_36 %>%
  table()

paat_rvet_36_aov <- aov(PAAT_RVET_36 ~ `Shenan (0=No BPD, 1= BPD, 4=dead, transfer, withdrew)`, data = paat_rvet_36)

summary(paat_rvet_36_aov)  
#wilberdingtrain
#again
#one more for the boys


