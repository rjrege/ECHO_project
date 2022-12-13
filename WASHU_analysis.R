
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

## PAAT RVET One

paat_rvet_One <- WASHU_all_weeks %>%
  select(`Shenan (0=No BPD, 1= BPD, 4=dead, transfer, withdrew)`, PAAT_RVET_One) %>%
  filter(!is.na(PAAT_RVET_One)) %>%
  filter(!(`Shenan (0=No BPD, 1= BPD, 4=dead, transfer, withdrew)` == 4))


paat_rvet_One %>%
  group_by(`Shenan (0=No BPD, 1= BPD, 4=dead, transfer, withdrew)`) %>%
  summarize(mean = mean(PAAT_RVET_One))

paat_rvet_One %>%
  table()

paat_rvet_One_aov <- aov(PAAT_RVET_One ~ `Shenan (0=No BPD, 1= BPD, 4=dead, transfer, withdrew)`, data = paat_rvet_One)

summary(paat_rvet_One_aov)  

## TAPSE 32

TAPSE_32 <- WASHU_all_weeks %>%
  select(`Shenan (0=No BPD, 1= BPD, 4=dead, transfer, withdrew)`, TAPSE_32) %>%
  filter(!is.na(TAPSE_32)) %>%
  filter(!(`Shenan (0=No BPD, 1= BPD, 4=dead, transfer, withdrew)` == 4))


TAPSE_32 %>%
  group_by(`Shenan (0=No BPD, 1= BPD, 4=dead, transfer, withdrew)`) %>%
  summarize(mean = mean(TAPSE_32))

TAPSE_32 %>%
  table()

TAPSE_32_aov <- aov(TAPSE_32 ~ `Shenan (0=No BPD, 1= BPD, 4=dead, transfer, withdrew)`, data = TAPSE_32)

summary(TAPSE_32_aov)  

## TAPSE 36

TAPSE_36 <- WASHU_all_weeks %>%
  select(`Shenan (0=No BPD, 1= BPD, 4=dead, transfer, withdrew)`, TAPSE_36) %>%
  filter(!is.na(TAPSE_36)) %>%
  filter(!(`Shenan (0=No BPD, 1= BPD, 4=dead, transfer, withdrew)` == 4))


TAPSE_36 %>%
  group_by(`Shenan (0=No BPD, 1= BPD, 4=dead, transfer, withdrew)`) %>%
  summarize(mean = mean(TAPSE_36))

TAPSE_36 %>%
  table()

TAPSE_36_aov <- aov(TAPSE_36 ~ `Shenan (0=No BPD, 1= BPD, 4=dead, transfer, withdrew)`, data = TAPSE_36)

summary(TAPSE_36_aov)

## TAPSE One

TAPSE_One <- WASHU_all_weeks %>%
  select(`Shenan (0=No BPD, 1= BPD, 4=dead, transfer, withdrew)`, TAPSE_One) %>%
  filter(!is.na(TAPSE_One)) %>%
  filter(!(`Shenan (0=No BPD, 1= BPD, 4=dead, transfer, withdrew)` == 4))


TAPSE_One %>%
  group_by(`Shenan (0=No BPD, 1= BPD, 4=dead, transfer, withdrew)`) %>%
  summarize(mean = mean(TAPSE_One))

TAPSE_One %>%
  table()

TAPSE_One_aov <- aov(TAPSE_One ~ `Shenan (0=No BPD, 1= BPD, 4=dead, transfer, withdrew)`, data = TAPSE_One)

summary(TAPSE_One_aov)
