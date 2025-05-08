# Title: Participant descriptive statistics
# Author: Abas Shkembi (ashkembi@umich.edu)
# Last updated: May 8, 2025

#------------------------------------------------------------
### Description
### Conducts descriptive statistics for all, video, and non-video participants
#------------------------------------------------------------

library(gbm)
library(randomForest)
library(tidyverse)
library(olsrr)
library(lme4)
library(glmnet)
library(mice)
library(ggh4x)

# according to 2020 ACGIH BEIs
bei <- tibble(
  metal = c("AdjBloodMn", "AdjBloodPb", "AdjUrinaryCu", "AdjUrinaryNi"),
  bei = c(NA, 20, NA, NA)
)

# descriptives - overall

Bv4 %>%
  select(AdjBloodPb, logBloodMn, AdjUrinaryNi, logUrinaryCu) %>%
  mutate(AdjBloodMn = exp(logBloodMn),
         AdjUrinaryCu = exp(logUrinaryCu)) %>%
  select(-logBloodMn, -logUrinaryCu) %>%
  gather("metal", "value") %>%
  left_join(bei) %>%
  mutate(above_bei = ifelse(value > bei, 1, 0)) %>%
  group_by(metal) %>%
  summarise(median = median(value, na.rm = TRUE),
            q25 = quantile(value, probs = 0.25, na.rm = TRUE),
            q75 = quantile(value, probs = 0.75, na.rm = TRUE),
            max = max(value, na.rm = TRUE),
            n_bei = sum(above_bei, na.rm = TRUE)) %>%
  mutate(perc_bei = n_bei/226*100)

# descriptives - video

Bv4 %>%
  filter(video_part == 1) %>%
  select(AdjBloodPb, logBloodMn, AdjUrinaryNi, logUrinaryCu) %>%
  mutate(AdjBloodMn = exp(logBloodMn),
         AdjUrinaryCu = exp(logUrinaryCu)) %>%
  select(-logBloodMn, -logUrinaryCu) %>%
  gather("metal", "value") %>%
  left_join(bei) %>%
  mutate(above_bei = ifelse(value > bei, 1, 0)) %>%
  group_by(metal) %>%
  summarise(median = median(value, na.rm = TRUE),
            q25 = quantile(value, probs = 0.25, na.rm = TRUE),
            q75 = quantile(value, probs = 0.75, na.rm = TRUE),
            max = max(value, na.rm = TRUE),
            n_bei = sum(above_bei, na.rm = TRUE)) %>%
  mutate(perc_bei = n_bei/41*100)

# descriptives - non-video

Bv4 %>%
  filter(video_part == 0) %>%
  select(AdjBloodPb, logBloodMn, AdjUrinaryNi, logUrinaryCu) %>%
  mutate(AdjBloodMn = exp(logBloodMn),
         AdjUrinaryCu = exp(logUrinaryCu)) %>%
  select(-logBloodMn, -logUrinaryCu) %>%
  gather("metal", "value") %>%
  left_join(bei) %>%
  mutate(above_bei = ifelse(value > bei, 1, 0)) %>%
  group_by(metal) %>%
  summarise(median = median(value, na.rm = TRUE),
            q25 = quantile(value, probs = 0.25, na.rm = TRUE),
            q75 = quantile(value, probs = 0.75, na.rm = TRUE),
            max = max(value, na.rm = TRUE),
            n_bei = sum(above_bei, na.rm = TRUE)) %>%
  mutate(perc_bei = n_bei/185*100)


