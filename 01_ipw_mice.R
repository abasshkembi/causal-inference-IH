# Title: IPW & MICE
# Author: Abas Shkembi (ashkembi@umich.edu)
# Last updated: May 8, 2025

#------------------------------------------------------------
### Description
### Conducts inverse probability weighting for selection bias
### Conducts multiple imputation for missing data
#------------------------------------------------------------

library(gbm)
library(randomForest)
library(tidyverse)
library(olsrr)
library(lme4)
library(glmnet)
library(mice)
library(ggh4x)

####### ipw of video participation

# logistic regression
ipw_mod <- Bv3 %>% glmer(video_part ~ Age + Numberfamsupported + OverMinWage + Education + (1|Country), family = "binomial", data = .)
summary(ipw_mod)

# assign ipw to participants
Bv4 <- Bv3 %>%
  mutate(index_ID = as.character(1:nrow(.))) %>%
  left_join(
    tibble(
      index_ID = names(predict(ipw_mod)),
      prob_video = predict(ipw_mod, type = "response")
    ),
    by = "index_ID"
  ) %>%
  mutate(
    ipw = ifelse(video_part == 0, 1/(1-prob_video), 1/prob_video)
  )

# get summary of ipws
summary(Bv4$ipw)
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#   1.113   1.166   1.219   1.993   1.393   8.628      17 













####### multiple imputation

tasks <- colnames(Bv4[-1])[which(!(colnames(Bv4[-1]) %in% 
                                     c("Age", "Numberfamsupported", "OverMinWage", "Education", "Country", 
                                       "Marital", "Employed", "Sex", "video_part",
                                       "index_ID", "prob_video", "ipw",
                                       "YearsEwaste", "BMI", "ResidenceTime", "Hoursworked",
                                       Bv4 %>% select(starts_with("log")) %>% colnames, 
                                       Bv4 %>% select(starts_with("Adj")) %>% colnames)))]

Bv5 <- Bv4 %>%
  filter(video_part == 1) %>%
  select(logBloodPb, logBloodMn, logUrinaryNi, logUrinaryCu, all_of(tasks), ipw)

sum(is.na(Bv5)) # number of missing cells in data
# 11
41*50 # number of cells in data
# 2050

11/2050*100
# 0.54% of data is missing


library(mice)
set.seed(722);Bv5_mice <- mice(Bv5)
Bv5_comp <- complete(Bv5_mice) %>% as_tibble()
