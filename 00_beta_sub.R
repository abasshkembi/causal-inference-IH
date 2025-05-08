# Title: Beta substitution for Pb
# Author: Abas Shkembi (ashkembi@umich.edu) & Jie He (jayhe@umich.edu)
# Last updated: May 8, 2025

#------------------------------------------------------------
### Description
### Conducts beta-substitution for values <LOD for blood Pb
#------------------------------------------------------------

library(gbm)
library(randomForest)
library(tidyverse)
library(olsrr)
library(lme4)
library(glmnet)
library(mice)
library(ggh4x)

# create function to perform beta-substitution

beta_substitution <- function(data, x_var, lod_var) {
  # Check for multiple LOD values in non-detects
  lod_values <- unique(data[[x_var]][data[[lod_var]] == 0])
  if (length(lod_values) > 1) {
    stop("Multiple LOD values detected. All non-detects must have the same LOD.")
  }
  LOD <- if (length(lod_values) == 0) NA else lod_values[1]
  
  n <- nrow(data)
  k <- sum(data[[lod_var]] == 0)
  
  if (k == 0) {
    return(list(
      beta_mean = NA_real_,
      beta_gm = NA_real_,
      substituted_mean = mean(data[[x_var]]),
      substituted_gm = exp(mean(log(data[[x_var]]))),
      substituted_gsd = NA_real_,
      substituted_x95 = NA_real_
    ))
  }
  
  if (n - k < 2) {
    stop("Insufficient detected values (need at least 2).")
  }
  
  # Step 1: Create a detects array
  detected <- data[[x_var]][data[[lod_var]] == 1]
  
  # Step 2: Calculate input and intermediate values
  y_bar <- mean(log(detected))
  z <- qnorm(k / n)
  fz <- dnorm(z) / (1 - pnorm(z))
  sy <- sqrt((y_bar - log(LOD))^2 / (fz - z)^2)
  f_sy_z <- (1 - pnorm(z - sy / n)) / (1 - pnorm(z))
  
  # Step 3: Calculate beta_mean
  beta_mean <- (n / k) * pnorm(z - sy) * exp(-sy * z + sy^2 / 2)
  
  # Step 4: Substitute each non-detect with beta_mean * LOD, and calculate the sample mean
  data_substituted_mean <- data[[x_var]]
  data_substituted_mean[data[[lod_var]] == 0] <- beta_mean * LOD
  mean_sub <- mean(data_substituted_mean)
  
  # Step 5: Calculate beta_gm
  beta_gm <- exp(- (n * (n - k) / k * log(f_sy_z) - sy * z - (n - k) / (2 * k * n) * sy^2))
  
  # Step 6: Substitute each non-detect with beta_gm * LOD, and calculate the sample GM         
  data_substituted_gm <- data[[x_var]]
  data_substituted_gm[data[[lod_var]] == 0] <- beta_gm * LOD
  gm_sub <- exp(mean(log(data_substituted_gm)))
  
  # Step 7: Recalculate sy and then calculate the sample GSD
  ratio <- mean_sub / gm_sub
  if (ratio > 1) {
    sy_new <- sqrt((2 * n / (n - 1)) * log(ratio))
  } else {
    sy_new <- 0
  }
  gsd_sub <- exp(sy_new)
  
  # Step 8: Calculate the sample 95th percentile
  x95 <- exp(log(gm_sub) - (sy_new^2) / (2 * n) + 1.645 * sy_new)
  
  list(
    beta_mean = beta_mean,
    beta_gm = beta_gm,
    substituted_mean = mean_sub,
    substituted_gm = gm_sub,
    substituted_gsd = gsd_sub,
    substituted_x95 = x95
  )
}


# AAS methods to analyze Blood Pb
## limit of detection for Pb : 1 ug/dL

bs_Pb <- Bio %>%
  arrange(PbBloodLOD) %>%
  mutate(PbBloodµgdL = as.numeric(PbBloodµgdL)) %>%
  mutate(newPbBloodLOD = ifelse(PbBloodLOD == 0 | PbBloodµgdL <= 1, 0, 1)) %>%
  select(PbBloodLOD, PbBloodµgdL, newPbBloodLOD) %>%
  filter(!is.na(PbBloodLOD)) %>%
  mutate(PbBloodµgdL = ifelse(newPbBloodLOD == 0, 1, PbBloodµgdL))

result_Pb <- beta_substitution(bs_Pb, "PbBloodµgdL", "newPbBloodLOD")
# $substituted_gm
# [1] 1.228665

#$substituted_mean
#[1] 2.67581
