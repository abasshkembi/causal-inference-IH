# Title: Sensitivity analysis
# Author: Abas Shkembi (ashkembi@umich.edu)
# Last updated: May 8, 2025

#------------------------------------------------------------
### Description
### Conducts sensitivity analysis of most important controls under four scenarios
### (1) using the best performing model along with IPW; 
### (2) using the best performing model but not IPW; 
### (3) using a more traditionally-utilized, forward selection linear regression model along with IPW; and 
### (4) using a forward selection linear regression model without IPW.
#------------------------------------------------------------



# begin comparison of methods
# 1. RF + IPW
# 2. RF - IPW
# 3. FWDSEL + IPW
# 4. FWDSEL-IPW

# 1. 
rf.mod.Pb # using the same model from before

# 2.
set.seed(722); rf.mod.Pb.noIPW <- randomForest(as.formula(paste("logBloodPb ~ ", paste(tasks, collapse = "+"))),
                                               data=Bv5_comp, 
                                               #weights = Bv5_comp$ipw, # commenting this out
                                               ntree = 500, 
                                               importance = TRUE)

# 3. 
lm.mod.Pb.IPW1 <- lm(as.formula(paste("logBloodPb ~ ", paste(tasks, collapse = "+"))),
                     data=Bv5_comp, 
                     weights = ipw)

lm.mod.Pb.IPW2 <- ols_step_forward_p(lm.mod.Pb.IPW1)


flm.mod.Pb.IPW <- lm(as.formula(paste("logBloodPb ~ ", paste(lm.mod.Pb.IPW2$predictors, collapse = "+"))),
                     data=Bv5_comp, 
                     weights = ipw)

summary(flm.mod.Pb.IPW)


# 4. 
lm.mod.Pb.noIPW1 <- lm(as.formula(paste("logBloodPb ~ ", paste(tasks, collapse = "+"))),
                       data=Bv5_comp#, 
                       #weights = ipw # commenting this out
)

lm.mod.Pb.noIPW2 <- ols_step_forward_p(lm.mod.Pb.noIPW1)


flm.mod.Pb.noIPW <- lm(as.formula(paste("logBloodPb ~ ", paste(lm.mod.Pb.noIPW2$predictors, collapse = "+"))),
                       data=Bv5_comp#, 
                       #weights = ipw # commenting this out
)

summary(flm.mod.Pb.noIPW)




# create function to get hypothetical reductions thru bootstrapping
ewaste_get_counterfactual_comparison <- function(model, boot_i) {
  
  setup_df <- tibble(var = tasks)
  boot_cf_est <- NULL
  diff <- c()
  perc_diff <- c()
  
  for(B in 1:boot_i) {
    set.seed(722+B); b_i <- sample(1:nrow(Bv5_comp), nrow(Bv5_comp), replace = TRUE)
    Bv5_boot <- Bv5_comp[b_i,]
    if(model == "rf.mod.Pb") pred_boot <- predict(rf.mod.Pb, newdata = Bv5_boot)
    if(model == "rf.mod.Pb.noIPW") pred_boot <- predict(rf.mod.Pb.noIPW, newdata = Bv5_boot)
    if(model == "flm.mod.Pb.IPW") pred_boot <- predict(flm.mod.Pb.IPW, newdata = Bv5_boot)
    if(model == "flm.mod.Pb.noIPW") pred_boot <- predict(flm.mod.Pb.noIPW, newdata = Bv5_boot)
    
    for(i in 1:length(tasks)) {
      Bv5_boot_cf <- Bv5_boot
      i_var <- tasks[i]
      if(str_detect(i_var, "^PPE")) {
        min_value <- max(Bv5_comp[[i_var]])
      } else {
        min_value <- min(Bv5_comp[[i_var]])
      }
      Bv5_boot_cf[[i_var]] <- min_value
      
      if(model == "rf.mod.Pb") pred_boot_cf <- predict(rf.mod.Pb, newdata = Bv5_boot_cf)
      if(model == "rf.mod.Pb.noIPW") pred_boot_cf <- predict(rf.mod.Pb.noIPW, newdata = Bv5_boot_cf)
      if(model == "flm.mod.Pb.IPW") pred_boot_cf <- predict(flm.mod.Pb.IPW, newdata = Bv5_boot_cf)
      if(model == "flm.mod.Pb.noIPW") pred_boot_cf <- predict(flm.mod.Pb.noIPW, newdata = Bv5_boot_cf)
      
      # average difference
      diff[i] <- mean(exp(pred_boot) - exp(pred_boot_cf))
      
      #percent difference
      perc_diff[i] <- mean((exp(pred_boot) - exp(pred_boot_cf))/exp(pred_boot))
    }
    
    temp_df <- setup_df %>% mutate(diff_est = diff, perc_diff_est = perc_diff) %>% mutate(iteration = B)
    boot_cf_est <- rbind(boot_cf_est, temp_df)
    
    if(B %% 100 == 0) print(paste0("Model: ", model, ". Bootstrap: ", B, " of ", boot_i))
  }
  
  boot_cf_est %>%
    group_by(var) %>%
    summarise(mean = round(mean(diff_est), 2),
              q2.5 = round(quantile(diff_est, probs = 0.025), 2),
              q97.5 = round(quantile(diff_est, probs = 0.975), 2),
              perc_mean = round(mean(perc_diff_est), 2),
              perc_q2.5 = round(quantile(perc_diff_est, probs = 0.025), 2),
              perc_q97.5 = round(quantile(perc_diff_est, probs = 0.975), 2)
    ) %>% 
    ungroup()
}

## extract estimates
cf_estimates_compare <- ewaste_get_counterfactual_comparison(model = "rf.mod.Pb", boot_i = 500) %>% mutate(model = "Random Forest & IPW") %>%
  rbind(ewaste_get_counterfactual_comparison(model = "rf.mod.Pb.noIPW", boot_i = 500) %>% mutate(model = "Random Forest & No IPW")) %>%
  rbind(ewaste_get_counterfactual_comparison(model = "flm.mod.Pb.IPW", boot_i = 500) %>% mutate(model = "Forward Selection & IPW")) %>%
  rbind(ewaste_get_counterfactual_comparison(model = "flm.mod.Pb.noIPW", boot_i = 500) %>% mutate(model = "Forward Selection & No IPW"))

rf_vars_for_plotting <- cf_estimates_compare %>%
  arrange(model, -mean) %>%
  filter(model == "Random Forest & IPW") %>%
  group_by(model) %>%
  mutate(first = row_number()) %>%
  ungroup() %>%
  filter(first %in% 1:5) %>% 
  .$var %>% 
  unique

rf_comparison_df <- cf_estimates_compare %>%
  arrange(model, -mean) %>%
  group_by(model) %>%
  mutate(first = row_number()) %>%
  ungroup() %>%
  filter(var %in% rf_vars_for_plotting) %>%
  mutate(first = ifelse(str_detect(model, "Forward") & (mean+q2.5+q97.5) == 0, 50, first)) %>%
  mutate(model = factor(model, levels = c("Random Forest & IPW", "Random Forest & No IPW", "Forward Selection & IPW", "Forward Selection & No IPW"))) %>%
  mutate(var = factor(var, levels = rf_vars_for_plotting))


fwd_vars_for_plotting <- cf_estimates_compare %>%
  arrange(model, -mean) %>%
  filter(model == "Forward Selection & No IPW") %>%
  group_by(model) %>%
  mutate(first = row_number()) %>%
  ungroup() %>%
  filter(first %in% 1:5) %>% 
  .$var %>% 
  unique

fwd_comparison_df <- cf_estimates_compare %>%
  arrange(model, -mean) %>%
  group_by(model) %>%
  mutate(first = row_number()) %>%
  ungroup() %>%
  filter(var %in% fwd_vars_for_plotting) %>%
  mutate(first = ifelse(str_detect(model, "Forward") & (mean+q2.5+q97.5) == 0, 50, first)) %>%
  mutate(model = factor(model, levels = c("Random Forest & IPW", "Random Forest & No IPW", "Forward Selection & IPW", "Forward Selection & No IPW"))) %>%
  mutate(var = factor(var, levels = fwd_vars_for_plotting))

