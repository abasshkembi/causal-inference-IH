# Title: Hypothetical reductions
# Author: Abas Shkembi (ashkembi@umich.edu)
# Last updated: May 8, 2025

#------------------------------------------------------------
### Description
### Estimates hypothetical reductions in heavy metal concentrations
### using random forest model
#------------------------------------------------------------

library(gbm)
library(randomForest)
library(tidyverse)
library(olsrr)
library(lme4)
library(glmnet)
library(mice)
library(ggh4x)

# begin hypothetical reductions

## build random forest models for each heavy metal
set.seed(722); rf.mod.Pb <- randomForest(as.formula(paste("logBloodPb ~ ", paste(tasks, collapse = "+"))),
                                         data=Bv5_comp, 
                                         weights = Bv5_comp$ipw,
                                         ntree = 500, 
                                         importance = TRUE)

set.seed(722); rf.mod.Mn <- randomForest(as.formula(paste("logBloodMn ~ ", paste(tasks, collapse = "+"))),
                                         data=Bv5_comp, 
                                         weights = Bv5_comp$ipw,
                                         ntree = 500, 
                                         importance = TRUE)

set.seed(722); rf.mod.Ni <- randomForest(as.formula(paste("logUrinaryNi ~ ", paste(tasks, collapse = "+"))),
                                         data=Bv5_comp, 
                                         weights = Bv5_comp$ipw,
                                         ntree = 500, 
                                         importance = TRUE)

set.seed(722); rf.mod.Cu <- randomForest(as.formula(paste("logUrinaryCu ~ ", paste(tasks, collapse = "+"))),
                                         data=Bv5_comp, 
                                         weights = Bv5_comp$ipw,
                                         ntree = 500, 
                                         importance = TRUE)



# create function to get hypothetical reductions thru bootstrapping
ewaste_get_counterfactual <- function(model, boot_i) {
  
  setup_df <- tibble(var = tasks)
  boot_cf_est <- NULL
  diff <- c()
  perc_diff <- c()
  
  for(B in 1:boot_i) {
    set.seed(722+B); b_i <- sample(1:nrow(Bv5_comp), nrow(Bv5_comp), replace = TRUE)
    Bv5_boot <- Bv5_comp[b_i,]
    if(model == "rf.mod.Pb") pred_boot <- predict(rf.mod.Pb, newdata = Bv5_boot)
    if(model == "rf.mod.Mn") pred_boot <- predict(rf.mod.Mn, newdata = Bv5_boot)
    if(model == "rf.mod.Ni") pred_boot <- predict(rf.mod.Ni, newdata = Bv5_boot)
    if(model == "rf.mod.Cu") pred_boot <- predict(rf.mod.Cu, newdata = Bv5_boot)
    
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
      if(model == "rf.mod.Mn") pred_boot_cf <- predict(rf.mod.Mn, newdata = Bv5_boot_cf)
      if(model == "rf.mod.Ni") pred_boot_cf <- predict(rf.mod.Ni, newdata = Bv5_boot_cf)
      if(model == "rf.mod.Cu") pred_boot_cf <- predict(rf.mod.Cu, newdata = Bv5_boot_cf)
      
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

# get hypothetical reductions by calling function
cf_estimates <- ewaste_get_counterfactual(model = "rf.mod.Pb", boot_i = 500) %>% mutate(metal = "Pb (Blood)") %>%
  rbind(ewaste_get_counterfactual(model = "rf.mod.Mn", boot_i = 500) %>% mutate(metal = "Mn (Blood)")) %>%
  rbind(ewaste_get_counterfactual(model = "rf.mod.Ni", boot_i = 500) %>% mutate(metal = "Ni (Urine)")) %>%
  rbind(ewaste_get_counterfactual(model = "rf.mod.Cu", boot_i = 500) %>% mutate(metal = "Cu (Urine)"))


## get clean variable crosswalk for figures
variable_crosswalk <- xlsx::read.xlsx("variables.xlsx", sheetIndex = 1) %>% as_tibble()
variable_crosswalk <- variable_crosswalk %>%
  mutate(hierarchy = factor(hierarchy, levels = c("Elimination", "Substitution", "Administrative Control", "PPE")))

## annotations to the figure
anno <- tibble(
  hierarchy = factor(c("Elimination", "Substitution", "Administrative Control", "PPE"), levels = c("Elimination", "Substitution", "Administrative Control", "PPE")),
  x1 = c(5.6, 5.6, 5.6, 5.6)-1,
  x2 = c(5.65, 5.65, 5.65, 5.65)-1,
  y1 = c(10, 11, 16, 7),
  y2 = c(1, 1, 1, 1),
  lab = c("Interpretation:\nHow much, on average, can\nwe reduce blood and urinary\nheavy metal concentrations if\ncertain e-waste dismantling was\neliminated, in comparison\nwith doing nothing?", 
          "Interpretation:\nHow much, on average, can\nwe reduce blood and urinary\nheavy metal concentrations if workers\nstopped using certain tools for\ndismantling e-waste, in comparison\nwith doing nothing?", 
          "Interpretation:\nHow much, on average, can\nwe reduce blood and urinary\nheavy metal concentrations if \nworkers were trained to stop certain\ndangerous work activities, in\ncomparison with doing nothing?", 
          "Interpretation:\nHow much, on average, can\nwe reduce blood and urinary\nheavy metal concentrations if\nevery worker wears certain PPE,\nin comparison with doing nothing?"),
  x_lab = c(5.7-1, 5.7-1, 5.7-1, 5.7-1),
  y_lab = c(6, 5.5, 8.5, 4)
)

## final figure
fig3_final <- cf_estimates %>% 
  mutate(metal = factor(metal, levels = c("Pb (Blood)", "Mn (Blood)", "Ni (Urine)", "Cu (Urine)"),
                        labels = c("Blood Pb\n(µg/dL)", "Blood Mn\n(µg/L)", "Urinary Ni\n(µg/g creatinine)", "Urinary Cu\n(µg/g creatinine)"))) %>%
  left_join(variable_crosswalk, by = "var") %>%
  arrange(hierarchy, -perc_mean) %>%
  mutate(id = row_number()) %>%
  group_by(clean_name) %>%
  mutate(id = min(id)) %>%
  ungroup() %>%
  mutate(perc_mean = perc_mean*100) %>%
  mutate(perc_color = ifelse(abs(perc_mean) >= 22, "white", "black")) %>%
  mutate(sig = ifelse(perc_mean < -1 & q97.5 < 0, "sig",
                      ifelse(perc_mean > 1 & q2.5 > 0, "sig",
                             "not_sig"))) %>%
  mutate(label = ifelse(sig == "not_sig", "", paste0(mean, " (", q2.5, ", ", q97.5, ")"))) %>%
  ggplot() + 
  geom_tile(aes(x=metal, y=fct_reorder(clean_name, -id), fill=perc_mean))  +
  geom_text(aes(x = metal, y = fct_reorder(clean_name, -id), label = label, color = perc_color)) +
  geom_segment(data = anno, aes(x = x2, xend = x2, y = y1, yend = y2), color = "grey70") +
  geom_segment(data = anno, aes(x = x1, xend = x2, y = y1, yend = y1), color = "grey70") +
  geom_segment(data = anno, aes(x = x1, xend = x2, y = y2, yend = y2), color = "grey70") +
  geom_text(data = anno, aes(x = x_lab, y = y_lab, label = lab), color = "grey45", hjust = 0, size = 3) +
  #annotate("segment", x = 6, xend = 6, y = 1, yend = 5) +
  facet_grid2(vars(hierarchy), 
              scales = "free_y", 
              space = "free",
              switch = "y",
              strip = strip_themed(
                background_y = list(element_rect(fill = "#5A82AD", colour = "#5A82AD"),
                                    element_rect(fill = "#91BD60", colour = "#91BD60"),
                                    element_rect(fill = "#DC6C39", colour = "#DC6C39"),
                                    element_rect(fill = "#D83730", colour = "#D83730")),
                text_y = list(element_text(colour = "white", face = "bold"), 
                              element_text(colour = "white", face = "bold"), 
                              element_text(colour = "white", face = "bold"), 
                              element_text(colour = "white", face = "bold")))) +
  scale_fill_gradientn(colours = c("#a9350f", "#df9a82", "white", "#9dafd4", "#2866a9"),
                       limits = c(-40, 40), breaks = c(-25, -15, 0, 15, 25),
                       labels = c("",
                                  expression("" %<-% "Increase in concentration"), 
                                  "", 
                                  expression("Reduction in concentration" %->% ""),
                                  ""),
                       name = "          ") +
  scale_color_manual(
    values = c("white" = "white",
               "black" = "black")
  ) + 
  theme_bw() +
  coord_cartesian(xlim = c(1, 5.5)) +
  theme(panel.spacing = unit(0.3, "mm"), 
        panel.border = element_rect(color = "white"),
        strip.placement = "outside",
        axis.ticks = element_blank(),
        panel.grid = element_blank(),
        #legend.text = element_blank(),
        legend.direction = "horizontal",
        legend.position = "bottom",
        legend.justification = "left",
        legend.key.width = unit(0.08, "npc"),
        legend.key.height = unit(0.01, "npc")) +
  labs(x = NULL, y = NULL,
       title = "Hypothetical Average Reductions in Heavy Metal Concentrations") +
  guides(color = "none")

fig3_final
# save as svg with 900 height, 1100 width






