# Title: Model comparison w/ LOOCV-RMSE
# Author: Abas Shkembi (ashkembi@umich.edu)
# Last updated: May 8, 2025

#------------------------------------------------------------
### Description
### Computes Leave One Out Cross Validation Root Mean Squared Errors (LOOCV-RMSE)
### To determine best model for causal inference
#------------------------------------------------------------



ewaste_get_rmse <- function(outcome, model) {
  
  outcome_i <- which(colnames(Bv5_comp) == outcome)
  remove <- which(!(c(1, 2, 3, 4, 5) %in% outcome_i))
  df_temp <- Bv5_comp %>% select(-all_of(remove))
  rmse <- c()
  
  if(model == "lm") {
    
    for(i in 1:nrow(df_temp)) {
      df_loocv <- df_temp[-i,]
      
      model <- lm(as.formula(paste(outcome, " ~ ", paste(tasks, collapse = "+"))),
                  data=df_loocv,
                  weights = ipw)
      pred.i <- predict(model, newdata = df_temp[i,])
      rmse[i] <- (as.numeric(df_temp[i,1]) - pred.i)^2
      
      if(i%%10 == 0) print(paste0("Iteration: ", i, " of ", nrow(df_temp)))
    }
    
  } else if(model == "forward2"){
    
    for(i in 1:nrow(df_temp)) {
      df_loocv2 <<- df_temp[-i,]
      
      temp_lm_model <- lm(as.formula(paste(outcome, " ~ ", paste(tasks, collapse = "+"))),
                          data=df_loocv2,
                          weights = ipw)
      flm_perform <- ols_step_forward_p(temp_lm_model)
      
      model <- lm(as.formula(paste(outcome, " ~ ", paste(flm_perform$predictors, collapse = "+"))),
                  data=df_loocv2,
                  weights = ipw)
      
      pred.i <- predict(model, newdata = df_temp[i,])
      rmse[i] <- (as.numeric(df_temp[i,1]) - pred.i)^2
      
      if(i%%10 == 0) print(paste0("Iteration: ", i, " of ", nrow(df_temp)))
    }
    
  } else if(model == "lasso"){
    
    ipw_i <- which(colnames(df_temp) %in% "ipw")
    x <- data.matrix(df_temp[c(-1, -ipw_i)]) # remove outcome and ipw columns
    
    for(i in 1:nrow(df_temp)) {
      df_loocv <- df_temp[-i,]
      
      x2 <- data.matrix(df_loocv[c(-1, -ipw_i)])
      y <- data.matrix(df_loocv[1])
      weights <- df_loocv$ipw
      
      model <- cv.glmnet(x2,y, weights = weights, family = "gaussian")
      
      pred.i <- predict(model, newx = x[i,], s = "lambda.1se")
      rmse[i] <- (as.numeric(df_temp[i,1]) - pred.i)^2
      
      if(i%%10 == 0) print(paste0("Iteration: ", i, " of ", nrow(df_temp)))
      
    }
    
  } else if(model == "gbm"){
    
    for(i in 1:nrow(df_temp)) {
      df_loocv <- df_temp[-i,]
      
      set.seed(722); model <- gbm(as.formula(paste(outcome, " ~ ", paste(tasks, collapse = "+"))),
                                  data=df_loocv, 
                                  weights = df_loocv$ipw,
                                  distribution='gaussian',
                                  cv.folds = 3,
                                  shrinkage = 0.1, 
                                  verbose=FALSE,
                                  n.trees = 500, 
                                  bag.fraction = 0.8,
                                  n.minobsinnode = 2,
                                  interaction.depth = 10)
      #ntrees_cv <- which.min(model$cv.error)
      pred.i <- predict(model, n.trees = 500, newdata = df_temp[i,])
      rmse[i] <- (as.numeric(df_temp[i,1]) - pred.i)^2
      
      if(i%%10 == 0) print(paste0("Iteration: ", i, " of ", nrow(df_temp)))
    }
    
  } else if(model == "rf"){
    
    for(i in 1:nrow(df_temp)) {
      df_loocv <- df_temp[-i,]
      
      set.seed(722); model <- randomForest(as.formula(paste(outcome, " ~ ", paste(tasks, collapse = "+"))),
                                           data=df_loocv, 
                                           weights = df_loocv$ipw,
                                           ntree = 500)
      pred.i <- predict(model, newdata = df_temp[i,])
      rmse[i] <- (as.numeric(df_temp[i,1]) - pred.i)^2
      
      if(i%%10 == 0) print(paste0("Iteration: ", i, " of ", nrow(df_temp)))
    }
    
    
  }
  
  return(sqrt(mean(rmse)))
}






### Lead

ewaste_get_rmse(outcome = "logBloodPb", model = "lm")
ewaste_get_rmse(outcome = "logBloodPb", model = "forward2")
ewaste_get_rmse(outcome = "logBloodPb", model = "lasso")
ewaste_get_rmse(outcome = "logBloodPb", model = "gbm")
ewaste_get_rmse(outcome = "logBloodPb", model = "rf")




### Mangenese

ewaste_get_rmse(outcome = "logBloodMn", model = "lm")
ewaste_get_rmse(outcome = "logBloodMn", model = "forward2")
ewaste_get_rmse(outcome = "logBloodMn", model = "lasso")
ewaste_get_rmse(outcome = "logBloodMn", model = "gbm")
ewaste_get_rmse(outcome = "logBloodMn", model = "rf")



### Nickel

ewaste_get_rmse(outcome = "logUrinaryNi", model = "lm")
ewaste_get_rmse(outcome = "logUrinaryNi", model = "forward2")
ewaste_get_rmse(outcome = "logUrinaryNi", model = "lasso")
ewaste_get_rmse(outcome = "logUrinaryNi", model = "gbm")
ewaste_get_rmse(outcome = "logUrinaryNi", model = "rf")



### Copper

ewaste_get_rmse(outcome = "logUrinaryCu", model = "lm")
ewaste_get_rmse(outcome = "logUrinaryCu", model = "forward2")
ewaste_get_rmse(outcome = "logUrinaryCu", model = "lasso")
ewaste_get_rmse(outcome = "logUrinaryCu", model = "gbm")
ewaste_get_rmse(outcome = "logUrinaryCu", model = "rf")


