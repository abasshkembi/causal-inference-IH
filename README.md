# causal-inference-IH
Using Casual Inference and Machine Learning with Exposure Determinant Modeling to Identify Important Workplace Controls

--------------------------------------------------------------------------------

*Code for*

**Shkembi A**, Virji MA, He J, Nambunmee K, Ruiz-Rudolph, Neitzel RL. Using Casual Inference and Machine Learning with Exposure Determinant Modeling to Identify Important Workplace Controls.

I (Abas Shkembi) conducted the analysis. For questions about this analysis, you can reach me at ashkembi@umich.edu. Jie He (jayhe@umich.edu) aided in the beta substitution for concentrations <LOD for blood Pb.

### Code

The R scripts in this repository document the code used to conduct the beta substitution, inverse probability weighting, multiple imputation, model comparison w/ LOOCV-RMSE, and the hypothetical reductions thru bootstrapping. The data is available for those interested upon reasonable request.

  * **00_beta_sub.R** - provides R code for beta substitution of values <LOD for blood Pb
  * **01_ipw_mice.R** - inverse probability weighting for selection bias and multiple imputation
  * **02_descriptives.R** - basic descriptive statistics of four heavy metal concentrations
  * **03_model_loocv_rmse.R** - comparison of potential exposure determinant models (i.e., linear regression, forward selection, LASSO, boosted regression tree, and random forest) using leave one out cross validation root mean squared error (LOOCV-RMSE)
  * **04_rf_hypothetical_reductions.R** - estimating the hypothetical reductions for each heavy metal concentrations under counterfactual workplace controls
  * **05_sens_analysis.R** - sensitivity analysis comparing model choice (random forest vs forward selection) and whether you account for selection bias (IPW vs no IPW)

If this approach is of interest to you, but you would like help integrating the code for your own data, please do not hesitate to contact me (Abas, ashkembi@umich.edu). I am more than happy to help and answer any questions.