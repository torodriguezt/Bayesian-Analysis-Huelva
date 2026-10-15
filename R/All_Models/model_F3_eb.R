################################################################################
# F3 with empirical Bayes, for the article's tables and figures.
# posterior_draws.R integrates over the hyperparameters by refitting at draws
# of theta; the CCD fit (model_F3.R) is kept for the robustness simulation.
################################################################################
source("R/All_Models/model_data.R")

fit_model(build_model(Y, scale = "both"), "F3_eb", lso = TRUE)
