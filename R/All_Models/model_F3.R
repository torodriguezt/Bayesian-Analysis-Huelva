################################################################################
# F3, main model: smooth field + patchy effect per signature.
# Hyperparameters integrated with CCD (needs ~8+ GB of RAM; run it alone).
################################################################################
source("R/All_Models/model_data.R")

fit_model(build_model(Y, scale = "both"), "F3", lso = TRUE, int_strategy = "ccd")
