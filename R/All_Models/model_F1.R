################################################################################
# F1: signatures with a patchy site effect only (phi = 0).
################################################################################
source("R/All_Models/model_data.R")

fit_model(build_model(Y, scale = "patchy"), "F1", lso = TRUE)
