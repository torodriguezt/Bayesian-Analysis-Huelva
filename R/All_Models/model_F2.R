################################################################################
# F2: signatures with a smooth Matern field only (phi = 1).
################################################################################
source("R/All_Models/model_data.R")

fit_model(build_model(Y, scale = "smooth"), "F2", lso = TRUE)
