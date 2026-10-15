################################################################################
# F0: independent elements, each with its own Matern field (separate maps).
################################################################################
source("R/All_Models/model_data.R")

fit_model(build_model(Y, shared = FALSE), "F0", lso = TRUE)
