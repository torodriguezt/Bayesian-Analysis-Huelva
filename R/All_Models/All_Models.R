################################################################################
######            Fitting every model of the article (repo root)           ######
################################################################################
# Each script saves its INLA fit to results/fits/.

rm(list = ls()); source("R/All_Models/model_F0.R")   # independent elements
rm(list = ls()); source("R/All_Models/model_F1.R")   # patchy only
rm(list = ls()); source("R/All_Models/model_F2.R")   # smooth only
rm(list = ls()); source("R/All_Models/model_F3_eb.R")  # main: smooth + patchy, tables and figures
rm(list = ls()); source("R/All_Models/model_F3.R")   # main with CCD, for the robustness simulation
