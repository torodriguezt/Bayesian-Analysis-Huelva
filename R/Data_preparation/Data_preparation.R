################################################################################
######                 Data preparation (run from repo root)              ######
################################################################################
# Needs the confidential survey files in data/raw/. Writes to results/data/.

rm(list = ls()); source("R/Data_preparation/1_samples.R")
rm(list = ls()); source("R/Data_preparation/2_covariates.R")
rm(list = ls()); source("R/Data_preparation/3_mesh.R")
