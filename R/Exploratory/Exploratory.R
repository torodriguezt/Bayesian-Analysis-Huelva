################################################################################
######          Exploratory checks before modelling (repo root)           ######
################################################################################
# Needs Data_preparation.R. Tables go to results/tables/exploratory/.

rm(list = ls()); source("R/Exploratory/signatures.R")          # element set and signature scores
rm(list = ls()); source("R/Exploratory/scale_by_signature.R")  # nugget fraction per signature
rm(list = ls()); source("R/Exploratory/site_contrast.R")       # coupling and co-exceedance by site
