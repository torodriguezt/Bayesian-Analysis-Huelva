################################################################################
######          Tables and figures of the article (repo root)             ######
################################################################################
# Requires All_Models.R (F0-F2 and F3_eb). Tables go to results/tables/ and
# figures (PNG and PDF) to figures/.

rm(list = ls()); source("R/Results/posterior_draws.R")         # joint draws of F3 (about 6 min)
rm(list = ls()); source("R/Results/Table_descriptive.R")       # Annex IV elements and NGR exceedance
rm(list = ls()); source("R/Results/Table_model_comparison.R")  # F0-F3
rm(list = ls()); source("R/Results/Table_loadings.R")          # source signatures in F3
rm(list = ls()); source("R/Results/Table_scale_fractions.R")   # mappable fraction
rm(list = ls()); source("R/Results/Table_regulatory.R")        # map vs inspection under the NGR
rm(list = ls()); source("R/Results/Figures.R")                 # all figures
