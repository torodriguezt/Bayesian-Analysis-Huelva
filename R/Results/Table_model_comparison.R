################################################################################
# F0-F3: WAIC, leave-site-out log score and fitting time.
################################################################################
source("R/All_Models/model_data.R")
dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)

MODELS <- c(F0 = "independent elements", F1 = "patchy signatures (phi = 0)",
            F2 = "smooth signatures (phi = 1)", F3 = "smooth + patchy signatures")
## all four with empirical Bayes, so the comparison is like for like
FITS <- c(F0 = "F0", F1 = "F1", F2 = "F2", F3 = "F3_eb")
cmp <- do.call(rbind, lapply(names(MODELS), function(nm) {
  f <- readRDS(file.path("results/fits", paste0(FITS[[nm]], ".rds")))
  assert_current_fit(f)
  if (!is.finite(f$huelva$waic_elements) || !is.finite(f$huelva$lso))
    stop(nm, " has incomplete WAIC/group-CV results; finish its fit before making the comparison table.")
  data.frame(model = nm, structure = MODELS[[nm]], waic_elements = f$huelva$waic_elements,
             lso_logscore = f$huelva$lso, cpu_s = f$cpu.used[["Total"]])
}))
write.csv(cmp, "results/tables/model_comparison.csv", row.names = FALSE)
print(cmp, digits = 4, row.names = FALSE)
