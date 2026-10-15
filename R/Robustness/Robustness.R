# Checks of the main model F3. Run after All_Models.R.

rm(list = ls()); source("R/Robustness/simulation_identifiability.R")   # is phi identifiable?

# The sequential driver aggregates explicitly after its worker completes.
# Parallel users run SIM_SUMMARY_ONLY=1 separately after ALL shards succeed.
if (!SIM$summary_only && SIM$shards == 1L && is.null(SIM$task)) (function() {
  previous_summary <- Sys.getenv("SIM_SUMMARY_ONLY", unset = NA_character_)
  on.exit(if (is.na(previous_summary)) Sys.unsetenv("SIM_SUMMARY_ONLY") else
    Sys.setenv(SIM_SUMMARY_ONLY = previous_summary))
  Sys.setenv(SIM_SUMMARY_ONLY = "1")
  source("R/Robustness/simulation_identifiability.R")
})()
