################################################################################
# F3 loadings and element noise SD.
################################################################################
source("R/All_Models/model_data.R")
dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)

fit <- readRDS("results/fits/F3_eb.rds")
assert_current_fit(fit)
hp <- fit$summary.hyperpar
load <- data.frame(element = ELEMENTS, signature = FACTOR_OF[ELEMENTS],
                   anchor = ELEMENTS %in% ANCHORS, loading = 1, q025 = NA, q975 = NA,
                   row.names = NULL)
for (j in which(!load$anchor)) {
  b <- sprintf("Beta for u_%s_%s", load$signature[j], load$element[j])
  load[j, c("loading", "q025", "q975")] <- hp[b, c("0.5quant", "0.025quant", "0.975quant")]
}
load$sd_noise <- 1 / sqrt(hp[lik_prec(seq_len(n_el)), "0.5quant"])
write.csv(load, "results/tables/loadings.csv", row.names = FALSE)
print(load, digits = 3, row.names = FALSE)
