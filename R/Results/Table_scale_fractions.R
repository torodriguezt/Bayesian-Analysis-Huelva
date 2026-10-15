################################################################################
# Mappable fraction of each signature in F3. Needs posterior_draws.R.
################################################################################
source("R/All_Models/model_data.R")
dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)

fit <- readRDS("results/fits/F3_eb.rds")
assert_current_fit(fit)
set.seed(2026)
hs <- inla.hyperpar.sample(4000, fit, improve.marginals = TRUE)
## realized fraction from the joint draws, integrated over the hyperparameters
draws <- readRDS("results/posterior/F3_draws.rds")
stopifnot(identical(draws$input_md5, MODEL_INPUT_MD5))
pr <- draws$phi_realized
summ <- function(x, signature, quantity)
  data.frame(signature = signature, quantity = quantity, median = median(x),
             q05 = quantile(x, 0.05, names = FALSE), q95 = quantile(x, 0.95, names = FALSE))

rows <- list()
for (k in seq_len(n_fac)) {
  f <- names(FACTORS)[k]
  sd_g <- hs[, paste0("Stdev for g_", f)]
  sd_e <- 1 / sqrt(hs[, lik_prec(n_el + k)])
  rows[[f]] <- rbind(summ(pr[, f], f, "phi realized (share across sites)"),
                     summ(sd_g^2 / (sd_g^2 + sd_e^2), f, "phi marginal-variance ratio"),
                     summ(sd_g, f, "sd smooth"), summ(sd_e, f, "sd patchy"),
                     summ(hs[, paste0("Range for g_", f)] / 1000, f, "range smooth (km)"))
}
contrast <- pr[, "sulphide"] - pr[, "lithogenic"]
tab <- rbind(do.call(rbind, rows),
             summ(contrast, "sulphide - lithogenic", "phi realized difference"),
             data.frame(signature = "sulphide - lithogenic",
                        quantity = "P(phi realized sulphide < lithogenic)",
                        median = mean(contrast < 0), q05 = NA, q95 = NA),
             summ(1 / sqrt(hs[, "Precision for d"]), "all elements", "sd dilution"))
write.csv(tab, "results/tables/scale_fractions.csv", row.names = FALSE)
print(tab, digits = 2, row.names = FALSE)
