################################################################################
# Map vs inspection: % of area above the NGR, with and without dependence,
# from 1,000 of the joint draws of posterior_draws.R; every grid cell, and
# the pointwise probabilities for the maps. V and Co are not modelled.
################################################################################
source("R/All_Models/model_data.R")
dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)
ngr  <- read.csv("results/data/ngr.csv")
REG  <- intersect(ngr$element, ELEMENTS)
USES <- c("industrial", "urban", "other")
stopifnot(all(FACTOR_OF[REG] == "sulphide"))
S <- 1000

grid <- read.csv("results/data/covariates_grid.csv")
stopifnot(nrow(grid) == nrow(geo$A_grid),
          nrow(grid) == nrow(geo$coords_grid),
          max(abs(as.matrix(grid[, c("x", "y")]) - geo$coords_grid)) < 1e-6,
          all(grid$ngr_use %in% USES))
set.seed(2026)
spacing <- median(sapply(sample(nrow(grid), 500), function(i) {
  d <- sqrt((grid$x - grid$x[i])^2 + (grid$y - grid$y[i])^2); min(d[d > 0]) }))
area_km2 <- nrow(grid) * spacing^2 / 1e6
ngr_mat  <- as.matrix(ngr[match(REG, ngr$element), USES])   # elements x land uses

N_GRID <- nrow(grid)
use   <- grid$ngr_use
lthr  <- log(t(ngr_mat)[match(use, USES), ])                # cells x elements

## posterior draws
draws <- readRDS("results/posterior/F3_draws.rds")
stopifnot(identical(draws$input_md5, MODEL_INPUT_MD5))
keep <- sort(sample.int(ncol(draws$g$sulphide), S))
g   <- as.matrix(geo$A_grid %*% draws$g$sulphide[, keep])
g_obs <- as.matrix(geo$A_obs %*% draws$g$sulphide[, keep])
mu  <- draws$intercept[REG, keep, drop = FALSE]
hyp <- draws$hyper[keep, ]
rm(draws)
lam <- sapply(REG, function(el) if (el == ANCHORS[["sulphide"]]) rep(1, S) else
  hyp[, sprintf("Beta for u_sulphide_%s", el)])
sd_patchy <- 1 / sqrt(hyp[, lik_prec(n_el + 1)])
sd_dil    <- 1 / sqrt(hyp[, "Precision for d"])
sd_noise  <- sapply(REG, function(el) 1 / sqrt(hyp[, lik_prec(match(el, ELEMENTS))]))

noise <- function(n, sd) sweep(matrix(rnorm(n * S), n, S), 2, sd, "*")
## exceedance per location and draw of each element in `els` (TRUE if any);
## gg: smooth part at the locations, lt: log NGR (locations x elements)
union_hits <- function(sample, joint, patchy = TRUE, gg = g, lt = lthr,
                       els = seq_along(REG)) {
  n <- nrow(gg)
  e <- noise(n, sd_patchy); d <- noise(n, sd_dil)
  hit <- matrix(FALSE, n, S)
  for (j in els) {
    p <- if (joint) seq_len(S) else sample.int(S)
    val <- sweep(sweep(gg[, p], 2, lam[p, j], "*"), 2, mu[j, p], "+")
    if (sample) {
      if (!joint) { e <- noise(n, sd_patchy); d <- noise(n, sd_dil) }
      val <- val + d[, p] + noise(n, sd_noise[p, j])
      if (patchy) val <- val + sweep(e[, p], 2, lam[p, j], "*")
    }
    hit <- hit | val > lt[, j]
  }
  hit
}

runs <- list(
  "map, joint"                        = union_hits(FALSE, TRUE),
  "map, independent"                  = union_hits(FALSE, FALSE),
  "inspection, joint"                 = union_hits(TRUE, TRUE),
  "inspection, independent"           = union_hits(TRUE, FALSE),
  "inspection, joint, no patchy part" = union_hits(TRUE, TRUE, patchy = FALSE))

## pointwise probabilities, for the maps
write.csv(data.frame(x = grid$x, y = grid$y, ngr_use = use,
                     p_map = rowMeans(runs[["map, joint"]]),
                     p_inspection = rowMeans(runs[["inspection, joint"]])),
          "results/posterior/compliance_grid.csv", row.names = FALSE)

scopes <- c(list(all = seq_len(N_GRID)), lapply(setNames(USES, USES), function(u) which(use == u)))
scope_km2 <- c(all = area_km2, sapply(setNames(USES, USES), function(u) area_km2 * mean(grid$ngr_use == u)))
tab <- do.call(rbind, lapply(names(runs), function(nm) do.call(rbind, lapply(names(scopes), function(sc) {
  idx <- scopes[[sc]]
  if (!length(idx)) return(NULL)
  pct <- 100 * colMeans(runs[[nm]][idx, , drop = FALSE])
  data.frame(quantity = nm, scope = sc, scope_km2 = scope_km2[[sc]], median_pct = median(pct),
             q05_pct = quantile(pct, 0.05, names = FALSE), q95_pct = quantile(pct, 0.95, names = FALSE),
             median_km2 = median(pct) / 100 * scope_km2[[sc]])
}))))

## observed samples
cov_obs <- read.csv("results/data/covariates_obs.csv")
stopifnot(nrow(cov_obs) == nrow(samples), all(cov_obs$id == samples$id),
          max(abs(as.matrix(cov_obs[, c("x", "y")]) - geo$coords_obs)) < 1e-6,
          all(cov_obs$ngr_use %in% USES))
use_obs <- cov_obs$ngr_use
x_obs <- as.matrix(samples[, REG])
obs_hit <- rowSums(x_obs > t(ngr_mat)[match(use_obs, USES), ]) > 0
## check: a new sample at each sampled site, against the observed samples
lthr_obs <- log(t(ngr_mat)[match(use_obs, USES), ])
site_hits <- union_hits(TRUE, TRUE, gg = g_obs, lt = lthr_obs)
scopes_obs <- c(list(all = seq_along(use_obs)),
                lapply(setNames(USES, USES), function(u) which(use_obs == u)))
tab <- rbind(tab, do.call(rbind, lapply(names(scopes_obs), function(sc) {
  idx <- scopes_obs[[sc]]
  if (!length(idx)) return(NULL)
  pct <- 100 * colMeans(site_hits[idx, , drop = FALSE])
  rbind(data.frame(quantity = "observed samples", scope = sc, scope_km2 = NA,
                   median_pct = 100 * mean(obs_hit[idx]), q05_pct = NA, q95_pct = NA,
                   median_km2 = NA),
        data.frame(quantity = "inspection at the sampled sites", scope = sc, scope_km2 = NA,
                   median_pct = median(pct), q05_pct = quantile(pct, 0.05, names = FALSE),
                   q95_pct = quantile(pct, 0.95, names = FALSE), median_km2 = NA))
})))
write.csv(tab, "results/tables/regulatory_union.csv", row.names = FALSE)
print(tab, digits = 3, row.names = FALSE)

## the same check element by element
check <- do.call(rbind, lapply(seq_along(REG), function(j) {
  hits <- union_hits(TRUE, TRUE, gg = g_obs, lt = lthr_obs, els = j)
  obs_j <- x_obs[, j] > ngr_mat[j, match(use_obs, USES)]
  do.call(rbind, lapply(names(scopes_obs), function(sc) {
    idx <- scopes_obs[[sc]]
    pct <- 100 * colMeans(hits[idx, , drop = FALSE])
    data.frame(element = REG[j], scope = sc, n = length(idx),
               observed_pct = 100 * mean(obs_j[idx]), predicted_pct = median(pct),
               q05_pct = quantile(pct, 0.05, names = FALSE),
               q95_pct = quantile(pct, 0.95, names = FALSE))
  }))
}))
write.csv(check, "results/tables/regulatory_check_elements.csv", row.names = FALSE)
print(check[check$scope == "all", ], digits = 3, row.names = FALSE)
