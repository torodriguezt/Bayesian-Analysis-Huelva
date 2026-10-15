################################################################################
# Spatial scale of each signature: semivariogram nugget fraction (1 = patchy).
################################################################################
source("R/Data_preparation/config.R")
suppressPackageStartupMessages(library(gstat))

N_BOOT <- 200; CUTOFF <- 8000; WIDTH <- 500
SCORES <- c("regulated", "witness", "sulphide", "lithogenic", "phosphogypsum")
dir_tab <- "results/tables/exploratory"

scores <- read_csv(file.path(DIR_DATA, "signature_scores.csv"), show_col_types = FALSE)
pts <- st_as_sf(scores, coords = c("x", "y"), crs = CRS_UTM)
D <- as.matrix(dist(scores[, c("x", "y")]))
n <- nrow(scores)

## Matern correlation, nu = 1
matern1 <- function(h, a) ifelse(h == 0, 1, (h / a) * besselK(pmax(h, 1e-9) / a, 1))

fit_vgm <- function(z) {
  pts$z <- z
  v <- variogram(z ~ 1, pts, cutoff = CUTOFF, width = WIDTH)
  init <- vgm(psill = var(z) / 2, model = "Mat", range = 1000, nugget = var(z) / 2, kappa = 1)
  f <- tryCatch(suppressWarnings(fit.variogram(v, init)), error = function(e) NULL)
  if (is.null(f) || isTRUE(attr(f, "singular")) || f$range[2] > 3 * CUTOFF)
    return(c(nugget = NA, psill = NA, a = NA))
  c(nugget = f$psill[1], psill = f$psill[2], a = f$range[2])
}

short_long <- function(z) {
  g <- function(lo, hi) {
    idx <- which(D > lo & D <= hi & upper.tri(D), arr.ind = TRUE)
    0.5 * mean((z[idx[, 1]] - z[idx[, 2]])^2)
  }
  g(0, 600) / g(3000, 8000)
}

set.seed(2026)
rows <- list()
for (s in SCORES) {
  z <- scores[[s]]
  p <- fit_vgm(z)
  nf <- p[["nugget"]] / (p[["nugget"]] + p[["psill"]])
  boot <- matrix(NA_real_, N_BOOT, 2)
  if (is.finite(nf)) {
    R <- chol(p[["psill"]] * matern1(D, p[["a"]]) + diag(p[["nugget"]] + 1e-8, n))
    for (b in seq_len(N_BOOT)) {
      zb <- drop(crossprod(R, rnorm(n)))
      pb <- fit_vgm(zb)
      boot[b, ] <- c(pb[["nugget"]] / (pb[["nugget"]] + pb[["psill"]]), short_long(zb))
    }
  }
  rows[[s]] <- data.frame(
    score = s, variance = var(z), nugget = p[["nugget"]], partial_sill = p[["psill"]],
    practical_range_m = sqrt(8) * p[["a"]], nugget_fraction = nf,
    nf_q05 = quantile(boot[, 1], 0.05, na.rm = TRUE),
    nf_q95 = quantile(boot[, 1], 0.95, na.rm = TRUE),
    ratio_short_long = short_long(z),
    ratio_q05 = quantile(boot[, 2], 0.05, na.rm = TRUE),
    ratio_q95 = quantile(boot[, 2], 0.95, na.rm = TRUE),
    boot_failed = sum(is.na(boot[, 1])), row.names = NULL)
}
tab <- do.call(rbind, rows)
write_csv(tab, file.path(dir_tab, "scale_by_signature.csv"))
print(tab[, c("score", "nugget_fraction", "nf_q05", "nf_q95", "practical_range_m",
              "ratio_short_long", "ratio_q05", "ratio_q95")], digits = 2)
