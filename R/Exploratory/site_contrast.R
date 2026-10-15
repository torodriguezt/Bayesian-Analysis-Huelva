################################################################################
# Huelva vs mining sites: coupling of regulated elements and NGR exceedance.
################################################################################
source("R/Data_preparation/config.R")
dir_tab <- "results/tables/exploratory"
N_COPULA <- 1e5

samples <- read_csv(file.path(DIR_DATA, "samples.csv"), show_col_types = FALSE)
use_huelva <- read_csv(file.path(DIR_DATA, "covariates_obs.csv"), show_col_types = FALSE)$ngr_use
el <- NGR$element; k <- length(el)
reg_sulphide <- intersect(SIGNATURES$sulphide, el)
witness      <- setdiff(SIGNATURES$sulphide, el)

poisson_binomial <- function(p) {
  out <- c(1, rep(0, length(p)))
  for (q in p) out <- c(out[1] * (1 - q), out[-1] * (1 - q) + head(out, -1) * q)
  out
}

## x, thr: samples x elements
coexceedance <- function(x, thr) {
  n <- nrow(x)
  p <- sapply(seq_len(k), function(j) sapply(thr[, j], function(t) mean(x[, j] > t)))
  observed    <- tabulate(rowSums(x > thr) + 1, k + 1)
  independent <- Reduce(`+`, lapply(seq_len(n), function(i) poisson_binomial(p[i, ])))
  ## Gaussian copula
  R  <- cor(qnorm(apply(log(x), 2, rank) / (n + 1)))
  zc <- pnorm(matrix(rnorm(N_COPULA * k), ncol = k) %*% chol(R + diag(1e-8, k)))
  xc <- sapply(seq_len(k), function(j) quantile(x[, j], zc[, j], type = 8, names = FALSE))
  key <- apply(thr, 1, paste, collapse = "|")
  copula <- Reduce(`+`, lapply(unique(key), function(kk) {
    t <- thr[match(kk, key), ]
    sum(key == kk) * tabulate(rowSums(sweep(xc, 2, t, ">")) + 1, k + 1) / N_COPULA
  }))
  list(union = c(observed = mean(rowSums(x > thr) > 0),
                 independent = mean(1 - apply(1 - p, 1, prod)),
                 comonotone = mean(apply(p, 1, max)),
                 gaussian_copula = 1 - copula[1] / n),
       contributing = paste(el[colSums(x > thr) > 0], collapse = " "),
       spectrum = data.frame(n_elements = 0:k, observed, independent, gaussian_copula = copula))
}

set.seed(2026)
coupling <- list(); unions <- list(); spectra <- list()
for (s in names(FILE_SITES)) {
  d <- as.data.frame(samples[samples$site == s, ])
  x <- as.matrix(d[, el])
  D <- as.matrix(dist(d[, c("x", "y")])); diag(D) <- Inf
  w <- intersect(witness, names(d)[colSums(is.na(d)) == 0])
  lx <- log(x[, reg_sulphide])
  coupling[[s]] <- data.frame(
    site = s, n = nrow(d), nn_median_m = median(apply(D, 1, min)),
    mean_cor_regulated_sulphide = mean(cor(lx)[upper.tri(cor(lx))]),
    cor_witness_regulated = cor(group_score(d, w), group_score(d, reg_sulphide)))
  uses <- c("industrial", "urban", "other", if (s == "Huelva") "applicable")
  for (u in uses) {
    thr <- if (u == "applicable") as.matrix(NGR[, use_huelva])
           else matrix(NGR[[u]], k, nrow(d))
    ce <- coexceedance(x, t(thr))
    unions[[paste(s, u)]] <- data.frame(site = s, ngr_use = u, t(ce$union),
                                        contributing = ce$contributing)
    spectra[[paste(s, u)]] <- data.frame(site = s, ngr_use = u, ce$spectrum)
  }
}
spec <- do.call(rbind, spectra)
spec <- spec[spec$n_elements <= max(5, spec$n_elements[spec$observed > 0]), ]
write_csv(do.call(rbind, coupling), file.path(dir_tab, "site_coupling.csv"))
write_csv(do.call(rbind, unions), file.path(dir_tab, "union_by_use.csv"))
write_csv(spec, file.path(dir_tab, "coexceedance.csv"))
print(do.call(rbind, unions)[, 1:6], digits = 2, row.names = FALSE)
