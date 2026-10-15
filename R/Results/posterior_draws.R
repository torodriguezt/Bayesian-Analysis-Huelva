################################################################################
# 2,000 joint posterior draws of F3 integrating over the hyperparameters:
# 40 draws of theta from the EB fit, F3 refitted at each theta, 50 latent
# draws per refit. Saves the smooth fields (mesh), the site scores, the
# intercepts, the hyperparameters of each draw and the realized fraction.
# Used by the tables and figures of the article (about 6 min).
################################################################################
source("R/All_Models/model_data.R")
dir.create("results/posterior", recursive = TRUE, showWarnings = FALSE)
N_THETA <- 40; N_PER <- 50

fit <- readRDS("results/fits/F3_eb.rds")
assert_current_fit(fit)
m <- build_model(Y, scale = "both")

set.seed(2026)
theta <- inla.hyperpar.sample(N_THETA, fit, intern = TRUE)
to_user <- function(th) setNames(vapply(seq_along(th), function(i)
  do.call(fit$misc$from.theta[[i]], list(th[[i]])), 0), rownames(fit$summary.hyperpar))
sel <- c(setNames(rep(list(seq_len(spde$n.spde)), n_fac), paste0("g_", names(FACTORS))),
         setNames(rep(list(seq_len(n_site)), n_fac), paste0("u_", names(FACTORS))),
         setNames(rep(list(1), n_el), paste0("int_", ELEMENTS)))

parts <- lapply(seq_len(N_THETA), function(t) {
  ft <- inla(m$formula, family = m$family,
             data = inla.stack.data(m$stack, spde = spde),
             control.family = m$control.family,
             control.predictor = list(A = inla.stack.A(m$stack)),
             control.compute = list(config = TRUE),
             control.fixed = list(prec = 1e-3),
             control.mode = list(theta = theta[t, ], fixed = TRUE))
  s <- inla.posterior.sample(N_PER, ft, selection = sel, seed = 2026 + t,
                             num.threads = "1:1", parallel.configs = FALSE)
  message(sprintf("  theta %d/%d", t, N_THETA))
  sapply(s, function(x) x$latent[, 1])
})
lat <- do.call(cbind, parts)
rows <- function(prefix) startsWith(rownames(lat), prefix)

draws <- list(
  g = lapply(setNames(names(FACTORS), names(FACTORS)),
             function(f) lat[rows(paste0("g_", f, ":")), ]),
  u = lapply(setNames(names(FACTORS), names(FACTORS)),
             function(f) lat[rows(paste0("u_", f, ":")), ]),
  intercept = lat[paste0("int_", ELEMENTS, ":1"), ],
  hyper = do.call(rbind, lapply(seq_len(N_THETA), function(t)
    matrix(to_user(theta[t, ]), N_PER, ncol(theta), byrow = TRUE,
           dimnames = list(NULL, rownames(fit$summary.hyperpar))))),
  theta_index = rep(seq_len(N_THETA), each = N_PER),
  input_md5 = MODEL_INPUT_MD5)
rownames(draws$intercept) <- ELEMENTS

## realized mappable fraction: share of the across-site variance of each score
draws$phi_realized <- sapply(names(FACTORS), function(f) {
  g_obs <- as.matrix(geo$A_obs %*% draws$g[[f]])
  v_smooth <- apply(g_obs, 2, var)
  v_smooth / (v_smooth + apply(draws$u[[f]] - g_obs, 2, var))
})
saveRDS(draws, "results/posterior/F3_draws.rds")
message(sprintf("Saved %d joint draws of F3 (%d theta x %d); realized fraction %s",
                ncol(lat), N_THETA, N_PER,
                paste(sprintf("%s %.2f", names(FACTORS), apply(draws$phi_realized, 2, median)),
                      collapse = ", ")))
