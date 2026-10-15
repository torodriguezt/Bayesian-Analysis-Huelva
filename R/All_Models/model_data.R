# Data, priors and model builder shared by every model.

suppressPackageStartupMessages({library(INLA); library(Matrix)})
inla.setOption(num.threads = "4:1")
dir.create("results/fits", recursive = TRUE, showWarnings = FALSE)

samples <- read.csv("results/data/samples.csv", check.names = FALSE)
samples <- samples[samples$site == "Huelva", ]
geo <- readRDS("results/data/mesh.rds")
stopifnot(nrow(samples) == nrow(geo$coords_obs),
          max(abs(as.matrix(samples[, c("x", "y")]) - geo$coords_obs)) < 1e-6)

## signatures; the first element is the anchor (loading 1)
FACTORS <- list(sulphide   = c("Pb", "As", "Cu", "Zn", "Cd", "Bi", "Sb", "Tl", "Hg"),
                lithogenic = c("Al", "Sc", "Th", "Ga"))
ELEMENTS  <- unlist(FACTORS, use.names = FALSE)
FACTOR_OF <- setNames(rep(names(FACTORS), lengths(FACTORS)), ELEMENTS)
ANCHORS   <- sapply(FACTORS, `[`, 1)

n_site <- nrow(samples); n_el <- length(ELEMENTS); n_fac <- length(FACTORS)
Y <- log(as.matrix(samples[, ELEMENTS]))
stopifnot(all(is.finite(Y)))

## Link saved fits to the inputs and builder that actually produced them.
model_input_md5 <- function() tools::md5sum(c(
  "results/data/samples.csv", "results/data/mesh.rds", "R/All_Models/model_data.R"))
MODEL_INPUT_MD5 <- model_input_md5()

assert_current_fit <- function(fit, require_integrated = FALSE) {
  expected <- model_input_md5()
  recorded <- fit$huelva$input_md5
  if (is.null(recorded) || !identical(recorded, expected))
    stop("Fit inputs/builder are missing or changed; refit the model before using these results.")
  if (!is.null(fit$fit_file_md5) &&
      (!file.exists(fit$fit_file) || !identical(unname(tools::md5sum(fit$fit_file)), fit$fit_file_md5)))
    stop("The lightweight reference does not match its saved fit; finish/refit F3 before continuing.")
  n_configs <- if (!is.null(fit$misc$configs$config)) length(fit$misc$configs$config) else
    if (!is.null(fit$huelva$n_configs)) fit$huelva$n_configs else 0L
  if (require_integrated && n_configs < 2L)
    stop("This result requires hyperparameter-integrated posterior draws; refit F3 with CCD first.")
  invisible(TRUE)
}

## priors; P(range < 1 km) = 0.05
PRIOR_RANGE <- c(1000, 0.05)
PRIOR_SIGMA <- c(2, 0.01)
PC_PATCHY   <- list(prec = list(prior = "pc.prec", param = c(2, 0.01), initial = 0))
PC_NOISE    <- list(prec = list(prior = "pc.prec", param = c(1, 0.01), initial = 2))
PC_DIL      <- list(prec = list(prior = "pc.prec", param = c(1, 0.01), initial = 3))
U_FREE      <- list(prec = list(initial = -10, fixed = TRUE))   # vague site scores
SMOOTH_ONLY <- list(prec = list(initial = 12, fixed = TRUE))    # patchy part switched off
LOAD_PRIOR  <- c(0, 1)
LOAD_INIT   <- 1
load_hyper  <- function() list(beta = list(prior = "normal", param = LOAD_PRIOR,
                                           initial = LOAD_INIT, fixed = FALSE))

spde <- inla.spde2.pcmatern(geo$mesh, prior.range = PRIOR_RANGE, prior.sigma = PRIOR_SIGMA)


## Model builder --------------------------------------------------------------
## site score per signature: u = g (smooth field) + e (patchy), via 0 = -u + g + e
## element = intercept + loading * u + dilution + noise
##   scale "both"   : g + e (main)
##   scale "smooth" : g only
##   scale "patchy" : e only;  shared = FALSE: independent elements
build_model <- function(Y, scale = c("both", "smooth", "patchy"), shared = TRUE,
                        dilution = shared, element_fields = !shared) {
  scale <- match.arg(scale)
  n <- nrow(Y); n_lik <- n_el + if (shared) n_fac else 0
  blocks <- list(); terms <- c("0", paste0("int_", ELEMENTS))
  if (shared) for (f in names(FACTORS)) {
    terms <- c(terms, sprintf("f(u_%s, model = 'iid', hyper = U_FREE)", f))
    if (scale != "patchy") terms <- c(terms, sprintf("f(g_%s, model = spde)", f))
  }

  ## one block per element
  for (j in seq_len(n_el)) {
    el <- ELEMENTS[j]; f <- FACTOR_OF[[el]]
    ymat <- matrix(NA_real_, n, n_lik); ymat[, j] <- Y[, j]
    eff <- setNames(list(rep(1, n)), paste0("int_", el))
    if (shared) {
      u <- if (el == ANCHORS[[f]]) paste0("u_", f) else paste0("u_", f, "_", el)
      eff[[u]] <- 1:n
      if (el != ANCHORS[[f]]) terms <- c(terms, sprintf(
        "f(%s, copy = 'u_%s', fixed = FALSE, hyper = load_hyper())", u, f))
    }
    if (dilution) eff$d <- 1:n
    A <- list(1); effects <- list(eff)
    if (element_fields) {
      A <- c(A, list(geo$A_obs))
      effects <- c(effects, list(setNames(list(1:spde$n.spde), paste0("w_", el))))
      terms <- c(terms, sprintf("f(w_%s, model = spde)", el))
    }
    blocks[[el]] <- inla.stack(data = list(y = ymat), A = A, effects = effects, tag = el)
  }

  ## one pseudo-likelihood block per signature
  if (shared) for (k in seq_len(n_fac)) {
    f <- names(FACTORS)[k]
    ymat <- matrix(NA_real_, n, n_lik); ymat[, n_el + k] <- 0
    A <- list(sparseMatrix(i = 1:n, j = 1:n, x = -1, dims = c(n, n)))
    effects <- list(setNames(list(1:n), paste0("u_", f)))
    if (scale != "patchy") {
      A <- c(A, list(geo$A_obs))
      effects <- c(effects, list(setNames(list(1:spde$n.spde), paste0("g_", f))))
    }
    blocks[[paste0("pseudo_", f)]] <- inla.stack(data = list(y = ymat), A = A,
                                                 effects = effects, tag = paste0("pseudo_", f))
  }
  if (dilution) terms <- c(terms, "f(d, model = 'iid', hyper = PC_DIL)")

  stk <- do.call(inla.stack, unname(blocks))
  rows <- lapply(setNames(names(blocks), names(blocks)),
                 function(b) inla.stack.index(stk, b)$data)
  list(stack = stk, n_rows = nrow(inla.stack.A(stk)),
       formula = as.formula(paste("y ~", paste(terms, collapse = " + "))),
       family = rep("gaussian", n_lik),
       control.family = c(rep(list(list(hyper = PC_NOISE)), n_el),
                          rep(list(list(hyper = if (scale == "smooth") SMOOTH_ONLY else PC_PATCHY)),
                              if (shared) n_fac else 0)),
       rows = list(element = rows[ELEMENTS], pseudo = rows[setdiff(names(rows), ELEMENTS)]),
       meta = list(scale = scale, shared = shared, dilution = dilution,
                   element_fields = element_fields, factors = FACTORS))
}


## Fitting --------------------------------------------------------------------
## leave-site-out log score
lso_score <- function(fit, rows, n_rows) {
  el_rows <- do.call(cbind, rows$element)              # sites x elements
  groups  <- as.list(seq_len(n_rows))
  for (i in seq_len(nrow(el_rows))) for (r in el_rows[i, ]) groups[[r]] <- el_rows[i, ]
  cv <- inla.group.cv(fit, groups = groups)
  predictive_density <- cv$cv[as.vector(el_rows)]
  if (any(!is.finite(predictive_density)) || any(predictive_density <= 0))
    stop("Group-CV returned invalid predictive densities; the saved fit remains a checkpoint.")
  -mean(log(predictive_density))
}

## EB by default; F3 requests CCD explicitly.
fit_model <- function(m, name = NULL, lso = FALSE, waic = TRUE, config = TRUE,
                      int_strategy = "eb") {
  if (!is.null(name)) message("Fitting ", name, " ...")
  fit <- inla(m$formula, family = m$family,
              data = inla.stack.data(m$stack, spde = spde),
              control.family = m$control.family,
              control.predictor = list(A = inla.stack.A(m$stack)),
              control.compute = list(config = config, waic = waic, dic = waic),
              control.fixed = list(prec = 1e-3),
              control.inla = list(int.strategy = int_strategy))
  el_rows <- unlist(m$rows$element)
  fit$huelva <- list(rows = m$rows, meta = m$meta, input_md5 = MODEL_INPUT_MD5,
                     n_configs = length(fit$misc$configs$config), int_strategy = int_strategy,
                     waic_elements = if (waic) sum(fit$waic$local.waic[el_rows]) else NA,
                     lso = NA_real_)
  ## Preserve an expensive fit even if the separate group-CV calculation fails.
  ## A checkpoint with lso = NA is deliberately rejected by the comparison table.
  if (!is.null(name)) saveRDS(fit, file.path("results/fits", paste0(name, ".rds")))
  if (lso) fit$huelva$lso <- lso_score(fit, m$rows, m$n_rows)
  if (!is.null(name)) {
    fit_file <- file.path("results/fits", paste0(name, ".rds"))
    saveRDS(fit, fit_file)
    ## Simulation workers need only hyperparameter summaries, not CCD configs.
    ## Record the full-fit hash so a stale reference cannot authorize a new checkpoint.
    if (identical(name, "F3")) saveRDS(
      list(summary.hyperpar = fit$summary.hyperpar, huelva = fit$huelva,
           fit_file = fit_file, fit_file_md5 = unname(tools::md5sum(fit_file))),
      "results/fits/F3_reference.rds")
    message(sprintf("  %s: WAIC(elements) %.1f, leave-site-out logS %.3f, %.0f s", name,
                    fit$huelva$waic_elements, fit$huelva$lso, fit$cpu.used[["Total"]]))
  }
  invisible(fit)
}

## INLA name of the precision of likelihood l
lik_prec <- function(l) paste0("Precision for the Gaussian observations",
                               ifelse(l == 1, "", sprintf("[%d]", l)))


## Mappable fraction ----------------------------------------------------------
## share of each signature's variance across sites carried by the smooth part
phi_realized <- function(fit, n = 1000, seed = 2026) {
  sel <- c(setNames(rep(list(seq_len(spde$n.spde)), n_fac), paste0("g_", names(FACTORS))),
           setNames(rep(list(seq_len(n_site)), n_fac), paste0("u_", names(FACTORS))))
  draws <- inla.posterior.sample(n, fit, selection = sel, seed = seed, num.threads = "1:1",
                                 parallel.configs = FALSE)
  lat <- sapply(draws, function(s) s$latent[, 1])
  sapply(names(FACTORS), function(f) {
    g <- as.matrix(geo$A_obs %*% lat[startsWith(rownames(lat), paste0("g_", f, ":")), ])
    u <- lat[startsWith(rownames(lat), paste0("u_", f, ":")), ]
    v_smooth <- apply(g, 2, var)
    v_smooth / (v_smooth + apply(u - g, 2, var))
  })
}

## the same, integrating over the hyperparameters (CCD needs ~8 GB here)
phi_realized_mixture <- function(fit, m, n_theta = 30, n_per = 50, seed = 2026) {
  set.seed(seed)
  theta <- inla.hyperpar.sample(n_theta, fit, intern = TRUE)
  do.call(rbind, lapply(seq_len(n_theta), function(t) {
    ft <- inla(m$formula, family = m$family,
               data = inla.stack.data(m$stack, spde = spde),
               control.family = m$control.family,
               control.predictor = list(A = inla.stack.A(m$stack)),
               control.compute = list(config = TRUE),
               control.fixed = list(prec = 1e-3),
               control.mode = list(theta = theta[t, ], fixed = TRUE))
    phi_realized(ft, n_per, seed = seed + t)
  }))
}

## ratio of marginal variances (secondary)
phi_draws <- function(fit, n = 4000) {
  hs <- inla.hyperpar.sample(n, fit, improve.marginals = TRUE)
  out <- sapply(seq_len(n_fac), function(k) {
    sd_g <- hs[, paste0("Stdev for g_", names(FACTORS)[k])]
    sd_g^2 / (sd_g^2 + 1 / hs[, lik_prec(n_el + k)])
  })
  colnames(out) <- names(FACTORS)
  out
}
