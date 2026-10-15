################################################################################
# Can F3 recover a known mappable fraction? Simulate at the Huelva sites, refit.
# Each replicate uses CCD. SIM_TASK=1..(9*SIM_REPS) runs a single scenario in
# a fresh R process; alternatively SIM_SHARDS=K, SIM_SHARD=1..K partitions them.
# Compatible completed replicates are reused. Workers never summarize; after
# all succeed, run SIM_SUMMARY_ONLY=1 with the same SIM_REPS and no SIM_TASK.
################################################################################
source("R/Robustness/simulation_cache.R")
SIM <- simulation_settings()
source("R/All_Models/model_data.R")

SIM_PHI    <- c(0.25, 0.5, 0.75)
SIM_RANGE  <- c(1500, 3000, 6000)                          # m; F3 estimates 6.6 km
SIM_REPS   <- SIM$reps
SIM_SHARDS <- SIM$shards
SIM_SHARD  <- SIM$shard
SIM_NLATENT <- 2000L
SIM_NMARGINAL <- 4000L

## The reference fit must belong to these inputs, not a previous data/model run.
# This lightweight reference contains hyperparameter summaries and fit
# provenance, avoiding four copies of the CCD latent configuration in RAM.
if (!file.exists("results/fits/F3_reference.rds"))
  stop("Missing F3_reference.rds: finish R/All_Models/model_F3.R before running the simulation.")
reference_fit <- readRDS("results/fits/F3_reference.rds")
assert_current_fit(reference_fit, require_integrated = TRUE)
manifest <- simulation_manifest(
  files = c("results/data/samples.csv", "results/data/mesh.rds", "results/fits/F3_reference.rds",
            "R/All_Models/model_data.R", "R/Robustness/simulation_identifiability.R",
            "R/Robustness/simulation_cache.R"),
  settings = list(phi = SIM_PHI, range = SIM_RANGE, reps = SIM_REPS,
                  int_strategy = "ccd", latent_draws = SIM_NLATENT,
                  marginal_draws = SIM_NMARGINAL,
                  simulation_seed_offset = 1000L),
  versions = list(R = R.version.string, INLA = as.character(packageVersion("INLA")),
                  Matrix = as.character(packageVersion("Matrix")),
                  fmesher = as.character(packageVersion("fmesher"))))
cache <- simulation_cache("results/tables/simulation", manifest, create = !SIM$summary_only)
message("Simulation run ", cache$run_id, " at ", cache$directory)

## other true values: F3 posterior medians
hp_f3 <- reference_fit$summary.hyperpar
rm(reference_fit)
med <- setNames(hp_f3[, "0.5quant"], rownames(hp_f3))
sd_smooth <- med[paste0("Stdev for g_", names(FACTORS))]
sd_patchy <- 1 / sqrt(med[lik_prec(n_el + seq_len(n_fac))])
TRUTH <- list(sd = setNames(sqrt(sd_smooth^2 + sd_patchy^2), names(FACTORS)),
              phi_lithogenic = sd_smooth[[2]]^2 / (sd_smooth[[2]]^2 + sd_patchy[[2]]^2),
              range_lithogenic = med[["Range for g_lithogenic"]],
              sd_dilution = 1 / sqrt(med[["Precision for d"]]),
              sd_noise = setNames(pmax(1 / sqrt(med[lik_prec(seq_len(n_el))]), 0.05), ELEMENTS))
LOADINGS <- setNames(rep(1, n_el), ELEMENTS)
for (el in setdiff(ELEMENTS, ANCHORS))
  LOADINGS[[el]] <- med[[sprintf("Beta for u_%s_%s", FACTOR_OF[[el]], el)]]

## Matern correlation, nu = 1
matern <- function(D, r) {
  k <- sqrt(8) / r
  ifelse(D == 0, 1, k * D * besselK(pmax(k * D, 1e-12), 1))
}
D <- as.matrix(dist(geo$coords_obs))

## simulated data; the true realized fraction is an attribute
simulate_Y <- function(phi, range, seed) {
  set.seed(seed)
  parts <- function(sd, phi, range) {
    g <- drop(crossprod(chol(matern(D, range) + diag(1e-8, n_site)), rnorm(n_site)))
    list(smooth = sd * sqrt(phi) * g, patchy = sd * sqrt(1 - phi) * rnorm(n_site))
  }
  p <- list(sulphide   = parts(TRUTH$sd[["sulphide"]], phi, range),
            lithogenic = parts(TRUTH$sd[["lithogenic"]], TRUTH$phi_lithogenic,
                               TRUTH$range_lithogenic))
  d <- rnorm(n_site, 0, TRUTH$sd_dilution)
  Y <- sapply(ELEMENTS, function(el) {
    q <- p[[FACTOR_OF[[el]]]]
    LOADINGS[[el]] * (q$smooth + q$patchy) + d + rnorm(n_site, 0, TRUTH$sd_noise[[el]])
  })
  attr(Y, "phi_realized") <- sapply(p, function(q) var(q$smooth) / (var(q$smooth) + var(q$patchy)))
  Y
}

scenarios <- expand.grid(phi = SIM_PHI, range = SIM_RANGE, rep = seq_len(SIM_REPS))
for (s in simulation_task_indices(SIM, nrow(scenarios))) {
  sc <- scenarios[s, ]
  file <- simulation_replicate_path(cache, sc)
  if (file.exists(file)) {
    simulation_read_replicate(file, sc, names(FACTORS))
    message("Reusing validated replicate ", basename(file))
    next
  }
  message("Starting CCD scenario ", s, "/", nrow(scenarios), ", shard ", SIM_SHARD,
          "/", SIM_SHARDS, ": ", basename(file))
  t0 <- Sys.time()
  Ysim <- simulate_Y(sc$phi, sc$range, seed = 1000 + s)
  m <- build_model(Ysim)
  # Do not hide failed fits: a worker must exit nonzero so the runner cannot
  # publish a summary from a selectively successful subset of replicates.
  posterior <- simulation_ccd_draws(m, seed = s,
    n_latent = SIM_NLATENT, n_marginal = SIM_NMARGINAL)

  ## realized and marginal fractions
  draws <- posterior$draws
  truth <- list(realized = attr(Ysim, "phi_realized"),
                marginal = c(sulphide = sc$phi, lithogenic = TRUTH$phi_lithogenic))
  out <- do.call(rbind, lapply(names(draws), function(e) data.frame(
    phi = sc$phi, range = sc$range, rep = sc$rep, estimand = e, signature = names(FACTORS),
    true = truth[[e]][names(FACTORS)],
    q05 = apply(draws[[e]], 2, quantile, 0.05), q50 = apply(draws[[e]], 2, quantile, 0.5),
    q95 = apply(draws[[e]], 2, quantile, 0.95))))
  out$range_hat <- posterior$summary.hyperpar[paste0("Range for g_", out$signature), "0.5quant"]
  out$seconds <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  simulation_validate_replicate(out, sc, names(FACTORS), basename(file))
  simulation_atomic_write(file, function(path) write.csv(out, path, row.names = FALSE))
  r <- out[out$estimand == "realized" & out$signature == "sulphide", ]
  message(sprintf("  %d/%d  phi %.2f range %d rep %d: realized sulphide %.2f, estimate %.2f [%.2f, %.2f], %.0f s",
                  s, nrow(scenarios), sc$phi, sc$range, sc$rep, r$true, r$q50, r$q05, r$q95, r$seconds))
  rm(posterior, draws, Ysim, m, truth, out, r)
  gc(verbose = FALSE)
}

## A single explicit aggregation, only after every expected replicate exists.
if (SIM$summary_only) {
  res <- simulation_collect(cache, scenarios, names(FACTORS))
  summ <- simulation_summarize(res, cache$run_id)
  simulation_atomic_write(file.path(cache$directory, "simulation_summary.csv"),
    function(path) write.csv(summ, path, row.names = FALSE), overwrite = TRUE)
  simulation_atomic_write("results/tables/simulation_summary.csv",
    function(path) write.csv(summ, path, row.names = FALSE), overwrite = TRUE)
  print(summ, digits = 2, row.names = FALSE)
} else {
  worker <- if (is.null(SIM$task)) paste0("Shard ", SIM_SHARD, "/", SIM_SHARDS) else
    paste0("Task ", SIM$task, "/", nrow(scenarios))
  message(worker,
          " finished. No shared summary was written; run SIM_SUMMARY_ONLY=1 after all workers succeed.")
}
