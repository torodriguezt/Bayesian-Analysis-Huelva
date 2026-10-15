# Execution, CCD dispatch and cache integrity for the identifiability simulation.

simulation_settings <- function(env = Sys.getenv) {
  integer_setting <- function(name, default) {
    value <- env(name, default)
    if (length(value) != 1L || !grepl("^[0-9]+$", value))
      stop(name, " must be a positive integer.")
    number <- suppressWarnings(as.numeric(value))
    if (!is.finite(number) || number < 1 || number > .Machine$integer.max)
      stop(name, " must be a positive integer within R's integer range.")
    as.integer(number)
  }
  out <- list(reps = integer_setting("SIM_REPS", "10"),
              shards = integer_setting("SIM_SHARDS", "1"),
              shard = integer_setting("SIM_SHARD", "1"))
  # SIM_NTHETA belonged to the former EB mixture approximation. It is ignored
  # by CCD, including when an older runner still exports the variable.
  task <- env("SIM_TASK", "")
  if (length(task) != 1L) stop("SIM_TASK must be empty or a positive scenario index.")
  out$task <- if (identical(task, "")) NULL else integer_setting("SIM_TASK", "")
  if (!is.null(out$task) && out$task > 9 * as.double(out$reps))
    stop("SIM_TASK must be between 1 and 9 * SIM_REPS.")
  if (out$shard > out$shards) stop("SIM_SHARD must be between 1 and SIM_SHARDS.")
  if (!is.null(out$task) && out$shards != 1L)
    stop("SIM_TASK requires SIM_SHARDS=1 and SIM_SHARD=1; do not combine tasks and shards.")
  summary_only <- env("SIM_SUMMARY_ONLY", "0")
  if (length(summary_only) != 1L || !summary_only %in% c("0", "1"))
    stop("SIM_SUMMARY_ONLY must be 0 or 1.")
  out$summary_only <- identical(summary_only, "1")
  if (out$summary_only && !is.null(out$task))
    stop("SIM_SUMMARY_ONLY cannot be combined with SIM_TASK.")
  if (out$summary_only && out$shards != 1L)
    stop("Summary-only mode requires SIM_SHARDS=1 and SIM_SHARD=1; run one aggregator.")
  out
}

simulation_task_indices <- function(settings, n_scenarios) {
  if (settings$summary_only) return(integer())
  if (!is.null(settings$task)) {
    if (settings$task > n_scenarios) stop("SIM_TASK exceeds the number of scenarios.")
    return(settings$task)
  }
  indices <- seq_len(n_scenarios)
  indices[(indices - 1L) %% settings$shards == settings$shard - 1L]
}

# This single dispatch point makes the integration strategy and both posterior
# sample counts explicit and testable without running any INLA fits.
simulation_ccd_draws <- function(model, seed, n_latent = 2000L, n_marginal = 4000L,
                                 fit_fun = fit_model, validate_fun = assert_current_fit,
                                 realized_fun = phi_realized, marginal_fun = phi_draws) {
  fit <- fit_fun(model, waic = FALSE, int_strategy = "ccd")
  validate_fun(fit, require_integrated = TRUE)
  draws <- list(realized = realized_fun(fit, n = n_latent, seed = seed),
                marginal = marginal_fun(fit, n = n_marginal))
  hyperpar <- fit$summary.hyperpar
  rm(fit)
  gc(verbose = FALSE)
  list(summary.hyperpar = hyperpar, draws = draws)
}

simulation_manifest <- function(files, settings, versions) {
  if (!length(files) || anyDuplicated(files) || any(!file.exists(files)))
    stop("Simulation provenance requires every declared input file to exist.")
  hashes <- tools::md5sum(files)
  if (anyNA(hashes)) stop("Could not hash simulation inputs.")
  list(schema = 1L, inputs = setNames(as.character(hashes), files),
       settings = settings, versions = versions)
}

simulation_fingerprint <- function(manifest) {
  path <- tempfile("huelva-simulation-manifest-", fileext = ".rds")
  on.exit(unlink(path), add = TRUE)
  saveRDS(manifest, path, version = 2, compress = FALSE)
  unname(as.character(tools::md5sum(path)))
}

# Same-directory rename exposes only a complete file. Existing completed
# replicates are never overwritten; only the explicit summary may be replaced.
simulation_atomic_write <- function(path, writer, overwrite = FALSE) {
  if (!dir.exists(dirname(path))) stop("Output directory does not exist: ", dirname(path))
  if (!overwrite && file.exists(path)) stop("Refusing to overwrite completed output: ", path)
  temporary <- tempfile(paste0(".", basename(path), "-"), tmpdir = dirname(path))
  on.exit(unlink(temporary), add = TRUE)
  writer(temporary)
  if (!file.exists(temporary)) stop("Output writer did not create its temporary file.")
  if (!overwrite && file.exists(path)) stop("Output appeared during this run: ", path)
  if (!file.rename(temporary, path)) stop("Could not atomically publish output: ", path)
  invisible(path)
}

simulation_cache <- function(root, manifest, create = TRUE) {
  run_id <- simulation_fingerprint(manifest)
  directory <- file.path(root, run_id)
  manifest_path <- file.path(directory, "manifest.rds")
  if (create) {
    dir.create(directory, recursive = TRUE, showWarnings = FALSE)
    if (!dir.exists(directory)) stop("Could not create simulation cache: ", directory)
    if (!file.exists(manifest_path)) {
      # mkdir elects a single initializer. Replacing even identical manifests
      # can interrupt concurrent readers on Windows, so never replace one.
      lock <- file.path(directory, "manifest.lock")
      if (dir.create(lock, showWarnings = FALSE)) {
        tryCatch({
          if (!file.exists(manifest_path)) simulation_atomic_write(manifest_path,
            function(path) saveRDS(manifest, path, version = 2))
        }, finally = unlink(lock, recursive = TRUE))
      } else {
        deadline <- Sys.time() + 30
        while (dir.exists(lock) && Sys.time() < deadline) Sys.sleep(0.05)
        if (dir.exists(lock) || !file.exists(manifest_path))
          stop("Simulation manifest initialization did not finish; inspect lock: ", lock)
      }
    }
  }
  if (!file.exists(manifest_path)) stop("No simulation cache for this configuration: ", directory)
  existing <- tryCatch(readRDS(manifest_path), error = function(e)
    stop("Unreadable simulation manifest: ", manifest_path, ": ", conditionMessage(e)))
  if (!identical(existing, manifest))
    stop("Simulation manifest differs from current inputs; refusing incompatible cache: ", directory)
  list(directory = directory, run_id = run_id, manifest = manifest)
}

simulation_replicate_path <- function(cache, scenario)
  file.path(cache$directory, sprintf("phi%.2f_range%d_rep%02d.csv",
                                    scenario$phi, scenario$range, scenario$rep))

simulation_validate_replicate <- function(result, scenario, signatures, label = "replicate") {
  columns <- c("phi", "range", "rep", "estimand", "signature", "true",
               "q05", "q50", "q95", "range_hat", "seconds")
  expected <- as.vector(outer(c("realized", "marginal"), signatures, paste, sep = ":"))
  fail <- function() stop("Invalid or incomplete simulation ", label, "; do not reuse this file.")
  if (!is.data.frame(result) || !all(columns %in% names(result)) ||
      nrow(result) != length(expected) || anyNA(result[, columns])) fail()
  keys <- paste(result$estimand, result$signature, sep = ":")
  if (anyDuplicated(keys) || !setequal(keys, expected)) fail()
  numeric_columns <- setdiff(columns, c("estimand", "signature"))
  if (!all(vapply(result[, numeric_columns], is.numeric, logical(1))) ||
      any(!is.finite(as.matrix(result[, numeric_columns])))) fail()
  if (any(abs(result$phi - scenario$phi) > 1e-12) ||
      any(result$range != scenario$range) || any(result$rep != scenario$rep) ||
      any(result$true < 0 | result$true > 1) ||
      any(result$q05 < 0 | result$q95 > 1) ||
      any(result$q05 > result$q50 | result$q50 > result$q95) ||
      any(result$range_hat <= 0 | result$seconds < 0)) fail()
  invisible(result)
}

simulation_read_replicate <- function(path, scenario, signatures) {
  result <- tryCatch(read.csv(path), error = function(e)
    stop("Cannot read simulation replicate ", path, ": ", conditionMessage(e)))
  simulation_validate_replicate(result, scenario, signatures, path)
  result
}

simulation_collect <- function(cache, scenarios, signatures) {
  paths <- vapply(seq_len(nrow(scenarios)), function(i)
    simulation_replicate_path(cache, scenarios[i, ]), character(1))
  missing <- paths[!file.exists(paths)]
  if (length(missing)) stop("Simulation is incomplete: ", length(missing), "/", length(paths),
                            " expected replicate(s) missing. First: ", basename(missing[1]))
  do.call(rbind, lapply(seq_along(paths), function(i)
    simulation_read_replicate(paths[i], scenarios[i, ], signatures)))
}

simulation_summarize <- function(result, run_id) {
  result$err <- result$q50 - result$true
  result$covered <- result$q05 <= result$true & result$true <= result$q95
  result$width <- result$q95 - result$q05
  summary <- do.call(rbind, lapply(split(result,
    list(result$estimand, result$signature, result$phi, result$range), drop = TRUE),
    function(d) data.frame(estimand = d$estimand[1], signature = d$signature[1],
      phi = d$phi[1], range = d$range[1], reps = nrow(d), bias = mean(d$err),
      rmse = sqrt(mean(d$err^2)), coverage90 = mean(d$covered),
      width90 = mean(d$width), seconds = mean(d$seconds), run_id = run_id)))
  summary[order(summary$estimand, summary$signature, summary$range, summary$phi), ]
}
