
rm(list = ls())
set.seed(123)

suppressPackageStartupMessages({
  library(here)
  library(sf)
  library(dplyr)
  library(tidyr)
  library(Matrix)
  library(fields)       
  library(INLA)
  library(spBayes)
  library(rstan)
  library(ggplot2)
  library(viridis)
  library(cowplot)
  library(plotly)
  library(tibble)
})

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

source(here::here("config", "paths.R"))
check_input_files()
create_output_dirs()

crs_utm <- 25829  

resolve_target_col <- function(nms) {
  cand <- c("Pb_ppm", "Pb__ppm_", "Pb (ppm)")
  hit <- cand[cand %in% nms]
  if (length(hit) == 0) stop("No encuentro columna de Pb (Pb_ppm / Pb__ppm_ / 'Pb (ppm)')")
  hit[1]
}

read_sf_utm <- function(path, epsg) {
  stopifnot(file.exists(path))
  x <- st_read(path, quiet = TRUE)
  if (is.na(st_crs(x))) stop(paste0("CRS faltante en ", path))
  if (st_crs(x)$epsg != epsg) x <- st_transform(x, epsg)
  x
}

min_dist_to <- function(points_sf, polys_sf) {
  as.numeric(st_distance(points_sf, st_union(polys_sf)))
}

fit_scale_log1p <- function(x_train) {
  z <- log1p(x_train)
  mu <- mean(z, na.rm = TRUE)
  sdv <- sd(z, na.rm = TRUE)
  if (!is.finite(sdv) || sdv == 0) sdv <- 1
  list(mu = mu, sd = sdv)
}
apply_scale_log1p <- function(x, scaler) {
  z <- log1p(x)
  as.numeric((z - scaler$mu) / scaler$sd)
}

chol_solve <- function(K, b) {
  U <- chol(K)                 # U' U = K
  y <- forwardsolve(t(U), b)
  x <- backsolve(U, y)
  as.numeric(x)
}

matern32_cov <- function(D, sigma2, rho) {
  a <- sqrt(3) * D / rho
  sigma2 * (1 + a) * exp(-a)
}

summ_vec <- function(v) {
  v <- as.numeric(v)
  v <- v[is.finite(v)]
  c(
    mean = mean(v),
    sd   = sd(v),
    q025 = unname(quantile(v, 0.025)),
    q50  = unname(quantile(v, 0.50)),
    q975 = unname(quantile(v, 0.975))
  )
}

path_bdry <- PATHS$boundaries
path_pts  <- PATHS$samples
path_grid <- PATHS$prediction_grid

path_muestras <- list(
  Agric     = PATHS$land_use_samples$agric,
  Urban     = PATHS$land_use_samples$urban,
  Industria = PATHS$land_use_samples$industria,
  Refineria = PATHS$land_use_samples$refineria,
  Phospho   = PATHS$land_use_samples$phospho,
  Marsh     = PATHS$land_use_samples$marsh,
  Bare      = PATHS$land_use_samples$bare,
  Park      = PATHS$land_use_samples$park
)

path_grid_cov <- list(
  Agric     = PATHS$land_use_covariates$agric,
  Urban     = PATHS$land_use_covariates$urban,
  Industria = PATHS$land_use_covariates$industria,
  Refineria = PATHS$land_use_covariates$refineria,
  Phospho   = PATHS$land_use_covariates$phospho,
  Marsh     = PATHS$land_use_covariates$marsh,
  Bare      = PATHS$land_use_covariates$bare,
  Park      = PATHS$land_use_covariates$park
)

andalucia_huelva <- read_sf_utm(path_bdry, crs_utm)
datos_huelva     <- read_sf_utm(path_pts,  crs_utm) |> st_zm()
grid_huelva      <- read_sf_utm(path_grid, crs_utm)

muestras_cov <- lapply(path_muestras, read_sf_utm, epsg = crs_utm)
grid_cov     <- lapply(path_grid_cov, read_sf_utm, epsg = crs_utm)

coords_obs  <- st_coordinates(datos_huelva)
coords_pred <- st_coordinates(grid_huelva)

for (nm in names(muestras_cov)) {
  datos_huelva[[paste0("Distancia_", nm)]] <- min_dist_to(datos_huelva, muestras_cov[[nm]])
}
for (nm in names(grid_cov)) {
  grid_huelva[[nm]] <- min_dist_to(grid_huelva, grid_cov[[nm]])
}

classes <- c("Industria","Marsh","Phospho","Urban","Agric","Refineria","Bare","Park")

expected_m <- paste0("Distancia_", classes)
present_m  <- intersect(expected_m, names(datos_huelva))
present_g  <- intersect(classes, names(grid_huelva))
if (!length(present_m) || !length(present_g)) stop("No hay columnas de distancia presentes.")

D_m <- as.data.frame(st_drop_geometry(datos_huelva[, present_m, drop = FALSE]))
classes_use_m <- sub("^Distancia_", "", present_m)
colnames(D_m) <- classes_use_m
idx_m <- max.col(-as.matrix(D_m), ties.method = "first")
datos_huelva$clasificacion <- factor(classes_use_m[idx_m], levels = classes)

D_g <- as.data.frame(st_drop_geometry(grid_huelva[, present_g, drop = FALSE]))
classes_use_g <- present_g
idx_g <- max.col(-as.matrix(D_g), ties.method = "first")
grid_huelva$clasificacion <- factor(classes_use_g[idx_g], levels = classes)

classes_cov <- intersect(classes_use_m, classes_use_g)

scalers <- list()
for (nm in classes_cov) {
  x_train <- st_drop_geometry(datos_huelva)[[paste0("Distancia_", nm)]]
  scalers[[nm]] <- fit_scale_log1p(x_train)
}

mk_X <- function(sfobj, is_grid, classes_cov, scalers) {
  df <- st_drop_geometry(sfobj)
  out <- data.frame(Intercept = rep(1, nrow(df)))
  for (nm in classes_cov) {
    colname <- if (is_grid) nm else paste0("Distancia_", nm)
    out[[paste0("d_", nm)]] <- apply_scale_log1p(df[[colname]], scalers[[nm]])
  }
  out
}

X_obs_df  <- mk_X(datos_huelva, is_grid = FALSE, classes_cov = classes_cov, scalers = scalers)
X_pred_df <- mk_X(grid_huelva,  is_grid = TRUE,  classes_cov = classes_cov, scalers = scalers)

X_obs  <- as.matrix(X_obs_df)
X_pred <- as.matrix(X_pred_df)

target_col <- resolve_target_col(names(st_drop_geometry(datos_huelva)))
y_obs <- log1p(st_drop_geometry(datos_huelva)[[target_col]])
stopifnot(all(is.finite(y_obs)))

N <- nrow(X_obs)
P <- ncol(X_obs)
Npred <- nrow(X_pred)

cat("N obs =", N, " | P =", P, " | N pred =", Npred, "\n")
cat("Target:", target_col, " | y = log1p(Pb)\n")
cat("Covariables:", paste(colnames(X_obs), collapse = ", "), "\n")

cat("\n========== INLA/SPDE ==========\n")

huelva.bdry <- andalucia_huelva |>
  as("Spatial") |>
  inla.sp2segment()

mesh_huelva <- inla.mesh.2d(
  boundary = huelva.bdry,
  max.edge = c(1000, 4000),
  offset   = c(1000, 4000),
  cutoff   = 300
)

A_obs  <- inla.spde.make.A(mesh_huelva, loc = coords_obs)
A_pred <- inla.spde.make.A(mesh_huelva, loc = coords_pred)

spde <- inla.spde2.pcmatern(
  mesh = mesh_huelva,
  prior.range = c(5000, 0.9),
  prior.sigma = c(1, 0.01)
)
s_index <- inla.spde.make.index("s", n.spde = spde$n.spde)

stk_obs <- inla.stack(
  data    = list(y = y_obs),
  A       = list(A_obs, 1),
  effects = list(s_index, X_obs_df),
  tag     = "obs"
)

stk_pred <- inla.stack(
  data    = list(y = NA),
  A       = list(A_pred, 1),
  effects = list(s_index, X_pred_df),
  tag     = "pred"
)

stk <- inla.stack(stk_obs, stk_pred)

covar_terms <- colnames(X_obs_df)
rhs <- paste(c("0", covar_terms, "f(s, model = spde)"), collapse = " + ")
formula_inla <- as.formula(paste("y ~", rhs))

t0 <- Sys.time()
fit_inla <- inla(
  formula_inla,
  data = inla.stack.data(stk),
  family = "gaussian",
  control.predictor = list(A = inla.stack.A(stk), compute = TRUE),
  control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE, config = TRUE, return.marginals.predictor = TRUE),
  control.fixed = list(mean = 0, prec = 1e-6)
)
t1 <- Sys.time()

time_inla <- as.numeric(difftime(t1, t0, units = "secs"))
cat("INLA time (s):", time_inla, "\n")

idx_pred_inla <- inla.stack.index(stk, "pred")$data
mu_pred_inla_mean <- fit_inla$summary.fitted.values[idx_pred_inla, "mean"]
mu_pred_inla_sd   <- fit_inla$summary.fitted.values[idx_pred_inla, "sd"]
mu_pred_inla_q025 <- fit_inla$summary.fitted.values[idx_pred_inla, "0.025quant"]
mu_pred_inla_q975 <- fit_inla$summary.fitted.values[idx_pred_inla, "0.975quant"]

pb_pred_inla_mean <- expm1(mu_pred_inla_mean)
pb_pred_inla_q025 <- expm1(mu_pred_inla_q025)
pb_pred_inla_q975 <- expm1(mu_pred_inla_q975)

grid_huelva$mu_inla_mean <- mu_pred_inla_mean
grid_huelva$mu_inla_sd   <- mu_pred_inla_sd
grid_huelva$pb_inla_mean <- pb_pred_inla_mean
grid_huelva$pb_inla_q025 <- pb_pred_inla_q025
grid_huelva$pb_inla_q975 <- pb_pred_inla_q975

cat("\nINLA fixed effects:\n")
print(fit_inla$summary.fixed)

cat("\nINLA hyperparameters:\n")
print(fit_inla$summary.hyperpar)

cat("\n========== MCMC (spBayes::spLM) ==========\n")

D_obs <- as.matrix(dist(coords_obs))
maxD <- max(D_obs)
minD <- quantile(D_obs[D_obs > 0], probs = 0.05)

phi_lower <- 3 / maxD
phi_upper <- 3 / minD
if (!is.finite(phi_upper) || phi_upper <= phi_lower) phi_upper <- phi_lower * 50

p <- P

priors_sp <- list(
  beta.Norm   = list(mean = rep(0, p), cov = diag(100, p)),
  sigma.sq.IG = c(2, 1),
  tau.sq.IG   = c(2, 1),
  phi.Unif    = c(phi_lower, phi_upper),
  nu.Unif     = c(0.5, 2.0)
)

starting_sp <- list(
  beta     = rep(0, p),
  sigma.sq = var(y_obs) * 0.5,
  tau.sq   = var(y_obs) * 0.5,
  phi      = (phi_lower + phi_upper) / 2,
  nu       = 1.0
)

tuning_sp <- list(
  sigma.sq = 0.1,
  tau.sq   = 0.1,
  phi      = 0.1,
  nu       = 0.05
)

n_samples <- 12000
burn_frac <- 0.50
thin      <- 6

dat_sp <- cbind.data.frame(y = y_obs, as.data.frame(X_obs))
rhs_sp <- paste(colnames(X_obs), collapse = " + ")
form_sp <- as.formula(paste0("y ~ 0 + ", rhs_sp))

t0 <- Sys.time()
fit_sp <- spLM(
  formula   = form_sp,
  coords    = coords_obs,
  data      = dat_sp,
  starting  = starting_sp,
  tuning    = tuning_sp,
  priors    = priors_sp,
  cov.model = "matern",
  n.samples = n_samples,
  verbose   = FALSE
)

burn_use <- floor(burn_frac * n_samples)
burn_use <- min(max(1, burn_use), n_samples - 1)
thin_use <- max(1, thin)

fit_sp_rec <- spRecover(
  fit_sp,
  start = burn_use,
  thin  = thin_use,
  get.beta = TRUE,
  get.w    = TRUE,
  verbose  = FALSE
)

pred_sp <- spPredict(
  fit_sp_rec,
  start       = burn_use,
  thin        = thin_use,
  pred.covars = X_pred,
  pred.coords = coords_pred,
  verbose     = FALSE
)
t1 <- Sys.time()

time_mcmc <- as.numeric(difftime(t1, t0, units = "secs"))
cat("MCMC time (s):", time_mcmc, "\n")

y_pred_mcmc_samps <- pred_sp$p.y.predictive.samples
mu_pred_mcmc_mean <- rowMeans(y_pred_mcmc_samps, na.rm = TRUE)
mu_pred_mcmc_sd   <- apply(y_pred_mcmc_samps, 1, sd, na.rm = TRUE)
mu_pred_mcmc_q025 <- apply(y_pred_mcmc_samps, 1, quantile, probs = 0.025, na.rm = TRUE)
mu_pred_mcmc_q975 <- apply(y_pred_mcmc_samps, 1, quantile, probs = 0.975, na.rm = TRUE)

grid_huelva$mu_mcmc_mean <- mu_pred_mcmc_mean
grid_huelva$mu_mcmc_sd   <- mu_pred_mcmc_sd
grid_huelva$pb_mcmc_mean <- expm1(mu_pred_mcmc_mean)
grid_huelva$pb_mcmc_q025 <- expm1(mu_pred_mcmc_q025)
grid_huelva$pb_mcmc_q975 <- expm1(mu_pred_mcmc_q975)

beta_samps  <- fit_sp_rec$p.beta.recover.samples
theta_samps <- fit_sp_rec$p.theta.recover.samples

cat("\nMCMC posterior beta (mean, sd, q025, q975):\n")
print(t(apply(beta_samps, 2, function(v) c(mean=mean(v), sd=sd(v), q025=quantile(v,0.025), q975=quantile(v,0.975)))))

cat("\nMCMC posterior theta (mean, sd, q025, q975):\n")
print(t(apply(theta_samps, 2, function(v) c(mean=mean(v), sd=sd(v), q025=quantile(v,0.025), q975=quantile(v,0.975)))))

cat("\n========== VB (Stan ADVI) ==========\n")

stan_code_vb <- "
data {
  int<lower=1> N;
  int<lower=1> P;
  matrix[N,P] X;
  vector[N] y;
  matrix[N,N] D;
}
parameters {
  vector[P] beta;
  real<lower=0> sigma;   // GP marginal sd
  real<lower=0> rho;     // range
  real<lower=0> tau;     // obs noise sd
  vector[N] w_raw;
}
transformed parameters {
  matrix[N,N] K;
  matrix[N,N] L_K;
  vector[N] w;
  vector[N] mu;

  for (i in 1:N) {
    for (j in i:N) {
      real d = D[i,j];
      real a = sqrt(3.0) * d / rho;
      real k = square(sigma) * (1 + a) * exp(-a);
      K[i,j] = k;
      K[j,i] = k;
    }
    K[i,i] = K[i,i] + 1e-6;
  }
  L_K = cholesky_decompose(K);
  w = L_K * w_raw;
  mu = X * beta + w;
}
model {
  beta  ~ normal(0, 10);
  sigma ~ normal(0, 5);
  tau   ~ normal(0, 5);
  rho   ~ lognormal(log(3000), 0.7);

  w_raw ~ normal(0, 1);
  y ~ normal(mu, tau);
}
"

stan_mod <- rstan::stan_model(model_code = stan_code_vb)

stan_data_vb <- list(
  N = N,
  P = P,
  X = X_obs,
  y = as.vector(y_obs),
  D = D_obs
)

t0 <- Sys.time()
fit_vb <- rstan::vb(
  object = stan_mod,
  data = stan_data_vb,
  algorithm = "fullrank",
  seed = 123,
  output_samples = 4000,
  iter = 40000,
  tol_rel_obj = 0.001
)
t1 <- Sys.time()

time_vb <- as.numeric(difftime(t1, t0, units = "secs"))
cat("VB time (s):", time_vb, "\n")

post <- rstan::extract(fit_vb, pars = c("beta", "sigma", "rho", "tau"))
beta_draws  <- post$beta
sigma_draws <- post$sigma
rho_draws   <- post$rho
tau_draws   <- post$tau
S <- length(sigma_draws)

cat("\nVB posterior summaries (mean, sd, q025, q975):\n")
print(rbind(
  sigma = summ_vec(sigma_draws),
  rho   = summ_vec(rho_draws),
  tau   = summ_vec(tau_draws)
))
cat("\nVB beta posterior (column-wise):\n")
print(t(apply(beta_draws, 2, summ_vec)))

D_pn <- fields::rdist(coords_pred, coords_obs)  # Npred x N

chunk_size <- 5000
chunks <- split(1:Npred, ceiling(seq_len(Npred) / chunk_size))

mu_pred_vb_mean <- numeric(Npred)
mu_pred_vb_sd   <- numeric(Npred)
mu_pred_vb_q025 <- numeric(Npred)
mu_pred_vb_q975 <- numeric(Npred)

for (cc in seq_along(chunks)) {
  idx <- chunks[[cc]]
  D_chunk <- D_pn[idx, , drop = FALSE]

  mu_draw_chunk <- matrix(NA_real_, nrow = S, ncol = length(idx))

  for (s in 1:S) {
    beta_s  <- beta_draws[s, ]
    sigma_s <- sigma_draws[s]
    rho_s   <- rho_draws[s]
    tau_s   <- tau_draws[s]

    K_nn   <- matern32_cov(D_obs, sigma2 = sigma_s^2, rho = rho_s) + diag(1e-6, N)
    K_pn_s <- matern32_cov(D_chunk, sigma2 = sigma_s^2, rho = rho_s)

    Ky <- K_nn + diag(tau_s^2 + 1e-6, N)
    r <- as.numeric(y_obs - X_obs %*% beta_s)

    alpha  <- chol_solve(Ky, r)
    w_pred <- as.numeric(K_pn_s %*% alpha)
    mu_pred <- as.numeric(X_pred[idx, , drop = FALSE] %*% beta_s + w_pred)

    mu_draw_chunk[s, ] <- mu_pred
  }

  mu_pred_vb_mean[idx] <- colMeans(mu_draw_chunk)
  mu_pred_vb_sd[idx]   <- apply(mu_draw_chunk, 2, sd)
  mu_pred_vb_q025[idx] <- apply(mu_draw_chunk, 2, quantile, probs = 0.025)
  mu_pred_vb_q975[idx] <- apply(mu_draw_chunk, 2, quantile, probs = 0.975)

  cat(sprintf("VB pred chunk %d/%d listo\n", cc, length(chunks)))
}

grid_huelva$mu_vb_mean <- mu_pred_vb_mean
grid_huelva$mu_vb_sd   <- mu_pred_vb_sd
grid_huelva$pb_vb_mean <- expm1(mu_pred_vb_mean)
grid_huelva$pb_vb_q025 <- expm1(mu_pred_vb_q025)
grid_huelva$pb_vb_q975 <- expm1(mu_pred_vb_q975)

summary_times <- data.frame(
  metodo = c("INLA", "MCMC", "VB"),
  tiempo_s = c(time_inla, time_mcmc, time_vb)
)
print(summary_times)

make_map <- function(sf_grid, var, title, subtitle) {
  ggplot(sf_grid, aes(x = X, y = Y, fill = .data[[var]])) +
    geom_tile() +
    scale_fill_viridis(na.value = "transparent", option = "A") +
    labs(title = title, subtitle = subtitle) +
    coord_equal() +
    theme(
      plot.title = element_text(face = "bold"),
      plot.subtitle = element_text(face = "italic"),
      axis.text = element_blank(),
      axis.title = element_blank(),
      panel.background = element_blank()
    )
}

if (!("X" %in% names(grid_huelva)) || !("Y" %in% names(grid_huelva))) {
  xy <- st_coordinates(grid_huelva)
  grid_huelva$X <- xy[,1]
  grid_huelva$Y <- xy[,2]
}

p_inla <- make_map(grid_huelva, "pb_inla_mean", "Pb (ppm) – Predicción", "INLA/SPDE")
p_mcmc <- make_map(grid_huelva, "pb_mcmc_mean", "Pb (ppm) – Predicción", "MCMC (spBayes)")
p_vb   <- make_map(grid_huelva, "pb_vb_mean",   "Pb (ppm) – Predicción", "VB (Stan)")

plot_grid(p_inla, p_mcmc, p_vb, ncol = 3)

cat("\nListo: grid_huelva contiene predicciones por método (mu_* y pb_*).\n")
cat("mu_* está en escala log (log1p(Pb)); pb_* está en escala Pb (ppm).\n")


df_inla_fixed <- fit_inla$summary.fixed %>%
  as.data.frame() %>%
  tibble::rownames_to_column("param") %>%
  transmute(
    metodo = "INLA",
    param,
    mean  = mean,
    sd    = sd,
    q025  = `0.025quant`,
    q50   = `0.5quant`,
    q975  = `0.975quant`
  )

beta_mcmc <- fit_sp_rec$p.beta.recover.samples
colnames(beta_mcmc) <- colnames(X_obs)

df_mcmc_fixed <- t(apply(beta_mcmc, 2, summ_vec)) %>%
  as.data.frame() %>%
  tibble::rownames_to_column("param") %>%
  mutate(metodo = "MCMC") %>%
  select(metodo, param, everything())

post_vb_beta <- rstan::extract(fit_vb, pars = "beta")
beta_vb <- post_vb_beta$beta
colnames(beta_vb) <- colnames(X_obs)

df_vb_fixed <- t(apply(beta_vb, 2, summ_vec)) %>%
  as.data.frame() %>%
  tibble::rownames_to_column("param") %>%
  mutate(metodo = "VB") %>%
  select(metodo, param, everything())

post_fixed_all <- bind_rows(df_inla_fixed, df_mcmc_fixed, df_vb_fixed) %>%
  arrange(param, metodo)

post_fixed_wide <- post_fixed_all %>%
  pivot_wider(
    names_from  = metodo,
    values_from = c(mean, sd, q025, q50, q975),
    names_sep   = "__"
  ) %>%
  arrange(param)

print(post_fixed_all, n = 200)
print(post_fixed_wide, n = 200)


inla_hyp <- fit_inla$summary.hyperpar %>%
  as.data.frame() %>%
  tibble::rownames_to_column("param") %>%
  transmute(
    metodo = "INLA",
    param,
    mean = mean,
    sd   = sd,
    q025 = `0.025quant`,
    q50  = `0.5quant`,
    q975 = `0.975quant`
  )

inla_hyp$tipo <- ifelse(
  grepl("Range|Stdev|Sigma|kappa|tau|s$", inla_hyp$param, ignore.case = TRUE),
  "spatial", "other"
)

theta_mcmc <- fit_sp_rec$p.theta.recover.samples
stopifnot(!is.null(theta_mcmc))

if (is.null(colnames(theta_mcmc))) {
  colnames(theta_mcmc) <- c("sigma.sq", "tau.sq", "phi", "nu")[1:ncol(theta_mcmc)]
}

df_mcmc_theta <- t(apply(theta_mcmc, 2, summ_vec)) %>%
  as.data.frame() %>%
  tibble::rownames_to_column("param") %>%
  mutate(metodo = "MCMC", tipo = "spatial") %>%
  select(metodo, tipo, param, mean, sd, q025, q50, q975)

post_vb_h <- rstan::extract(fit_vb, pars = c("sigma","rho","tau"))

df_vb_hyp <- bind_rows(
  data.frame(t(summ_vec(post_vb_h$sigma)), metodo="VB", tipo="spatial", param="sigma"),
  data.frame(t(summ_vec(post_vb_h$rho)),   metodo="VB", tipo="spatial", param="rho"),
  data.frame(t(summ_vec(post_vb_h$tau)),   metodo="VB", tipo="spatial", param="tau")
) %>%
  select(metodo, tipo, param, mean, sd, q025, q50, q975)

post_hyper_all <- bind_rows(inla_hyp, df_mcmc_theta, df_vb_hyp) %>%
  arrange(param, metodo)

post_hyper_wide <- post_hyper_all %>%
  select(metodo, param, mean, sd, q025, q50, q975) %>%
  pivot_wider(
    names_from  = metodo,
    values_from = c(mean, sd, q025, q50, q975),
    names_sep   = "__"
  ) %>%
  arrange(param)

print(post_hyper_all, n = 200)
print(post_hyper_wide, n = 200)


make_map_points_log <- function(sf_grid, var, title, subtitle, legend_label) {
  ggplot(sf_grid) +
    geom_sf(aes(color = .data[[var]]), size = 0.6) +
    scale_color_viridis(option = "A", na.value = "transparent") +
    labs(
      title = title,
      subtitle = subtitle,
      color = legend_label
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold"),
      plot.subtitle = element_text(face = "italic"),
      axis.title = element_blank(),
      axis.text  = element_blank(),
      panel.grid = element_blank()
    )
}

p_mean_log <- plot_grid(
  make_map_points_log(grid_huelva, "mu_inla_mean",
                      "Pb – Posterior mean (log scale)", "INLA / SPDE", "Mean log1p(Pb)"),
  make_map_points_log(grid_huelva, "mu_mcmc_mean",
                      "Pb – Posterior mean (log scale)", "MCMC (spBayes)", "Mean log1p(Pb)"),
  make_map_points_log(grid_huelva, "mu_vb_mean",
                      "Pb – Posterior mean (log scale)", "VB (Stan)", "Mean log1p(Pb)"),
  ncol = 3
)
print(p_mean_log)

ggsave(
  filename = "pb_posterior_mean_logscale.png",
  plot     = p_mean_log,
  width    = 14,
  height   = 5,
  dpi      = 300
)

p_sd_log <- plot_grid(
  make_map_points_log(grid_huelva, "mu_inla_sd",
                      "Pb – Posterior SD (log scale)", "INLA / SPDE", "SD log1p(Pb)"),
  make_map_points_log(grid_huelva, "mu_mcmc_sd",
                      "Pb – Posterior SD (log scale)", "MCMC (spBayes)", "SD log1p(Pb)"),
  make_map_points_log(grid_huelva, "mu_vb_sd",
                      "Pb – Posterior SD (log scale)", "VB (Stan)", "SD log1p(Pb)"),
  ncol = 3
)
print(p_sd_log)

ggsave(
  filename = "pb_posterior_sd_logscale.png",
  plot     = p_sd_log,
  width    = 14,
  height   = 5,
  dpi      = 300
)
