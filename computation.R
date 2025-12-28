
rm(list = ls())
set.seed(123)

library(spBayes)
library(INLA)
library(dplyr)
library(ggplot2)
library(MASS)
library(Matrix)
library(rstan)
library(writexl)
library(readxl)

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())


matern32_cov <- function(D, sigma2, rho) {
  a <- sqrt(3) * D / rho
  sigma2 * (1 + a) * exp(-a)
}

exp_cov <- function(D, sigma2, rho) {
  sigma2 * exp(-D / rho)
}

clip01 <- function(x, lo = 0, hi = 10000) pmin(hi, pmax(lo, x))

gen_coords <- function(n, design = c("uniform", "band"), domain = c(0, 10000)) {
  design <- match.arg(design)
  lo <- domain[1]; hi <- domain[2]

  if (design == "uniform") {
    coords <- cbind(runif(n, lo, hi), runif(n, lo, hi))
    return(coords)
  }

  n_band <- round(0.7 * n)
  n_un   <- n - n_band

  x_band <- runif(n_band, lo, hi)
  sd_band <- 400
  y_band <- x_band + rnorm(n_band, 0, sd_band)

  coords_band <- cbind(x_band, clip01(y_band, lo, hi))
  coords_un <- cbind(runif(n_un, lo, hi), runif(n_un, lo, hi))

  rbind(coords_band, coords_un)
}

train_test_split <- function(n, test_frac = 0.2) {
  n_test <- max(1, floor(test_frac * n))
  idx_test <- sample.int(n, n_test)
  idx_train <- setdiff(seq_len(n), idx_test)
  list(train = idx_train, test = idx_test)
}

scale_train_apply <- function(x_train, x_test) {
  mu <- mean(x_train); sdv <- sd(x_train)
  if (!is.finite(sdv) || sdv == 0) sdv <- 1
  list(train = (x_train - mu) / sdv, test = (x_test - mu) / sdv, mu = mu, sd = sdv)
}

chol_solve <- function(K, b) {
  U <- chol(K)
  y <- forwardsolve(t(U), b)
  x <- backsolve(U, y)
  as.numeric(x)
}


escenarios <- list(
  spatial_wins = list(
    n = 150,
    sigma_eps = 0.3, sigma_w = 1.5, rho_w = 3500,
    beta = c(1.0, 0.6, -0.4),
    design = "uniform",
    misspec = "none",
    n_iter = 50,
    mcmc = list(n.samples = 8000, burn_frac = 0.40, thin = 4),
    inla = list(max_edge = 1800, cutoff = 300),
    vb   = list(output_samples = 1000)
  ),
  nugget_short_band = list(
    n = 300,
    sigma_eps = 1.2, sigma_w = 0.7, rho_w = 800,
    beta = c(1.0, 0.4, -0.2),
    design = "band",
    misspec = "none",
    n_iter = 50,
    mcmc = list(n.samples = 10000, burn_frac = 0.50, thin = 5),
    inla = list(max_edge = 700, cutoff = 150),
    vb   = list(output_samples = 800)
  ),
  misspecified = list(
    n = 350,
    sigma_eps = 0.8, sigma_w = 1.0, rho_w = 2000,
    beta = c(0.7, 0.3, -0.25),
    design = "uniform",
    misspec = "rho_by_region",  # "rho_by_region" o "exp_kernel"
    n_iter = 50,
    mcmc = list(n.samples = 12000, burn_frac = 0.50, thin = 6),
    inla = list(max_edge = 1000, cutoff = 200),
    vb   = list(output_samples = 600)
  )
)


stan_code <- "
data {
  int<lower=1> N;
  int<lower=1> P;
  matrix[N,P] X;
  vector[N] y;
  matrix[N,N] D;
}
parameters {
  vector[P] beta;
  real<lower=0> sigma;
  real<lower=0> rho;
  real<lower=0> tau;
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
  beta ~ normal(0, 10);
  sigma ~ normal(0, 5);
  tau   ~ normal(0, 5);
  rho   ~ lognormal(log(3000), 0.7);

  w_raw ~ normal(0, 1);
  y ~ normal(mu, tau);
}
"

stan_mod <- rstan::stan_model(model_code = stan_code)


results <- data.frame(
  escenario = character(),
  metodo    = character(),
  MAE       = numeric(),
  tiempo    = numeric(),
  stringsAsFactors = FALSE
)


for (esc in names(escenarios)) {
  pars <- escenarios[[esc]]

  for (it in 1:pars$n_iter) {
    cat("Scenario:", esc, " Iteration:", it, "\n")
    set.seed(10000 + it)

    n <- pars$n
    coords_all <- gen_coords(n, design = pars$design, domain = c(0, 10000))

    z1_raw <- runif(n, 0, 10)
    z2_raw <- runif(n, 0, 5)

    split <- train_test_split(n, test_frac = 0.2)
    idx_tr <- split$train
    idx_te <- split$test

    s1 <- scale_train_apply(z1_raw[idx_tr], z1_raw[idx_te])
    s2 <- scale_train_apply(z2_raw[idx_tr], z2_raw[idx_te])

    X_tr <- cbind(Intercept = 1, z1 = as.numeric(s1$train), z2 = as.numeric(s2$train))
    X_te <- cbind(Intercept = 1, z1 = as.numeric(s1$test),  z2 = as.numeric(s2$test))

    coords_tr <- coords_all[idx_tr, , drop = FALSE]
    coords_te <- coords_all[idx_te, , drop = FALSE]

    D_all <- as.matrix(dist(coords_all))

    if (is.null(pars$misspec) || pars$misspec == "none") {
      Sigma_w_all <- matern32_cov(D_all, sigma2 = pars$sigma_w^2, rho = pars$rho_w) + diag(1e-6, n)
      w_all <- as.numeric(MASS::mvrnorm(1, mu = rep(0, n), Sigma = Sigma_w_all))

    } else if (pars$misspec == "exp_kernel") {
      Sigma_w_all <- exp_cov(D_all, sigma2 = pars$sigma_w^2, rho = pars$rho_w) + diag(1e-6, n)
      w_all <- as.numeric(MASS::mvrnorm(1, mu = rep(0, n), Sigma = Sigma_w_all))

    } else if (pars$misspec == "rho_by_region") {
      left <- coords_all[,1] < 5000
      rho1 <- max(500, round(1.4 * pars$rho_w))
      rho2 <- max(300, round(0.7 * pars$rho_w))

      Sigma1 <- matern32_cov(D_all, sigma2 = pars$sigma_w^2, rho = rho1) + diag(1e-6, n)
      Sigma2 <- matern32_cov(D_all, sigma2 = pars$sigma_w^2, rho = rho2) + diag(1e-6, n)

      w1 <- as.numeric(MASS::mvrnorm(1, mu = rep(0, n), Sigma = Sigma1))
      w2 <- as.numeric(MASS::mvrnorm(1, mu = rep(0, n), Sigma = Sigma2))
      w_all <- ifelse(left, w1, w2)

    } else {
      stop("misspec no reconocida: ", pars$misspec)
    }

    z1_all_sc <- (z1_raw - s1$mu) / s1$sd
    z2_all_sc <- (z2_raw - s2$mu) / s2$sd
    X_all <- cbind(Intercept = 1, z1 = z1_all_sc, z2 = z2_all_sc)

    y_all <- as.numeric(X_all %*% pars$beta + w_all + rnorm(n, 0, pars$sigma_eps))

    y_tr <- y_all[idx_tr]
    y_te <- y_all[idx_te]

    dat_tr <- data.frame(
      y = y_tr,
      Intercept = X_tr[, "Intercept"],
      z1 = X_tr[, "z1"],
      z2 = X_tr[, "z2"]
    )

    p <- ncol(X_tr)

    priors <- list(
      beta.Norm   = list(mean = rep(0, p), cov = diag(100, p)),
      sigma.sq.IG = c(2, 1),
      tau.sq.IG   = c(2, 1),
      phi.Unif    = c(0.0005, 0.005),
      nu.Unif     = c(0.5, 2.0)
    )

    starting <- list(
      beta     = rep(0, p),
      sigma.sq = var(y_tr) * 0.5,
      tau.sq   = var(y_tr) * 0.5,
      phi      = 0.001,
      nu       = 1
    )

    tuning <- list(sigma.sq = 0.1, tau.sq = 0.1, phi = 0.1, nu = 0.1)

    t0 <- Sys.time()
    m_fit <- spLM(
      formula   = y ~ Intercept + z1 + z2 - 1,
      coords    = coords_tr,
      data      = dat_tr,
      starting  = starting,
      tuning    = tuning,
      priors    = priors,
      cov.model = "matern",
      n.samples = pars$mcmc$n.samples,
      verbose   = FALSE
    )

    burn_in_use <- floor(pars$mcmc$burn_frac * pars$mcmc$n.samples)
    burn_use <- min(max(1, burn_in_use), m_fit$n.samples - 1)
    thin_use <- max(1, pars$mcmc$thin)

    m_fit <- spRecover(m_fit, start = burn_use, thin = thin_use, get.beta = TRUE, get.w = TRUE, verbose = FALSE)

    m_pred <- spPredict(
      m_fit,
      start       = burn_use,
      thin        = thin_use,
      pred.covars = X_te,
      pred.coords = coords_te,
      verbose     = FALSE
    )

    yhat_mcmc <- rowMeans(m_pred$p.y.predictive.samples, na.rm = TRUE)
    t1 <- Sys.time()

    results <- rbind(results, data.frame(
      escenario = esc,
      metodo    = "MCMC",
      MAE       = mean(abs(y_te - yhat_mcmc)),
      tiempo    = as.numeric(difftime(t1, t0, units = "secs"))
    ))

    mesh <- inla.mesh.2d(
      loc = rbind(coords_tr, coords_te),
      max.edge = pars$inla$max_edge,
      cutoff = pars$inla$cutoff
    )
    spde <- inla.spde2.matern(mesh, alpha = 2)

    A_tr <- inla.spde.make.A(spde$mesh, loc = coords_tr)
    A_te <- inla.spde.make.A(spde$mesh, loc = coords_te)

    df_tr <- data.frame(y = y_tr, Intercept = 1, z1 = X_tr[, "z1"], z2 = X_tr[, "z2"])
    df_te <- data.frame(y = NA,  Intercept = 1, z1 = X_te[, "z1"], z2 = X_te[, "z2"])

    stk_tr <- inla.stack(
      data    = list(y = df_tr$y),
      A       = list(1, A_tr),
      effects = list(df_tr[, c("Intercept","z1","z2")], s = 1:spde$n.spde),
      tag     = "train"
    )

    stk_te <- inla.stack(
      data    = list(y = df_te$y),
      A       = list(1, A_te),
      effects = list(df_te[, c("Intercept","z1","z2")], s = 1:spde$n.spde),
      tag     = "test"
    )

    stk <- inla.stack(stk_tr, stk_te)
    formula_spde <- y ~ 0 + Intercept + z1 + z2 + f(s, model = spde)

    t0 <- Sys.time()
    fit_inla <- inla(
      formula_spde,
      data = inla.stack.data(stk),
      family = "gaussian",
      control.predictor = list(A = inla.stack.A(stk), compute = TRUE),
      control.compute   = list(config = TRUE)
    )
    t1 <- Sys.time()

    idx_te_inla <- inla.stack.index(stk, "test")$data
    yhat_inla <- fit_inla$summary.fitted.values[idx_te_inla, "mean"]

    results <- rbind(results, data.frame(
      escenario = esc,
      metodo    = "INLA",
      MAE       = mean(abs(y_te - yhat_inla)),
      tiempo    = as.numeric(difftime(t1, t0, units = "secs"))
    ))

    D_tr <- as.matrix(dist(coords_tr))
    stan_data <- list(
      N = nrow(coords_tr),
      P = p,
      X = X_tr,
      y = as.vector(y_tr),
      D = D_tr
    )

    t0 <- Sys.time()
    vb_fit <- rstan::vb(
      object = stan_mod,
      data = stan_data,
      algorithm = "meanfield",
      seed = 2000 + it,
      output_samples = pars$vb$output_samples,
      iter = 20000,
      tol_rel_obj = 0.01
    )
    t1 <- Sys.time()

    post <- rstan::extract(vb_fit, pars = c("beta", "sigma", "rho", "tau"))
    beta_hat  <- colMeans(post$beta)
    sigma_hat <- mean(post$sigma)
    rho_hat   <- mean(post$rho)
    tau_hat   <- mean(post$tau)

    n_tr <- nrow(coords_tr)
    n_te <- nrow(coords_te)

    D_join <- as.matrix(dist(rbind(coords_te, coords_tr)))
    D_tn  <- D_join[1:n_te, (n_te + 1):(n_te + n_tr)]

    K_nn <- matern32_cov(D_tr, sigma2 = sigma_hat^2, rho = rho_hat)
    K_tn <- matern32_cov(D_tn, sigma2 = sigma_hat^2, rho = rho_hat)

    K_nn_y <- K_nn + diag(tau_hat^2 + 1e-6, n_tr)
    r_tr <- as.numeric(y_tr - X_tr %*% beta_hat)

    alpha <- chol_solve(K_nn_y, r_tr)
    w_te_hat <- as.numeric(K_tn %*% alpha)

    yhat_vb <- as.numeric(X_te %*% beta_hat + w_te_hat)

    results <- rbind(results, data.frame(
      escenario = esc,
      metodo    = "VB",
      MAE       = mean(abs(y_te - yhat_vb)),
      tiempo    = as.numeric(difftime(t1, t0, units = "secs"))
    ))
  }
}


results <- results %>%
  mutate(
    escenario = case_when(
      escenario == "spatial_wins" ~ "Spatial wins",
      escenario == "nugget_short_band" ~ "Nugget + short-range + band",
      escenario == "misspecified" ~ "Mis-specified",
      TRUE ~ escenario
    ),
    escenario = factor(escenario, levels = c("Spatial wins", "Nugget + short-range + band", "Mis-specified")),
    metodo = factor(metodo, levels = c("MCMC", "INLA", "VB"))
  )

plot_mae <- ggplot(results, aes(x = metodo, y = MAE, fill = metodo)) +
  geom_boxplot(alpha = 0.7, outlier.shape = 21, outlier.fill = "red") +
  facet_wrap(~ escenario, scales = "free_y", ncol = 3) +
  labs(
    title = "Out-of-sample MAE by method across scenarios",
    x = "Method",
    y = "MAE (test)"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title   = element_text(face = "bold", hjust = 0.5),
    legend.position = "none",
    strip.text   = element_text(face = "bold")
  )

plot_time <- ggplot(results, aes(x = metodo, y = tiempo, fill = metodo)) +
  geom_boxplot(alpha = 0.7, outlier.shape = 21, outlier.fill = "red") +
  facet_wrap(~ escenario, scales = "free_y", ncol = 3) +
  labs(
    title = "Computation time by method across scenarios",
    x = "Method",
    y = "Time (seconds)"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title   = element_text(face = "bold", hjust = 0.5),
    legend.position = "none",
    strip.text   = element_text(face = "bold")
  )

ggsave("plot_mae.pdf",  plot_mae,  width = 12, height = 6, units = "in")
ggsave("plot_time.pdf", plot_time, width = 12, height = 6, units = "in")

writexl::write_xlsx(results, "mi_lista_vb.xlsx")


df <- read_excel("mi_lista_vb.xlsx")

df <- df %>%
  mutate(
    escenario = factor(
      escenario,
      levels = c(
        "Spatial wins",
        "Mis-specified",
        "Nugget + short-range + band"
      ),
      labels = c(
        "Escenario 1",
        "Escenario 2",
        "Escenario 3"
      )
    ),
    metodo = factor(metodo, levels = c("MCMC", "INLA", "VB"))
  )


library(ggplot2)
library(dplyr)

set.seed(123)

clip01 <- function(x, lo = 0, hi = 10000) pmin(hi, pmax(lo, x))

gen_coords <- function(n, design = c("uniform", "band"), domain = c(0, 10000)) {
  design <- match.arg(design)
  lo <- domain[1]; hi <- domain[2]

  if (design == "uniform") {
    return(data.frame(
      x = runif(n, lo, hi),
      y = runif(n, lo, hi)
    ))
  }

  n_band <- round(0.7 * n)
  n_un   <- n - n_band

  x_band <- runif(n_band, lo, hi)
  sd_band <- 400
  y_band <- clip01(x_band + rnorm(n_band, 0, sd_band), lo, hi)

  coords_band <- data.frame(x = x_band, y = y_band)
  coords_un   <- data.frame(
    x = runif(n_un, lo, hi),
    y = runif(n_un, lo, hi)
  )

  rbind(coords_band, coords_un)
}

coords_s1 <- gen_coords(150, design = "uniform") %>%
  mutate(Scenario = "Scenario 1")

coords_s2 <- gen_coords(300, design = "band") %>%
  mutate(Scenario = "Scenario 2")

coords_s3 <- gen_coords(350, design = "uniform") %>%
  mutate(Scenario = "Scenario 3")

coords_all <- bind_rows(coords_s1, coords_s2, coords_s3)

regions_plot <- ggplot(coords_all, aes(x = x, y = y)) +
  geom_point(alpha = 0.7, size = 1) +
  facet_wrap(~ Scenario, ncol = 3) +
  coord_equal() +
  labs(
    title = "Spatial sampling designs across simulation scenarios",
    x = "x-coordinate",
    y = "y-coordinate"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    strip.text = element_text(face = "bold"),
    plot.title = element_text(face = "bold", hjust = 0.5)
  )

ggsave("plot_regions.pdf",  regions_plot,  width = 12, height = 6, units = "in")


plot_mae <- ggplot(df, aes(x = metodo, y = MAE, fill = metodo)) +
  geom_boxplot(alpha = 0.7, outlier.shape = 21, outlier.fill = "red") +
  facet_wrap(~ escenario, scales = "free_y", ncol = 3) +
  labs(
    title = "Out-of-sample MAE by method across scenarios",
    x = "Method",
    y = "MAE (test)"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title   = element_text(face = "bold", hjust = 0.5),
    legend.position = "none",
    strip.text   = element_text(face = "bold")
  )

plot_time <- ggplot(df, aes(x = metodo, y = tiempo, fill = metodo)) +
  geom_boxplot(alpha = 0.7, outlier.shape = 21, outlier.fill = "red") +
  facet_wrap(~ escenario, scales = "free_y", ncol = 3) +
  labs(
    title = "Computation time by method across scenarios",
    x = "Method",
    y = "Time (seconds)"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title   = element_text(face = "bold", hjust = 0.5),
    legend.position = "none",
    strip.text   = element_text(face = "bold")
  )

ggsave("plot_mae.pdf",  plot_mae,  width = 12, height = 6, units = "in")
ggsave("plot_time.pdf", plot_time, width = 12, height = 6, units = "in")
