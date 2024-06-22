library(INLA)
library(leaflet)
library(viridis)
library(ggplot2)
library(cowplot)
library(plotly)
library(leaflet)
library(tidyverse)

d <- as.data.frame(coordenadas_muestras)
coo <- cbind(d$X, d$Y)
mesh <- mesh_huelva
spde <- inla.spde2.matern(mesh = mesh, alpha = 2, constr = TRUE)

indexs <- inla.spde.make.index("s", spde$n.spde)
A <- inla.spde.make.A(mesh = mesh, loc = coo)

coop <- cbind(grid_huelva$X,grid_huelva$Y)
Ap <- inla.spde.make.A(mesh = mesh, loc = coop)

stk.e <- inla.stack(tag = "est", data = list(y = datos_huelva$Pb__ppm_)
                    , A = list(1, A),
                    effects = list(data.frame(b0 = 1, 
                                              ind = datos_huelva$Distancia_Industria,
                                              mar = datos_huelva$Distancia_Marsh,
                                              agric = datos_huelva$Distancia_Agric,
                                              refi = datos_huelva$Distancia_Refineria,
                                              urban = datos_huelva$Distancia_Urbano,
                                              phospho = datos_huelva$Distancia_Phospho,
                                              park = datos_huelva$Distancia_Park,
                                              bare = datos_huelva$Distancia_Bare
                    ),
                    s = indexs))

stk.p <- inla.stack(tag = "pred", data = list(y = NA),
                    A = list(1, Ap), 
                    effects = list(data.frame(b0 = 1, 
                                              ind = grid_huelva$Industria,
                                              mar = grid_huelva$Marsh,
                                              agric = grid_huelva$Agric,
                                              refi = grid_huelva$Refineria,
                                              urban = grid_huelva$Urbano,
                                              phospho = grid_huelva$Phospho,
                                              park = grid_huelva$Park,
                                              bare = grid_huelva$Bare
                    ), s = indexs))


stk.full <- inla.stack(stk.e, stk.p)

formula <- y ~ -1 + b0 + mar + agric + f(s, model = spde)

res <- inla(formula, family = "lognormal", 
            data = inla.stack.data(stk.full),
            control.predictor = list(A = inla.stack.A(stk.full), compute = TRUE),
            control.compute = list(cpo = TRUE, dic = TRUE, waic = TRUE,
                                   return.marginals.predictor = TRUE,
                                   config = TRUE),
            control.fixed=list(mean = 0,prec = 0.01,
                               mean.intercept = 0,prec.intercept = 0.01))
summary(res)


index <- inla.stack.index(stack = stk.full, tag = "pred")$data

excprob_pb_ind <- sapply(res$marginals.fitted.values[index],
                         FUN = function(marg){1-inla.pmarginal(q = log(2750), marginal = marg)})

excprob_pb_urb <- sapply(res$marginals.fitted.values[index],
                         FUN = function(marg){1-inla.pmarginal(q = log(275), marginal = marg)})

excprob_pb_agric <- sapply(res$marginals.fitted.values[index],
                           FUN = function(marg){1-inla.pmarginal(q = log(100), marginal = marg)})

excprob_pb_otro <- sapply(res$marginals.fitted.values[index],
                          FUN = function(marg){1-inla.pmarginal(q = log(275), marginal = marg)})

rang <- apply(mesh$loc[, c(1, 2)], 2, range)
proj <- inla.mesh.projector(mesh, xlim = rang[, 1], ylim = rang[, 2], dims = c(500, 500))
mean_s <- inla.mesh.project(proj, res$summary.random$s$mean)
sd_s <- inla.mesh.project(proj, res$summary.random$s$sd)
df <- expand.grid(x = proj$x, y = proj$y)
df$mean_s <- as.vector(mean_s)
df$sd_s <- as.vector(sd_s)

gmean <- ggplot(df, aes(x = x, y = y, fill = mean_s)) + geom_raster() +
  scale_fill_viridis(name = "Mean Spatial Effect", na.value = "transparent", 
                     option = "inferno") + 
  coord_fixed(ratio = 1) + theme_bw()

gsd <- ggplot(df, aes(x = x, y = y, fill = sd_s)) + geom_raster() +
  scale_fill_viridis(name = "SD Spatial Effect", na.value = "transparent",
                     option = "inferno") + 
  coord_fixed(ratio = 1) + theme_bw()

plot_grid(gmean, gsd)

index.pred <- inla.stack.index(stk.full, "pred")$data
grid_huelva$pb.spde <- res$summary.fitted.values[index.pred, "mean"]
grid_huelva$pb.spde.sd <- res$summary.fitted.values[index.pred, "sd"]

p2 <- ggplot(grid_huelva, aes(x = X, y = Y, fill = pb.spde)) +
  geom_tile() +
  scale_fill_viridis(name = "Hg", na.value = 'transparent', option = "D") +
  labs(title = 'Concentración de Plomo (ppm)',
       subtitle = 'INLA') +
  coord_equal() +
  theme(plot.title = element_text(margin = margin(b = 2), size = 12,
                                  hjust = 0.0, color = 'black',
                                  face = "bold"),
        plot.subtitle = element_text(margin = margin(b = 2), size = 10,
                                     hjust = 0.0, color = 'black',
                                     face = "italic"),
        line = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        panel.background = element_blank())

ggplotly(p2)

num_negativos <- grid_huelva %>%
  summarize(n_negativos = sum(pb.spde < 0)) %>%
  pull(n_negativos)

num_negativos

fun2 <- function()
  return(
    b0 +
      mar * grid_huelva$Marsh +
      agric * grid_huelva$Agric +
      Ap %*% s
  )

n_samples <- 5000

samples <- inla.posterior.sample(n_samples, res)

prediccion_posterior_2 <- inla.posterior.sample.eval(fun2, samples)

prediccion_posterior_2[[1]]@x

simulaciones_pb <- list()

for (i in 1:n_samples){
  simulaciones_pb[[i]] <- prediccion_posterior_2[[i]]@x
}

simulaciones_pb <- as.data.frame(simulaciones_pb)

names(simulaciones_pb) <- paste("Simulation", 1:n_samples, sep = " ")

simulaciones_pb$clasificacion <- grid_huelva_clasificacion$clasificacion

simulaciones_pb$mean <- rowMeans(simulaciones_pb[1:5000])

simulaciones_pb$X <- grid_huelva$X
simulaciones_pb$Y <- grid_huelva$Y