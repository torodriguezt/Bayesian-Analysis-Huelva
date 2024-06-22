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
spde <- inla.spde2.pcmatern(
  mesh = mesh_huelva,
  prior.range = c(1, 0.2),
  prior.sigma = c(7, 0.75))
spde <- inla.spde2.matern(mesh = mesh, alpha = 2, constr = TRUE)
indexs <- inla.spde.make.index("s", spde$n.spde)
A <- inla.spde.make.A(mesh = mesh, loc = coo)

coop <- cbind(grid_huelva$X,grid_huelva$Y)
Ap <- inla.spde.make.A(mesh = mesh, loc = coop)

stk.e <- inla.stack(tag = "est", data = list(y = datos_huelva$Hg__ppb_)
                    , A = list(1, A),
                    effects = list(data.frame(b0 = 1, 
                                              ind = datos_huelva$Distancia_Industria,
                                              mar = datos_huelva$Distancia_Marsh,
                                              agric = datos_huelva$Distancia_Agric,
                                              refi = datos_huelva$Distancia_Refineria,
                                              urban = datos_huelva$Distancia_Urbano,
                                              phospho = datos_huelva$Distancia_Phospho,
                                              park = datos_huelva$Distancia_Park,
                                              bare = datos_huelva$Distancia_Bare),
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

formula <- y ~ -1 + b0 + ind + bare + f(s, model = spde)

res <- inla(formula, family = "normal", 
            data = inla.stack.data(stk.full),
            control.predictor = list(A = inla.stack.A(stk.full), compute = TRUE),
            control.compute = list(cpo = TRUE, dic = TRUE, waic = TRUE,
                                   return.marginals.predictor = TRUE,
                                   config = TRUE),
            control.fixed=list(mean = 0,prec = 0.01,
                               mean.intercept = 0,prec.intercept = 0.01))
summary(res)


r <- inla.group.cv(res, num.level.sets = 1)

-mean(log(r$cv[1:138]))