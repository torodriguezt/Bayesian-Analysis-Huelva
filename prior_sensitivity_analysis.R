library(scales)

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

formula <- y ~ -1 + b0 + mar + agric + ind + f(s, model = spde)

res <- inla(formula, family = "lognormal", 
            data = inla.stack.data(stk.full),
            control.predictor = list(A = inla.stack.A(stk.full), compute = TRUE),
            control.compute = list(cpo = TRUE, dic = TRUE, waic = TRUE,
                                   return.marginals.predictor = TRUE,
                                   config = TRUE),
            control.fixed=list(mean = 0,prec = 0.001,
                               mean.intercept = 0,prec.intercept = 0.001))

sexta_b3 <- inla.smarginal(res$marginals.fixed$ind)

###Intercepto---------------------------------------------------------------------------------------------

data <- bind_rows(primera = primera_b0, segunda = segunda_b0, tercera = tercera_b0, 
                  cuarta = cuarta_b0, quinta = quinta_b0, sexta = sexta_b0,
                  .id = "grupo")

g <- ggplot(data, aes(x = x, y = y, color = grupo)) +
  geom_line() +
  ggtitle(expression("Prior Analysis for" ~ beta[0])) +
  scale_color_manual(
    values = c("red", "blue", "green", "black", "deeppink", "skyblue"),
    labels = c(bquote(beta[0] ~ .(" ~ N(0, 100)")),
               bquote(beta[0] ~ .(" ~ N(0, 10)")),
               bquote(beta[0] ~ .(" ~ N(0, 1)")),
               bquote(beta[0] ~ .(" ~ N(0, 0.1)")),
               bquote(beta[0] ~ .(" ~ N(0, 0.01)")),
               bquote(beta[0] ~ .(" ~ N(0, 0.001)"))),
    breaks = c("primera", "segunda", "tercera", "cuarta", "quinta", "sexta")
  ) +
  scale_x_continuous(labels = label_number()) +  # Formato normal para el eje x
  scale_y_continuous(labels = label_number()) +  # Formato normal para el eje y
  labs(color = "Distribution") +
  theme_minimal()

# Mostrar el gráfico
print(g)


###Beta1---------------------------------------------------------------------------------------------

data <- bind_rows(primera = primera_b1, segunda = segunda_b1, tercera = tercera_b1, 
                  cuarta = cuarta_b1, quinta = quinta_b1, sexta = sexta_b1,
                  .id = "grupo")

g <- ggplot(data, aes(x = x, y = y, color = grupo)) +
  geom_line() +
  ggtitle(expression("Prior Analysis for" ~ beta[1])) +
  scale_color_manual(
    values = c("red", "blue", "green", "black", "deeppink", "skyblue"),
    labels = c(bquote(beta[1] ~ .(" ~ N(0, 100)")),
               bquote(beta[1] ~ .(" ~ N(0, 10)")),
               bquote(beta[1] ~ .(" ~ N(0, 1)")),
               bquote(beta[1] ~ .(" ~ N(0, 0.1)")),
               bquote(beta[1] ~ .(" ~ N(0, 0.01)")),
               bquote(beta[1] ~ .(" ~ N(0, 0.001)"))),
    breaks = c("primera", "segunda", "tercera", "cuarta", "quinta", "sexta")
  ) +
  scale_x_continuous(labels = label_number()) +  # Formato normal para el eje x
  scale_y_continuous(labels = label_number()) +  # Formato normal para el eje y
  labs(color = "Distribution") +
  theme_minimal()

# Mostrar el gráfico
print(g)


###Beta2---------------------------------------------------------------------------------------------

data <- bind_rows(primera = primera_b2, segunda = segunda_b2, tercera = tercera_b2, 
                  cuarta = cuarta_b2, quinta = quinta_b2, sexta = sexta_b2,
                  .id = "grupo")

g <- ggplot(data, aes(x = x, y = y, color = grupo)) +
  geom_line() +
  ggtitle(expression("Prior Analysis for" ~ beta[2])) +
  scale_color_manual(
    values = c("red", "blue", "green", "black", "deeppink", "skyblue"),
    labels = c(bquote(beta[2] ~ .(" ~ N(0, 100)")),
               bquote(beta[2] ~ .(" ~ N(0, 10)")),
               bquote(beta[2] ~ .(" ~ N(0, 1)")),
               bquote(beta[2] ~ .(" ~ N(0, 0.1)")),
               bquote(beta[2] ~ .(" ~ N(0, 0.01)")),
               bquote(beta[2] ~ .(" ~ N(0, 0.001)"))),
    breaks = c("primera", "segunda", "tercera", "cuarta", "quinta", "sexta")
  ) +
  scale_x_continuous(labels = label_number()) +  # Formato normal para el eje x
  scale_y_continuous(labels = label_number()) +  # Formato normal para el eje y
  labs(color = "Distribution") +
  theme_minimal()

# Mostrar el gráfico
print(g)


###Beta_3----------------------------------------------------------------------------------------

data <- bind_rows(primera = primera_b3, segunda = segunda_b3, tercera = tercera_b3, 
                  cuarta = cuarta_b3, quinta = quinta_b3, sexta = sexta_b3,
                  .id = "grupo")

g <- ggplot(data, aes(x = x, y = y, color = grupo)) +
  geom_line() +
  ggtitle(expression("Prior Analysis for" ~ beta[3])) +
  scale_color_manual(
    values = c("red", "blue", "green", "black", "deeppink", "skyblue"),
    labels = c(bquote(beta[2] ~ .(" ~ N(0, 100)")),
               bquote(beta[2] ~ .(" ~ N(0, 10)")),
               bquote(beta[2] ~ .(" ~ N(0, 1)")),
               bquote(beta[2] ~ .(" ~ N(0, 0.1)")),
               bquote(beta[2] ~ .(" ~ N(0, 0.01)")),
               bquote(beta[2] ~ .(" ~ N(0, 0.001)"))),
    breaks = c("primera", "segunda", "tercera", "cuarta", "quinta", "sexta")
  ) +
  scale_x_continuous(labels = label_number()) +  # Formato normal para el eje x
  scale_y_continuous(labels = label_number()) +  # Formato normal para el eje y
  labs(color = "Distribution") +
  theme_minimal()

# Mostrar el gráfico
print(g)


###Beta_4---------------------------------------------------------------------------------------------


data <- bind_rows(primera = primera_b4, segunda = segunda_b4, tercera = tercera_b4, 
                  cuarta = cuarta_b4, quinta = quinta_b4, sexta = sexta_b4,
                  .id = "grupo")

g <- ggplot(data, aes(x = x, y = y, color = grupo)) +
  geom_line() +
  ggtitle(expression("Prior Analysis for" ~ beta[2])) +
  scale_color_manual(
    values = c("red", "blue", "green", "black", "deeppink", "skyblue"),
    labels = c(bquote(beta[2] ~ .(" ~ N(0, 100)")),
               bquote(beta[2] ~ .(" ~ N(0, 10)")),
               bquote(beta[2] ~ .(" ~ N(0, 1)")),
               bquote(beta[2] ~ .(" ~ N(0, 0.1)")),
               bquote(beta[2] ~ .(" ~ N(0, 0.01)")),
               bquote(beta[2] ~ .(" ~ N(0, 0.001)"))),
    breaks = c("primera", "segunda", "tercera", "cuarta", "quinta", "sexta")
  ) +
  scale_x_continuous(labels = label_number()) +  # Formato normal para el eje x
  scale_y_continuous(labels = label_number()) +  # Formato normal para el eje y
  labs(color = "Distribution") +
  theme_minimal()

# Mostrar el gráfico
print(g)
