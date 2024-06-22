res <- inla(formula, family = "lognormal", 
            data = inla.stack.data(stk.full),
            control.predictor = list(A = inla.stack.A(stk.full), compute = TRUE),
            control.compute = list(cpo = TRUE, dic = TRUE, waic = TRUE,
                                   return.marginals.predictor = TRUE,
                                   config = TRUE),
            control.fixed=list(mean = 0,prec = 0.01,
                               mean.intercept = 0,prec.intercept = 0.01))


sim <- inla.posterior.sample(500, res)

df <- data.frame(matrix(ncol = 0, nrow = 46722))

for (i in 1:length(sim)) {
  df <- cbind(df, sim[[i]]$latent)
}

df <- exp(df)

df$media <- rowMeans(df)

data_selected <- df[1:138,]

# Calcular la densidad para el conjunto de datos principal de Huelva
dens_huelva <- density(datos_huelva$Cu__ppm_)
df_huelva <- data.frame(x = dens_huelva$x, y = dens_huelva$y, group = "Observed data")

# Calcular las densidades para cada una de las 500 columnas
densities <- lapply(1:500, function(i) {
  dens <- density(unlist(data_selected[, i]))
  data.frame(x = dens$x, y = dens$y, group = "Simulated sample")
})

# Unir todas las densidades en un solo data frame
densities_df <- do.call(rbind, densities)

# Calcular la densidad para la columna 'media'
dens_media <- density(data_selected$media)
df_media <- data.frame(x = dens_media$x, y = dens_media$y, group = "Mean simulated samples")

# Unir todos los data frames
all_densities <- bind_rows(df_huelva, densities_df, df_media)

# Crear el gráfico con ggplot2
p <- ggplot() +
  geom_line(data = subset(all_densities, group == "Simulated sample"), aes(x = x, y = y, group = group, color = group), size = 1, alpha = 0.5) +
  geom_line(data = subset(all_densities, group == "Observed data"), aes(x = x, y = y, group = group, color = group), size = 2) +
  geom_line(data = subset(all_densities, group == "Mean simulated samples"), aes(x = x, y = y, group = group, color = group), size = 2) +
  labs(title = "Density Plots of Pb Concentrations and Media",
       x = "Density",
       y = "Frequency") +
  scale_color_manual(values = c("Observed data" = "blue", 
                                "Mean simulated samples" = "red", 
                                "Simulated sample" = "grey"),
                     labels = c("Observed data", "Mean simulated samples", "Simulated sample")) +
  theme_minimal()

print(p)