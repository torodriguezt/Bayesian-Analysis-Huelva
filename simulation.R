simulaciones_cu_1 <- simulaciones_cu %>%
  filter(clasificacion == "Industria")

simulaciones_hg_1 <- simulaciones_hg %>%
  filter(clasificacion == "Industria")

simulaciones_pb_1 <- simulaciones_pb %>%
  filter(clasificacion == "Industria")

simulaciones_pb_1$Media <- rowMeans(simulaciones_pb_1[1:5000])
simulaciones_cu_1$Media <- rowMeans(simulaciones_cu_1[1:5000])
simulaciones_hg_1$Media <- rowMeans(simulaciones_hg_1[1:5000])


cor(simulaciones_pb_1$Media, simulaciones_cu_1$Media)
cor(simulaciones_cu_1$Media, simulaciones_hg_1$Media)
cor(simulaciones_hg_1$Media, simulaciones_pb_1$Media)


simulaciones_cu_2 <- simulaciones_cu %>%
  filter(clasificacion == "Urbano")

simulaciones_hg_2 <- simulaciones_hg %>%
  filter(clasificacion == "Urbano")

simulaciones_pb_2 <- simulaciones_pb %>%
  filter(clasificacion == "Urbano")

simulaciones_pb_2$Media <- rowMeans(simulaciones_pb_2[1:5000])
simulaciones_cu_2$Media <- rowMeans(simulaciones_cu_2[1:5000])
simulaciones_hg_2$Media <- rowMeans(simulaciones_hg_2[1:5000])


cor(simulaciones_pb_2$Media, simulaciones_cu_2$Media)
cor(simulaciones_cu_2$Media, simulaciones_hg_2$Media)
cor(simulaciones_hg_2$Media, simulaciones_pb_2$Media)

simulaciones_cu_3 <- simulaciones_cu %>%
  filter(clasificacion == "Agric")

simulaciones_hg_3 <- simulaciones_hg %>%
  filter(clasificacion == "Agric")

simulaciones_pb_3 <- simulaciones_pb %>%
  filter(clasificacion == "Agric")

simulaciones_pb_3$Media <- rowMeans(simulaciones_pb_3[1:5000])
simulaciones_cu_3$Media <- rowMeans(simulaciones_cu_3[1:5000])
simulaciones_hg_3$Media <- rowMeans(simulaciones_hg_3[1:5000])


cor(simulaciones_pb_3$Media, simulaciones_cu_3$Media)
cor(simulaciones_cu_3$Media, simulaciones_hg_3$Media)
cor(simulaciones_hg_3$Media, simulaciones_pb_3$Media)

simulaciones_cu_4 <- simulaciones_cu %>%
  filter(clasificacion == "Marsh")

simulaciones_hg_4 <- simulaciones_hg %>%
  filter(clasificacion == "Marsh")

simulaciones_pb_4 <- simulaciones_pb %>%
  filter(clasificacion == "Marsh")

simulaciones_pb_4$Media <- rowMeans(simulaciones_pb_4[1:5000])
simulaciones_cu_4$Media <- rowMeans(simulaciones_cu_4[1:5000])
simulaciones_hg_4$Media <- rowMeans(simulaciones_hg_4[1:5000])


cor(simulaciones_pb_4$Media, simulaciones_cu_4$Media)
cor(simulaciones_cu_4$Media, simulaciones_hg_4$Media)
cor(simulaciones_hg_4$Media, simulaciones_pb_4$Media)

simulaciones_cu_5 <- simulaciones_cu %>%
  filter(clasificacion == "Phospho")

simulaciones_hg_5 <- simulaciones_hg %>%
  filter(clasificacion == "Phospho")

simulaciones_pb_5 <- simulaciones_pb %>%
  filter(clasificacion == "Phospho")

simulaciones_pb_5$Media <- rowMeans(simulaciones_pb_5[1:5000])
simulaciones_cu_5$Media <- rowMeans(simulaciones_cu_5[1:5000])
simulaciones_hg_5$Media <- rowMeans(simulaciones_hg_5[1:5000])


cor(simulaciones_pb_5$Media, simulaciones_cu_5$Media)
cor(simulaciones_cu_5$Media, simulaciones_hg_5$Media)
cor(simulaciones_hg_5$Media, simulaciones_pb_5$Media)

simulaciones_cu_6 <- simulaciones_cu %>%
  filter(clasificacion == "Bare")

simulaciones_hg_6 <- simulaciones_hg %>%
  filter(clasificacion == "Bare")

simulaciones_pb_6 <- simulaciones_pb %>%
  filter(clasificacion == "Bare")

simulaciones_pb_6$Media <- rowMeans(simulaciones_pb_6[1:5000])
simulaciones_cu_6$Media <- rowMeans(simulaciones_cu_6[1:5000])
simulaciones_hg_6$Media <- rowMeans(simulaciones_hg_6[1:5000])


cor(simulaciones_pb_6$Media, simulaciones_cu_6$Media)
cor(simulaciones_cu_6$Media, simulaciones_hg_6$Media)
cor(simulaciones_hg_6$Media, simulaciones_pb_6$Media)



simulaciones_pb_1$Media <- rowMeans(simulaciones_pb_1[1:5000])
simulaciones_cu_1$Media <- rowMeans(simulaciones_cu_1[1:5000])

cor(simulaciones_pb_1$Media, simulaciones_cu_1$Media)

correlaciones_pb_hg <- vector()
spearman_pb_hg <- vector()

for (i in 1:22426){
  pb <- (unlist(simulaciones_pb[i,1:5000]))
  hg <- (unlist(simulaciones_hg[i,1:5000]))
  
  correlaciones_pb_hg[i] <- cor(pb, hg)
  spearman_pb_hg[i] <- cor(pb, hg, method = "spearman")
}

correlaciones_pb_hg <- as.data.frame(correlaciones_pb_hg)

correlaciones_pb_hg$X <- grid_huelva$X
correlaciones_pb_hg$Y <- grid_huelva$Y


correlaciones_hg_cu <- vector()
spearman_hg_cu <- vector()

for (i in 1:22426){
  hg <- (unlist(simulaciones_hg[i,1:5000]))
  cu <- (unlist(simulaciones_cu[i,1:5000]))
  
  correlaciones_hg_cu[i] <- cor(hg, cu)
  spearman_hg_cu[i] <- cor(hg, cu, method = "spearman")
}

correlaciones_hg_cu <- as.data.frame(correlaciones_hg_cu)

correlaciones_hg_cu$X <- grid_huelva$X
correlaciones_hg_cu$Y <- grid_huelva$Y


correlaciones_pb_cu <- vector()
spearman_pb_cu <- vector()

for (i in 1:22426){
  pb <- (unlist(simulaciones_pb[i,1:5000]))
  cu <- (unlist(simulaciones_cu[i,1:5000]))
  
  correlaciones_pb_cu[i] <- cor(pb, cu)
  spearman_pb_cu[i] <- cor(pb, cu, method = "spearman")
}

correlaciones_pb_cu <- as.data.frame(correlaciones_pb_cu)

correlaciones_pb_cu$X <- grid_huelva$X
correlaciones_pb_cu$Y <- grid_huelva$Y



number_values_pb_hg <- setNames(
  c(0.610, 0.851, 0.540, 0.03, 0.740, 0.844),
  c('Agric', 'Marsh', 'Urbano', 'Phospho', 'Industria', "Bare")
)

correlaciones_pb_hg$Clasificacion <- simulaciones_pb$clasificacion
correlaciones_pb_hg$Corr_Zona <- number_values_pb_hg[correlaciones_pb_hg$Clasificacion]

number_values_hg_cu <- setNames(
  c(0.504, 0.646, 0.523, 0.03, 0.735, 0.676),
  c('Agric', 'Marsh', 'Urbano', 'Phospho', 'Industria', "Bare")
)

correlaciones_hg_cu$Clasificacion <- simulaciones_pb$clasificacion
correlaciones_hg_cu$Corr_Zona <- number_values_hg_cu[correlaciones_hg_cu$Clasificacion]

number_values_pb_cu <- setNames(
  c(0.863, 0.761, 0.938, 0.709, 0.922, 0.825),
  c('Agric', 'Marsh', 'Urbano', 'Phospho', 'Industria', "Bare")
)

correlaciones_pb_cu$Clasificacion <- simulaciones_pb$clasificacion
correlaciones_pb_cu$Corr_Zona <- number_values_pb_cu[correlaciones_pb_cu$Clasificacion]


p2 <- ggplot(correlaciones_pb_hg, aes(x = X, y = Y, fill = (Corr_Zona))) +
  geom_tile() +
  scale_fill_viridis(name = "Hg", na.value = 'transparent', option = "D") +
  labs(title = 'Correlacion',
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

datos <- correlaciones_pb_hg


library(leaflet)
library(viridis)  # Para usar la paleta de colores viridis

# Supongamos que tienes algunos rangos específicos que quieres en la leyenda
breaks <- c(min(datos$Corr_Zona), 0.54, 0.61, 0.74, 0.844, 0.851, 1)
labels <- c("0.03 Phospho", "0.54 Urban", "0.61 Agric", "0.74 Industry", "0.84 Bare", "0.85 Marsh")

# Colores para cada uno de los rangos
colors <- viridis_pal(option = "D", direction = 1)(length(breaks) - 1)

mapa <- leaflet(data = datos) %>%
  addTiles() %>%  # Añadir las baldosas básicas del mapa
  addCircleMarkers(
    lng = ~X, lat = ~Y,  # Longitud y latitud
    color = ~colors,  # Color según Corr_Zona
    radius = 1,  # Tamaño del radio de los círculos
    fillOpacity = 1,
    stroke = FALSE
  ) %>%
  addLegend(
    position = "bottomright",  # Posición de la leyenda
    title = "Correlation Hg(ppb) and Pb(ppm)",  # Título de la leyenda
    colors = colors,  # Colores de la leyenda
    labels = labels,  # Textos de la leyenda
    opacity = 1
  )

# Mostrar el mapa
mapa