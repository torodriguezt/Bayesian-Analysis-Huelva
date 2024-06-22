library(leaflet)
library(dplyr)
library(RColorBrewer)

# Asumiendo que correlaciones_pb_hg es tu dataframe y tiene las columnas X, Y, Corr_Zona
correlaciones_pb_hg <- correlaciones_pb_hg %>%
  mutate(
    Label = cut(Corr_Zona,
                breaks = c(0.501, 0.607, 0.711, 0.767, 0.787, 0.809, 0.947),
                labels = c("0.607 Marsh", "0.711 Urban", "0.767 Phospho", "0.787 Agric", "0.809 Industry","0.942 Bare"),
                include.lowest = TRUE)
  )

colors <- colorFactor(palette = "viridis", domain = correlaciones_pb_hg$Label)

# Crear el mapa en Leaflet
leaflet(data = correlaciones_pb_hg) %>%
  addTiles() %>%
  addCircles(
    ~X, ~Y,
    opacity = 0.35,
    fillColor = ~colors(Label),
    color = "#FFFFFF", # Color de borde de cada punto
    fillOpacity = 0.8, weight = 1,
    popup = ~paste("Correlación: ", Corr_Zona)  # Popup para mostrar la correlación al hacer clic en un punto
  ) %>%
  addLegend(pos = "bottomright", pal = colors, values = ~Label,
            title = '<span style="color: black; font-weight: bold;">Correlation Hg(ppb) and Pb(ppm)</span>',
            opacity = 0.7,
            labels = unique(correlaciones_pb_hg$Label)) %>%
  addScaleBar(position = "bottomleft")


library(leaflet)
library(dplyr)
library(RColorBrewer)

# Asumiendo que correlaciones_pb_hg es tu dataframe y tiene las columnas X, Y, Corr_Zona
correlaciones_hg_cu <- correlaciones_hg_cu %>%
  mutate(
    Label = cut(Corr_Zona,
                breaks = c(0.4, 0.45, 0.708, 0.786, 0.794 ,0.823, 0.947),
                labels = c("0.45 Marsh", "0.708 Urban", "0.786 Phospho","0.794 Bare" ,"0.823 Agric", "0.894 Industry"),
                include.lowest = TRUE)
  )

colors <- colorFactor(palette = "viridis", domain = correlaciones_hg_cu$Label)

# Crear el mapa en Leaflet
leaflet(data = correlaciones_hg_cu) %>%
  addTiles() %>%
  addCircles(
    ~X, ~Y,
    opacity = 0.35,
    fillColor = ~colors(Label),
    color = "#FFFFFF", # Color de borde de cada punto
    fillOpacity = 0.8, weight = 1,
    popup = ~paste("Correlación: ", Corr_Zona)  # Popup para mostrar la correlación al hacer clic en un punto
  ) %>%
  addLegend(pos = "bottomright", pal = colors, values = ~Label,
            title = '<span style="color: black; font-weight: bold;">Correlation Hg(ppb) and Cu(ppm)</span>',
            opacity = 0.7,
            labels = unique(correlaciones_hg_cu$Label)) %>%
  addScaleBar(position = "bottomleft")



correlaciones_pb_cu <- correlaciones_pb_cu %>%
  mutate(
    Label = cut(Corr_Zona,
                breaks = c(0.7, 0.705, 0.829, 0.860, 0.891, 0.9355 ,0.9510),
                labels = c("0.705 Marsh", "0.829 Urban", "0.860 Industry","0.891 Agric" ,"0.9355 Bare", "0.9510 Phospho"),
                include.lowest = TRUE)
  )

colors <- colorFactor(palette = "viridis", domain = correlaciones_pb_cu$Label)

# Crear el mapa en Leaflet
leaflet(data = correlaciones_pb_cu) %>%
  addTiles() %>%
  addCircles(
    ~X, ~Y,
    opacity = 0.35,
    fillColor = ~colors(Label),
    color = "#FFFFFF", # Color de borde de cada punto
    fillOpacity = 0.8, weight = 1,
    popup = ~paste("Correlación: ", Corr_Zona)  # Popup para mostrar la correlación al hacer clic en un punto
  ) %>%
  addLegend(pos = "bottomright", pal = colors, values = ~Label,
            title = '<span style="color: black; font-weight: bold;">Correlation Pb(ppm) and Cu(ppm)</span>',
            opacity = 0.7,
            labels = unique(correlaciones_pb_cu$Label)) %>%
  addScaleBar(position = "bottomleft")
