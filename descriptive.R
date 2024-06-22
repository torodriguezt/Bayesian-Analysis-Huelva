create_plot <- function(data, x, y, log_transform = FALSE, title) {
  y_label <- if(log_transform) {
    "Log Hg Concentration (ppb)"
  } else {
    "Hg Concentration (ppb)"
  }
  
  p1 <- ggplot(data, aes_string(x = x, y = y)) +
    geom_point() +
    geom_smooth(method = "lm", col = "blue") +
    labs(title = paste("Scatter -", title), x = "Distance (m)", y = y_label) +
    theme(plot.title = element_text(hjust = 0.5))
  
  if (log_transform) {
    p1 <- p1 + scale_y_log10()
  } else {
    p1 <- p1 + ylim(0, 1000)
  }
  
  p2 <- ggplot(data, aes_string(x = y)) +
    geom_density(fill = "#1E90FF", alpha = 0.6) +
    labs(title = paste("Density Hg"), x = y_label) +
    theme(plot.title = element_text(hjust = 0.6, size = 30))
  
  combined <- grid.arrange(p1, p2, ncol = 2)
  
  return(list(p1, p2))
}

g1 <- create_plot(datos_huelva, "Distancia_Refineria", "log(Hg__ppb_)", TRUE, "Distance to Refinery")
g2 <- create_plot(datos_huelva, "Distancia_Refineria", "Hg__ppb_", FALSE, "Distance to Refinery")

grid.arrange(grobs = c(g2, g1), ncol = 2, nrow = 2)

grid.arrange(grobs = c(g2[2], g1[2]), ncol = 2, nrow = 1)

g3 <- create_plot(datos_huelva, "Distancia_Industria", "log(Hg__ppb_)", TRUE, "Distance to Industry")
g4 <- create_plot(datos_huelva, "Distancia_Industria", "Hg__ppb_", FALSE, "Distance to Industry")

grid.arrange(grobs = c(g4, g3), ncol = 2, nrow = 2)


g5 <- create_plot(datos_huelva, "Distancia_Agric", "log(Hg__ppb_)", TRUE, "Distance to Agricultural Areas")
g6 <- create_plot(datos_huelva, "Distancia_Agric", "Hg__ppb_", FALSE, "Distance to Agricultural Areas")

grid.arrange(grobs = c(g6, g5), ncol = 2, nrow = 2)


g7 <- create_plot(datos_huelva, "Distancia_Marsh", "log(Hg__ppb_)", TRUE, "Distance to Salt Marsh")
g8 <- create_plot(datos_huelva, "Distancia_Marsh", "Hg__ppb_", FALSE, "Distance to Salt Marsh")

grid.arrange(grobs = c(g8, g7), ncol = 2, nrow = 2)

g9 <- create_plot(datos_huelva, "Distancia_Phospho", "log(Hg__ppb_)", TRUE, "Distance to Phosphogypsum Piles")
g10 <- create_plot(datos_huelva, "Distancia_Phospho", "Hg__ppb_", FALSE, "Distance to Phosphogypsum Piles")

grid.arrange(grobs = c(g10, g9), ncol = 2, nrow = 2)

g11 <- create_plot(datos_huelva, "Distancia_Urbano", "log(Hg__ppb_)", TRUE, "Distance to Urban Areas")
g12 <- create_plot(datos_huelva, "Distancia_Urbano", "Hg__ppb_", FALSE, "Distance to Urban Areas")

grid.arrange(grobs = c(g12, g11), ncol = 2, nrow = 2)
