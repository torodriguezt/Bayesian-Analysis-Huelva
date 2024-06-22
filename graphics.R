addLegend_decreasing <- function (map, position = c("topright", "bottomright", "bottomleft","topleft"),
                                  pal, values, na.label = "NA", bins = 7, colors, 
                                  opacity = 0.5, labels = NULL, labFormat = labelFormat(), 
                                  title = NULL, className = "info legend", layerId = NULL, 
                                  group = NULL, data = getMapData(map), decreasing = FALSE) {
  
  position <- match.arg(position)
  type <- "unknown"
  na.color <- NULL
  extra <- NULL
  if (!missing(pal)) {
    if (!missing(colors)) 
      stop("You must provide either 'pal' or 'colors' (not both)")
    if (missing(title) && inherits(values, "formula")) 
      title <- deparse(values[[2]])
    values <- evalFormula(values, data)
    type <- attr(pal, "colorType", exact = TRUE)
    args <- attr(pal, "colorArgs", exact = TRUE)
    na.color <- args$na.color
    if (!is.null(na.color) && col2rgb(na.color, alpha = TRUE)[[4]] == 
        0) {
      na.color <- NULL
    }
    if (type != "numeric" && !missing(bins)) 
      warning("'bins' is ignored because the palette type is not numeric")
    if (type == "numeric") {
      cuts <- if (length(bins) == 1) 
        pretty(values, bins)
      else bins   
      if (length(bins) > 2) 
        if (!all(abs(diff(bins, differences = 2)) <= 
                 sqrt(.Machine$double.eps))) 
          stop("The vector of breaks 'bins' must be equally spaced")
      n <- length(cuts)
      r <- range(values, na.rm = TRUE)
      cuts <- cuts[cuts >= r[1] & cuts <= r[2]]
      n <- length(cuts)
      p <- (cuts - r[1])/(r[2] - r[1])
      extra <- list(p_1 = p[1], p_n = p[n])
      p <- c("", paste0(100 * p, "%"), "")
      if (decreasing == TRUE){
        colors <- pal(rev(c(r[1], cuts, r[2])))
        labels <- rev(labFormat(type = "numeric", cuts))
      }else{
        colors <- pal(c(r[1], cuts, r[2]))
        labels <- rev(labFormat(type = "numeric", cuts))
      }
      colors <- paste(colors, p, sep = " ", collapse = ", ")
    }
    else if (type == "bin") {
      cuts <- args$bins
      n <- length(cuts)
      mids <- (cuts[-1] + cuts[-n])/2
      if (decreasing == TRUE){
        colors <- pal(rev(mids))
        labels <- rev(labFormat(type = "bin", cuts))
      }else{
        colors <- pal(mids)
        labels <- labFormat(type = "bin", cuts)
      }
    }
    else if (type == "quantile") {
      p <- args$probs
      n <- length(p)
      cuts <- quantile(values, probs = p, na.rm = TRUE)
      mids <- quantile(values, probs = (p[-1] + p[-n])/2, na.rm = TRUE)
      if (decreasing == TRUE){
        colors <- pal(rev(mids))
        labels <- rev(labFormat(type = "quantile", cuts, p))
      }else{
        colors <- pal(mids)
        labels <- labFormat(type = "quantile", cuts, p)
      }
    }
    else if (type == "factor") {
      v <- sort(unique(na.omit(values)))
      colors <- pal(v)
      labels <- labFormat(type = "factor", v)
      if (decreasing == TRUE){
        colors <- pal(rev(v))
        labels <- rev(labFormat(type = "factor", v))
      }else{
        colors <- pal(v)
        labels <- labFormat(type = "factor", v)
      }
    }
    else stop("Palette function not supported")
    if (!any(is.na(values))) 
      na.color <- NULL
  }
  else {
    if (length(colors) != length(labels)) 
      stop("'colors' and 'labels' must be of the same length")
  }
  legend <- list(colors = I(unname(colors)), labels = I(unname(labels)), 
                 na_color = na.color, na_label = na.label, opacity = opacity, 
                 position = position, type = type, title = title, extra = extra, 
                 layerId = layerId, className = className, group = group)
  invokeMethod(map, data, "addLegend", legend)
}

#####Mercurio---------------------------------------------------------------------------------
library(leaflet)
library(scales)  # para la función alpha

# Crear una paleta de colores con transparencia ajustada
pal <- colorNumeric("viridis", range(grid_huelva$hg.spde), na.color = "transparent", alpha = 0.35)  # nivel muy bajo de alpha para la paleta

leaflet() %>% addTiles() %>%
  addCircles(
    lng = grid_huelva$X, 
    lat = grid_huelva$Y,
    color = pal(grid_huelva$hg.spde),
    opacity = 0.35,  
    fillOpacity = 0.35  
  ) %>%
  addLegend_decreasing("bottomright", pal = pal, values = grid_huelva$hg.spde, 
                       title = "log Hg (ppb)", decreasing = TRUE) %>%
  addScaleBar(position = "bottomleft")

pal <- colorNumeric("viridis", range(grid_huelva$hg.spde.sd), na.color = "transparent", alpha = 0.35)

leaflet() %>% addTiles() %>%
  addCircles(lng = grid_huelva$X, lat = grid_huelva$Y, 
             color = pal(grid_huelva$hg.spde.sd),
             opacity = 0.35,
             fillOpacity = 0.35) %>%
  addLegend_decreasing("bottomright", pal = pal, values = grid_huelva$hg.spde.sd, 
                       title = "Mercury SD", decreasing = TRUE) %>%
  addScaleBar(position = c("bottomleft"))


#####Plomo------------------------------------------------------------------------------------

pal <- colorNumeric("viridis", range(grid_huelva$pb.spde), na.color = "transparent", alpha = 0.35)

leaflet() %>% addTiles() %>%
  addCircles(lng = grid_huelva$X, lat = grid_huelva$Y, 
             color = pal(grid_huelva$pb.spde),
             opacity = 0.35,
             fillOpacity = 0.35) %>%
  addLegend_decreasing("bottomright", pal = pal, values = grid_huelva$pb.spde, 
                       title = "log Pb (ppm)", decreasing = TRUE) %>%
  addScaleBar(position = c("bottomleft"))


pal <- colorNumeric("viridis", range(grid_huelva$pb.spde.sd), na.color = "transparent", alpha = 0.35)

leaflet() %>% addTiles() %>%
  addCircles(lng = grid_huelva$X, lat = grid_huelva$Y, 
             color = pal(grid_huelva$pb.spde.sd),
             opacity = 0.35,
             fillOpacity = 0.35) %>%
  addLegend_decreasing("bottomright", pal = pal, values = grid_huelva$pb.spde.sd, 
                       title = "Pb SD", decreasing = TRUE) %>%
  addScaleBar(position = c("bottomleft"))

###Cobre--------------------------------------------------------------------------------------
pal <- colorNumeric("viridis", range(grid_huelva$cu.spde), na.color = "transparent", alpha = 0.35)

leaflet() %>% addTiles() %>%
  addCircles(lng = grid_huelva$X, lat = grid_huelva$Y, 
             color = pal(grid_huelva$cu.spde),
             opacity = 0.35,
             fillOpacity = 0.35) %>%
  addLegend_decreasing("bottomright", pal = pal, values = grid_huelva$cu.spde, 
                       title = "log Cu (ppm)", decreasing = TRUE) %>%
  addScaleBar(position = c("bottomleft"))


pal <- colorNumeric("viridis", range(grid_huelva$cu.spde.sd), na.color = "transparent", alpha = 0.35)

leaflet() %>% addTiles() %>%
  addCircles(lng = grid_huelva$X, lat = grid_huelva$Y, 
             color = pal(grid_huelva$cu.spde.sd),
             opacity = 0.35,
             fillOpacity = 0.35) %>%
  addLegend_decreasing("bottomright", pal = pal, values = grid_huelva$cu.spde.sd, 
                       title = "Cu SD", decreasing = TRUE) %>%
  addScaleBar(position = c("bottomleft"))
###Probabilidad de excedencia----------------------------------------------------------------

##Mercurio

pal <- colorNumeric("viridis", c(0, 1), na.color = "transparent", alpha = 0.35)

leaflet() %>% addTiles() %>%
  addCircles(lng = grid_huelva$X, lat = grid_huelva$Y, 
             color = pal(excprob_hg_umbral),
             opacity = 0.35,
             fillOpacity = 0.35) %>%
  addLegend("bottomright", pal = pal, values = excprob_hg_umbral, title = "P(concen>log(3000))") %>%
  addScaleBar(position = c("bottomleft"))%>%
  addControl('<div style="background-color: rgba(255, 255, 255, 0.7); border-radius: 5px; padding: 8px; font-family: Arial, sans-serif; text-align: center;">
                <h2 style="color: #333; margin: 0; font-size: 16px;">log Hg(ppb) Threshold</h2>
              </div>',
             position = "topright")

pal <- colorNumeric("viridis", c(0, 1), na.color = "transparent", alpha = 0.35)

leaflet() %>% addTiles() %>%
  addCircles(lng = grid_huelva$X, lat = grid_huelva$Y, 
             color = pal(excprob_hg_agric),
             opacity = 0.35,
             fillOpacity = 0.35) %>%
  addLegend("bottomright", pal = pal, values = excprob_hg_agric, title = "P(concen>log(1000))") %>%
  addScaleBar(position = c("bottomleft"))%>%
  addControl('<div style="background-color: rgba(255, 255, 255, 0.7); border-radius: 5px; padding: 8px; font-family: Arial, sans-serif; text-align: center;">
                <h2 style="color: #333; margin: 0; font-size: 16px;">log Hg(ppb) Agriculture</h2>
              </div>',
             position = "topright")

#Plomo
pal <- colorNumeric("viridis", c(0,1), na.color = "transparent", alpha = 0.35)

leaflet() %>% addTiles() %>%
  addCircles(lng = grid_huelva$X, lat = grid_huelva$Y, 
             color = pal(excprob_pb_urb),
             opacity = 0.35,
             fillOpacity = 0.35) %>%
  addLegend("bottomright", pal = pal, values = excprob_pb_urb, title = "P(concen>log(275))") %>%
  addScaleBar(position = c("bottomleft"))%>%
  addControl('<div style="background-color: rgba(255, 255, 255, 0.7); border-radius: 5px; padding: 8px; font-family: Arial, sans-serif; text-align: center;">
                <h2 style="color: #333; margin: 0; font-size: 16px;">log Pb(ppm) Urban</h2>
              </div>',
             position = "topright")

pal <- colorNumeric("viridis", c(0, 1), na.color = "transparent", alpha = 0.35)

leaflet() %>% addTiles() %>%
  addCircles(lng = grid_huelva$X, lat = grid_huelva$Y, 
             color = pal(excprob_pb_agric),
             opacity = 0.35,
             fillOpacity = 0.35) %>%
  addLegend("bottomright", pal = pal, values = excprob_pb_agric, title = "P(concen>log(100))") %>%
  addScaleBar(position = c("bottomleft"))%>%
  addControl('<div style="background-color: rgba(255, 255, 255, 0.7); border-radius: 5px; padding: 8px; font-family: Arial, sans-serif; text-align: center;">
                <h2 style="color: #333; margin: 0; font-size: 16px;">log Pb(ppm) Agric</h2>
              </div>',
             position = "topright")



#Cobre

pal <- colorNumeric("viridis", c(0,1), na.color = "transparent", alpha = 0.35)

leaflet() %>% addTiles() %>%
  addCircles(lng = grid_huelva$X, lat = grid_huelva$Y, 
             color = pal(excprob_cu_urb),
             opacity = 0.35,
             fillOpacity = 0.35) %>%
  addLegend("bottomright", pal = pal, values = excprob_cu_urb, title = "P(concen>3130)") %>%
  addScaleBar(position = c("bottomleft"))%>%
  addControl('<div style="background-color: rgba(255, 255, 255, 0.7); border-radius: 5px; padding: 8px; font-family: Arial, sans-serif; text-align: center;">
                <h2 style="color: #333; margin: 0; font-size: 16px;">Cu(ppm) Urban</h2>
              </div>',
             position = "topright")


pal <- colorNumeric("viridis", c(0,1), na.color = "transparent", alpha = 0.35)

leaflet() %>% addTiles() %>%
  addCircles(lng = grid_huelva$X, lat = grid_huelva$Y, 
             color = pal(excprob_cu_otro),
             opacity = 0.35,
             fillOpacity = 0.35) %>%
  addLegend("bottomright", pal = pal, values = excprob_cu_otro, title = "P(log(concen)>595)") %>%
  addScaleBar(position = c("bottomleft"))%>%
  addControl('<div style="background-color: rgba(255, 255, 255, 0.7); border-radius: 5px; padding: 8px; font-family: Arial, sans-serif; text-align: center;">
                <h2 style="color: #333; margin: 0; font-size: 16px;">log Cu(ppm) Other</h2>
              </div>',
             position = "topright")

pal <- colorNumeric("viridis", c(0,1), na.color = "transparent", alpha = 0.35)

leaflet() %>% addTiles() %>%
  addCircles(lng = grid_huelva$X, lat = grid_huelva$Y, 
             color = pal(excprob_cu_agric),
             opacity = 0.35,
             fillOpacity = 0.35) %>%
  addLegend("bottomright", pal = pal, values = excprob_cu_agric, title = "P(log(concen)>50)") %>%
  addScaleBar(position = c("bottomleft"))%>%
  addControl('<div style="background-color: rgba(255, 255, 255, 0.7); border-radius: 5px; padding: 8px; font-family: Arial, sans-serif; text-align: center;">
                <h2 style="color: #333; margin: 0; font-size: 16px;">log Cu(ppm) Agric</h2>
              </div>',
             position = "topright")



###Simulacion Posteriori Predictiva-------------------------------------------------------------------

pal <- colorNumeric("viridis", range(correlaciones_pb_hg$correlaciones_pb_hg), na.color = "transparent", alpha = 0.35)

leaflet() %>% 
  addTiles() %>%
  addCircles(lng = grid_huelva$X, lat = grid_huelva$Y, 
             color = pal(correlaciones_pb_hg$correlaciones_pb_hg),
             opacity = 0.35,
             fillOpacity = 0.35) %>%
  addLegend_decreasing("bottomright", pal = pal, values = correlaciones_pb_hg$correlaciones_pb_hg, 
                       title = "Correlation", decreasing = TRUE) %>%
  addScaleBar(position = "bottomleft") %>%
  addControl('<div style="background-color: rgba(255, 255, 255, 0.7); border-radius: 5px; padding: 8px; font-family: Arial, sans-serif; text-align: center;">
                <h2 style="color: #333; margin: 0; font-size: 16px;">Pb and Hg</h2>
              </div>',
             position = "topright")


pal <- colorNumeric("viridis", range(correlaciones_hg_cu$correlaciones_hg_cu), na.color = "transparent", alpha = 0.35)

leaflet() %>% 
  addTiles() %>%
  addCircles(lng = grid_huelva$X, lat = grid_huelva$Y, 
             color = pal(correlaciones_hg_cu$correlaciones_hg_cu),
             opacity = 0.35,
             fillOpacity = 0.35) %>%
  addLegend_decreasing("bottomright", pal = pal, values = correlaciones_hg_cu$correlaciones_hg_cu, 
                       title = "Correlation", decreasing = TRUE) %>%
  addScaleBar(position = "bottomleft") %>%
  addControl('<div style="background-color: rgba(255, 255, 255, 0.7); border-radius: 5px; padding: 8px; font-family: Arial, sans-serif; text-align: center;">
                <h2 style="color: #333; margin: 0; font-size: 16px;">Hg and Cu</h2>
              </div>',
             position = "topright")



pal <- colorNumeric("viridis", range(correlaciones_pb_cu$correlaciones_pb_cu), na.color = "transparent", alpha = 0.35)

leaflet() %>% 
  addTiles() %>%
  addCircles(lng = grid_huelva$X, lat = grid_huelva$Y, 
             color = pal(correlaciones_pb_cu$correlaciones_pb_cu),
             opacity = 0.35,
             fillOpacity = 0.35) %>%
  addLegend_decreasing("bottomright", pal = pal, values = correlaciones_pb_cu$correlaciones_pb_cu, 
                       title = "Correlation", decreasing = TRUE) %>%
  addScaleBar(position = "bottomleft") %>%
  addControl('<div style="background-color: rgba(255, 255, 255, 0.7); border-radius: 5px; padding: 8px; font-family: Arial, sans-serif; text-align: center;">
                <h2 style="color: #333; margin: 0; font-size: 16px;">Pb and Cu</h2>
              </div>',
             position = "topright")






pal <- colorNumeric(
  palette = viridis(100, option = "D"),
  domain = correlaciones_pb_hg$Corr_Zona
)

leaflet(correlaciones_pb_hg) %>%
  addProviderTiles(providers$CartoDB.Positron) %>%
  addCircles(
    ~X, ~Y,
    color = ~pal(Corr_Zona),
    radius = 1,  # Ajusta el radio según sea necesario
    fillOpacity = 0.7,
    popup = ~paste("Correlación: ", Corr_Zona)
  ) %>%
  addLegend(
    "bottomright", 
    pal = pal, 
    values = ~Corr_Zona,
    title = "Hg",
    opacity = 1
  ) %>%
  addControl(
    tags$div(
      tags$h3("Correlación"),
      tags$h4("INLA")
    ), 
    position = "topleft"
  )










