grid_huelva <- st_read("SinAgua/capa_sin_agua.shp")

grid_puntos <- st_coordinates(grid_huelva)

plot(grid_puntos, type = "p", cex = 0.1)

zona <- raster::raster("C:/Users/Tomas/Documents/INLA_Part2/Raster/raster_huelva.tif")

