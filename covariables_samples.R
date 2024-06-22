agric <- st_read("Muestra_part2/Agric_2/Agric.shp")
urban <- st_read("Muestra_part2/Urban/Urban.shp")
industria <- st_read("Muestra_part2/Industry/Industry.shp")
refineria <- st_read("Muestra_part2/Refinery/Refinery.shp")
phospho <- st_read("Muestra_part2/Phospho/Phospho.shp")
marsh <- st_read("Muestra_part2/Marsh/Marsh.shp")
bare <- st_read("Muestra_part2/Bare/Bare.shp")
park <- st_read("Muestra_part2/Park/Park.shp")

datos_huelva$Hg__ppb_ <- agric$Hg__ppb_
datos_huelva$Pb__ppm_ <- agric$Pb__ppm_
datos_huelva$Distancia_Industria <- industria$distance
datos_huelva$Distancia_Marsh <- marsh$distance
datos_huelva$Distancia_Phospho <- phospho$distance
datos_huelva$Distancia_Refineria <- refineria$distance
datos_huelva$Distancia_Urbano <- urban$distance
datos_huelva$Distancia_Agric <- agric$distance
datos_huelva$Distancia_Bare <- bare$distance
datos_huelva$Distancia_Park <- park$distance

coordenadas_muestras <- as.data.frame(st_coordinates(datos_huelva))
datos_huelva$X <- coordenadas_muestras$X
datos_huelva$Y <- coordenadas_muestras$Y
