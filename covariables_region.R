Agric <- st_read("Covariables_muestra_part2/Agric/Agric.shp")
Industria <- st_read("Covariables_muestra_part2/Industria/Industry.shp")
Marsh <- st_read("Covariables_muestra_part2/Marsh/Marsh.shp")
Bare <- st_read("Covariables_muestra_part2/Bare/Bare.shp")
Refineria <- st_read("Covariables_muestra_part2/Refinery/Refineria.shp")
Park <- st_read("Covariables_muestra_part2/Natural Park/Park.shp")
Urban <- st_read("Covariables_muestra_part2/Urban/Urban.shp")
Phospho <- st_read("Covariables_muestra_part2/Phospho/Phospho.shp")



grid_huelva$Industria <- Industria$distance
grid_huelva$Marsh <- Marsh$distance
grid_huelva$Phospho <- Phospho$distance
grid_huelva$Urbano <- Urban$distance
grid_huelva$Agric <- Agric$distance
grid_huelva$Refineria <- Refineria$distance
grid_huelva$Bare <- Bare$distance
grid_huelva$Park <- Park$distance


coordenadas_region <- as.data.frame(st_coordinates(grid_huelva))
grid_huelva$X <- coordenadas_region$X
grid_huelva$Y <- coordenadas_region$Y

grid_huelva_clasificacion <- select(grid_huelva, c("Industria", "Marsh", "Agric", "Urbano", "Phospho", "Bare"))
grid_huelva_clasificacion <- st_drop_geometry(grid_huelva_clasificacion)

check_zeros <- function(row) {
  zero_cols <- which(row == 0)
  if (length(zero_cols) > 0) {
    return(names(row)[zero_cols[1]])
  } else {
    return("Phospho")
  }
}

grid_huelva_clasificacion$clasificacion <- apply(grid_huelva_clasificacion, 1, check_zeros)
