
BASE_DIR <- here::here()
DATA_DIR <- file.path(BASE_DIR, "data", "raw")
RESULTS_DIR <- file.path(BASE_DIR, "results")
SRC_DIR <- file.path(BASE_DIR, "src")

PATHS <- list(
  boundaries = file.path(DATA_DIR, "boundaries", "Malla_no_water_2", "malla.shp"),
  
  samples = file.path(DATA_DIR, "samples", "Datos", "HuelvaResultados.shp"),
  
  prediction_grid = file.path(DATA_DIR, "prediction_grid", "SinAgua", "capa_sin_agua.shp"),
  
  land_use_samples = list(
    agric = file.path(DATA_DIR, "land_use", "samples", "Muestra_part2", "Agric_2", "Agric.shp"),
    urban = file.path(DATA_DIR, "land_use", "samples", "Muestra_part2", "Urban", "Urban.shp"),
    industria = file.path(DATA_DIR, "land_use", "samples", "Muestra_part2", "Industry", "Industry.shp"),
    refineria = file.path(DATA_DIR, "land_use", "samples", "Muestra_part2", "Refinery", "Refinery.shp"),
    phospho = file.path(DATA_DIR, "land_use", "samples", "Muestra_part2", "Phospho", "Phospho.shp"),
    marsh = file.path(DATA_DIR, "land_use", "samples", "Muestra_part2", "Marsh", "Marsh.shp"),
    bare = file.path(DATA_DIR, "land_use", "samples", "Muestra_part2", "Bare", "Bare.shp"),
    park = file.path(DATA_DIR, "land_use", "samples", "Muestra_part2", "Park", "Park.shp")
  ),
  
  land_use_covariates = list(
    agric = file.path(DATA_DIR, "land_use", "covariates", "Covariables_muestra_part2", "Agric", "Agric.shp"),
    urban = file.path(DATA_DIR, "land_use", "covariates", "Covariables_muestra_part2", "Urban", "Urban.shp"),
    industria = file.path(DATA_DIR, "land_use", "covariates", "Covariables_muestra_part2", "Industria", "Industry.shp"),
    refineria = file.path(DATA_DIR, "land_use", "covariates", "Covariables_muestra_part2", "Refinery", "Refineria.shp"),
    phospho = file.path(DATA_DIR, "land_use", "covariates", "Covariables_muestra_part2", "Phospho", "Phospho.shp"),
    marsh = file.path(DATA_DIR, "land_use", "covariates", "Covariables_muestra_part2", "Marsh", "Marsh.shp"),
    bare = file.path(DATA_DIR, "land_use", "covariates", "Covariables_muestra_part2", "Bare", "Bare.shp"),
    park = file.path(DATA_DIR, "land_use", "covariates", "Covariables_muestra_part2", "Natural Park", "Park.shp")
  ),
  
  raster = file.path(DATA_DIR, "raster", "Raster", "tercero.tif")
)

OUTPUT_PATHS <- list(
  figures = file.path(RESULTS_DIR, "figures"),
  reports = file.path(RESULTS_DIR, "reports"),
  models = file.path(RESULTS_DIR, "models"),
  predictions = file.path(RESULTS_DIR, "predictions")
)

create_output_dirs <- function() {
  lapply(OUTPUT_PATHS, function(path) {
    if (!dir.exists(path)) {
      dir.create(path, recursive = TRUE)
      cat("Creado directorio:", path, "\n")
    }
  })
  invisible(TRUE)
}

check_input_files <- function() {
  missing_files <- c()
  
  individual_files <- list(
    PATHS$boundaries,
    PATHS$samples,
    PATHS$prediction_grid,
    PATHS$raster
  )
  
  for (file in individual_files) {
    if (!file.exists(file)) {
      missing_files <- c(missing_files, file)
    }
  }
  
  for (file in PATHS$land_use_samples) {
    if (!file.exists(file)) {
      missing_files <- c(missing_files, file)
    }
  }
  
  for (file in PATHS$land_use_covariates) {
    if (!file.exists(file)) {
      missing_files <- c(missing_files, file)
    }
  }
  
  if (length(missing_files) > 0) {
    cat("ADVERTENCIA: Los siguientes archivos no fueron encontrados:\n")
    for (file in missing_files) {
      cat("  -", file, "\n")
    }
    return(FALSE)
  } else {
    cat("Todos los archivos de entrada están disponibles.\n")
    return(TRUE)
  }
}

cat("Configuración de rutas cargada exitosamente.\n")
cat("Use create_output_dirs() para crear directorios de salida.\n")
cat("Use check_input_files() para verificar archivos de entrada.\n")