# Paths, settings and helpers for data preparation and exploratory scripts.

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(readr); library(sf); library(INLA)
})

## inputs (data/raw/ is confidential)
FILE_SITES <- c(Huelva   = "data/raw/samples/Datos/HuelvaResultados.shp",
                Aroche   = "data/raw/samples/Datos/ArocheResultados.shp",
                Campillo = "data/raw/samples/Datos/CampilloResultados.shp")
FILE_GRID   <- "data/raw/prediction_grid/SinAgua/capa_sin_agua.shp"
DIR_LU_OBS  <- "data/raw/land_use/samples/Muestra_part2"
DIR_LU_GRID <- "data/raw/land_use/covariates/Covariables_muestra_part2"

## land-use classes and the NGR use of each
LAND_USE <- data.frame(
  class   = c("Industria", "Marsh", "Phospho", "Urban", "Agric", "Refineria", "Bare", "Park"),
  ngr_use = c("industrial", "other", "industrial", "urban", "other", "industrial",
              "other", "other"),
  obs     = c("Industry/Industry.shp", "Marsh/Marsh.shp", "Phospho/Phospho.shp",
              "Urban/Urban.shp", "Agric_2/Agric.shp", "Refinery/Refinery.shp",
              "Bare/Bare.shp", "Park/Park.shp"),
  grid    = c("Industria/Industry.shp", "Marsh/Marsh.shp", "Phospho/Phospho.shp",
              "Urban/Urban.shp", "Agric/Agric.shp", "Refinery/Refineria.shp",
              "Bare/Bare.shp", "Natural Park/Park.shp"))

## outputs
DIR_DATA <- "results/data"
dir.create(DIR_DATA, recursive = TRUE, showWarnings = FALSE)

CRS_UTM  <- 25829    # ETRS89 / UTM 29N
N_HUELVA <- 138

## NGR (mg/kg), Annex IV of Decreto 18/2015; Hg at the elemental level, no Cr
NGR <- data.frame(
  element    = c("As", "Pb", "Cu", "Zn", "Cd", "Tl", "Sb", "Hg", "Co", "V",
                 "Ni", "Mo", "Se", "Be", "Ba", "Sn"),
  industrial = c(40, 2750, 10000, 10000, 750, 2.3, 4700, 15, 250, 3650,
                 10000, 3900, 3900, 1450, 10000, 10000),
  urban      = c(36, 275, 3130, 10000, 75, 0.23, 470, 6, 25, 365,
                 1530, 390, 390, 145, 10000, 10000),
  other      = c(36, 275, 595, 10000, 25, 0.23, 90, 6, 24, 50,
                 1530, 115, 370, 145, 10000, 10000))

## candidate signatures (Ce left out: r = 0.985 with La)
SIGNATURES <- list(
  sulphide      = c("As", "Pb", "Cu", "Zn", "Cd", "Bi", "Sb", "Tl", "Hg", "Ag", "Sn"),
  lithogenic    = c("Th", "La", "Sc", "Al", "Zr", "Ti", "Rb", "Ga", "Nb"),
  phosphogypsum = c("P", "U", "Sr"))
MAX_PCT_AT_MIN <- 10   # max % of values at the detection floor

## SPDE mesh
MESH <- list(max.edge = c(350, 3000), cutoff = 75, hull = c(1500, 6000))

## sample columns are "symbol__unit_"
META_COLS <- c("MUESTRA_or", "Sample", "Analyte", "Campo", "X", "Y", "Z",
               "POINT_X", "POINT_Y")

element_unit <- function(col)
  ifelse(grepl("ppm", col), "ppm", ifelse(grepl("ppb", col), "ppb", "pct"))
to_ppm <- function(v, unit) switch(unit, ppm = v, ppb = v / 1000, pct = v * 1e4)

read_utm <- function(path) st_transform(st_zm(st_read(path, quiet = TRUE)), CRS_UTM)

## mean standardised log-concentration of a group of elements
group_score <- function(d, elements) rowMeans(scale(log(as.matrix(d[, elements]))))
