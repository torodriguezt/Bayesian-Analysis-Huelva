################################################################################
# Land use of each sample and grid point: distances, class and NGR use.
################################################################################
source("R/Data_preparation/config.R")
source("R/utils/land_use.R")
source("R/utils/land_use_assignment.R")

samples <- read_csv(file.path(DIR_DATA, "samples.csv"), show_col_types = FALSE) %>%
  filter(site == "Huelva")
obs <- st_as_sf(samples[, c("id", "x", "y")], coords = c("x", "y"), crs = CRS_UTM,
                remove = FALSE)

grid <- read_utm(FILE_GRID)
if (!all(st_geometry_type(grid) == "POINT")) grid <- st_centroid(grid)
xy_grid <- st_coordinates(grid)
grid <- st_as_sf(data.frame(x = xy_grid[, 1], y = xy_grid[, 2]), coords = c("x", "y"),
                 crs = CRS_UTM, remove = FALSE)

dist_to <- function(pts, file)
  land_use_distances(pts, read_utm(file))

## distances (m), points x classes
raw_obs  <- sapply(LAND_USE$obs,  function(f) dist_to(obs,  file.path(DIR_LU_OBS,  f)))
raw_grid <- sapply(LAND_USE$grid, function(f) dist_to(grid, file.path(DIR_LU_GRID, f)))
colnames(raw_obs) <- colnames(raw_grid) <- LAND_USE$class

cov_obs  <- st_drop_geometry(obs)
cov_grid <- st_drop_geometry(grid)
for (cl in LAND_USE$class) {
  l_obs <- log1p(raw_obs[, cl]); l_grid <- log1p(raw_grid[, cl])
  scale_sd <- sd(l_obs)
  if (!is.finite(scale_sd) || scale_sd == 0) scale_sd <- 1
  cov_obs[[paste0("distance_", cl, "_m")]] <- raw_obs[, cl]
  cov_grid[[paste0("distance_", cl, "_m")]] <- raw_grid[, cl]
  cov_obs[[paste0("d_", cl)]]  <- (l_obs  - mean(l_obs)) / scale_sd
  cov_grid[[paste0("d_", cl)]] <- (l_grid - mean(l_obs)) / scale_sd
}

## Distances identify a zero-distance class or a nearest-use proxy. Preserve
## the evidence and reject ties between different NGR categories. Proximity
## alone does not verify the legal land use at a location.
ngr_by_class <- setNames(LAND_USE$ngr_use, LAND_USE$class)
cov_obs <- cbind(cov_obs, assign_land_use(raw_obs, ngr_by_class))
cov_grid <- cbind(cov_grid, assign_land_use(raw_grid, ngr_by_class))
stopifnot(!anyNA(cov_obs), !anyNA(cov_grid))

write_csv(cov_obs,  file.path(DIR_DATA, "covariates_obs.csv"))
write_csv(cov_grid, file.path(DIR_DATA, "covariates_grid.csv"))
message(sprintf("Land use: %d samples (%s) and %d grid points",
                nrow(cov_obs),
                paste(names(table(cov_obs$ngr_use)), table(cov_obs$ngr_use), collapse = ", "),
                nrow(cov_grid)))
message(sprintf("Nearest-use proxies: %d samples, %d grid points (see land_use_method)",
                sum(cov_obs$land_use_method == "nearest_proxy"),
                sum(cov_grid$land_use_method == "nearest_proxy")))
