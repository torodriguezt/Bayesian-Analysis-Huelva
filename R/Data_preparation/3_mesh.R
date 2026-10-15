################################################################################
# SPDE mesh and projector matrices for samples and grid.
################################################################################
source("R/Data_preparation/config.R")

samples <- read_csv(file.path(DIR_DATA, "samples.csv"), show_col_types = FALSE) %>%
  filter(site == "Huelva")
grid <- read_csv(file.path(DIR_DATA, "covariates_grid.csv"), show_col_types = FALSE)

coords_obs  <- as.matrix(samples[, c("x", "y")])
coords_grid <- as.matrix(grid[, c("x", "y")])

hull <- lapply(MESH$hull, function(d)
  fmesher::fm_nonconvex_hull(coords_obs, convex = d, format = "fm"))
mesh <- fmesher::fm_mesh_2d_inla(loc = coords_obs, boundary = hull,
                                 max.edge = MESH$max.edge, cutoff = MESH$cutoff)
A_obs  <- fmesher::fm_basis(mesh, loc = coords_obs)
A_grid <- fmesher::fm_basis(mesh, loc = coords_grid)
stopifnot(all(abs(Matrix::rowSums(A_obs) - 1) < 1e-8),
          all(abs(Matrix::rowSums(A_grid) - 1) < 1e-8))

## edge lengths near the samples
tv <- mesh$graph$tv
edge <- function(a, b) sqrt(rowSums((mesh$loc[tv[, a], 1:2] - mesh$loc[tv[, b], 1:2])^2))
tri_obs <- unique(fmesher::fm_bary(mesh, coords_obs)$index)
edges_obs <- c(edge(1, 2)[tri_obs], edge(2, 3)[tri_obs], edge(1, 3)[tri_obs])

saveRDS(list(mesh = mesh, A_obs = A_obs, A_grid = A_grid,
             coords_obs = coords_obs, coords_grid = coords_grid, settings = MESH),
        file.path(DIR_DATA, "mesh.rds"))
message(sprintf("Mesh: %d vertices; edges at the samples %.0f-%.0f m (median %.0f)",
                mesh$n, min(edges_obs), max(edges_obs), median(edges_obs)))
