# Distances in metres from points to a land-use class.
#
# Polygon layers are measured directly. The Huelva export instead contains one
# LINESTRING per source point: its first vertex is that point, its last vertex
# is the nearest point of the land-use class, and `distance` is the distance in
# ETRS89 / UTM 29N metres. These lines are connectors, not land-use boundaries.
# Coordinate matching deliberately ignores row order and non-unique source IDs.
land_use_distances <- function(points_sf, features_sf, tolerance_m = 0.01,
                               metric_crs = 25829) {
  if (!inherits(points_sf, "sf") || !inherits(features_sf, "sf"))
    stop("points_sf and features_sf must be sf objects.")
  if (is.na(sf::st_crs(points_sf)) || is.na(sf::st_crs(features_sf)))
    stop("Both input layers must have a declared CRS.")
  if (length(tolerance_m) != 1L || !is.finite(tolerance_m) || tolerance_m <= 0)
    stop("tolerance_m must be a positive finite distance in metres.")
  target_crs <- sf::st_crs(metric_crs)
  if (is.na(target_crs) || isTRUE(sf::st_is_longlat(target_crs)) ||
      !tolower(target_crs$units_gdal) %in% c("metre", "meter", "metres", "meters", "m"))
    stop("metric_crs must be a projected CRS with metre units.")
  if (!nrow(features_sf)) stop("The land-use layer is empty.")
  if (any(sf::st_is_empty(features_sf)) || any(sf::st_is_empty(points_sf)))
    stop("Empty geometries cannot be used for land-use distances.")
  if (!all(sf::st_geometry_type(points_sf) == "POINT"))
    stop("The source layer must contain POINT geometries.")

  kind <- as.character(sf::st_geometry_type(features_sf))
  polygons <- all(kind %in% c("POLYGON", "MULTIPOLYGON"))
  connectors <- all(kind == "LINESTRING")
  if (!polygons && !connectors)
    stop("Land-use features must be polygons or LINESTRING connectors with distance.")

  points <- sf::st_transform(sf::st_zm(sf::st_geometry(points_sf)), target_crs)
  features <- sf::st_transform(sf::st_zm(sf::st_geometry(features_sf)), target_crs)
  if (!nrow(points_sf)) return(numeric())

  if (polygons) {
    boundary <- sf::st_union(sf::st_make_valid(features))
    return(as.numeric(sf::st_distance(points, boundary)))
  }

  if (!"distance" %in% names(features_sf))
    stop("LINESTRING connectors require a distance attribute in metres.")
  distance <- features_sf$distance
  if (inherits(distance, "units")) distance <- units::set_units(distance, "m")
  if (!is.numeric(distance) || any(!is.finite(distance)) || any(distance < 0))
    stop("Connector distances must be finite, non-negative numeric metres.")
  distance <- as.numeric(distance)
  connector_length <- as.numeric(sf::st_length(features))
  if (any(abs(distance - connector_length) > pmax(tolerance_m, distance * 1e-6)))
    stop("Connector distance does not match its length in metric_crs; check units and CRS.")

  origin_xy <- t(vapply(features, function(g) g[1L, 1:2], numeric(2)))
  if (any(!is.finite(origin_xy))) stop("Connector origins must have finite coordinates.")
  if (anyDuplicated(as.data.frame(origin_xy)))
    stop("Duplicate connector origins make coordinate matching ambiguous.")
  origins <- sf::st_sfc(lapply(seq_len(nrow(origin_xy)), function(i)
    sf::st_point(origin_xy[i, ])), crs = target_crs)

  # GEOS indexes intersections, avoiding an all-pairs distance matrix for the
  # 22,426 grid points. The diamond buffer contains the tolerance circle; exact
  # distances then discard the extra corners before checking uniqueness.
  nearby <- sf::st_intersects(
    sf::st_buffer(points, dist = tolerance_m * sqrt(2), nQuadSegs = 1), origins)
  point_index <- rep(seq_along(nearby), lengths(nearby))
  origin_index <- unlist(nearby, use.names = FALSE)
  if (length(origin_index)) {
    gap <- as.numeric(sf::st_distance(points[point_index], origins[origin_index],
                                      by_element = TRUE))
    keep <- gap <= tolerance_m
    point_index <- point_index[keep]
    origin_index <- origin_index[keep]
  }
  n_matches <- tabulate(point_index, nbins = length(points))
  if (any(n_matches == 0L))
    stop(sprintf("Missing connector origins for %d source point(s).", sum(n_matches == 0L)))
  if (any(n_matches > 1L))
    stop("Multiple connector origins match a source point within tolerance_m.")
  if (anyDuplicated(origin_index))
    stop("Duplicate source points match the same connector origin.")

  matched <- integer(length(points))
  matched[point_index] <- origin_index
  distance[matched]
}
