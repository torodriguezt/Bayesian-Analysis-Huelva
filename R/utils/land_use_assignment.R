# Infer a land-use category from validated distances. This records proximity
# evidence, not verification of the legal land use at a location.
assign_land_use <- function(distances, class_ngr_use, tolerance_m = 1e-5) {
  distances <- as.matrix(distances)
  if (!is.numeric(distances) || !nrow(distances) || !ncol(distances) ||
      any(!is.finite(distances)) || any(distances < 0))
    stop("Land-use distances must be a nonempty, finite, nonnegative matrix.")
  if (is.null(colnames(distances)) || anyDuplicated(colnames(distances)) ||
      is.null(names(class_ngr_use)) || anyDuplicated(names(class_ngr_use)) ||
      !all(colnames(distances) %in% names(class_ngr_use)))
    stop("Every distance column must have a unique named NGR-use mapping.")
  if (length(tolerance_m) != 1L || !is.finite(tolerance_m) || tolerance_m < 0)
    stop("tolerance_m must be a finite nonnegative number.")
  categories <- unname(class_ngr_use[colnames(distances)])
  if (anyNA(categories) || !all(categories %in% c("industrial", "urban", "other")))
    stop("Invalid NGR-use mapping.")

  nearest <- max.col(-distances, ties.method = "first")
  minimum <- distances[cbind(seq_len(nrow(distances)), nearest)]
  tied <- abs(distances - minimum) <= tolerance_m
  n_ties <- rowSums(tied)
  ambiguous <- which(n_ties > 1L)
  conflicting <- ambiguous[vapply(ambiguous, function(i)
    length(unique(categories[tied[i, ]])) > 1L, logical(1))]
  if (length(conflicting)) stop(
    "Equally close land-use classes imply different NGR uses in ",
    length(conflicting), " location(s); resolve their land use before applying thresholds.")

  method <- ifelse(minimum <= tolerance_m, "zero_distance", "nearest_proxy")
  data.frame(land_use = colnames(distances)[nearest],
             ngr_use = categories[nearest],
             land_use_distance_m = minimum,
             land_use_method = method,
             land_use_n_ties = n_ties,
             ngr_use_basis = ifelse(method == "zero_distance",
                                    "land_use_at_zero_distance", "nearest_land_use_proxy"),
             stringsAsFactors = FALSE)
}
