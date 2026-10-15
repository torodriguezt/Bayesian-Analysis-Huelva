################################################################################
# Samples of the three surveys in mg/kg, element quality and NGR table.
################################################################################
source("R/Data_preparation/config.R")

read_site <- function(site) {
  x  <- read_utm(FILE_SITES[[site]])
  df <- st_drop_geometry(x)
  cols <- setdiff(names(df), META_COLS)
  cols <- cols[grepl("^[A-Z][a-z]?_", cols)]
  vals <- lapply(cols, function(cl) {
    v <- suppressWarnings(as.numeric(df[[cl]]))
    v[v <= 0] <- NA                                # missing or below detection
    to_ppm(v, element_unit(cl))
  })
  names(vals) <- sub("_.*$", "", cols)
  xy <- st_coordinates(x)
  data.frame(site = site, id = seq_len(nrow(x)), x = xy[, 1], y = xy[, 2],
             vals, check.names = FALSE)
}

samples <- bind_rows(lapply(names(FILE_SITES), read_site))
elements <- setdiff(names(samples), c("site", "id", "x", "y"))

huelva <- samples[samples$site == "Huelva", ]
stopifnot(nrow(huelva) == N_HUELVA, all(NGR$element %in% elements),
          !anyNA(huelva[, NGR$element]))

## ties at the minimum flag a detection floor
quality <- samples %>%
  pivot_longer(all_of(elements), names_to = "element", values_to = "ppm") %>%
  group_by(site, element) %>%
  summarise(n = n(), n_missing = sum(is.na(ppm)),
            min_ppm = suppressWarnings(min(ppm, na.rm = TRUE)),
            median_ppm = median(ppm, na.rm = TRUE),
            max_ppm = suppressWarnings(max(ppm, na.rm = TRUE)),
            n_at_min = sum(ppm == min_ppm, na.rm = TRUE),
            n_at_max = sum(ppm == max_ppm, na.rm = TRUE), .groups = "drop") %>%
  mutate(pct_at_min = 100 * n_at_min / n)

write_csv(samples, file.path(DIR_DATA, "samples.csv"))
write_csv(quality, file.path(DIR_DATA, "element_quality.csv"))
write_csv(NGR, file.path(DIR_DATA, "ngr.csv"))
message(sprintf("Samples: %s; %d elements",
                paste(sprintf("%s %d", names(table(samples$site)), table(samples$site)),
                      collapse = ", "), length(elements)))
