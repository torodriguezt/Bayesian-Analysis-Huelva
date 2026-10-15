################################################################################
# Annex IV elements: statistics, NGR and % of samples above it; land use.
################################################################################
dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)

samples <- read.csv("results/data/samples.csv", check.names = FALSE)
samples <- samples[samples$site == "Huelva", ]
cov  <- read.csv("results/data/covariates_obs.csv")
grid <- read.csv("results/data/covariates_grid.csv")
ngr  <- read.csv("results/data/ngr.csv")
stopifnot(nrow(samples) == nrow(cov), all(samples$id == cov$id))
uses <- c("industrial", "urban", "other")

x <- as.matrix(samples[, ngr$element])
## NGR per sample: its own land use, or one use for all
thr_applicable <- t(as.matrix(ngr[, cov$ngr_use]))
thr_use <- function(u) matrix(ngr[[u]], nrow(x), ncol(x), byrow = TRUE)
pct <- function(t) 100 * colMeans(x > t)
pct_any <- function(t) 100 * mean(rowSums(x > t) > 0)

tab <- data.frame(element = ngr$element, median = apply(x, 2, median), mean = colMeans(x),
                  sd = apply(x, 2, sd), cv = apply(x, 2, sd) / colMeans(x),
                  max = apply(x, 2, max),
                  setNames(ngr[, uses], paste0("ngr_", uses)),
                  pct_applicable = pct(thr_applicable),
                  setNames(as.data.frame(lapply(uses, function(u) pct(thr_use(u)))),
                           paste0("pct_", uses)),
                  row.names = NULL)
any_row <- tab[1, ]
any_row[] <- NA
any_row$element <- "any"
any_row$pct_applicable <- pct_any(thr_applicable)
for (u in uses) any_row[[paste0("pct_", u)]] <- pct_any(thr_use(u))
tab <- rbind(tab, any_row)
write.csv(tab, "results/tables/descriptive.csv", row.names = FALSE)

lu <- merge(as.data.frame(table(land_use = cov$land_use), responseName = "n_samples"),
            as.data.frame(table(land_use = grid$land_use), responseName = "n_grid"),
            all = TRUE)
lu[is.na(lu)] <- 0
lu$pct_grid <- 100 * lu$n_grid / nrow(grid)
lu <- merge(unique(grid[, c("land_use", "ngr_use")]), lu, all.y = TRUE)
write.csv(lu[, c("land_use", "ngr_use", "n_samples", "pct_grid")],
          "results/tables/land_use.csv", row.names = FALSE)
print(tab[, c("element", "max", "ngr_urban", "pct_applicable", paste0("pct_", uses))],
      digits = 3, row.names = FALSE)
