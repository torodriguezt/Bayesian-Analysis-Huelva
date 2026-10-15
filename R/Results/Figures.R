################################################################################
# Figures of the article. PNG and PDF in figures/.
# Needs Exploratory.R, posterior_draws.R and the result tables.
################################################################################
source("R/All_Models/model_data.R")
suppressPackageStartupMessages({
  library(patchwork); library(sf); library(gstat); library(ggspatial)
})
source("R/Results/figure_style.R")

grid    <- read.csv("results/data/covariates_grid.csv")
cov_obs <- read.csv("results/data/covariates_obs.csv")
ngr     <- read.csv("results/data/ngr.csv")
draws   <- readRDS("results/posterior/F3_draws.rds")
stopifnot(all(cov_obs$id == samples$id), identical(draws$input_md5, MODEL_INPUT_MD5))
set.seed(2026)
SPACING <- grid_spacing(grid$x, grid$y)
UTM <- st_crs(25829)
map_coord <- function() coord_sf(crs = UTM, datum = NA, expand = FALSE)
scale_bar <- function() annotation_scale(location = "bl", width_hint = 0.3,
  style = "ticks", line_width = 0.5, text_cex = 0.6, height = unit(0.12, "cm"))
sig_levels <- function(x) factor(unname(LAB_SIG[x]), levels = LAB_SIG)


## Fig 1: land use (NGR category), samples and SPDE mesh --------------------
grid$use <- factor(grid$ngr_use, levels = c("industrial", "urban", "other"))
p_use <- ggplot() +
  geom_cells(grid, aes(x, y, fill = use), SPACING) +
  geom_point(data = cov_obs, aes(x, y), size = 0.7, color = "black") +
  scale_fill_manual(values = COL_USE, labels = LAB_USE, name = "Land use (NGR)") +
  map_coord() + scale_bar() +
  labs(title = sprintf("Land use and the %d samples", nrow(cov_obs))) + theme_map()

tv <- geo$mesh$graph$tv; loc <- geo$mesh$loc
e  <- unique(t(apply(rbind(tv[, 1:2], tv[, 2:3], tv[, c(1, 3)]), 1, sort)))
edges <- data.frame(x = loc[e[, 1], 1], y = loc[e[, 1], 2],
                    xend = loc[e[, 2], 1], yend = loc[e[, 2], 2])
p_mesh <- ggplot() +
  geom_cells(grid, aes(x, y), SPACING, fill = COL_USE[["other"]]) +
  geom_segment(data = edges, aes(x, y, xend = xend, yend = yend),
               linewidth = 0.08, color = "grey35") +
  geom_point(data = cov_obs, aes(x, y), size = 0.5, color = DIV_HIGH) +
  map_coord() + scale_bar() +
  labs(title = sprintf("SPDE mesh (%s vertices)", format(geo$mesh$n, big.mark = ","))) +
  theme_map()
save_fig("fig01_study_area.png", (p_use | p_mesh) + plot_annotation(tag_levels = "a"),
         7, 4.7)


## Fig 2: element signatures (exploratory factor analysis) -------------------
fa <- read.csv("results/tables/exploratory/fa_loadings.csv")
facs <- c("sulphide", "lithogenic", "phosphogypsum")
fa$label <- paste0(fa$element, ifelse(fa$element %in% ELEMENTS, "*", ""))
own <- as.matrix(fa[, facs])[cbind(seq_len(nrow(fa)), match(fa$signature, facs))]
lev <- rev(fa$label[order(match(fa$signature, facs), -own)])
long <- do.call(rbind, lapply(facs, function(f) data.frame(
  label = factor(fa$label, levels = lev), signature = sig_levels(fa$signature),
  factor = sig_levels(f), loading = fa[[f]])))
p <- ggplot(long, aes(factor, label, fill = loading)) +
  geom_tile(color = "white", linewidth = 0.7) +
  geom_text(data = subset(long, abs(loading) >= 0.4),
            aes(label = sprintf("%.2f", loading), color = abs(loading) > 0.75),
            size = 2.3, show.legend = FALSE) +
  scale_color_manual(values = c(`TRUE` = "white", `FALSE` = "black")) +
  scale_fill_gradient2(low = DIV_LOW, mid = DIV_MID, high = DIV_HIGH, midpoint = 0,
                       limits = c(-1.1, 1.1), breaks = c(-1, -0.5, 0, 0.5, 1),
                       name = "Loading") +
  scale_x_discrete(position = "top") +
  facet_grid(signature ~ ., scales = "free_y", space = "free_y") +
  labs(x = NULL, y = NULL) + theme_paper() +
  theme(legend.position = "right", legend.title = element_text(size = 8),
        axis.line = element_blank(), axis.ticks = element_blank(),
        strip.text.y = element_text(angle = 0, hjust = 0))
save_fig("fig02_signatures.png", p, 4.6, 5.2)


## Fig 3: semivariograms of the signature scores -----------------------------
scores <- read.csv("results/data/signature_scores.csv")
vpar   <- read.csv("results/tables/exploratory/scale_by_signature.csv")
pts <- st_as_sf(scores, coords = c("x", "y"), crs = UTM)
matern1 <- function(h, a) ifelse(h == 0, 1, (h / a) * besselK(pmax(h, 1e-9) / a, 1))
emp <- do.call(rbind, lapply(facs, function(s) {
  pts$z <- scores[[s]]
  v <- variogram(z ~ 1, pts, cutoff = 8000, width = 500)
  data.frame(score = s, dist = v$dist / 1000, gamma = v$gamma / var(scores[[s]]))
}))
h <- seq(5, 8000, length.out = 300)
fitted <- do.call(rbind, lapply(facs, function(s) {
  r <- vpar[vpar$score == s, ]
  data.frame(score = s, dist = h / 1000, gamma = (r$nugget + r$partial_sill *
    (1 - matern1(h, r$practical_range_m / sqrt(8)))) / r$variance)
}))
lab <- vpar[match(facs, vpar$score), ]
lab$text <- sprintf("nugget fraction %.2f\n90%% interval [%.2f, %.2f]",
                    lab$nugget_fraction, lab$nf_q05, lab$nf_q95)
for (d in c("emp", "fitted", "lab")) assign(d, transform(get(d), score = sig_levels(score)))
col_score <- setNames(c(COL_SIG, COL_NEUTRAL), LAB_SIG)
p <- ggplot(emp, aes(dist, gamma, color = score)) +
  geom_hline(yintercept = 1, color = "grey80", linewidth = 0.35) +
  geom_line(data = fitted, linewidth = LINE_W) +
  geom_point(size = PT_SZ) +
  geom_text(data = lab, aes(x = 8, y = 0.03, label = text), hjust = 1, vjust = 0,
            size = 2.4, color = "grey20", lineheight = 0.95) +
  facet_wrap(~score) +
  scale_color_manual(values = col_score, guide = "none") +
  scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.05))) +
  labs(x = "Distance (km)", y = "Semivariance / variance") + theme_paper()
save_fig("fig03_variograms.png", p, 7, 2.6)


## Fig 4: co-exceedance of the NGR at the samples ----------------------------
ce <- read.csv("results/tables/exploratory/coexceedance.csv")
d  <- subset(ce, site == "Huelva" & ngr_use == "applicable")
d  <- d[d$n_elements <= max(d$n_elements[d$observed > 0 | d$independent >= 0.5 |
                                           d$gaussian_copula >= 0.5]), ]
ref <- rbind(data.frame(n = d$n_elements, count = d$independent, ref = "Independent elements"),
             data.frame(n = d$n_elements, count = d$gaussian_copula, ref = "Gaussian copula"))
p <- ggplot() +
  geom_col(data = d, aes(n_elements, observed, fill = "Observed"), width = 0.62) +
  geom_point(data = ref, aes(n, count, shape = ref), size = 1.9,
             position = position_dodge(width = 0.32)) +
  scale_fill_manual(values = c(Observed = COL_USE[["urban"]])) +
  scale_shape_manual(values = c("Independent elements" = 1, "Gaussian copula" = 17)) +
  scale_x_continuous(breaks = d$n_elements) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  labs(x = "Number of elements above their NGR", y = "Number of samples") +
  theme_paper() + theme(legend.position.inside = c(0.98, 0.98),
                        legend.justification = c(1, 1))
save_fig("fig04_coexceedance.png", p, 3.6, 2.8)


## Fig 5: loadings of F3 --------------------------------------------------------
ld <- read.csv("results/tables/loadings.csv")
ld$label <- factor(ld$element, levels = rev(ld$element[order(match(ld$signature, facs),
                                                             -ld$loading)]))
ld$signature <- sig_levels(ld$signature)
ld$kind <- ifelse(ld$anchor, "Anchor (fixed at 1)", "Estimated, 95% interval")
p <- ggplot(ld, aes(loading, label, color = signature)) +
  geom_vline(xintercept = 1, color = "grey80", linewidth = 0.35) +
  geom_errorbar(aes(xmin = q025, xmax = q975), width = 0, linewidth = 0.5,
                orientation = "y", na.rm = TRUE) +
  geom_point(aes(shape = kind), size = 1.9, fill = "white") +
  facet_grid(signature ~ ., scales = "free_y", space = "free_y") +
  scale_color_manual(values = setNames(COL_SIG, LAB_SIG[names(COL_SIG)]), guide = "none") +
  scale_shape_manual(values = c("Anchor (fixed at 1)" = 21, "Estimated, 95% interval" = 16)) +
  labs(x = "Loading on the signature score", y = NULL) +
  theme_paper(legend = "top") + theme(strip.text.y = element_text(angle = 0, hjust = 0))
save_fig("fig05_loadings.png", p, 3.8, 4.2)


## Fig 6: realized mappable fraction ---------------------------------------
pr <- draws$phi_realized
dens <- do.call(rbind, lapply(colnames(pr), function(f) {
  k <- density(pr[, f], from = 0, to = 1, adjust = 1.2)
  data.frame(signature = f, x = k$x, y = k$y)
}))
dens$signature <- sig_levels(dens$signature)
peaks <- do.call(rbind, lapply(split(dens, dens$signature), function(s) s[which.max(s$y), ]))
peaks$hjust <- ifelse(peaks$x < 0.7, 1, 0)   # labels away from each other
peaks$x <- peaks$x + ifelse(peaks$x < 0.7, -0.02, 0.02)
diff_sl <- pr[, "sulphide"] - pr[, "lithogenic"]
kd <- density(diff_sl, adjust = 1.2)
p_phi <- ggplot(dens, aes(x, y, color = signature, fill = signature)) +
  geom_area(alpha = 0.10, linewidth = 0, position = "identity") +
  geom_line(linewidth = 0.6) +
  geom_text(data = peaks, aes(label = signature, hjust = hjust), vjust = 0, size = 2.8,
            color = "black") +
  scale_color_manual(values = setNames(COL_SIG, LAB_SIG[names(COL_SIG)]), guide = "none") +
  scale_fill_manual(values = setNames(COL_SIG, LAB_SIG[names(COL_SIG)]), guide = "none") +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
  labs(x = "Realized mappable fraction", y = "Posterior density") + theme_paper()
p_diff <- ggplot(data.frame(x = kd$x, y = kd$y), aes(x, y)) +
  geom_vline(xintercept = 0, color = "grey70", linewidth = 0.35) +
  geom_area(fill = "grey20", alpha = 0.10) + geom_line(linewidth = 0.6) +
  annotate("text", x = min(kd$x), y = 1.22 * max(kd$y), hjust = 0, vjust = 1, size = 2.6,
           label = sprintf("P(sulphide < lithogenic) = %.3f", mean(diff_sl < 0))) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.02))) +
  labs(x = "Difference, sulphide - lithogenic", y = NULL) + theme_paper()
save_fig("fig06_mappable_fraction.png",
         (p_phi | p_diff) + plot_layout(widths = c(1.4, 1)) +
           plot_annotation(tag_levels = "a"), 6.6, 2.7)


## Fig 7: smooth (mappable) and patchy parts of each signature --------------
maps <- lapply(names(FACTORS), function(f) {
  smooth <- as.vector(geo$A_grid %*% rowMeans(draws$g[[f]]))
  patchy <- rowMeans(draws$u[[f]] - as.matrix(geo$A_obs %*% draws$g[[f]]))
  lim <- max(abs(c(smooth, patchy)))
  fill <- scale_fill_gradient2(low = DIV_LOW, mid = DIV_MID, high = DIV_HIGH, midpoint = 0,
                               limits = c(-lim, lim), name = "Log scale")
  p_s <- ggplot() + geom_cells(transform(grid, v = smooth), aes(x, y, fill = v), SPACING) +
    fill + map_coord() + scale_bar() +
    labs(title = sprintf("%s, smooth part", LAB_SIG[[f]])) + theme_map()
  p_p <- ggplot() +
    geom_cells(grid, aes(x, y), SPACING, fill = COL_USE[["other"]]) +
    geom_point(data = transform(cov_obs, v = patchy), aes(x, y, fill = v), shape = 21,
               size = 1.9, stroke = 0.25, color = "white") +
    fill + map_coord() +
    labs(title = sprintf("%s, patchy part (samples)", LAB_SIG[[f]])) + theme_map()
  (p_s | p_p) + plot_layout(guides = "collect")
})
save_fig("fig07_scale_maps.png", (maps[[1]] / maps[[2]]) +
           plot_annotation(tag_levels = "a"), 7, 9.6)


## Fig 8: probability of exceeding any NGR, map vs inspection ---------------
cg  <- read.csv("results/posterior/compliance_grid.csv")
reg <- intersect(ngr$element, ELEMENTS)
thr <- t(as.matrix(ngr[match(reg, ngr$element), cov_obs$ngr_use]))
cov_obs$exceeds <- rowSums(as.matrix(samples[, reg]) > thr) > 0
long <- rbind(data.frame(cg[, c("x", "y")], p = cg$p_map, what = "Map (smooth part)"),
              data.frame(cg[, c("x", "y")], p = cg$p_inspection,
                         what = "Inspection (new sample)"))
long$what <- factor(long$what, levels = c("Map (smooth part)", "Inspection (new sample)"))
p <- ggplot() +
  geom_cells(long, aes(x, y, fill = p), SPACING) +
  geom_point(data = subset(cov_obs, exceeds), aes(x, y, shape = "Sample above an NGR"),
             size = 1.3, stroke = 0.45, color = "black") +
  facet_wrap(~what) +
  scale_fill_gradientn(colours = SEQ_BLUE, limits = c(0, 1), name = "Probability") +
  scale_shape_manual(values = c("Sample above an NGR" = 1), name = NULL) +
  map_coord() + scale_bar() + theme_map() +
  theme(legend.position = "right")
save_fig("fig08_compliance_maps.png", p, 7, 5.2)

message("Figures written to figures/")
