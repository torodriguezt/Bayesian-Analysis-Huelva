# single visual system, sourced at the start of every figure script
library(ggplot2)

## signatures: validated pair (deutan dE 10.1, normal dE 22.3, >= 3:1 on white)
COL_SIG  <- c(sulphide = "#c8664a",   # terracotta
              lithogenic = "#2a9d8f") # teal
COL_NEUTRAL <- "grey45"               # signatures outside the model
LAB_SIG  <- c(sulphide = "Sulphide", lithogenic = "Lithogenic",
              phosphogypsum = "Phosphogypsum")

## NGR land use: neutral greys, darker = more permissive level
COL_USE <- c(other = "#e1e0d9", urban = "#c3c2b7", industrial = "#898781")
LAB_USE <- c(industrial = "Industrial", urban = "Urban", other = "Other uses")

## sequential (probability) and diverging (log-scale deviations) ramps
SEQ_BLUE <- c("#cde2fb", "#9ec5f4", "#6da7ec", "#3987e5", "#256abf", "#184f95",
              "#0d366b")
DIV_LOW <- "#2a78d6"; DIV_MID <- "#f0efec"; DIV_HIGH <- "#e34948"

## marks: thin lines and small points, drawn at roughly their printed size
LINE_W <- 0.45
PT_SZ  <- 1.5

## theme for axis-based plots: close to base R graphics, no chrome
theme_paper <- function(base = 10, legend = c("inside", "top", "none")) {
  legend <- match.arg(legend)
  th <- theme_classic(base_size = base) +
    theme(
      text = element_text(color = "black"),
      plot.title = element_text(size = base + 1, color = "black", face = "bold",
                                hjust = 0.5, margin = margin(b = 7)),
      axis.line  = element_line(color = "black", linewidth = 0.4),
      axis.ticks = element_line(color = "black", linewidth = 0.4),
      axis.ticks.length = unit(3.5, "pt"),
      axis.title = element_text(size = base - 1, color = "black"),
      axis.text  = element_text(size = base - 2, color = "black"),
      panel.grid = element_blank(),
      legend.background = element_blank(),
      legend.key        = element_blank(),
      legend.title      = element_blank(),
      legend.text       = element_text(size = base - 2, color = "black"),
      legend.key.width  = unit(0.5, "cm"),
      legend.key.height = unit(0.42, "cm"),
      legend.margin     = margin(0, 0, 0, 0),
      strip.background = element_blank(),
      strip.text = element_text(size = base - 1, color = "black", face = "bold",
                                margin = margin(2, 2, 6, 2)),
      panel.spacing = unit(11, "pt"),
      plot.background = element_rect(fill = "white", color = NA),
      plot.margin = margin(10, 14, 8, 10)
    )
  th + switch(
    legend,
    inside = theme(legend.position = "inside",
                   legend.position.inside = c(0.01, 1),
                   legend.justification = c(0, 1)),
    top    = theme(legend.position = "top", legend.justification = "left"),
    none   = theme(legend.position = "none")
  )
}

## theme for maps
theme_map <- function(base = 10) {
  theme_void(base_size = base) +
    theme(
      plot.title = element_text(size = base + 1, color = "black", face = "bold",
                                hjust = 0.5, margin = margin(b = 5)),
      legend.position = "right",
      legend.title = element_text(size = base - 1, color = "grey20"),
      legend.text  = element_text(size = base - 2, color = "grey20"),
      legend.key.height = unit(0.8, "cm"),
      legend.key.width  = unit(0.3, "cm"),
      strip.text = element_text(size = base, color = "black", face = "bold",
                                margin = margin(3, 3, 5, 3)),
      plot.background = element_rect(fill = "white", color = NA),
      plot.margin = margin(6, 6, 6, 6)
    )
}

## grid cells as tiles: the prediction grid is not axis-aligned in UTM, so
## tiles are drawn slightly larger than the cell spacing to leave no gaps
grid_spacing <- function(x, y, n = 500) {
  i <- sample.int(length(x), min(n, length(x)))
  median(sapply(i, function(k) {
    d <- sqrt((x - x[k])^2 + (y - y[k])^2); min(d[d > 0]) }))
}
geom_cells <- function(data, mapping, spacing, ...)
  geom_tile(data = data, mapping = mapping, width = 1.12 * spacing,
            height = 1.12 * spacing, ...)

## saves a PNG (preview) and a vector PDF (article); cairo_pdf embeds unicode
save_fig <- function(name, plot, w, h, dir_fig = "figures") {
  dir.create(dir_fig, showWarnings = FALSE)
  ggsave(file.path(dir_fig, name), plot, width = w, height = h,
         dpi = 320, bg = "white")
  ggsave(file.path(dir_fig, sub("\\.png$", ".pdf", name)), plot, width = w,
         height = h, bg = "white", device = cairo_pdf)
}
