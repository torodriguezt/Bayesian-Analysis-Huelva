################################################################################
# Element signatures in Huelva: factor analysis and a score per signature.
################################################################################
source("R/Data_preparation/config.R")
dir_tab <- "results/tables/exploratory"
dir.create(dir_tab, recursive = TRUE, showWarnings = FALSE)

samples <- read_csv(file.path(DIR_DATA, "samples.csv"), show_col_types = FALSE) %>%
  filter(site == "Huelva")
quality <- read_csv(file.path(DIR_DATA, "element_quality.csv"), show_col_types = FALSE) %>%
  filter(site == "Huelva")

## drop elements with missing values or many at the detection floor
ok <- quality$element[quality$n_missing == 0 & quality$pct_at_min <= MAX_PCT_AT_MIN]
sig <- lapply(SIGNATURES, intersect, ok)
message("Screened out (detection floor or missing): ",
        paste(setdiff(unlist(SIGNATURES), ok), collapse = ", "))

## factor analysis
el <- unlist(sig, use.names = FALSE)
Z  <- scale(log(as.matrix(samples[, el])))
fa <- factanal(Z, factors = length(sig), rotation = "promax")
L  <- unclass(fa$loadings)
colnames(L) <- make.unique(apply(L, 2, function(l)
  names(sig)[which.max(sapply(sig, function(s) mean(l[s])))]))
fa_tab <- data.frame(element = el, signature = rep(names(sig), lengths(sig)),
                     round(L, 3), uniqueness = round(fa$uniquenesses, 3),
                     check.names = FALSE)
write_csv(fa_tab, file.path(dir_tab, "fa_loadings.csv"))

## scores; witness = sulphide elements without NGR (Bi, Ag)
groups <- c(sig, list(regulated = intersect(sig$sulphide, NGR$element),
                      witness   = setdiff(sig$sulphide, NGR$element)))
scores <- data.frame(samples[, c("id", "x", "y")],
                     sapply(groups, function(s) group_score(samples, s)))
write_csv(scores, file.path(DIR_DATA, "signature_scores.csv"))
write_csv(data.frame(score = names(groups), round(cor(scores[, names(groups)]), 3)),
          file.path(dir_tab, "score_correlation.csv"))

## anchor candidates
own <- sapply(names(sig), function(s) {
  col <- grep(paste0("^", s), colnames(L))[1]
  if (is.na(col)) NA_character_ else sig[[s]][which.max(L[sig[[s]], col])]
})
message(sprintf("Factor analysis: %d elements, %d factors; anchor candidates %s",
                length(el), ncol(L), paste(names(own), own, sep = " = ", collapse = ", ")))
