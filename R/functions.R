# Shared data loading and functions for the reproducibility scripts.
# Source this file from the repository root: source("R/functions.R")

suppressPackageStartupMessages({
  library(tidyverse)
  library(vegan)
})

# Data ---------------------------------------------------------------------------------------------

# Sample metadata: one row per sample (sampling unit x date).
# `unit` and `location` use the names of the manuscript (Fig. 1). The factor levels keep the order
# of the original analysis, because the permutation tests draw random numbers in that order: with
# a different order the p-values would differ within Monte Carlo error.
load_samples <- function() {
  read_csv("data/samples.csv", show_col_types = FALSE) %>%
    mutate(date = as.Date(date),
           unit = factor(unit, levels = read_lines("data/unit_levels.txt")))
}

# The five long-term locations (sampled before the refilling), in the order of the original analysis
long_term_locations <- c("CLM", "TLE", "TLR", "TSP", "DLM")

# Dissimilarity matrices between samples, one per community descriptor
descriptor_names <- c("traits_cells", "traits_biovolume", "traits_cells_dom", "traits_cells_klais",
                      "traits_cells_rimet", "taxa_genus", "taxa_genus_biovol", "taxa_class")
main_descriptors <- c("traits_cells", "traits_biovolume", "taxa_genus", "taxa_genus_biovol", "taxa_class")

load_descriptors <- function() {
  set_names(descriptor_names) %>%
    map(~ read.csv(file.path("data/distances", paste0(.x, ".csv")), row.names = 1, check.names = FALSE) %>%
          as.matrix() %>% as.dist())
}

# Community-weighted means of traits (cells or biovolume weighting) and trait types
load_cwm <- function(weighting) {
  read_csv(file.path("data", paste0("cwm_", weighting, ".csv")), show_col_types = FALSE) %>% select(-sample)
}
trait_types <- read_lines("data/traits/trait_types.txt") %>% str_split_fixed(": ", 2)
cont_traits <- str_split(trait_types[trait_types[, 1] == "continuous", 2], ",")[[1]]
bin_traits  <- str_split(trait_types[trait_types[, 1] == "binary", 2], ",")[[1]]

# Distances in trait space ---------------------------------------------------------------------------

# Gower distance on community-weighted means. All CWMs are numeric (for binary traits they are the
# weighted share of the community that has the trait), so vegdist(method = "gower") applies: a
# range-scaled Manhattan distance in which every trait weighs the same. Continuous traits are
# standardized first, as in the original analysis; this does not change Gower distances.
gower_dist <- function(cwm_tbl) {
  cwm_tbl %>%
    mutate(across(all_of(cont_traits), ~ as.vector(scale(.x)))) %>%
    vegdist(method = "gower", na.rm = TRUE)
}

# Alternative metrics for the sensitivity analysis (Table S12)
metric_set <- function(cwm_tbl) {
  m <- cwm_tbl %>% select(where(~ sd(.x, na.rm = TRUE) > 0))
  list(gower      = gower_dist(cwm_tbl),
       eucl_z     = dist(scale(m)),
       eucl_range = dist(apply(m, 2, \(x) (x - min(x)) / diff(range(x)))),
       manh_z     = dist(scale(m), method = "manhattan"))
}

# The statistic: PERMANOVA pseudo-F for a two-group split ---------------------------------------------

# D2: matrix of squared distances; group: logical or two-level vector, one value per sample
pseudo_F <- function(D2, group) {
  a   <- group == group[1]
  N   <- length(group)
  SST <- sum(D2) / (2 * N)
  SSW <- sum(D2[a, a]) / (2 * sum(a)) + sum(D2[!a, !a]) / (2 * sum(!a))
  c(F = (SST - SSW) / (SSW / (N - 2)), R2 = (SST - SSW) / SST)
}
sq <- function(D) as.matrix(D)^2

# Change-point search ----------------------------------------------------------------------------------

# Pseudo-F of every split "before = dates 1..k, after = dates k+1..n", with at least min_seg dates per side
split_F <- function(D2, d_idx, min_seg = 3) {
  ks <- min_seg:(max(d_idx) - min_seg)
  map_dbl(ks, ~ pseudo_F(D2, d_idx <= .x)[["F"]])
}
scan_splits <- function(D2, d_idx, min_seg = 3) {
  tibble(k = min_seg:(max(d_idx) - min_seg), F = split_F(D2, d_idx, min_seg))
}

# Break search in a subset of samples with dates `d`. The p-value comes from permuting whole dates and
# repeating the full search in every permutation, so it accounts for having chosen the best split.
# n_perm = 0 returns only the estimate (used by the bootstrap and the jackknife).
find_break <- function(D2, d, n_perm, min_seg = NULL) {
  u  <- sort(unique(d)); di <- match(d, u); n <- length(u)
  ms <- min_seg %||% max(2, floor(0.15 * n))
  best <- slice_max(scan_splits(D2, di, ms), F, n = 1, with_ties = FALSE)
  p <- if (n_perm > 0) {
    null_F <- map_dbl(seq_len(n_perm), ~ { q <- sample(n); max(split_F(D2, q[di], ms)) })
    (sum(null_F >= best$F) + 1) / (n_perm + 1)
  } else NA_real_
  tibble(last_pre = u[best$k], first_post = u[best$k + 1], F = best$F,
         R2 = pseudo_F(D2, di <= best$k)[["R2"]], p = p, n_dates = n)
}

# A priori detection window: a significant break counts as a detection if the change happened within
# `months` of the onset of the refilling (1 December 2017). The refilling is used only to evaluate
# the procedure; the procedure itself uses community data alone.
refilling_onset <- as.Date("2017-12-01")
in_window <- function(last_pre, first_post, months)
  first_post > refilling_onset & last_pre < refilling_onset %m+% months(months)

# Periods used in the descriptive analyses
phase_of <- function(d) case_when(d <= as.Date("2018-02-08") ~ "pre",
                                  d >= as.Date("2018-08-13") ~ "post",
                                  TRUE ~ "transition")

# Permutation test of a difference in means (post minus pre) among dates
perm_diff <- function(values, phase, nperm = 999) {
  obs  <- mean(values[phase == "post"]) - mean(values[phase == "pre"])
  null <- map_dbl(seq_len(nperm), ~ { p <- sample(phase); mean(values[p == "post"]) - mean(values[p == "pre"]) })
  c(diff = obs, p = (sum(abs(null) >= abs(obs)) + 1) / (nperm + 1))
}

# adonis2 R2 of a model; the formula is built here, next to D, so that every vegan version finds D
r2_model <- function(D, rhs, d) {
  f <- as.formula(paste("D ~", rhs))
  environment(f) <- environment()
  adonis2(f, data = d, by = NULL, permutations = 0)["Model", "R2"]
}

dir.create("output", showWarnings = FALSE)
