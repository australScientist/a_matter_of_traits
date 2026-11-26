# Shared data loading and functions for the reproducibility scripts.
# Source this file from the repository root: source("R/functions.R")
#
# A note on random numbers. Several functions below draw random numbers (permutations, resampling).
# The scripts call them in the same order as the original analysis, so that every p-value of the
# paper is reproduced exactly. If you reorder calls, results change within Monte Carlo error.

suppressPackageStartupMessages({
  library(tidyverse)
  library(vegan)
})

# Data ---------------------------------------------------------------------------------------------

# Sample metadata: one row per sample (sampling unit x date).
# `unit` and `location` use the names of the manuscript (Fig. S1). The factor levels keep the order
# of the original analysis, because the permutation tests draw random numbers in that order: with
# a different order the p-values would differ within Monte Carlo error.
load_samples <- function() {
  read_csv("data/samples.csv", show_col_types = FALSE) %>%
    mutate(date = as.Date(date),
           unit = factor(unit, levels = read_lines("data/unit_levels.txt")))
}

# The five long-term locations (sampled before the refilling), in the order of the original analysis
long_term_locations <- c("CLM", "TLE", "TLR", "TSP", "DLM")

# Community descriptors: one dissimilarity matrix between samples for each
descriptor_names <- c("traits_cells", "traits_biovolume", "traits_cells_dom", "traits_cells_klais",
                      "traits_cells_rimet", "taxa_genus", "taxa_genus_biovol", "taxa_class")
main_descriptors <- c("traits_cells", "traits_biovolume", "taxa_genus", "taxa_genus_biovol", "taxa_class")

# Read one distributed dissimilarity matrix (a square CSV with sample names) as a `dist` object
read_distance <- function(descriptor) {
  file.path("data/distances", paste0(descriptor, ".csv")) %>%
    read.csv(row.names = 1, check.names = FALSE) %>%
    as.matrix() %>%
    as.dist()
}

# All descriptors as a named list: descriptors$traits_cells, descriptors$taxa_genus, ...
load_descriptors <- function() {
  descriptor_names %>%
    set_names() %>%
    map(read_distance)
}

# Community-weighted means of traits (one row per sample, one column per trait)
load_cwm <- function(weighting) {
  file.path("data", paste0("cwm_", weighting, ".csv")) %>%
    read_csv(show_col_types = FALSE) %>%
    select(-sample)
}

# Which traits are continuous and which binary (the file has two lines, "continuous: a,b,c" and
# "binary: d,e,f")
trait_types <- read_lines("data/traits/trait_types.txt") %>%
  str_split_fixed(": ", 2)
traits_of_type <- function(type) {
  trait_types[trait_types[, 1] == type, 2] %>%
    str_split(",") %>%
    pluck(1)
}
cont_traits <- traits_of_type("continuous")
bin_traits  <- traits_of_type("binary")

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

# Alternative metrics for the sensitivity analysis (Table S12). Traits without variation are dropped,
# because they cannot be standardized.
metric_set <- function(cwm_tbl) {
  varying <- cwm_tbl %>%
    select(where(~ sd(.x, na.rm = TRUE) > 0))
  range_scaled <- varying %>%
    mutate(across(everything(), ~ (.x - min(.x)) / diff(range(.x))))
  list(gower      = gower_dist(cwm_tbl),
       eucl_z     = dist(scale(varying)),
       eucl_range = dist(as.matrix(range_scaled)),
       manh_z     = dist(scale(varying), method = "manhattan"))
}

# The statistic: PERMANOVA pseudo-F for a two-group split ---------------------------------------------
#
# These two functions run hundreds of thousands of times inside the permutation tests, so they work on
# matrices instead of tables.
#
# In PERMANOVA, a sum of squares can be computed from distances alone: the sum of squares of a group of
# samples is the sum of all squared distances among them divided by twice their number (Anderson 2001).
# `sum(D2)` counts each pair twice (i-j and j-i), hence the division by 2 * n.

# D2: matrix of squared distances; group: logical or two-level vector, one value per sample
pseudo_F <- function(D2, group) {
  in_first <- group == group[1]                  # TRUE for the samples of the first group
  n_total  <- length(group)

  ss_total   <- sum(D2) / (2 * n_total)
  ss_within  <- sum(D2[in_first, in_first])   / (2 * sum(in_first)) +
                sum(D2[!in_first, !in_first]) / (2 * sum(!in_first))
  ss_between <- ss_total - ss_within

  # two groups: 1 degree of freedom between groups, n - 2 within
  c(F  = ss_between / (ss_within / (n_total - 2)),
    R2 = ss_between / ss_total)
}

# Squared distances as a full matrix (what pseudo_F needs)
sq <- function(D) {
  as.matrix(D)^2
}

# Change-point search ----------------------------------------------------------------------------------
#
# `date_rank` is the rank of the date of each sample (1 = first date). The split "after date k" puts the
# samples of dates 1..k before the break and the rest after it. Every split leaves at least `min_seg`
# dates on each side.

# Pseudo-F of every candidate split
split_F <- function(D2, date_rank, min_seg = 3) {
  candidate_k <- min_seg:(max(date_rank) - min_seg)
  map_dbl(candidate_k, \(k) pseudo_F(D2, date_rank <= k)[["F"]])
}

# The same, as a table with the split (k) and its pseudo-F
scan_splits <- function(D2, date_rank, min_seg = 3) {
  tibble(k = min_seg:(max(date_rank) - min_seg),
         F = split_F(D2, date_rank, min_seg))
}

# Maximum pseudo-F after shuffling the dates: the samples of each date move together to a new position
# in time, and the full search is repeated. Its distribution under many shuffles is the null of the
# change-point test.
max_F_with_permuted_dates <- function(D2, date_rank, n_dates, min_seg = 3) {
  new_rank_of_each_date <- sample(n_dates)
  split_F(D2, new_rank_of_each_date[date_rank], min_seg) %>%
    max()
}

# Permutation p-value: share of permuted statistics at least as large as the observed one, adding the
# observed value to both counts so that p is never zero (Phipson and Smyth 2010)
perm_p_value <- function(null, observed) {
  (sum(null >= observed) + 1) / (length(null) + 1)
}

# Break search in a subset of samples with dates `dates`. The p-value comes from permuting whole dates
# and repeating the full search in every permutation, so it accounts for having chosen the best split
# (Andrews 1993). n_perm = 0 returns only the estimate (used by the bootstrap and the jackknife).
find_break <- function(D2, dates, n_perm, min_seg = NULL) {
  unique_dates <- sort(unique(dates))
  date_rank    <- match(dates, unique_dates)
  n_dates      <- length(unique_dates)
  min_seg      <- min_seg %||% max(2, floor(0.15 * n_dates))   # default: 15% of dates, at least 2

  best_split <- scan_splits(D2, date_rank, min_seg) %>%
    slice_max(F, n = 1, with_ties = FALSE)

  p_value <- NA_real_
  if (n_perm > 0) {
    null_max_F <- map_dbl(seq_len(n_perm),
                          \(i) max_F_with_permuted_dates(D2, date_rank, n_dates, min_seg))
    p_value <- perm_p_value(null_max_F, best_split$F)
  }

  tibble(last_pre   = unique_dates[best_split$k],        # last date before the break
         first_post = unique_dates[best_split$k + 1],    # first date after it
         F          = best_split$F,
         R2         = pseudo_F(D2, date_rank <= best_split$k)[["R2"]],
         p          = p_value,
         n_dates    = n_dates)
}

# Run find_break() on the samples `rows` for every descriptor of `D2_list` (a named list of squared
# distance matrices) and stack the results, one row per descriptor. All descriptors use the same
# samples, so their results can be compared directly.
search_descriptors <- function(D2_list, rows, n_perm, min_seg = NULL, sample_dates = samples$date) {
  find_in_one <- function(D2, descriptor) {
    find_break(D2[rows, rows], sample_dates[rows], n_perm, min_seg) %>%
      mutate(descriptor = descriptor)
  }
  D2_list %>%
    imap(find_in_one) %>%
    list_rbind()
}

# A priori detection window: a significant break counts as a detection if the change happened within
# `months` of the onset of the refilling (1 December 2017). The refilling is used only to evaluate
# the procedure; the procedure itself uses community data alone.
refilling_onset <- as.Date("2017-12-01")
in_window <- function(last_pre, first_post, months) {
  first_post > refilling_onset & last_pre < refilling_onset %m+% months(months)
}

# Periods used in the descriptive analyses
phase_of <- function(date) {
  case_when(date <= as.Date("2018-02-08") ~ "pre",
            date >= as.Date("2018-08-13") ~ "post",
            TRUE                          ~ "transition")
}

# Permutation test of a difference in means (after minus before) among dates: the period labels are
# shuffled among dates. Two-sided p-value.
perm_diff <- function(values, phase, nperm = 999) {
  mean_difference <- function(labels) {
    mean(values[labels == "post"]) - mean(values[labels == "pre"])
  }
  observed <- mean_difference(phase)
  null     <- map_dbl(seq_len(nperm), \(i) mean_difference(sample(phase)))
  c(diff = observed,
    p    = (sum(abs(null) >= abs(observed)) + 1) / (nperm + 1))
}

# adonis2 R2 of a model; the formula is built here, next to D, so that every vegan version finds D
r2_model <- function(D, rhs, d) {
  f <- as.formula(paste("D ~", rhs))
  environment(f) <- environment()
  adonis2(f, data = d, by = NULL, permutations = 0)["Model", "R2"]
}

dir.create("output", showWarnings = FALSE)
