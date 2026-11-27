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

# Tests that permute whole dates ------------------------------------------------------------------------
#
# The samples of one date share variation (the whole reservoir changes from one date to the next), so
# they cannot be permuted one by one: a permutation of samples, even within sampling units, treats
# samples of the same date as independent and gives p-values that are far too small (in records
# simulated without a change, about half of the steps came out significant; traitshift 0.2.0, NEWS).
# The tests below move all the samples of a date together, as the change-point search does.

# Run `code` with its own seed and then put the random number generator back where it was, so that a
# test added here does not change the random numbers, and hence the results, of later sections
with_own_seed <- function(seed, code) {
  saved <- get(".Random.seed", envir = globalenv())
  on.exit(assign(".Random.seed", saved, envir = globalenv()))
  set.seed(seed)
  code
}

# New position in time of each date, for one permutation (a number per date)
permuted_dates <- function(n_dates) sample(n_dates)

# Step versus trend. Distance-based linear model (McArdle and Anderson 2001) with sampling unit, a time
# trend and a step after `last_pre`; the step is tested by its marginal contribution (the same F and R2
# as adonis2(by = "margin")). In each permutation every date takes the time, and the side of the
# break, of another date, with all its samples. trend: "linear" or "spline" (natural cubic spline,
# 3 df). Samples must be in the order of the rows of D.
step_test <- function(D, samples, last_pre, trend = "linear", n_perm = 999) {
  stopifnot(identical(as.character(samples$sample), labels(D)))
  G <- gower_centred(D)
  ss_total <- sum(diag(G))

  dates      <- sort(unique(samples$date))
  date_index <- match(samples$date, dates)
  years      <- as.numeric(dates - min(dates)) / 365         # one value per date
  trend_of_dates <- if (trend == "linear") {
    matrix(years)
  } else {
    basis <- splines::ns(as.numeric(samples$date - min(dates)) / 365, df = 3)   # knots from the samples
    predict(basis, years)
  }
  step_of_dates <- as.numeric(dates > as.Date(last_pre))
  unit_part     <- model.matrix(~ unit, samples)

  # F and R2 of the step when date i takes the time of date order[i]
  step_F <- function(order) {
    assigned  <- order[date_index]
    X_reduced <- cbind(unit_part, trend_of_dates[assigned, , drop = FALSE])
    X_full    <- cbind(X_reduced, step_of_dates[assigned])
    reduced   <- explained_ss(G, X_reduced)
    full      <- explained_ss(G, X_full)
    c(F  = (full$ss - reduced$ss) / ((ss_total - full$ss) / (nrow(G) - full$rank)),
      R2 = (full$ss - reduced$ss) / ss_total)
  }

  observed <- step_F(seq_along(dates))
  null_F   <- map_dbl(seq_len(n_perm), \(i) step_F(permuted_dates(length(dates)))[["F"]])
  tibble(R2_step = observed[["R2"]],
         F_step  = observed[["F"]],
         p_step  = perm_p_value(null_F, observed[["F"]]))
}

# Gower-centred matrix of a distance matrix: G = -1/2 J D2 J (Gower 1966; McArdle and Anderson 2001).
# The sum of squares a design explains is the trace of G projected on the columns of the design.
gower_centred <- function(D) {
  n <- attr(D, "Size")
  J <- diag(n) - 1 / n
  -0.5 * J %*% sq(D) %*% J
}
explained_ss <- function(G, X) {
  q <- qr(X)
  Q <- qr.Q(q)[, seq_len(q$rank), drop = FALSE]
  list(ss = sum(Q * (G %*% Q)), rank = q$rank)
}

# Two periods (before and after a date): PERMANOVA pseudo-F with whole dates permuted
period_test <- function(D, samples, last_pre, n_perm = 9999) {
  D2 <- sq(D)
  dates      <- sort(unique(samples$date))
  date_index <- match(samples$date, dates)
  before     <- dates <= as.Date(last_pre)                    # one value per date
  observed   <- pseudo_F(D2, before[date_index])
  null_F     <- map_dbl(seq_len(n_perm), \(i) pseudo_F(D2, before[permuted_dates(length(dates))][date_index])[["F"]])
  perm_p_value(null_F, observed[["F"]])
}

# Homogeneity of dispersion between two periods (betadisper) with whole dates permuted: the distances of
# the samples to their period centroid are kept and the period of each date is shuffled
dispersion_test <- function(D, samples, last_pre, n_perm = 999) {
  dates      <- sort(unique(samples$date))
  date_index <- match(samples$date, dates)
  before     <- dates <= as.Date(last_pre)
  to_centroid <- betadisper(D, factor(before[date_index]))$distances
  F_of <- function(side) {
    fit <- anova(lm(to_centroid ~ factor(side)))
    fit[1, "F value"]
  }
  observed <- F_of(before[date_index])
  null_F   <- map_dbl(seq_len(n_perm), \(i) F_of(before[permuted_dates(length(dates))][date_index]))
  perm_p_value(null_F, observed)
}

# adonis2 R2 of a model; the formula is built here, next to D, so that every vegan version finds D
r2_model <- function(D, rhs, d) {
  f <- as.formula(paste("D ~", rhs))
  environment(f) <- environment()
  adonis2(f, data = d, by = NULL, permutations = 0)["Model", "R2"]
}

dir.create("output", showWarnings = FALSE)

# Temporal memory and a change-point test that allows for it -------------------------------------------
#
# The date permutation assumes that, without a change, dates are exchangeable. Seasonality and slow
# community dynamics make dates a few weeks apart alike, so permuted series are less structured than
# real ones and the p-value is too small. These functions (the same algorithms as traitshift 0.2.0,
# temporal_memory() and false_alarm_rate()) estimate the memory of the record and build a null made of
# records with that memory and no change.
#
# Everything works on dates: S[a, b] is the sum of squared distances between the samples of date a and
# those of date b, from which the sums of squares of every split, and the distances between date
# centroids, follow without going back to the samples.

date_blocks <- function(D2, date_index, n_dates) {
  S <- rowsum(t(rowsum(D2, date_index, reorder = TRUE)), date_index, reorder = TRUE)
  list(S = S, m = tabulate(date_index, n_dates))
}

# Pseudo-F of every split of the dates (at least min_seg dates per side); `position` is the place of each
# date in time. Same values as split_F() on the samples.
split_F_blocks <- function(blocks, position, min_seg = 3) {
  n <- length(blocks$m)
  o <- order(position)
  S <- blocks$S[o, o]
  m <- blocks$m[o]
  N <- sum(m)
  cum2 <- t(apply(apply(S, 2, cumsum), 1, cumsum))        # cum2[i, j] = sum of S[1:i, 1:j]
  k <- min_seg:(n - min_seg)
  n_first    <- cumsum(m)[k]
  ss_total   <- cum2[n, n] / (2 * N)
  ss_first   <- cum2[cbind(k, k)] / (2 * n_first)
  ss_rest    <- (cum2[n, n] - 2 * cum2[cbind(k, n)] + cum2[cbind(k, k)]) / (2 * (N - n_first))
  ss_between <- ss_total - ss_first - ss_rest
  ss_between / ((ss_first + ss_rest) / (N - 2))
}

# Memory of a record (traitshift::temporal_memory). Semivariogram of the squared distances between date
# centroids, after removing the scatter of samples around them, fitted by
# gamma(lag) = nugget + sill (1 - exp(-lag / tau)); only pairs of dates on the same side of `break_after`.
# r_lag1: correlation between consecutive dates; p_memory: permutation test of r_lag1 (dates shuffled
# within each side).
temporal_memory <- function(D, samples, break_after, n_perm = 999, max_dim = 20) {
  dates  <- sort(unique(samples$date))
  blocks <- date_blocks(sq(D), match(samples$date, dates), length(dates))
  m      <- blocks$m
  within_var <- sum(diag(blocks$S) / (2 * m)) / sum(m - 1)        # scatter of samples around their date

  self <- diag(blocks$S) / (2 * m^2)
  C <- blocks$S / outer(m, m) - outer(self, self, "+") - within_var * outer(1 / m, 1 / m, "+")

  side  <- as.integer(dates > as.Date(break_after))
  pairs <- which(upper.tri(C), arr.ind = TRUE)
  pairs <- pairs[side[pairs[, 1]] == side[pairs[, 2]], , drop = FALSE]
  lag   <- as.numeric(dates[pairs[, 2]] - dates[pairs[, 1]])
  gamma <- C[pairs] / 2

  # Exponential fit with the plateau fixed at the mean semivariance beyond the median lag (so a flat,
  # noisy variogram is not fitted by a huge sill with a huge tau); tau on a log grid
  total <- max(mean(gamma[lag > median(lag)]), .Machine$double.eps)
  fit_tau <- function(tau) {
    x    <- exp(-lag / tau)
    sill <- min(max(sum(x * (total - gamma)) / sum(x^2), 0), total)
    list(tau = tau, nugget = total - sill, sill = sill, rss = sum((gamma - total + sill * x)^2))
  }
  fits <- map(exp(seq(log(1), log(10 * as.numeric(diff(range(dates)))), length.out = 400)), fit_tau)
  fit  <- fits[[which.min(map_dbl(fits, "rss"))]]

  # Model-free test: are consecutive dates more alike than any two dates on the same side?
  G <- C / 2
  same_side <- outer(side, side, "==") & upper.tri(G)
  consecutive_mean <- function(order) {
    a <- order[-length(order)]
    b <- order[-1]
    keep <- side[a] == side[b]
    mean(G[cbind(a[keep], b[keep])])
  }
  observed <- consecutive_mean(seq_along(side))
  shuffle_within_sides <- function(i) {
    order <- seq_along(side)
    for (s in unique(side)) {
      at <- which(side == s)
      order[at] <- at[sample.int(length(at))]
    }
    consecutive_mean(order)
  }
  null <- map_dbl(seq_len(n_perm), shuffle_within_sides)

  # Share of the date-level variance on each axis (classical scaling of the date centroids)
  C_all  <- C + within_var * outer(1 / m, 1 / m, "+")
  J      <- diag(nrow(C_all)) - 1 / nrow(C_all)
  values <- eigen(-0.5 * J %*% C_all %*% J, symmetric = TRUE, only.values = TRUE)$values
  values <- head(values[values > 1e-10 * max(values)], max_dim)

  gap <- as.numeric(median(diff(dates)))
  list(r_lag1 = 1 - observed / mean(G[same_side]),
       p_memory = (sum(null <= observed) + 1) / (n_perm + 1),
       tau_days = fit$tau, nugget = fit$nugget, sill = fit$sill, within_var = within_var,
       phi = fit$sill / (fit$sill + fit$nugget) * exp(-gap / fit$tau), median_gap_days = gap,
       spectrum = values / sum(values))
}

# One record with no change and the memory `memory`, sampled like the real one (its dates and samples per
# date): an Ornstein-Uhlenbeck process at the date level, a nugget, and independent scatter of samples.
# Returns its date blocks.
simulate_memory_record <- function(memory, dates, n_samples) {
  t <- as.numeric(dates - dates[1])
  n <- length(t)
  w <- memory$spectrum
  k <- length(w)
  sd_signal <- sqrt(memory$sill * w)
  x <- matrix(0, n, k)
  x[1, ] <- rnorm(k, 0, sd_signal)
  for (i in seq_len(n)[-1]) {
    rho    <- exp(-(t[i] - t[i - 1]) / memory$tau_days)
    x[i, ] <- rho * x[i - 1, ] + rnorm(k, 0, sd_signal * sqrt(1 - rho^2))
  }
  x <- x + matrix(rnorm(n * k, 0, sqrt(memory$nugget * w)), n, k, byrow = TRUE)
  date_index <- rep(seq_len(n), n_samples)
  y <- x[date_index, , drop = FALSE] +
    matrix(rnorm(length(date_index) * k, 0, sqrt(memory$within_var / k)), ncol = k)
  date_blocks(as.matrix(dist(y))^2, date_index, n)
}

# Strength of the best split of a record, on a scale that does not depend on the descriptor: how many
# standard deviations the maximum pseudo-F lies above the maxima of the date-permuted series (z), with
# the date-permutation p-value and the last date before the best split (its rank k)
break_strength <- function(blocks, n_perm = 999, min_seg = 3) {
  n       <- length(blocks$m)
  profile <- split_F_blocks(blocks, seq_len(n), min_seg)
  null    <- map_dbl(seq_len(n_perm), \(i) max(split_F_blocks(blocks, permuted_dates(n), min_seg)))
  c(z = (max(profile) - mean(null)) / sd(null),
    p = perm_p_value(null, max(profile)),
    k = which.max(profile) + min_seg - 1)
}
