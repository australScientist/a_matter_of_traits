# Main analysis: baseline, change-point search, step versus trend, metric sensitivity and spatial mosaic.
# Reproduces Table 1, Tables S4 and S12, and the data behind Figs. 3a, 4, S5 and S10a.
# Run from the repository root:  Rscript 01_main_analysis.R   (about 20 to 40 minutes)
#
# The random seed and the order of the calculations are those of the original analysis, so the
# p-values are reproduced exactly. Two steps of the original analysis need the phytoplankton counts
# and are not run here: the null model that shuffles traits among genera (its results are in
# data/precomputed/) and the leave-one-genus-out analysis. Because the null model drew random numbers
# before the change-point search, the random stream is advanced by the same number of draws (see below).
#
# Revision of 27 Nov 2026: the p-values of the baseline PERMANOVA, of the dispersion test and of the step
# versus trend now come from permutations of whole dates (step_test(), period_test() and
# dispersion_test() in R/functions.R). The adonis2() and permutest() calls are kept, because they give the
# R2 and F and because their permutations are part of the random stream on which the later sections
# depend; their p-values, which permute single samples, are no longer reported. The new tests run with
# their own seeds (with_own_seed()), so every other result of this script is unchanged.

source("R/functions.R")
set.seed(2026)

samples     <- load_samples()
descriptors <- load_descriptors()
cwm_cells   <- load_cwm("cells")
cwm_biov    <- load_cwm("biovolume")

# Check: the Gower distances recomputed from the CWM tables equal the distributed matrices
max_abs_difference <- function(d1, d2) max(abs(as.matrix(d1) - as.matrix(d2)))
stopifnot(max_abs_difference(gower_dist(cwm_cells), descriptors$traits_cells)     < 1e-8,
          max_abs_difference(gower_dist(cwm_biov),  descriptors$traits_biovolume) < 1e-8)

# Two periods split in mid-June 2018, used only for the baseline description
samples <- samples %>%
  mutate(period = if_else(date < as.Date("2018-06-15"), "pre", "post"))

# Check: our own pseudo-F (used in the permutation loops) equals the one of vegan::adonis2
own_F     <- pseudo_F(sq(descriptors$traits_cells), samples$period)
adonis_F  <- adonis2(descriptors$traits_cells ~ period, data = samples, permutations = 0)
stopifnot(abs(own_F[["F"]]  - adonis_F$F[1])  < 1e-8,
          abs(own_F[["R2"]] - adonis_F$R2[1]) < 1e-8)

# 1. Baseline: PERMANOVA, dispersion and PCoA ------------------------------------------------------------
# For each descriptor: how much of the variation the two periods explain (PERMANOVA), whether the
# spread of samples differs between periods (betadisper), and how much variation the first two
# principal coordinates carry. p and betadisper_p permute whole dates (see the note at the top).
baseline_of <- function(D, descriptor) {
  permanova  <- adonis2(D ~ period, data = samples, permutations = 9999)
  dispersion <- permutest(betadisper(D, samples$period), permutations = 999)
  eigenvalues <- cmdscale(D, k = 2, eig = TRUE)$eig %>%
    keep(\(x) x > 0)                                   # negative eigenvalues carry no variance
  share_of <- function(axis) round(100 * eigenvalues[axis] / sum(eigenvalues), 1)

  tibble(descriptor   = descriptor,
         R2           = round(permanova$R2[1], 3),
         F            = round(permanova$F[1], 1),
         p            = with_own_seed(3001, period_test(D, samples, "2018-06-14", n_perm = 9999)),
         betadisper_p = with_own_seed(3002, dispersion_test(D, samples, "2018-06-14", n_perm = 999)),
         PCoA1_pct    = share_of(1),
         PCoA2_pct    = share_of(2))
}
baseline <- descriptors %>%
  imap(baseline_of) %>%
  list_rbind()
write_csv(baseline, "output/1_baseline.csv")

# 2. Null model (needs the counts; results in data/precomputed/2_null_model.csv) -------------------------
# The original analysis permuted the 38 genera 999 times for each weighting. The same number of draws
# is taken here so that the following tests use the same random numbers as the original analysis.
for (i in seq_len(2 * 999)) invisible(sample(38))

# 3. Change-point search with date permutation (Table 1, Table S4, Fig. 3a) ------------------------------
# For each descriptor:
#   1. try every split between two consecutive dates and keep the one with the largest pseudo-F;
#   2. shuffle the dates 999 times and, in each shuffle, repeat the full search. Repeating the search
#      matters: some split always has a high F by chance, so the null must include "choosing the best";
#   3. the p-value is the share of shuffles whose maximum is at least as large as the observed one.
dates      <- sort(unique(samples$date))
date_index <- match(samples$date, dates)               # rank of the date of each sample
n_perm_cp  <- 999

test_change_point <- function(D, descriptor) {
  D2 <- sq(D)

  # Step 1: observed break
  best_split <- scan_splits(D2, date_index) %>%
    slice_max(F, n = 1)

  # Step 2: null distribution of the maximum pseudo-F
  null_max_F <- map_dbl(seq_len(n_perm_cp),
                        \(i) max_F_with_permuted_dates(D2, date_index, length(dates)))

  # Step 3: p-value and effect size at the break
  R2_at_break <- pseudo_F(D2, date_index <= best_split$k)[["R2"]]

  tibble(descriptor  = descriptor,
         last_pre    = dates[best_split$k],             # last date before the break
         first_post  = dates[best_split$k + 1],         # first date after it
         maxF        = round(best_split$F, 1),
         R2          = round(R2_at_break, 3),
         p_date_perm = perm_p_value(null_max_F, best_split$F))
}
change_point <- descriptors %>%
  imap(test_change_point) %>%
  list_rbind()
write_csv(change_point, "output/3_change_point.csv")

# Pseudo-F of every split, for the profiles of Fig. 3a
profile_of <- function(D, descriptor) {
  scan_splits(sq(D), date_index) %>%
    mutate(descriptor = descriptor,
           last_pre   = dates[k])
}
F_profiles <- descriptors %>%
  imap(profile_of) %>%
  list_rbind()
write_csv(F_profiles, "output/3_F_profiles.csv")

# 4. Step versus gradual trend (Table 1, Fig. S5) --------------------------------------------------------
# Does a step explain variation beyond a linear time trend? Model: unit + years (+ step), with the
# step tested by its marginal contribution. Steps on two fixed dates (July and March 2018). R2 come from
# adonis2; p from step_test(), which permutes whole dates (see the note at the top).
samples <- samples %>%
  mutate(years    = as.numeric(date - min(date)) / 365,       # time since the first sample
         step_jul = factor(date > as.Date("2018-07-01")),
         step_mar = factor(date > as.Date("2018-03-01")))
within_units <- how(blocks = samples$unit, nperm = 999)

step_vs_trend_of <- function(D, descriptor) {
  trend_only <- adonis2(D ~ unit + years,            data = samples, by = "margin", permutations = within_units)
  step_july  <- adonis2(D ~ unit + years + step_jul, data = samples, by = "margin", permutations = within_units)
  step_march <- adonis2(D ~ unit + years + step_mar, data = samples, by = "margin", permutations = within_units)

  tibble(descriptor              = descriptor,
         R2_trend                = round(trend_only["years", "R2"], 3),
         R2_stepJul_beyond_trend = round(step_july["step_jul", "R2"], 3),
         p_Jul                   = with_own_seed(3003, step_test(D, samples, "2018-07-01", "linear"))$p_step,
         R2_stepMar_beyond_trend = round(step_march["step_mar", "R2"], 3),
         p_Mar                   = with_own_seed(3004, step_test(D, samples, "2018-03-01", "linear"))$p_step)
}
step_trend <- descriptors %>%
  imap(step_vs_trend_of) %>%
  list_rbind()
write_csv(step_trend, "output/4_step_vs_trend.csv")

# R2 with and without a step at the detected break, beyond a linear and a flexible trend (Table 1).
# The flexible trend is a natural cubic spline with 3 degrees of freedom.
R2_with_and_without_step <- function(descriptor, last_pre) {
  with_step <- samples %>%
    mutate(step = factor(date > last_pre))
  D <- descriptors[[descriptor]]

  tibble(descriptor,
         last_pre,
         R2_linear      = r2_model(D, "unit + years", with_step),
         R2_linear_step = r2_model(D, "unit + years + step", with_step),
         R2_spline      = r2_model(D, "unit + splines::ns(years, df = 3)", with_step),
         R2_spline_step = r2_model(D, "unit + splines::ns(years, df = 3) + step", with_step)) %>%
    mutate(dR2_linear = R2_linear_step - R2_linear,      # what the step adds to each trend
           dR2_spline = R2_spline_step - R2_spline)
}
step_trend_r2 <- change_point %>%
  filter(descriptor %in% main_descriptors) %>%
  select(descriptor, last_pre) %>%
  pmap(R2_with_and_without_step) %>%
  list_rbind() %>%
  mutate(across(where(is.double), ~ round(.x, 3)))
write_csv(step_trend_r2, "output/19_step_trend_model_R2.csv")

# 5. Preliminary thinning (not reported) -------------------------------------------------------------------
# A first version of the data-requirement analysis, with a detection window derived from the data. It is
# superseded by the a priori window of 02_robustness.R but is kept because it draws random numbers that
# the spatial-mosaic tests below depend on.

# TRUE if a break is significant (99 date permutations) and falls between 8 Feb and 13 Aug 2018
detect_break <- function(D2, sample_dates) {
  unique_dates <- sort(unique(sample_dates))
  date_rank    <- match(sample_dates, unique_dates)
  n_dates      <- length(unique_dates)
  min_seg      <- max(2, floor(0.15 * n_dates))

  best_split <- scan_splits(D2, date_rank, min_seg) %>%
    slice_max(F, n = 1, with_ties = FALSE)
  null_max_F <- map_dbl(1:99, \(i) max_F_with_permuted_dates(D2, date_rank, n_dates, min_seg))
  p_value    <- (sum(null_max_F >= best_split$F) + 1) / 100

  p_value < 0.05 &&
    unique_dates[best_split$k]     >= as.Date("2018-02-08") &&
    unique_dates[best_split$k + 1] <= as.Date("2018-08-13")
}

scenarios <- tribble(
  ~scenario,               ~frac_dates, ~n_loc,
  "all data",              1,           NA,
  "2/3 dates",             2/3,         NA,
  "1/2 dates",             1/2,         NA,
  "1/3 dates",             1/3,         NA,
  "3 locations",           1,           3,
  "1 location",            1,           1,
  "1 location, 1/2 dates", 1/2,         1)
D2_list <- map(descriptors[main_descriptors], sq)

# One thinned replicate: keep a random subset of dates and of locations, then test every descriptor
preliminary_replicate <- function(frac_dates, n_loc) {
  kept_dates     <- sample(dates, round(frac_dates * length(dates)))
  kept_locations <- if (is.na(n_loc)) unique(samples$location) else sample(long_term_locations, n_loc)
  rows <- which(samples$date %in% kept_dates & samples$location %in% kept_locations)
  D2_list %>%
    map_lgl(\(D2) detect_break(D2[rows, rows], samples$date[rows])) %>%
    enframe("descriptor", "hit")
}
preliminary <- scenarios %>%
  mutate(reps = if_else(scenario == "all data", 1, 150)) %>%   # the full data are run once
  uncount(reps, .id = "rep") %>%
  mutate(hits = map2(frac_dates, n_loc, preliminary_replicate)) %>%
  unnest(hits) %>%
  group_by(scenario, descriptor) %>%
  summarise(detect_rate = round(mean(hit), 2), .groups = "drop")
write_csv(preliminary, "output/6_preliminary_thinning_not_reported.csv")

# 6. Sensitivity to the distance metric in trait space (Table S12) -----------------------------------------
# For each weighting and each metric: R2 between periods, break dates, and rank correlation of the
# distances with the Gower distances of the main analysis.
metric_results <- function(cwm_tbl, weighting) {
  distances <- metric_set(cwm_tbl)

  one_metric <- function(d, metric) {
    best_split <- scan_splits(sq(d), date_index) %>%
      slice_max(F, n = 1)
    tibble(weights      = weighting,
           metric       = metric,
           R2_period    = round(pseudo_F(sq(d), samples$period)[["R2"]], 3),
           last_pre     = dates[best_split$k],
           first_post   = dates[best_split$k + 1],
           rho_vs_gower = round(cor(as.vector(d), as.vector(distances$gower), method = "spearman"), 3))
  }
  distances %>%
    imap(one_metric) %>%
    list_rbind()
}
metric_sens <- list(cells = cwm_cells, biovolume = cwm_biov) %>%
  imap(metric_results) %>%
  list_rbind()
write_csv(metric_sens, "output/18_metric_sensitivity.csv")

# 7. Spatial and temporal heterogeneity, and dominance (Fig. S10) --------------------------------------------
# Was the community more patchy in space, or more variable in time, before or after the shift?
# Only the long-term sampling units, and only samples outside the transition.
samples <- samples %>%
  mutate(phase = phase_of(date))
long_term_units <- samples %>%
  filter(location %in% long_term_locations) %>%
  distinct(unit) %>%
  pull(unit)
lt <- samples %>%
  filter(unit %in% long_term_units, phase != "transition")

heterogeneity_of <- function(D, descriptor) {
  M <- as.matrix(D)

  # every pair of samples taken on the same date, with their dissimilarity
  pairs <- lt %>%
    select(sample, date, phase, location, unit) %>%
    inner_join(lt %>% select(sample2 = sample, date, location2 = location),
               by = "date", relationship = "many-to-many") %>%
    filter(sample < sample2) %>%                          # each pair once
    mutate(d = M[cbind(sample, sample2)])

  # spatial heterogeneity per date: pairs of different locations, dates with at least 3 locations
  spatial <- pairs %>%
    filter(location != location2) %>%
    group_by(date, phase) %>%
    filter(n_distinct(c(location, location2)) >= 3) %>%
    summarise(het = mean(d), .groups = "drop")
  spatial_test <- perm_diff(spatial$het, spatial$phase)

  # difference between the two depths of the same location and date
  depth_pairs <- pairs %>%
    filter(location == location2)

  # temporal turnover per sampling unit: dissimilarity between consecutive samples of the same phase
  turnover <- lt %>%
    arrange(unit, date) %>%
    group_by(unit) %>%
    mutate(prev = lag(sample), prev_phase = lag(phase)) %>%
    filter(!is.na(prev), prev_phase == phase) %>%
    mutate(d = map2_dbl(sample, prev, \(a, b) M[a, b])) %>%
    ungroup()

  # variation partitioned between date and location within each phase
  partition_in_phase <- function(ph) {
    phase_samples <- lt %>% filter(phase == ph)
    fit <- adonis2(as.dist(M[phase_samples$sample, phase_samples$sample]) ~ factor(date) + location,
                   data = phase_samples, by = "margin", permutations = 999)
    tibble(phase = ph, R2_date = fit["factor(date)", "R2"], R2_location = fit["location", "R2"])
  }
  partition <- map(c("pre", "post"), partition_in_phase) %>%
    list_rbind()

  mean_in <- function(x, phases, which_phase) mean(x[phases == which_phase])
  tibble(descriptor       = descriptor,
         spatial_pre      = mean_in(spatial$het, spatial$phase, "pre"),
         spatial_post     = mean_in(spatial$het, spatial$phase, "post"),
         spatial_p        = spatial_test[["p"]],
         depth_pre        = mean_in(depth_pairs$d, depth_pairs$phase, "pre"),
         depth_post       = mean_in(depth_pairs$d, depth_pairs$phase, "post"),
         turnover_pre     = mean_in(turnover$d, turnover$phase, "pre"),
         turnover_post    = mean_in(turnover$d, turnover$phase, "post"),
         R2_date_pre      = partition$R2_date[1],
         R2_location_pre  = partition$R2_location[1],
         R2_date_post     = partition$R2_date[2],
         R2_location_post = partition$R2_location[2])
}
mosaic <- descriptors[c("traits_cells", "traits_biovolume", "taxa_genus")] %>%
  imap(heterogeneity_of) %>%
  list_rbind() %>%
  mutate(across(where(is.numeric), ~ round(.x, 3)))

# Berger-Parker dominance (share of the most abundant genus), by cells and by biovolume, averaged by date
dom_dates <- lt %>%
  group_by(date, phase) %>%
  summarise(bp_cells = mean(berger_parker_cells),
            bp_biov  = mean(berger_parker_biovolume),
            .groups = "drop")
mean_by_phase <- function(x, which_phase) mean(x[dom_dates$phase == which_phase])
p_cells <- perm_diff(dom_dates$bp_cells, dom_dates$phase)[["p"]]   # cells first, as in the original
p_biov  <- perm_diff(dom_dates$bp_biov,  dom_dates$phase)[["p"]]
dominance <- tibble(
  weighting = c("cells", "biovolume"),
  BP_pre    = c(mean_by_phase(dom_dates$bp_cells, "pre"),  mean_by_phase(dom_dates$bp_biov, "pre")),
  BP_post   = c(mean_by_phase(dom_dates$bp_cells, "post"), mean_by_phase(dom_dates$bp_biov, "post")),
  p         = c(p_cells, p_biov)) %>%
  mutate(across(where(is.numeric), ~ round(.x, 3)))

write_csv(mosaic, "output/8_mosaic_heterogeneity.csv")
write_csv(dominance, "output/8_mosaic_dominance.csv")

saveRDS(sessionInfo(), "output/sessionInfo_01.rds")
message("Done: results in output/")
