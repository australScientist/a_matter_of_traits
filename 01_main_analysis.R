# Main analysis: baseline, change-point search, step versus trend, metric sensitivity and spatial mosaic.
# Reproduces Table 1, Tables S4 and S12, and the data behind Figs. 3a, 4, S4 and S9a.
# Run from the repository root:  Rscript 01_main_analysis.R   (about 20 to 40 minutes)
#
# The random seed and the order of the calculations are those of the original analysis, so the
# p-values are reproduced exactly. Two steps of the original analysis need the phytoplankton counts
# and are not run here: the null model that shuffles traits among genera (its results are in
# data/precomputed/) and the leave-one-genus-out analysis. Because the null model drew random numbers
# before the change-point search, the random stream is advanced by the same number of draws (see below).

source("R/functions.R")
set.seed(2026)

samples     <- load_samples()
descriptors <- load_descriptors()
cwm_cells   <- load_cwm("cells")
cwm_biov    <- load_cwm("biovolume")

# Check: the Gower distances recomputed from the CWM tables equal the distributed matrices
stopifnot(max(abs(as.matrix(gower_dist(cwm_cells)) - as.matrix(descriptors$traits_cells))) < 1e-8,
          max(abs(as.matrix(gower_dist(cwm_biov))  - as.matrix(descriptors$traits_biovolume))) < 1e-8)

# Two periods split in mid-June 2018, used only for the baseline description
samples <- samples %>% mutate(period = if_else(date < as.Date("2018-06-15"), "pre", "post"))

# Check: own pseudo-F equals adonis2
own <- pseudo_F(sq(descriptors$traits_cells), samples$period)
ref <- adonis2(descriptors$traits_cells ~ period, data = samples, permutations = 0)
stopifnot(abs(own[["F"]] - ref$F[1]) < 1e-8, abs(own[["R2"]] - ref$R2[1]) < 1e-8)

# 1. Baseline: PERMANOVA, dispersion and PCoA ------------------------------------------------------------
baseline <- imap_dfr(descriptors, ~ {
  a   <- adonis2(.x ~ period, data = samples, permutations = 9999)
  bd  <- permutest(betadisper(.x, samples$period), permutations = 999)
  eig <- cmdscale(.x, k = 2, eig = TRUE)$eig %>% keep(~ .x > 0)
  tibble(descriptor = .y, R2 = round(a$R2[1], 3), F = round(a$F[1], 1), p = a$`Pr(>F)`[1],
         betadisper_p = bd$tab$`Pr(>F)`[1],
         PCoA1_pct = round(100 * eig[1] / sum(eig), 1),
         PCoA2_pct = round(100 * eig[2] / sum(eig), 1))
})
write_csv(baseline, "output/1_baseline.csv")

# 2. Null model (needs the counts; results in data/precomputed/2_null_model.csv) -------------------------
# The original analysis permuted the 38 genera 999 times for each weighting. The same number of draws
# is taken here so that the following tests use the same random numbers as the original analysis.
for (i in seq_len(2 * 999)) invisible(sample(38))

# 3. Change-point search with date permutation (Table 1, Table S4, Fig. 3a) ------------------------------
dates      <- sort(unique(samples$date))
date_index <- match(samples$date, dates)
n_perm_cp  <- 999

change_point <- imap_dfr(descriptors, ~ {
  D2   <- sq(.x)
  best <- slice_max(scan_splits(D2, date_index), F, n = 1)
  null_F <- map_dbl(seq_len(n_perm_cp), ~ {
    p <- sample(length(dates))                     # new order of the dates
    max(split_F(D2, p[date_index]))                # full search in the permuted series
  })
  tibble(descriptor = .y, last_pre = dates[best$k], first_post = dates[best$k + 1],
         maxF = round(best$F, 1), R2 = round(pseudo_F(D2, date_index <= best$k)[["R2"]], 3),
         p_date_perm = (sum(null_F >= best$F) + 1) / (n_perm_cp + 1))
})
write_csv(change_point, "output/3_change_point.csv")

F_profiles <- imap_dfr(descriptors, ~ scan_splits(sq(.x), date_index) %>%
                         mutate(descriptor = .y, last_pre = dates[k]))
write_csv(F_profiles, "output/3_F_profiles.csv")

# 4. Step versus gradual trend (Table 1, Fig. S4) --------------------------------------------------------
samples <- samples %>%
  mutate(years    = as.numeric(date - min(date)) / 365,
         step_jul = factor(date > as.Date("2018-07-01")),
         step_mar = factor(date > as.Date("2018-03-01")))
ctrl <- how(blocks = samples$unit, nperm = 999)        # permutations within sampling units

step_trend <- imap_dfr(descriptors, ~ {
  trend <- adonis2(.x ~ unit + years,            data = samples, by = "margin", permutations = ctrl)
  jul   <- adonis2(.x ~ unit + years + step_jul, data = samples, by = "margin", permutations = ctrl)
  mar   <- adonis2(.x ~ unit + years + step_mar, data = samples, by = "margin", permutations = ctrl)
  tibble(descriptor = .y, R2_trend = round(trend["years", "R2"], 3),
         R2_stepJul_beyond_trend = round(jul["step_jul", "R2"], 3), p_Jul = jul["step_jul", "Pr(>F)"],
         R2_stepMar_beyond_trend = round(mar["step_mar", "R2"], 3), p_Mar = mar["step_mar", "Pr(>F)"])
})
write_csv(step_trend, "output/4_step_vs_trend.csv")

# R2 with and without a step at the detected break, beyond a linear and a flexible trend (Table 1)
step_trend_r2 <- change_point %>%
  filter(descriptor %in% main_descriptors) %>%
  select(descriptor, last_pre) %>%
  pmap_dfr(function(descriptor, last_pre) {
    d <- samples %>% mutate(step = factor(date > last_pre))
    D <- descriptors[[descriptor]]
    tibble(descriptor, last_pre,
           R2_linear      = r2_model(D, "unit + years", d),
           R2_linear_step = r2_model(D, "unit + years + step", d),
           R2_spline      = r2_model(D, "unit + splines::ns(years, df = 3)", d),
           R2_spline_step = r2_model(D, "unit + splines::ns(years, df = 3) + step", d)) %>%
      mutate(dR2_linear = R2_linear_step - R2_linear, dR2_spline = R2_spline_step - R2_spline)
  }) %>%
  mutate(across(where(is.double), ~ round(.x, 3)))
write_csv(step_trend_r2, "output/19_step_trend_model_R2.csv")

# 5. Preliminary thinning (not reported) -------------------------------------------------------------------
# A first version of the data-requirement analysis, with a detection window derived from the data. It is
# superseded by the a priori window of 02_robustness.R but is kept because it draws random numbers that
# the spatial-mosaic tests below depend on.
detect_break <- function(D2, d) {
  u  <- sort(unique(d)); di <- match(d, u); n <- length(u)
  ms <- max(2, floor(0.15 * n))
  best   <- slice_max(scan_splits(D2, di, ms), F, n = 1, with_ties = FALSE)
  null_F <- map_dbl(1:99, ~ { p <- sample(n); max(split_F(D2, p[di], ms)) })
  p_val  <- (sum(null_F >= best$F) + 1) / 100
  p_val < 0.05 && u[best$k] >= as.Date("2018-02-08") && u[best$k + 1] <= as.Date("2018-08-13")
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
preliminary <- scenarios %>%
  mutate(reps = if_else(scenario == "all data", 1, 150)) %>%
  uncount(reps, .id = "rep") %>%
  mutate(hits = pmap(list(frac_dates, n_loc), function(frac_dates, n_loc) {
    keep_dates <- sample(dates, round(frac_dates * length(dates)))
    keep_loc <- if (is.na(n_loc)) unique(samples$location) else sample(long_term_locations, n_loc)
    i <- which(samples$date %in% keep_dates & samples$location %in% keep_loc)
    map_lgl(D2_list, ~ detect_break(.x[i, i], samples$date[i])) %>% enframe("descriptor", "hit")
  })) %>%
  unnest(hits) %>%
  group_by(scenario, descriptor) %>%
  summarise(detect_rate = round(mean(hit), 2), .groups = "drop")
write_csv(preliminary, "output/6_preliminary_thinning_not_reported.csv")

# 6. Sensitivity to the distance metric in trait space (Table S12) -----------------------------------------
metric_sens <- list(cells = cwm_cells, biovolume = cwm_biov) %>%
  imap_dfr(function(cw, nm) {
    ds <- metric_set(cw)
    imap_dfr(ds, function(d, metric) {
      best <- slice_max(scan_splits(sq(d), date_index), F, n = 1)
      tibble(weights = nm, metric = metric,
             R2_period = round(pseudo_F(sq(d), samples$period)[["R2"]], 3),
             last_pre = dates[best$k], first_post = dates[best$k + 1],
             rho_vs_gower = round(cor(as.vector(d), as.vector(ds$gower), method = "spearman"), 3))
    })
  })
write_csv(metric_sens, "output/18_metric_sensitivity.csv")

# 7. Spatial and temporal heterogeneity, and dominance (Fig. S9) --------------------------------------------
samples <- samples %>% mutate(phase = phase_of(date))
long_term_units <- samples %>% filter(location %in% long_term_locations) %>% distinct(unit) %>% pull(unit)
lt <- samples %>% filter(unit %in% long_term_units, phase != "transition")

mosaic <- imap_dfr(descriptors[c("traits_cells", "traits_biovolume", "taxa_genus")], ~ {
  M <- as.matrix(.x)
  pairs <- lt %>% select(sample, date, phase, location, unit) %>%
    inner_join(lt %>% select(sample2 = sample, date, location2 = location), by = "date",
               relationship = "many-to-many") %>%
    filter(sample < sample2) %>%
    mutate(d = M[cbind(sample, sample2)])
  # spatial heterogeneity per date: pairs of different locations, dates with at least 3 locations
  spatial <- pairs %>%
    filter(location != location2) %>%
    group_by(date, phase) %>%
    filter(n_distinct(c(location, location2)) >= 3) %>%
    summarise(het = mean(d), .groups = "drop")
  sp_test <- perm_diff(spatial$het, spatial$phase)
  # difference between the two depths of the same location and date
  depth_d <- pairs %>% filter(location == location2)
  # temporal turnover per sampling unit (consecutive samples in the same phase)
  turnover <- lt %>%
    arrange(unit, date) %>%
    group_by(unit) %>%
    mutate(prev = lag(sample), prev_phase = lag(phase)) %>%
    filter(!is.na(prev), prev_phase == phase) %>%
    mutate(d = map2_dbl(sample, prev, ~ M[.x, .y])) %>%
    ungroup()
  # variation partitioned between date and location within each phase
  part <- map_dfr(c("pre", "post"), function(ph) {
    s <- lt %>% filter(phase == ph)
    a <- adonis2(as.dist(M[s$sample, s$sample]) ~ factor(date) + location, data = s,
                 by = "margin", permutations = 999)
    tibble(phase = ph, R2_date = a["factor(date)", "R2"], R2_location = a["location", "R2"])
  })
  tibble(descriptor = .y,
         spatial_pre  = mean(spatial$het[spatial$phase == "pre"]),
         spatial_post = mean(spatial$het[spatial$phase == "post"]),
         spatial_p    = sp_test[["p"]],
         depth_pre    = mean(depth_d$d[depth_d$phase == "pre"]),
         depth_post   = mean(depth_d$d[depth_d$phase == "post"]),
         turnover_pre  = mean(turnover$d[turnover$phase == "pre"]),
         turnover_post = mean(turnover$d[turnover$phase == "post"]),
         R2_date_pre = part$R2_date[1], R2_location_pre = part$R2_location[1],
         R2_date_post = part$R2_date[2], R2_location_post = part$R2_location[2])
}) %>% mutate(across(where(is.numeric), ~ round(.x, 3)))

# Berger-Parker dominance (share of the most abundant genus), by cells and by biovolume
dom_dates <- lt %>%
  group_by(date, phase) %>%
  summarise(bp_cells = mean(berger_parker_cells), bp_biov = mean(berger_parker_biovolume), .groups = "drop")
dominance <- tibble(
  weighting = c("cells", "biovolume"),
  BP_pre  = c(mean(dom_dates$bp_cells[dom_dates$phase == "pre"]),  mean(dom_dates$bp_biov[dom_dates$phase == "pre"])),
  BP_post = c(mean(dom_dates$bp_cells[dom_dates$phase == "post"]), mean(dom_dates$bp_biov[dom_dates$phase == "post"])),
  p = c(perm_diff(dom_dates$bp_cells, dom_dates$phase)[["p"]], perm_diff(dom_dates$bp_biov, dom_dates$phase)[["p"]])) %>%
  mutate(across(where(is.numeric), ~ round(.x, 3)))

write_csv(mosaic, "output/8_mosaic_heterogeneity.csv")
write_csv(dominance, "output/8_mosaic_dominance.csv")

saveRDS(sessionInfo(), "output/sessionInfo_01.rds")
message("Done: results in output/")
