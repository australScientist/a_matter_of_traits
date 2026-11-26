# Robustness: data requirements with an a priori detection window, false alarms in stretches without a
# shift, uncertainty of the break date, and the step test beyond a flexible trend.
# Reproduces Fig. 3c, Tables S7 to S10, and the data behind Figs. S2, S3 and S4.
# Run from the repository root:  Rscript 02_robustness.R   (a few hours)
#
# Each section sets its own seed, as in the original analysis, so the results are reproduced exactly.
# Two analyses of the original document need the phytoplankton counts and are not run here: the
# known-truth tests with planted shifts (Fig. 5, results in data/precomputed/) and the sensitivity to
# cell volumes (Fig. S7).

source("R/functions.R")

samples     <- load_samples()
descriptors <- load_descriptors()
D2_list     <- map(descriptors[main_descriptors], sq)
dates       <- sort(unique(samples$date))

n_perm_cp  <- 999    # permutations for tests on the full data
n_perm_sub <- 999    # permutations for each thinned replicate
n_reps_sub <- 150    # replicates per thinning scenario
n_boot     <- 999    # bootstrap resamples of locations

# Check: the full-data breaks equal those of the main analysis (01_main_analysis.R)
full_breaks <- imap_dfr(D2_list, ~ find_break(.x, samples$date, 0, 3) %>% mutate(descriptor = .y))
print(full_breaks %>% select(descriptor, last_pre, first_post))

# 1. Data requirements with the a priori window (Fig. 3c, Fig. S3, Table S7) -------------------------------
set.seed(2027)
scenarios <- tribble(
  ~scenario,               ~frac_dates, ~n_loc,
  "all data",              1,           NA,
  "2/3 dates",             2/3,         NA,
  "1/2 dates",             1/2,         NA,
  "1/3 dates",             1/3,         NA,
  "3 locations",           1,           3,
  "1 location",            1,           1,
  "1 location, 1/2 dates", 1/2,         1)

# One replicate: keep a random subset of dates and locations from `pool_dates` and search for a break
thin_once <- function(frac_dates, n_loc, pool_dates, n_perm) {
  keep_dates <- sample(pool_dates, round(frac_dates * length(pool_dates)))
  keep_loc   <- if (is.na(n_loc)) unique(samples$location) else sample(long_term_locations, n_loc)
  i <- which(samples$date %in% keep_dates & samples$location %in% keep_loc)
  if (n_distinct(samples$date[i]) < 4) return(NULL)          # too few dates to search
  imap_dfr(D2_list, ~ find_break(.x[i, i], samples$date[i], n_perm) %>% mutate(descriptor = .y))
}
run_thinning <- function(pool_dates) {
  scenarios %>%
    mutate(reps = if_else(scenario == "all data", 1, n_reps_sub)) %>%
    uncount(reps, .id = "rep") %>%
    mutate(res = pmap(list(frac_dates, n_loc), thin_once, pool_dates = pool_dates, n_perm = n_perm_sub)) %>%
    unnest(res)
}

thin_reps <- run_thinning(dates)
detection <- thin_reps %>%
  mutate(sig = p < 0.05,
         data_window = sig & last_pre >= as.Date("2018-02-08") & first_post <= as.Date("2018-08-13"),
         W6  = sig & in_window(last_pre, first_post, 6),
         W12 = sig & in_window(last_pre, first_post, 12),
         W18 = sig & in_window(last_pre, first_post, 18)) %>%
  group_by(scenario, descriptor) %>%
  summarise(n_reps = n(), any_significant = mean(sig), old_window = mean(data_window),
            W6 = mean(W6), W12 = mean(W12), W18 = mean(W18), .groups = "drop") %>%
  mutate(scenario = factor(scenario, levels = scenarios$scenario)) %>%
  arrange(scenario, descriptor)
write_csv(detection, "output/11_detection_apriori_window.csv")
saveRDS(thin_reps, "output/11_thinning_replicates.rds")

# 2. False alarms in stretches without a regime shift (Table S8, Fig. S2b) ---------------------------------
set.seed(2028)
stable <- list(
  post = list(dates = dates[dates >= as.Date("2018-08-13")], min_seg = 3),
  pre  = list(dates = dates[dates <  refilling_onset],        min_seg = 2))

main_R2 <- read_csv("output/3_change_point.csv", show_col_types = FALSE) %>%
  select(descriptor, R2_main_break = R2)
stable_full <- imap_dfr(stable, function(st, nm) {
  i <- which(samples$date %in% st$dates)
  imap_dfr(D2_list, ~ find_break(.x[i, i], samples$date[i], n_perm_cp, st$min_seg) %>%
             mutate(descriptor = .y)) %>%
    mutate(stretch = nm, min_seg = st$min_seg)
}) %>%
  left_join(main_R2, by = "descriptor") %>%
  select(stretch, descriptor, n_dates, min_seg, last_pre, first_post, R2, R2_main_break, p) %>%
  mutate(across(c(R2, R2_main_break), ~ round(.x, 3)))
write_csv(stable_full, "output/12_false_alarm_stable_stretches.csv")

set.seed(2029)
false_alarm <- run_thinning(stable$post$dates) %>%
  group_by(scenario, descriptor) %>%
  summarise(n_reps = n(), false_alarm_rate = mean(p < 0.05), .groups = "drop") %>%
  mutate(scenario = factor(scenario, levels = scenarios$scenario)) %>%
  arrange(scenario, descriptor)
write_csv(false_alarm, "output/12_false_alarm_rates.csv")

# What the procedure found in the post-shift stretch: long-term locations together and one by one
set.seed(2031)
post_dates <- stable$post$dates
i_lt <- which(samples$date %in% post_dates & samples$location %in% long_term_locations)
post_lt <- imap_dfr(D2_list, ~ find_break(.x[i_lt, i_lt], samples$date[i_lt], n_perm_cp, 3) %>%
                      mutate(descriptor = .y)) %>%
  select(descriptor, last_pre, first_post, R2, p) %>% mutate(R2 = round(R2, 3))
post_by_location <- map_dfr(long_term_locations, function(l) {
  i <- which(samples$date %in% post_dates & samples$location == l)
  imap_dfr(D2_list[c("traits_biovolume", "taxa_genus_biovol")],
           ~ find_break(.x[i, i], samples$date[i], n_perm_cp, 3) %>% mutate(descriptor = .y, location = l))
}) %>% select(location, descriptor, last_pre, first_post, R2, p) %>% mutate(R2 = round(R2, 3))
write_csv(post_lt, "output/12_post_stretch_long_term.csv")
write_csv(post_by_location, "output/12_post_stretch_by_location.csv")

# Post hoc: the stretch after the late-2018 succession
set.seed(2032)
late_dates <- dates[dates >= as.Date("2018-12-10")]
i_late <- which(samples$date %in% late_dates)
late_full <- imap_dfr(D2_list, ~ find_break(.x[i_late, i_late], samples$date[i_late], n_perm_cp, 3) %>%
                        mutate(descriptor = .y)) %>%
  select(descriptor, n_dates, last_pre, first_post, R2, p) %>% mutate(R2 = round(R2, 3))
late_alarm <- run_thinning(late_dates) %>%
  group_by(scenario, descriptor) %>%
  summarise(n_reps = n(), false_alarm_rate = mean(p < 0.05), .groups = "drop") %>%
  mutate(scenario = factor(scenario, levels = scenarios$scenario)) %>%
  arrange(scenario, descriptor)
write_csv(late_full, "output/12_false_alarm_posthoc_full.csv")
write_csv(late_alarm, "output/12_false_alarm_posthoc_rates.csv")

# 3. Uncertainty of the break date (Table S10, Fig. S2c) ----------------------------------------------------
set.seed(2030)
lt <- which(samples$location %in% long_term_locations)
reference_lt <- imap_dfr(D2_list, ~ find_break(.x[lt, lt], samples$date[lt], 0, 3) %>%
                           mutate(descriptor = .y)) %>%
  select(descriptor, ref_last_pre = last_pre, ref_first_post = first_post)

# Bootstrap: resample the five long-term locations with replacement
boot <- map_dfr(seq_len(n_boot), function(b) {
  locs <- sample(long_term_locations, replace = TRUE)
  i    <- unlist(map(locs, ~ which(samples$location == .x)))
  imap_dfr(D2_list, ~ find_break(.x[i, i], samples$date[i], 0, 3) %>% mutate(descriptor = .y)) %>%
    mutate(boot = b, n_distinct_locations = n_distinct(locs))
})
saveRDS(boot, "output/13_bootstrap_replicates.rds")

boot_freq <- boot %>% count(descriptor, last_pre, first_post) %>%
  group_by(descriptor) %>% mutate(prop = n / sum(n)) %>% arrange(descriptor, desc(prop)) %>%
  mutate(cum = cumsum(prop)) %>% ungroup()
boot_summary <- boot_freq %>%
  group_by(descriptor) %>%
  summarise(modal_last_pre = first(last_pre), modal_first_post = first(first_post),
            prop_modal = round(first(prop), 3),
            n_breaks_95 = which(cum >= 0.95)[1],
            set95_from = min(last_pre[seq_len(which(cum >= 0.95)[1])]),
            set95_to   = max(first_post[seq_len(which(cum >= 0.95)[1])]),
            .groups = "drop") %>%
  left_join(boot %>% group_by(descriptor) %>%
              summarise(prop_in_W12 = round(mean(in_window(last_pre, first_post, 12)), 3)),
            by = "descriptor") %>%
  left_join(reference_lt, by = "descriptor")
write_csv(boot_summary, "output/13_break_date_bootstrap.csv")
write_csv(boot_freq, "output/13_break_date_bootstrap_freq.csv")

# Jackknife: drop one date at a time
full_last <- full_breaks %>% select(descriptor, full_last_pre = last_pre)
jack <- map_dfr(dates, function(dd) {
  i <- which(samples$date != dd)
  imap_dfr(D2_list, ~ find_break(.x[i, i], samples$date[i], 0, 3) %>% mutate(descriptor = .y)) %>%
    mutate(dropped = dd)
}) %>%
  left_join(full_last, by = "descriptor") %>%
  mutate(expected_last_pre = map2_dbl(dropped, full_last_pre,
                                      ~ max(dates[dates != .x & dates <= .y])) %>% as.Date(),
         same_break = last_pre == expected_last_pre)
jack_summary <- jack %>%
  group_by(descriptor) %>%
  summarise(n_dropped = n(), n_same = sum(same_break),
            alternatives = paste(unique(format(last_pre[!same_break], "%Y-%m-%d")), collapse = "; "),
            dropped_when_moved = paste(format(dropped[!same_break], "%Y-%m-%d"), collapse = "; "),
            .groups = "drop")
write_csv(jack_summary, "output/13_break_date_jackknife.csv")
write_csv(jack, "output/13_break_date_jackknife_all.csv")

# 4. Step at the detected break beyond a linear and a flexible trend (Fig. S4) -----------------------------
set.seed(2034)
samples_t <- samples %>% mutate(years = as.numeric(date - min(date)) / 365)
breaks_main <- read_csv("output/3_change_point.csv", show_col_types = FALSE) %>%
  filter(descriptor %in% main_descriptors) %>% select(descriptor, last_pre)
step_main <- breaks_main %>% pmap_dfr(function(descriptor, last_pre) {
  d <- samples_t %>% mutate(step = factor(date > last_pre))
  D <- descriptors[[descriptor]]
  ctrl <- how(blocks = d$unit, nperm = 999)
  lin <- adonis2(D ~ unit + years + step, data = d, by = "margin", permutations = ctrl)
  spl <- adonis2(D ~ unit + splines::ns(years, df = 3) + step, data = d, by = "margin", permutations = ctrl)
  tibble(descriptor, last_pre,
         R2_step_linear = round(lin["step", "R2"], 3), p_linear = lin["step", "Pr(>F)"],
         R2_step_spline = round(spl["step", "R2"], 3), p_spline = spl["step", "Pr(>F)"])
})
write_csv(step_main, "output/14_step_beyond_spline_main.csv")

saveRDS(sessionInfo(), "output/sessionInfo_02.rds")
message("Done: results in output/")
