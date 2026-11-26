# Robustness: data requirements with an a priori detection window, false alarms in stretches without a
# shift, uncertainty of the break date, and the step test beyond a flexible trend.
# Reproduces Fig. 3d, Tables S7 to S10, and the data behind Figs. S3, S4 and S5.
# Run from the repository root:  Rscript 02_robustness.R   (a few hours)
#
# Each section sets its own seed, as in the original analysis, so the results are reproduced exactly.
# Two analyses of the original document need the phytoplankton counts and are not run here: the
# known-truth tests with planted shifts (Fig. 5, results in data/precomputed/) and the sensitivity to
# cell volumes (Fig. S8).

source("R/functions.R")

samples     <- load_samples()
descriptors <- load_descriptors()
D2_list     <- map(descriptors[main_descriptors], sq)   # squared distances of the five main descriptors
dates       <- sort(unique(samples$date))
all_rows    <- seq_len(nrow(samples))

n_perm_cp  <- 999    # permutations for tests on the full data
n_perm_sub <- 999    # permutations for each thinned replicate
n_reps_sub <- 150    # replicates per thinning scenario
n_boot     <- 999    # bootstrap resamples of locations

# Check: the full-data breaks equal those of the main analysis (01_main_analysis.R)
full_breaks <- search_descriptors(D2_list, all_rows, n_perm = 0, min_seg = 3)
print(full_breaks %>% select(descriptor, last_pre, first_post))

# 1. Data requirements with the a priori window (Fig. 3d, Fig. S4, Table S7) -------------------------------
# Thin the record (fewer dates, fewer locations, or both) 150 times per scenario, repeat the full test
# in each replicate, and count how often a significant break falls in a window set before looking at
# the community data: within 6, 12 or 18 months of the onset of the refilling.
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

# One replicate: keep a random subset of dates (from `pool_dates`) and of locations, and search for a
# break in every descriptor. Locations are dropped whole: both depths of a location go together.
thin_once <- function(frac_dates, n_loc, pool_dates, n_perm) {
  kept_dates     <- sample(pool_dates, round(frac_dates * length(pool_dates)))
  kept_locations <- if (is.na(n_loc)) unique(samples$location) else sample(long_term_locations, n_loc)
  rows <- which(samples$date %in% kept_dates & samples$location %in% kept_locations)
  if (n_distinct(samples$date[rows]) < 4) {
    return(NULL)                                       # too few dates to search
  }
  search_descriptors(D2_list, rows, n_perm)
}

# All replicates of all scenarios (the full data are run once)
run_thinning <- function(pool_dates) {
  scenarios %>%
    mutate(reps = if_else(scenario == "all data", 1, n_reps_sub)) %>%
    uncount(reps, .id = "rep") %>%
    mutate(res = map2(frac_dates, n_loc,
                      \(f, n) thin_once(f, n, pool_dates = pool_dates, n_perm = n_perm_sub))) %>%
    unnest(res)
}

# Share of replicates of each scenario and descriptor, in the original scenario order
rate_by_scenario <- function(replicates, ...) {
  replicates %>%
    group_by(scenario, descriptor) %>%
    summarise(n_reps = n(), ..., .groups = "drop") %>%
    mutate(scenario = factor(scenario, levels = scenarios$scenario)) %>%
    arrange(scenario, descriptor)
}

thin_reps <- run_thinning(dates)
detection <- thin_reps %>%
  mutate(sig = p < 0.05,
         # window derived from the full data (used in earlier versions; circular)
         data_window = sig & last_pre >= as.Date("2018-02-08") & first_post <= as.Date("2018-08-13"),
         # windows defined a priori from the refilling
         W6  = sig & in_window(last_pre, first_post, 6),
         W12 = sig & in_window(last_pre, first_post, 12),
         W18 = sig & in_window(last_pre, first_post, 18)) %>%
  rate_by_scenario(any_significant = mean(sig), old_window = mean(data_window),
                   W6 = mean(W6), W12 = mean(W12), W18 = mean(W18))
write_csv(detection, "output/11_detection_apriori_window.csv")
saveRDS(thin_reps, "output/11_thinning_replicates.rds")

# 2. False alarms in stretches without a regime shift (Table S8, Fig. S3b) ---------------------------------
# Does the procedure stay quiet where no regime shift is expected? Two stretches: after the shift
# (from 13 Aug 2018) and before the refilling (before 1 Dec 2017, short, so 2 dates per side minimum).
set.seed(2028)
stable <- list(
  post = list(dates = dates[dates >= as.Date("2018-08-13")], min_seg = 3),
  pre  = list(dates = dates[dates <  refilling_onset],        min_seg = 2))

main_R2 <- read_csv("output/3_change_point.csv", show_col_types = FALSE) %>%
  select(descriptor, R2_main_break = R2)

search_stretch <- function(stretch, stretch_name) {
  rows <- which(samples$date %in% stretch$dates)
  search_descriptors(D2_list, rows, n_perm_cp, stretch$min_seg) %>%
    mutate(stretch = stretch_name, min_seg = stretch$min_seg)
}
stable_full <- stable %>%
  imap(search_stretch) %>%
  list_rbind() %>%
  left_join(main_R2, by = "descriptor") %>%
  select(stretch, descriptor, n_dates, min_seg, last_pre, first_post, R2, R2_main_break, p) %>%
  mutate(across(c(R2, R2_main_break), ~ round(.x, 3)))
write_csv(stable_full, "output/12_false_alarm_stable_stretches.csv")

# Thinned replicates of the post-shift stretch: any significant break is a false alarm
set.seed(2029)
false_alarm <- run_thinning(stable$post$dates) %>%
  rate_by_scenario(false_alarm_rate = mean(p < 0.05))
write_csv(false_alarm, "output/12_false_alarm_rates.csv")

# What the procedure found in the post-shift stretch: long-term locations together and one by one
set.seed(2031)
post_dates <- stable$post$dates
rows_post_lt <- which(samples$date %in% post_dates & samples$location %in% long_term_locations)
post_lt <- search_descriptors(D2_list, rows_post_lt, n_perm_cp, 3) %>%
  select(descriptor, last_pre, first_post, R2, p) %>%
  mutate(R2 = round(R2, 3))

post_in_location <- function(loc) {
  rows <- which(samples$date %in% post_dates & samples$location == loc)
  search_descriptors(D2_list[c("traits_biovolume", "taxa_genus_biovol")], rows, n_perm_cp, 3) %>%
    mutate(location = loc)
}
post_by_location <- long_term_locations %>%
  map(post_in_location) %>%
  list_rbind() %>%
  select(location, descriptor, last_pre, first_post, R2, p) %>%
  mutate(R2 = round(R2, 3))
write_csv(post_lt, "output/12_post_stretch_long_term.csv")
write_csv(post_by_location, "output/12_post_stretch_by_location.csv")

# Post hoc: the stretch after the late-2018 succession (from 10 Dec 2018)
set.seed(2032)
late_dates <- dates[dates >= as.Date("2018-12-10")]
rows_late  <- which(samples$date %in% late_dates)
late_full <- search_descriptors(D2_list, rows_late, n_perm_cp, 3) %>%
  select(descriptor, n_dates, last_pre, first_post, R2, p) %>%
  mutate(R2 = round(R2, 3))
late_alarm <- run_thinning(late_dates) %>%
  rate_by_scenario(false_alarm_rate = mean(p < 0.05))
write_csv(late_full, "output/12_false_alarm_posthoc_full.csv")
write_csv(late_alarm, "output/12_false_alarm_posthoc_rates.csv")

# 3. Uncertainty of the break date (Table S10, Fig. S3c) ----------------------------------------------------
set.seed(2030)
rows_lt <- which(samples$location %in% long_term_locations)
reference_lt <- search_descriptors(D2_list, rows_lt, n_perm = 0, min_seg = 3) %>%
  select(descriptor, ref_last_pre = last_pre, ref_first_post = first_post)

# Bootstrap: resample the five long-term locations with replacement (a location drawn twice enters
# twice) and search for the break again
one_bootstrap <- function(b) {
  drawn <- sample(long_term_locations, replace = TRUE)
  rows  <- unlist(map(drawn, \(loc) which(samples$location == loc)))
  search_descriptors(D2_list, rows, n_perm = 0, min_seg = 3) %>%
    mutate(boot = b, n_distinct_locations = n_distinct(drawn))
}
boot <- map(seq_len(n_boot), one_bootstrap) %>%
  list_rbind()
saveRDS(boot, "output/13_bootstrap_replicates.rds")

# How often each break was chosen, from the most to the least frequent
boot_freq <- boot %>%
  count(descriptor, last_pre, first_post) %>%
  group_by(descriptor) %>%
  mutate(prop = n / sum(n)) %>%
  arrange(descriptor, desc(prop)) %>%
  mutate(cum = cumsum(prop)) %>%
  ungroup()

# Most frequent break, number of breaks needed to cover 95% of resamples and the dates they span
share_in_window <- boot %>%
  group_by(descriptor) %>%
  summarise(prop_in_W12 = round(mean(in_window(last_pre, first_post, 12)), 3))
boot_summary <- boot_freq %>%
  group_by(descriptor) %>%
  summarise(modal_last_pre   = first(last_pre),
            modal_first_post = first(first_post),
            prop_modal       = round(first(prop), 3),
            n_breaks_95      = which(cum >= 0.95)[1],
            set95_from       = min(last_pre[seq_len(which(cum >= 0.95)[1])]),
            set95_to         = max(first_post[seq_len(which(cum >= 0.95)[1])]),
            .groups = "drop") %>%
  left_join(share_in_window, by = "descriptor") %>%
  left_join(reference_lt, by = "descriptor")
write_csv(boot_summary, "output/13_break_date_bootstrap.csv")
write_csv(boot_freq, "output/13_break_date_bootstrap_freq.csv")

# Jackknife: drop one date at a time and check whether the break still separates the same dates.
# When the dropped date is the last one before the break, the "same" break is the previous date.
full_last <- full_breaks %>%
  select(descriptor, full_last_pre = last_pre)

without_date <- function(dropped_date) {
  rows <- which(samples$date != dropped_date)
  search_descriptors(D2_list, rows, n_perm = 0, min_seg = 3) %>%
    mutate(dropped = dropped_date)
}
expected_break <- function(dropped, full_last_pre) {
  max(dates[dates != dropped & dates <= full_last_pre])
}
jack <- map(dates, without_date) %>%
  list_rbind() %>%
  left_join(full_last, by = "descriptor") %>%
  mutate(expected_last_pre = as.Date(map2_dbl(dropped, full_last_pre, expected_break)),
         same_break        = last_pre == expected_last_pre)

jack_summary <- jack %>%
  group_by(descriptor) %>%
  summarise(n_dropped          = n(),
            n_same             = sum(same_break),
            alternatives       = paste(unique(format(last_pre[!same_break], "%Y-%m-%d")), collapse = "; "),
            dropped_when_moved = paste(format(dropped[!same_break], "%Y-%m-%d"), collapse = "; "),
            .groups = "drop")
write_csv(jack_summary, "output/13_break_date_jackknife.csv")
write_csv(jack, "output/13_break_date_jackknife_all.csv")

# 4. Step at the detected break beyond a linear and a flexible trend (Fig. S5) -----------------------------
# Model: unit + trend + step, with the step tested by its marginal contribution and permutations within
# sampling units. The flexible trend is a natural cubic spline with 3 degrees of freedom.
set.seed(2034)
samples_t <- samples %>%
  mutate(years = as.numeric(date - min(date)) / 365)
breaks_main <- read_csv("output/3_change_point.csv", show_col_types = FALSE) %>%
  filter(descriptor %in% main_descriptors) %>%
  select(descriptor, last_pre)

step_beyond_trends <- function(descriptor, last_pre) {
  with_step <- samples_t %>%
    mutate(step = factor(date > last_pre))
  D <- descriptors[[descriptor]]
  within_units <- how(blocks = with_step$unit, nperm = 999)

  linear <- adonis2(D ~ unit + years + step, data = with_step, by = "margin",
                    permutations = within_units)
  spline <- adonis2(D ~ unit + splines::ns(years, df = 3) + step, data = with_step, by = "margin",
                    permutations = within_units)

  tibble(descriptor, last_pre,
         R2_step_linear = round(linear["step", "R2"], 3), p_linear = linear["step", "Pr(>F)"],
         R2_step_spline = round(spline["step", "R2"], 3), p_spline = spline["step", "Pr(>F)"])
}
step_main <- breaks_main %>%
  pmap(step_beyond_trends) %>%
  list_rbind()
write_csv(step_main, "output/14_step_beyond_spline_main.csv")

saveRDS(sessionInfo(), "output/sessionInfo_02.rds")
message("Done: results in output/")
