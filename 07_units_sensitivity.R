# Sensitivity of the data-requirement and break-date-uncertainty results to the resampling unit (Table S15).
# Run from the repository root:  Rscript 07_units_sensitivity.R   (about 30 minutes)
#
# In 02_robustness.R the units that are subsampled (thinning) and resampled (bootstrap) are the five
# long-term locations; the two depths of the reservoir locations (CLM1 and CLM2, DLM1 and DLM2) are kept
# together. Here the seven long-term sampling units are treated as separate units instead. Part 1
# measures how similar and how coupled the units are, which is the reason to consider this alternative:
# two depths of one location are about as coupled as the three inflows, which are treated as separate
# locations. Part 2 repeats the thinning by locations and the bootstrap with units.

source("R/functions.R")

samples     <- load_samples()
descriptors <- load_descriptors()
D2_list     <- map(descriptors[main_descriptors], sq)
lt_units    <- c("CLM1", "CLM2", "DLM1", "DLM2", "TLE", "TLR", "TSP")
reservoir   <- c("CLM", "DLM")

# 1. Similarity and coupling of the long-term sampling units ----------------------------------------------
pair_type <- function(a, b) {
  la <- substr(a, 1, 3); lb <- substr(b, 1, 3)
  case_when(la == lb & la %in% reservoir ~ "same location, other depth",
            la %in% reservoir & lb %in% reservoir & substr(a, 4, 4) == substr(b, 4, 4) ~ "other location, same depth",
            la %in% reservoir & lb %in% reservoir ~ "other location, other depth",
            !(la %in% reservoir) & !(lb %in% reservoir) ~ "between inflows",
            TRUE ~ "reservoir and inflow")
}
lt <- samples %>% mutate(i = row_number()) %>% filter(unit %in% lt_units) %>% mutate(unit = as.character(unit))
same_date_pairs <- lt %>% inner_join(lt, by = "date", suffix = c("_a", "_b"), relationship = "many-to-many") %>%
  filter(unit_a < unit_b) %>% mutate(type = pair_type(unit_a, unit_b))
next_date_pairs <- lt %>% arrange(unit, date) %>% group_by(unit) %>% mutate(j = lead(i)) %>% ungroup() %>%
  filter(!is.na(j))

unit_similarity <- imap_dfr(descriptors[c("traits_cells", "traits_biovolume", "taxa_genus")], function(D, dsc) {
  M <- as.matrix(D)
  # coupling: cosine between the deviations of two samples from the centroid of their date (PCoA space);
  # under independence its mean is about -1/(k - 1), with k units per date, because deviations sum to zero
  X   <- cmdscale(D, k = 20)
  dev <- X - apply(X, 2, function(col) ave(col, samples$date))
  cosine <- function(a, b) sum(dev[a, ] * dev[b, ]) / sqrt(sum(dev[a, ]^2) * sum(dev[b, ]^2))
  bind_rows(
    same_date_pairs %>% transmute(type, d = M[cbind(i_a, i_b)], cos = map2_dbl(i_a, i_b, cosine)),
    next_date_pairs %>% transmute(type = "same unit, next date", d = M[cbind(i, j)], cos = NA_real_)) %>%
    group_by(type) %>%
    summarise(n_pairs = n(), median_dissimilarity = round(median(d), 3), mean_coupling = round(mean(cos), 2),
              .groups = "drop") %>%
    mutate(descriptor = dsc)
})
k_per_date <- samples %>% count(date) %>% pull(n) %>% mean()
unit_similarity <- unit_similarity %>% mutate(coupling_if_independent = round(-1 / (k_per_date - 1), 2))
write_csv(unit_similarity, "output/23_unit_similarity.csv")
print(unit_similarity, n = Inf, width = Inf)

# 2. Thinning and bootstrap with sampling units as resampling units -----------------------------------------
dates <- sort(unique(samples$date))
set.seed(2038)
scen_units <- tribble(
  ~scenario,          ~frac_dates, ~n_units,
  "3 units",          1,           3,
  "1 unit",           1,           1,
  "1 unit, 1/2 dates", 1/2,        1)
thin_units <- function(frac_dates, n_units) {
  keep_dates <- sample(dates, round(frac_dates * length(dates)))
  keep_units <- sample(lt_units, n_units)
  i <- which(samples$date %in% keep_dates & samples$unit %in% keep_units)
  if (n_distinct(samples$date[i]) < 4) return(NULL)
  imap_dfr(D2_list, ~ find_break(.x[i, i], samples$date[i], 999) %>% mutate(descriptor = .y))
}
thin_reps <- scen_units %>% uncount(150, .id = "rep") %>%
  mutate(res = map2(frac_dates, n_units, thin_units)) %>% unnest(res)
detection_units <- thin_reps %>%
  mutate(sig = p < 0.05, W12 = sig & in_window(last_pre, first_post, 12)) %>%
  group_by(scenario, descriptor) %>%
  summarise(n_reps = n(), any_significant = mean(sig), W12 = mean(W12), .groups = "drop")
write_csv(detection_units, "output/23_detection_units.csv")

set.seed(2039)
lt_idx <- which(samples$unit %in% lt_units)
boot_units <- map_dfr(seq_len(999), function(b) {
  u <- sample(lt_units, replace = TRUE)
  i <- unlist(map(u, ~ which(samples$unit == .x)))
  imap_dfr(D2_list, ~ find_break(.x[i, i], samples$date[i], 0, 3) %>% mutate(descriptor = .y)) %>% mutate(boot = b)
})
boot_units_summary <- boot_units %>% count(descriptor, last_pre, first_post) %>%
  group_by(descriptor) %>% mutate(prop = n / sum(n)) %>% arrange(descriptor, desc(prop)) %>%
  mutate(cum = cumsum(prop)) %>%
  summarise(modal_last_pre = first(last_pre), modal_first_post = first(first_post), prop_modal = round(first(prop), 3),
            n_breaks_95 = which(cum >= 0.95)[1], .groups = "drop") %>%
  left_join(boot_units %>% group_by(descriptor) %>%
              summarise(prop_in_W12 = round(mean(in_window(last_pre, first_post, 12)), 3)), by = "descriptor")
write_csv(boot_units_summary, "output/23_bootstrap_units.csv")
print(detection_units, n = Inf); print(boot_units_summary, width = Inf)
saveRDS(sessionInfo(), "output/sessionInfo_07.rds")
