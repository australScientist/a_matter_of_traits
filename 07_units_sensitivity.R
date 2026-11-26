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
# Type of a pair of sampling units. Unit names are location (3 letters) plus depth (1 or 2) for the
# reservoir locations, so CLM1 and CLM2 are the two depths of CLM.
pair_type <- function(a, b) {
  location_a <- substr(a, 1, 3)
  location_b <- substr(b, 1, 3)
  depth_a    <- substr(a, 4, 4)
  depth_b    <- substr(b, 4, 4)
  both_in_reservoir <- location_a %in% reservoir & location_b %in% reservoir
  both_inflows      <- !(location_a %in% reservoir) & !(location_b %in% reservoir)

  case_when(location_a == location_b & location_a %in% reservoir ~ "same location, other depth",
            both_in_reservoir & depth_a == depth_b              ~ "other location, same depth",
            both_in_reservoir                                    ~ "other location, other depth",
            both_inflows                                         ~ "between inflows",
            TRUE                                                 ~ "reservoir and inflow")
}
# Long-term samples, with their row number (i) in the distance matrices
lt <- samples %>%
  mutate(i = row_number()) %>%
  filter(unit %in% lt_units) %>%
  mutate(unit = as.character(unit))

# Every pair of samples of the same date (each pair once), labelled by type of pair
same_date_pairs <- lt %>%
  inner_join(lt, by = "date", suffix = c("_a", "_b"), relationship = "many-to-many") %>%
  filter(unit_a < unit_b) %>%
  mutate(type = pair_type(unit_a, unit_b))

# Every pair of consecutive samples of the same unit (j is the row of the next sample)
next_date_pairs <- lt %>%
  arrange(unit, date) %>%
  group_by(unit) %>%
  mutate(j = lead(i)) %>%
  ungroup() %>%
  filter(!is.na(j))

# Coupling: two samples of the same date can resemble each other because both reflect the date, or
# because they also deviate together from it. To measure the second, take the deviation of each sample
# from the centroid of its date (in a 20-axis PCoA space) and compute the cosine between the deviations
# of the two samples. Without coupling its mean is about -1/(k - 1), with k samples per date, because
# the deviations of a date sum to zero.
similarity_of <- function(D, descriptor) {
  M <- as.matrix(D)
  coords     <- cmdscale(D, k = 20)
  date_means <- apply(coords, 2, \(axis) ave(axis, samples$date))   # centroid of each sample's date
  deviation  <- coords - date_means
  cosine <- function(a, b) {
    sum(deviation[a, ] * deviation[b, ]) / sqrt(sum(deviation[a, ]^2) * sum(deviation[b, ]^2))
  }

  same_date <- same_date_pairs %>%
    transmute(type, d = M[cbind(i_a, i_b)], cos = map2_dbl(i_a, i_b, cosine))
  next_date <- next_date_pairs %>%
    transmute(type = "same unit, next date", d = M[cbind(i, j)], cos = NA_real_)

  bind_rows(same_date, next_date) %>%
    group_by(type) %>%
    summarise(n_pairs              = n(),
              median_dissimilarity = round(median(d), 3),
              mean_coupling        = round(mean(cos), 2),
              .groups = "drop") %>%
    mutate(descriptor = descriptor)
}
unit_similarity <- descriptors[c("traits_cells", "traits_biovolume", "taxa_genus")] %>%
  imap(similarity_of) %>%
  list_rbind()

k_per_date <- samples %>%
  count(date) %>%
  pull(n) %>%
  mean()
unit_similarity <- unit_similarity %>%
  mutate(coupling_if_independent = round(-1 / (k_per_date - 1), 2))
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
# One replicate: keep a random subset of dates and of sampling units (not whole locations)
thin_units <- function(frac_dates, n_units) {
  kept_dates <- sample(dates, round(frac_dates * length(dates)))
  kept_units <- sample(lt_units, n_units)
  rows <- which(samples$date %in% kept_dates & samples$unit %in% kept_units)
  if (n_distinct(samples$date[rows]) < 4) {
    return(NULL)
  }
  search_descriptors(D2_list, rows, n_perm = 999)
}
thin_reps <- scen_units %>%
  uncount(150, .id = "rep") %>%
  mutate(res = map2(frac_dates, n_units, thin_units)) %>%
  unnest(res)

# Detection: a significant break within 12 months of the onset of the refilling
detection_units <- thin_reps %>%
  mutate(sig = p < 0.05, W12 = sig & in_window(last_pre, first_post, 12)) %>%
  group_by(scenario, descriptor) %>%
  summarise(n_reps          = n(),
            any_significant = mean(sig),
            W12             = mean(W12),
            .groups = "drop")
write_csv(detection_units, "output/23_detection_units.csv")

set.seed(2039)
# Bootstrap: resample the seven sampling units with replacement and search for the break again
one_bootstrap <- function(b) {
  drawn <- sample(lt_units, replace = TRUE)
  rows  <- unlist(map(drawn, \(u) which(samples$unit == u)))
  search_descriptors(D2_list, rows, n_perm = 0, min_seg = 3) %>%
    mutate(boot = b)
}
boot_units <- map(seq_len(999), one_bootstrap) %>%
  list_rbind()

# Most frequent break and how many breaks are needed to cover 95% of resamples
boot_units_summary <- boot_units %>%
  count(descriptor, last_pre, first_post) %>%
  group_by(descriptor) %>%
  mutate(prop = n / sum(n)) %>%
  arrange(descriptor, desc(prop)) %>%
  mutate(cum = cumsum(prop)) %>%
  summarise(modal_last_pre   = first(last_pre),
            modal_first_post = first(first_post),
            prop_modal       = round(first(prop), 3),
            n_breaks_95      = which(cum >= 0.95)[1],
            .groups = "drop") %>%
  left_join(boot_units %>%
              group_by(descriptor) %>%
              summarise(prop_in_W12 = round(mean(in_window(last_pre, first_post, 12)), 3)),
            by = "descriptor")
write_csv(boot_units_summary, "output/23_bootstrap_units.csv")
print(detection_units, n = Inf)
print(boot_units_summary, width = Inf)
saveRDS(sessionInfo(), "output/sessionInfo_07.rds")
