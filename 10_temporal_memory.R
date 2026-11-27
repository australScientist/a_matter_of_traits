# Temporal memory of the record and a change-point test that allows for it (Methods, "Assumptions of the
# test"; Table S16). Run from the repository root:  Rscript 10_temporal_memory.R   (about 1 to 2 hours)
#
# The date permutation of the change-point search assumes that, without a change, sampling dates would
# be exchangeable. Seasonality and slow community dynamics make dates a few weeks apart alike, and a
# record with such memory but no change often gives a "significant" break under that permutation.
#
# 1. Memory. For each descriptor, the semivariogram of the distances between date centroids is fitted
#    with an exponential model, using only pairs of dates within the two stable stretches (before the
#    refilling, up to 8 Feb 2018, and after the transition, from 13 Aug 2018), so that neither the shift
#    nor the transition between the two states is counted as memory. r_lag1 is the correlation between
#    consecutive dates and p_memory its permutation test.
# 2. Calibrated test. 1000 records with that memory and no change are simulated with the real sampling
#    design (27 dates, the real number of samples per date) and searched like the real record. The
#    strength of the best split is measured as z, the number of standard deviations by which its maximum
#    pseudo-F exceeds the maxima of the date-permuted series; z does not have the floor of a permutation
#    p-value. p_calibrated is the share of simulated records with a z at least as large as the observed
#    one; p_calibrated_window also requires the simulated break to fall within 12 months of the onset of
#    the refilling, the a priori detection window of the paper. false_alarm_05 is the share of simulated
#    records that the date permutation calls significant (p < 0.05).
#
# The simulated records are Gaussian and their memory is a single exponential decay, a simplification
# of seasonality. The calibration is therefore an approximation, not an exact test.

source("R/functions.R")

samples     <- load_samples()
descriptors <- load_descriptors()
dates       <- sort(unique(samples$date))
date_index  <- match(samples$date, dates)
n_samples   <- tabulate(date_index, length(dates))

refilling_window_end <- refilling_onset %m+% months(12)
in_detection_window  <- function(k) dates[k + 1] > refilling_onset & dates[k] < refilling_window_end

stable_rows <- which(samples$date <= as.Date("2018-02-08") | samples$date >= as.Date("2018-08-13"))
stable_samples <- samples[stable_rows, ]
stable_distance <- function(D) as.dist(as.matrix(D)[stable_rows, stable_rows])

# 1. Memory ---------------------------------------------------------------------------------------------------
set.seed(2040)
memory <- descriptors[main_descriptors] %>%
  map(\(D) temporal_memory(stable_distance(D), stable_samples, break_after = "2018-02-08"))

memory_table <- memory %>%
  imap(\(m, descriptor) tibble(descriptor,
                               r_lag1          = round(m$r_lag1, 2),
                               p_memory        = m$p_memory,
                               tau_days        = round(m$tau_days),
                               median_gap_days = m$median_gap_days,
                               phi             = round(m$phi, 2))) %>%
  list_rbind()
write_csv(memory_table, "output/25_temporal_memory.csv")
print(memory_table)

# 2. Calibrated change-point test ---------------------------------------------------------------------------
n_sim <- as.integer(Sys.getenv("N_SIM", "1000"))   # N_SIM=10 for a quick trial
set.seed(2041)
calibrate <- function(D, descriptor) {
  observed  <- break_strength(date_blocks(sq(D), date_index, length(dates)))
  simulated <- map(seq_len(n_sim), \(i) break_strength(simulate_memory_record(memory[[descriptor]], dates, n_samples))) %>%
    do.call(rbind, .)
  stronger  <- simulated[, "z"] >= observed[["z"]]
  in_window <- in_detection_window(simulated[, "k"])
  tibble(descriptor,
         z_observed          = round(observed[["z"]], 2),
         p_date_permutation  = observed[["p"]],
         observed_in_window  = in_detection_window(observed[["k"]]),
         false_alarm_05      = mean(simulated[, "p"] < 0.05),
         p_calibrated        = (sum(stronger) + 1) / (n_sim + 1),
         p_calibrated_window = (sum(stronger & in_window) + 1) / (n_sim + 1))
}
calibrated <- descriptors[main_descriptors] %>%
  imap(calibrate) %>%
  list_rbind()
write_csv(calibrated, "output/25_calibrated_change_point.csv")
print(calibrated, width = Inf)
saveRDS(sessionInfo(), "output/sessionInfo_10.rds")
