# Permutation null distribution of the maximum pseudo-F (Fig. 3b).
# Run from the repository root:  Rscript 08_null_maxF.R   (about 10 to 20 minutes)
#
# Repeats sections 1 to 3 of 01_main_analysis.R with the same seed and the same order of calculations,
# so that the permuted maxima are exactly those behind the p-values of Table 1, and stores them.
# The script stops if the p-values differ from output/3_change_point.csv.

source("R/functions.R")
set.seed(2026)

samples     <- load_samples()
descriptors <- load_descriptors()
samples <- samples %>%
  mutate(period = if_else(date < as.Date("2018-06-15"), "pre", "post"))

# Advance the random stream exactly as sections 1 and 2 of 01_main_analysis.R do -------------------------
# section 1: baseline PERMANOVA and dispersion test of every descriptor
for (D in descriptors) {
  adonis2(D ~ period, data = samples, permutations = 9999)
  permutest(betadisper(D, samples$period), permutations = 999)
}
# section 2: the trait null model (2 weightings x 999 permutations of 38 genera)
for (i in seq_len(2 * 999)) invisible(sample(38))

# Section 3 of 01_main_analysis.R, keeping every permuted maximum ------------------------------------------
dates      <- sort(unique(samples$date))
date_index <- match(samples$date, dates)
n_perm_cp  <- 999

null_of <- function(D, descriptor) {
  D2 <- sq(D)
  best_split <- scan_splits(D2, date_index) %>%
    slice_max(F, n = 1)
  null_max_F <- map_dbl(seq_len(n_perm_cp),
                        \(i) max_F_with_permuted_dates(D2, date_index, length(dates)))
  tibble(descriptor = descriptor,
         perm       = seq_len(n_perm_cp),
         null_F     = null_max_F,
         observed_F = best_split$F)
}
null_maxF <- descriptors %>%
  imap(null_of) %>%
  list_rbind()

# Check: these permutations give exactly the p-values of Table 1
check <- null_maxF %>%
  group_by(descriptor) %>%
  summarise(p = perm_p_value(null_F, first(observed_F))) %>%
  left_join(read_csv("output/3_change_point.csv", show_col_types = FALSE), by = "descriptor")
print(check %>% select(descriptor, p, p_date_perm))
stopifnot(isTRUE(all.equal(check$p, check$p_date_perm)))

write_csv(null_maxF, "output/24_null_maxF.csv")
