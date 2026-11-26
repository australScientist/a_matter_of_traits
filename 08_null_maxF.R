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
samples <- samples %>% mutate(period = if_else(date < as.Date("2018-06-15"), "pre", "post"))

# random numbers used by section 1 of 01_main_analysis.R (baseline PERMANOVA and dispersion)
walk(descriptors, ~ {
  adonis2(.x ~ period, data = samples, permutations = 9999)
  permutest(betadisper(.x, samples$period), permutations = 999)
})
# random numbers used by the trait null model (section 2)
for (i in seq_len(2 * 999)) invisible(sample(38))

dates      <- sort(unique(samples$date))
date_index <- match(samples$date, dates)
n_perm_cp  <- 999

null_maxF <- imap_dfr(descriptors, ~ {
  D2   <- sq(.x)
  best <- slice_max(scan_splits(D2, date_index), F, n = 1)
  null_F <- map_dbl(seq_len(n_perm_cp), ~ {
    p <- sample(length(dates))
    max(split_F(D2, p[date_index]))
  })
  tibble(descriptor = .y, perm = seq_len(n_perm_cp), null_F = null_F, observed_F = best$F)
})

check <- null_maxF %>% group_by(descriptor) %>%
  summarise(p = (sum(null_F >= first(observed_F)) + 1) / (n() + 1)) %>%
  left_join(read_csv("output/3_change_point.csv", show_col_types = FALSE), by = "descriptor")
print(check %>% select(descriptor, p, p_date_perm))
stopifnot(isTRUE(all.equal(check$p, check$p_date_perm)))

write_csv(null_maxF, "output/24_null_maxF.csv")
