# Two checks of the assumptions behind the change-point test (Methods, "Change-point search"; Table S14).
# Run from the repository root:  Rscript 06_dependence_dispersion.R   (about 10 minutes)
#
# 1. Temporal autocorrelation. The date permutation assumes that, without a change, sampling dates would
#    be exchangeable. Autocorrelation violates that assumption in part and could make the p-value too
#    small. Here the null distribution is rebuilt by permuting blocks of 2, 3 or 4 consecutive dates
#    (a block permutation), which keeps the dependence between neighbouring dates within each block.
#    If short-range autocorrelation inflated the original p-values, they would increase markedly here.
# 2. Heterogeneous dispersion. The pseudo-F responds to differences in location (centroids) and also
#    to differences in dispersion, especially with unbalanced groups (Anderson and Walsh 2013), and
#    many candidate splits are unbalanced. Here the whole search is repeated with a statistic that is
#    robust to unequal dispersion, the two-group Welch statistic on distances (Alekseyenko 2016):
#    T2w = d(c1, c2)^2 / (s1^2/n1 + s2^2/n2), where d(c1, c2) is the distance between group centroids
#    and s_i^2 the within-group variance. If the break were mainly a change in dispersion, T2w would
#    not locate it, or it would not be significant.
#
# Both checks use their own seed and run after the main analysis, so they do not change any other result.

source("R/functions.R")

samples     <- load_samples()
descriptors <- load_descriptors()
D2_list     <- map(descriptors[main_descriptors], sq)
dates       <- sort(unique(samples$date))
date_index  <- match(samples$date, dates)
n_dates     <- length(dates)
n_perm      <- 999

# Statistics for a two-group split (group: logical, one value per sample) --------------------------------

# Welch statistic computed from squared distances, like the pseudo-F (see pseudo_F in R/functions.R):
# sums of squares come from squared distances, and the squared distance between the two group centroids
# is N / (n1 n2) times the between-group sum of squares.
welch_T2 <- function(D2, group) {
  in_first <- group == group[1]
  n1 <- sum(in_first)
  n2 <- sum(!in_first)
  N  <- n1 + n2

  ss_total  <- sum(D2) / (2 * N)
  ss_first  <- sum(D2[in_first, in_first])   / (2 * n1)
  ss_second <- sum(D2[!in_first, !in_first]) / (2 * n2)

  d2_centroids <- N / (n1 * n2) * (ss_total - ss_first - ss_second)
  variance_first  <- ss_first  / (n1 - 1)
  variance_second <- ss_second / (n2 - 1)
  d2_centroids / (variance_first / n1 + variance_second / n2)
}
stat_fun <- list(pseudo_F = function(D2, g) pseudo_F(D2, g)[["F"]],
                 welch    = welch_T2)

# Statistic of every split "before = dates 1..k, after = dates k+1..n" (at least 3 dates per side)
split_stat <- function(D2, date_rank, stat, min_seg = 3) {
  candidate_k <- min_seg:(max(date_rank) - min_seg)
  map_dbl(candidate_k, \(k) stat_fun[[stat]](D2, date_rank <= k))
}

# A new order of the dates: free permutation, or permutation of blocks of L consecutive dates. The block
# boundaries start at a random offset in every permutation, so that no fixed partition of the record
# is favoured; the order of the dates within each block is kept.
permute_dates <- function(n, L = 1) {
  if (L == 1) {
    return(sample(n))
  }
  offset    <- sample(0:(L - 1), 1)
  block_id  <- (seq_len(n) + offset - 1) %/% L          # e.g. L = 3, offset 1: 0 0 1 1 1 2 2 2 ...
  blocks    <- split(seq_len(n), block_id)
  new_order <- blocks[sample(length(blocks))]            # shuffle whole blocks
  unlist(new_order, use.names = FALSE)
}

# Full search with a given statistic and null; the permuted series is searched exactly like the real one
search_break <- function(D2, stat, L, n_perm) {
  observed <- split_stat(D2, date_index, stat)
  best_k   <- which.max(observed) + 2                    # the first admissible split is k = 3

  max_after_permuting <- function(i) {
    new_order     <- permute_dates(n_dates, L)
    position_of   <- match(seq_len(n_dates), new_order)  # position of each date in the permuted series
    max(split_stat(D2, position_of[date_index], stat))
  }
  null <- map_dbl(seq_len(n_perm), max_after_permuting)

  tibble(last_pre   = dates[best_k],
         first_post = dates[best_k + 1],
         stat_value = max(observed),
         p          = perm_p_value(null, max(observed)))
}

# 1. Block permutation of dates, pseudo-F (as in the main analysis) -----------------------------------------
set.seed(2035)
block_res <- expand_grid(descriptor = main_descriptors, L = c(1, 2, 3, 4)) %>%
  mutate(res = map2(descriptor, L, \(dsc, L) search_break(D2_list[[dsc]], "pseudo_F", L, n_perm))) %>%
  unnest(res) %>%
  mutate(null = if_else(L == 1, "free permutation of dates", paste("blocks of", L, "dates")))

# 2. Dispersion-robust statistic (Welch), free permutation of dates ------------------------------------------
set.seed(2036)
welch_res <- tibble(descriptor = main_descriptors) %>%
  mutate(res = map(descriptor, \(dsc) search_break(D2_list[[dsc]], "welch", 1, n_perm))) %>%
  unnest(res)

# Dispersion test at the detected break (betadisper, 999 permutations), to state where dispersion differs
set.seed(2037)
dispersion_p <- function(descriptor, last_pre) {
  side <- factor(samples$date > last_pre, labels = c("before", "after"))
  permutest(betadisper(descriptors[[descriptor]], side), permutations = 999)$tab$`Pr(>F)`[1]
}
breaks <- block_res %>%
  filter(L == 1) %>%
  select(descriptor, last_pre)
disp_res <- breaks %>%
  mutate(betadisper_p = map2_dbl(descriptor, last_pre, dispersion_p))

# Table S14: one row per descriptor
table_S14 <- block_res %>%
  select(descriptor, null, p) %>%
  pivot_wider(names_from = null, values_from = p) %>%
  left_join(welch_res %>% transmute(descriptor, welch_last_pre = last_pre,
                                    welch_first_post = first_post, welch_p = p),
            by = "descriptor") %>%
  left_join(disp_res %>% select(descriptor, betadisper_p_at_break = betadisper_p),
            by = "descriptor")

write_csv(block_res, "output/22_block_permutation.csv")
write_csv(welch_res, "output/22_welch_statistic.csv")
write_csv(table_S14, "output/22_table_S14.csv")
print(table_S14, width = Inf)
saveRDS(sessionInfo(), "output/sessionInfo_06.rds")
