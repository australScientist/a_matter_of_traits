# Descriptive hierarchical clustering (Table S13, Fig. S12). Clusters are defined from the same
# dissimilarities they describe, so this is a description, not a test.
# Run from the repository root:  Rscript 05_clusters.R   (a few seconds)

source("R/functions.R")
source("R/palette.R")
dir.create("output/figures", showWarnings = FALSE)

samples <- load_samples() %>%
  mutate(phase = factor(phase_of(date), levels = names(phase_cols)))
descriptors <- load_descriptors()

# Cluster membership of every sample when the UPGMA tree is cut into k groups
hc_groups <- function(D, k) {
  D %>%
    hclust(method = "average") %>%
    cutree(k = k)
}

# Table S13: two-group clusters against the periods before and after the transition ---------------------
cluster_descriptors <- list(
  `Traits, cells`     = descriptors$traits_cells,
  `Traits, biovolume` = descriptors$traits_biovolume,
  `Genera, cells`     = descriptors$taxa_genus,
  `Genera, biovolume` = descriptors$taxa_genus_biovol)
in_period <- samples$phase %in% c("pre", "post")        # samples of the transition are left out

agreement_of <- function(D, descriptor) {
  cluster <- hc_groups(D, 2)[in_period]
  period  <- droplevels(samples$phase[in_period])
  crosstab <- table(cluster, period)
  # agreement: share of samples that fall in the majority period of their cluster
  tibble(descriptor                   = descriptor,
         `samples (before + after)`   = sum(crosstab),
         `agreement with periods (%)` = round(100 * sum(apply(crosstab, 1, max)) / sum(crosstab)),
         `smaller cluster (n)`        = min(rowSums(crosstab)))
}
cluster_agreement <- cluster_descriptors %>%
  imap(agreement_of) %>%
  list_rbind()
write_csv(cluster_agreement, "output/21_cluster_agreement.csv")
print(cluster_agreement)

# Fig. S12: cluster of every sample by sampling unit and date ----------------------------------------------
# Each label is a genus cluster (6 groups); its fill is the trait cluster (2 groups, biovolume weighting).
theme_paper <- theme_minimal(base_size = 9) +
  theme(panel.grid.minor  = element_blank(),
        panel.grid.major  = element_line(colour = grid_col, linewidth = 0.3),
        axis.text         = element_text(colour = ink2),
        axis.title        = element_text(colour = ink),
        legend.title      = element_text(colour = ink2, size = 8),
        legend.text       = element_text(colour = ink, size = 8),
        plot.background   = element_rect(fill = "white", colour = NA))

samples_cl <- samples %>%
  mutate(T = paste0("T", hc_groups(descriptors$traits_biovolume, 2)),
         G = paste0("G", hc_groups(descriptors$taxa_genus, 6)))
unit_order <- c("CLM1", "CLM2", "DLM1", "DLM2", "TLE", "TLR", "TSP", "WILE", "WISP")

# each trait cluster takes the light colour of the period in which it has most samples
t_major <- samples_cl %>%
  filter(phase != "transition") %>%
  count(T, phase) %>%
  group_by(T) %>%
  slice_max(n, n = 1) %>%
  ungroup() %>%
  mutate(lab = paste0(T, " (mostly ", if_else(phase == "pre", "before", "after"), ")"))

# coloured bar along the bottom showing the periods
phase_bar <- samples %>%
  group_by(phase) %>%
  summarise(start = min(date), end = max(date))

fig_s12 <- ggplot(samples_cl, aes(date, factor(as.character(unit), levels = unit_order))) +
  geom_label(aes(label = str_remove(G, "G"), fill = T), size = 2.6, label.size = 0.15,
             label.padding = unit(0.12, "lines"), colour = ink, key_glyph = draw_key_rect) +
  scale_fill_manual(values = setNames(phase_light[as.character(t_major$phase)], t_major$T),
                    labels = setNames(t_major$lab, t_major$T), name = "Trait cluster (biovolume)") +
  geom_segment(data = phase_bar, aes(x = start, xend = end, y = 0.3, yend = 0.3, colour = phase),
               inherit.aes = FALSE, linewidth = 3, show.legend = FALSE) +
  scale_colour_manual(values = phase_cols) +
  scale_x_date(date_breaks = "6 months", date_labels = "%b %Y") +
  labs(x = NULL, y = "Sampling unit") +
  theme_paper +
  theme(panel.grid.major.y = element_blank(), legend.position = "bottom")

# File name kept from the earlier numbering (the figure is Fig. S12 in the manuscript)
ggsave("output/figures/Fig_S11_clusters.png", fig_s12, width = 180, height = 95, units = "mm", dpi = 300, bg = "white")
ggsave("output/figures/Fig_S11_clusters.pdf", fig_s12, width = 180, height = 95, units = "mm", bg = "white")
message("Done")
