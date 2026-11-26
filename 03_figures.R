# Figures that can be drawn from the derived data: Figs. 3 and 4 of the manuscript and
# Figs. S2, S3, S4 and S9 of the Supporting Information.
# Run after 01_main_analysis.R and 02_robustness.R:  Rscript 03_figures.R   (about 1 minute)
# Fig. 3b needs output/24_null_maxF.csv (08_null_maxF.R). Fig. 3c uses the trait null model, which needs the counts; its distributions are read from
# data/precomputed/2_null_distributions.rds.

source("R/functions.R")
source("R/palette.R")
suppressPackageStartupMessages({
  library(patchwork)
  library(ggrepel)
})
dir.create("output/figures", showWarnings = FALSE)

samples     <- load_samples() %>% mutate(phase = factor(phase_of(date), levels = names(phase_cols)))
descriptors <- load_descriptors()

theme_base <- theme_minimal(base_size = 9) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_line(colour = grid_col, linewidth = 0.3),
        axis.text = element_text(colour = ink2), axis.title = element_text(colour = ink),
        plot.title = element_text(colour = ink, face = "bold", size = 9),
        strip.text = element_text(colour = ink, hjust = 0),
        plot.background = element_rect(fill = "white", colour = NA))
theme_paper    <- theme_base + theme(legend.title = element_text(colour = ink2, size = 8),
                                     legend.text = element_text(colour = ink, size = 8))   # Figs. 4 and 5
theme_nolegend <- theme_base + theme(legend.position = "none")
save_fig <- function(fig, name, w, h) {
  ggsave(file.path("output/figures", paste0(name, ".png")), fig, width = w, height = h, units = "mm", dpi = 300, bg = "white")
  ggsave(file.path("output/figures", paste0(name, ".pdf")), fig, width = w, height = h, units = "mm", bg = "white")
}
descriptor_info <- desc_tbl
pal <- deframe(select(descriptor_info, descriptor, colour))
lty <- deframe(select(descriptor_info, descriptor, linetype))
scenario_labels <- c("all data" = "All data", "2/3 dates" = "2/3 of dates", "1/2 dates" = "1/2 of dates",
                     "1/3 dates" = "1/3 of dates", "3 locations" = "3 of 5 locations",
                     "1 location" = "1 location", "1 location, 1/2 dates" = "1 location, 1/2 dates")

# Fig. 3: break detection, trait null model and data requirements --------------------------------------
F_profiles <- read_csv("output/3_F_profiles.csv", show_col_types = FALSE)
nulls      <- readRDS("data/precomputed/2_null_distributions.rds")
detection  <- read_csv("output/11_detection_apriori_window.csv", show_col_types = FALSE)

profiles <- F_profiles %>%
  inner_join(descriptor_info, by = "descriptor") %>%
  group_by(descriptor) %>% mutate(F_scaled = F / max(F)) %>% ungroup()
peaks <- profiles %>% group_by(descriptor) %>% slice_max(F_scaled, n = 1) %>% ungroup()
end_labels <- profiles %>%
  group_by(descriptor) %>% slice_max(last_pre, n = 1) %>% ungroup() %>%
  arrange(F_scaled) %>% mutate(y_label = seq(0.12, 0.62, length.out = n()))

fig3a <- ggplot(profiles, aes(last_pre, F_scaled, colour = descriptor, linetype = descriptor)) +
  annotate("rect", xmin = as.Date("2017-12-01"), xmax = as.Date("2018-02-28"),
           ymin = -Inf, ymax = Inf, fill = flood_fill, alpha = 0.9) +
  annotate("text", x = as.Date("2018-01-15"), y = 0.04, label = "refilling", size = 2.6, colour = flood_ink, fontface = "bold") +
  geom_line(linewidth = 0.7) +
  geom_point(data = peaks, size = 2.2, stroke = 0.6, fill = "white", shape = 21) +
  geom_segment(data = end_labels, aes(x = last_pre, xend = last_pre + 25, y = F_scaled, yend = y_label),
               linewidth = 0.3, linetype = "solid") +
  geom_text(data = end_labels, aes(x = last_pre + 30, y = y_label, label = label),
            hjust = 0, size = 2.6, colour = ink) +
  scale_colour_manual(values = pal) + scale_linetype_manual(values = lty) +
  scale_x_date(date_breaks = "6 months", date_labels = "%b %Y") +
  scale_y_continuous(limits = c(0, 1.02)) +
  coord_cartesian(clip = "off") +
  labs(x = "Last sampling date before the split", y = "Pseudo-F / maximum",
       title = "a  Where the community splits best") +
  theme_nolegend

# b: permutation null of the maximum pseudo-F (whole search repeated in every permutation of dates)
null_maxF <- read_csv("output/24_null_maxF.csv", show_col_types = FALSE) %>%
  inner_join(descriptor_info, by = "descriptor") %>%
  mutate(ratio = null_F / observed_F, label = factor(label, levels = rev(descriptor_info$label)))
p_maxF <- null_maxF %>% group_by(descriptor, label) %>%
  summarise(p = (sum(null_F >= first(observed_F)) + 1) / (n() + 1), .groups = "drop") %>%
  mutate(p_lab = sprintf("p = %.3f", p))
x_lab <- 1.08
fig3b <- ggplot(null_maxF, aes(ratio, label)) +
  geom_violin(fill = "#c9c7c1", colour = NA, scale = "width", width = 0.8) +
  geom_point(data = p_maxF, aes(x = 1, colour = descriptor), shape = 18, size = 3) +
  geom_label(data = p_maxF, aes(x = x_lab, label = p_lab), hjust = 0, size = 2.4, colour = ink, fill = "white",
             label.size = 0, label.padding = unit(0.5, "mm")) +
  scale_colour_manual(values = pal) +
  scale_x_continuous(breaks = c(0, 0.5, 1)) +
  coord_cartesian(xlim = c(0, 1.5), clip = "on") +   # 3 permuted maxima of cell-weighted genera lie beyond 1.5
  labs(x = "Maximum pseudo-F / observed", y = NULL, title = "b  Significance of the break") +
  theme_nolegend + theme(panel.grid.major.y = element_blank())

null_tbl <- bind_rows(
  tibble(weighting = "Cell-weighted",      R2 = nulls$cells$null,     obs = nulls$cells$obs),
  tibble(weighting = "Biovolume-weighted", R2 = nulls$biovolume$null, obs = nulls$biovolume$obs)) %>%
  mutate(weighting = fct_inorder(weighting))
observed <- null_tbl %>%
  group_by(weighting) %>%
  summarise(obs = first(obs), p = (sum(R2 >= obs) + 1) / (n() + 1), .groups = "drop") %>%
  mutate(label = sprintf("observed\np = %.3f", p))
fig3c <- ggplot(null_tbl, aes(R2)) +
  geom_histogram(bins = 40, fill = "#c9c7c1", colour = "white", linewidth = 0.2) +
  geom_vline(data = observed, aes(xintercept = obs), colour = ink, linewidth = 0.8) +
  geom_text(data = observed, aes(x = obs, y = Inf, label = label),
            hjust = 1.08, vjust = 1.3, size = 2.6, colour = ink, lineheight = 0.9) +
  facet_wrap(~weighting, ncol = 1) +
  labs(x = expression(R^2~"pre vs post"), y = "Random trait assignments", title = "c  Trait null model") +
  theme_nolegend

heat <- function(df, value, title, legend_name, key_width = 8) {
  df %>%
    inner_join(descriptor_info, by = "descriptor") %>%
    mutate(scenario = factor(scenario_labels[as.character(scenario)], levels = scenario_labels),
           label = factor(label, levels = rev(descriptor_info$label)),
           v = {{ value }}) %>%
    ggplot(aes(scenario, label, fill = v)) +
    geom_tile(colour = "white", linewidth = 1) +
    geom_text(aes(label = sprintf("%.2f", v), colour = v > 0.55), size = 2.4) +
    scale_colour_manual(values = c(`TRUE` = "white", `FALSE` = ink), guide = "none") +
    scale_fill_gradient(low = mag_low, high = mag_high, limits = c(0, 1), name = legend_name) +
    scale_x_discrete(drop = FALSE) +
    labs(x = NULL, y = NULL, title = title) +
    theme_nolegend +
    theme(legend.position = "bottom", legend.key.height = unit(2.5, "mm"), legend.key.width = unit(key_width, "mm"),
          legend.title = element_text(size = 8, colour = ink2, vjust = 0.8),
          panel.grid = element_blank(), axis.text.x = element_text(angle = 40, hjust = 1, vjust = 1))
}
fig3d <- heat(detection, W12, "d  Detection with thinned data", "Detection rate", key_width = 12)

fig3 <- (fig3a + theme(plot.margin = margin(6, 105, 6, 6))) /
  ((fig3b | fig3c) + plot_layout(widths = c(1.25, 1))) /
  (fig3d + theme(plot.margin = margin(6, 60, 6, 6))) +
  plot_layout(heights = c(1, 0.85, 0.8))
save_fig(fig3, "Figure_3", 180, 215)

# Fig. 4a-c: similarity between samples ordered by date ------------------------------------------------------
ord <- samples %>% arrange(date, unit) %>% mutate(pos = row_number())
similarity_long <- function(D, name) {
  M <- as.matrix(D)[ord$sample, ord$sample]
  S <- 1 - M / max(M)
  expand_grid(row = seq_len(nrow(S)), col = seq_len(ncol(S))) %>%
    mutate(similarity = S[cbind(row, col)], panel = name)
}
panels <- c(traits_cells = "a  Traits, cell-weighted",
            traits_biovolume = "b  Traits, biovolume-weighted",
            taxa_genus = "c  Genera, cell-weighted")
sim <- imap_dfr(descriptors[names(panels)], ~ similarity_long(.x, panels[[.y]])) %>%
  mutate(panel = factor(panel, levels = panels))
break_pos <- function(last_pre) max(ord$pos[ord$date <= last_pre]) + 0.5
breaks_tbl <- tibble(pos = c(break_pos(as.Date("2018-06-26")), break_pos(as.Date("2018-02-08"))),
                     type = c("cells", "biovolume"))
year_ticks <- ord %>% group_by(year = format(date, "%Y")) %>% summarise(pos = min(pos)) %>%
  filter(year != "2016")                     # 2016 has a single date; its label would overlap
fig4 <- ggplot(sim, aes(col, row, fill = similarity)) +
  geom_raster() +
  geom_vline(data = breaks_tbl, aes(xintercept = pos, linetype = type), colour = "white", linewidth = 0.5) +
  geom_hline(data = breaks_tbl, aes(yintercept = pos, linetype = type), colour = "white", linewidth = 0.5) +
  scale_linetype_manual(values = c(cells = "solid", biovolume = "22"), guide = "none") +
  scale_fill_gradient(low = mag_low, high = mag_high, limits = c(0, 1), name = "Similarity") +
  scale_x_continuous(breaks = year_ticks$pos, labels = year_ticks$year, expand = c(0, 0)) +
  scale_y_continuous(breaks = year_ticks$pos, labels = year_ticks$year, expand = c(0, 0)) +
  coord_equal() + facet_wrap(~panel) +
  labs(x = "Samples ordered by date", y = "Samples ordered by date") +
  theme_paper +
  theme(panel.grid = element_blank(), strip.text = element_text(colour = ink, face = "bold", hjust = 0, size = 8),
        legend.position = "right", legend.key.height = unit(10, "mm"), legend.key.width = unit(2.5, "mm"))

# Fig. 4d-f: trajectory of date centroids in PCoA space ------------------------------------------------------
pcoa_scores <- function(D) {
  p <- cmdscale(D, k = 2, eig = TRUE)
  pos_eig <- p$eig[p$eig > 0]
  list(scores = as_tibble(p$points, .name_repair = ~ c("Axis1", "Axis2")) %>% mutate(sample = samples$sample),
       pct = round(100 * p$eig[1:2] / sum(pos_eig), 1))
}
pc <- map(descriptors[c("traits_cells", "traits_biovolume", "taxa_genus")], pcoa_scores)
first_post <- c(traits_cells = as.Date("2018-08-13"), traits_biovolume = as.Date("2018-05-16"),
                taxa_genus = as.Date("2018-08-13"))
trajectory_plot <- function(name, title) {
  s <- pc[[name]]$scores %>% left_join(samples, by = "sample")
  cent <- s %>% group_by(date, phase) %>% summarise(Axis1 = mean(Axis1), Axis2 = mean(Axis2), .groups = "drop") %>%
    arrange(date) %>% mutate(xend = lead(Axis1), yend = lead(Axis2), seg_phase = lead(phase))
  # hue = period of the arrival date; lightness = date within the period (light = early, dark = late)
  cent <- cent %>% group_by(seg_phase) %>%
    mutate(rel = if (n() > 1) (row_number() - 1) / (n() - 1) else 1,
           col = case_when(seg_phase == "pre"  ~ colorRampPalette(phase_ramp$pre)(101)[round(rel * 100) + 1],
                           seg_phase == "post" ~ colorRampPalette(phase_ramp$post)(101)[round(rel * 100) + 1],
                           TRUE ~ phase_cols[["transition"]]),
           lw = if_else(seg_phase == "transition", 1.1, 0.6)) %>% ungroup()
  ring <- filter(cent, date == first_post[[name]])
  ggplot(s, aes(Axis1, Axis2)) +
    stat_density_2d(aes(alpha = after_stat(level)), geom = "polygon", fill = "#c9c7c1", show.legend = FALSE) +
    scale_alpha(range = c(0.1, 0.6)) +
    geom_point(size = 0.8, colour = ink2, alpha = 0.35) +
    geom_segment(data = drop_na(cent, xend), aes(xend = xend, yend = yend, colour = col, linewidth = lw),
                 arrow = arrow(type = "closed", length = unit(1.5, "mm"))) +
    scale_colour_identity() + scale_linewidth_identity() +
    geom_point(data = ring, shape = 21, size = 3.4, stroke = 1.1, colour = ink, fill = "white", alpha = 0.9) +
    geom_label_repel(data = ring, aes(label = format(date, "%d %b %Y")), size = 2.4, colour = ink,
                     label.size = 0.15, label.padding = unit(0.12, "lines"), fill = "white",
                     min.segment.length = 0, nudge_y = 0.06 * diff(range(s$Axis2)), seed = 2) +
    geom_point(data = tibble(Axis1 = NA_real_, Axis2 = NA_real_, phase = factor(names(phase_cols), levels = names(phase_cols))),
               aes(fill = phase), shape = 22, size = 0, na.rm = TRUE) +
    scale_fill_manual(values = phase_cols, labels = phase_labels, name = NULL, drop = FALSE,
                      guide = guide_legend(override.aes = list(size = 3.5))) +
    labs(x = sprintf("PCoA 1 (%.1f%%)", pc[[name]]$pct[1]), y = sprintf("PCoA 2 (%.1f%%)", pc[[name]]$pct[2]),
         title = title) +
    theme_paper
}
fig_traj <- trajectory_plot("traits_cells", "d  Traits, cell-weighted") +
  trajectory_plot("traits_biovolume", "e  Traits, biovolume-weighted") +
  trajectory_plot("taxa_genus", "f  Genera, cell-weighted") +
  plot_layout(guides = "collect") & theme(legend.position = "bottom")
# similarity on top, ordination below, one column per descriptor
fig4_all <- wrap_elements(full = fig4) / wrap_elements(full = fig_traj) + plot_layout(heights = c(1, 1.08))
save_fig(fig4_all, "Figure_4", 180, 150)

# Fig. S2: robustness ----------------------------------------------------------------------------------------
false_alarm <- read_csv("output/12_false_alarm_rates.csv", show_col_types = FALSE)
boot_freq   <- read_csv("output/13_break_date_bootstrap_freq.csv", show_col_types = FALSE)
figS2a <- heat(detection, W12, "a  Detection (12-month window)", "Rate")
figS2b <- heat(false_alarm, false_alarm_rate, "b  Significant breaks after the shift", "Rate") +
  theme(axis.text.y = element_blank())
boot_plot <- boot_freq %>% inner_join(descriptor_info, by = "descriptor") %>%
  mutate(label = factor(label, levels = rev(descriptor_info$label)))
figS2c <- ggplot(boot_plot, aes(last_pre, label, size = prop, colour = label)) +
  annotate("rect", xmin = refilling_onset, xmax = as.Date("2018-02-28"), ymin = -Inf, ymax = Inf,
           fill = flood_fill, alpha = 0.9) +
  geom_point(alpha = 0.85) +
  scale_colour_manual(values = desc_cols, guide = "none") +
  scale_size_area(max_size = 7, name = "Share of resamples", breaks = c(0.1, 0.5, 0.9)) +
  scale_x_date(date_breaks = "2 months", date_labels = "%b\n%Y") +
  labs(x = "Last sampling date before the break", y = NULL,
       title = "c  Break date when locations are resampled (bootstrap)") +
  theme_base + theme(legend.position = "bottom", legend.title = element_text(size = 8, colour = ink2))
figS2 <- ((figS2a | figS2b) + plot_layout(guides = "collect") & theme(legend.position = "bottom")) / figS2c +
  plot_layout(heights = c(1.15, 0.85))
save_fig(figS2, "Fig_S2_robustness", 180, 190)

# Fig. S3: detection with different definitions of a detection --------------------------------------------
window_cols <- c("Data-defined window" = "#a3a29c", setNames(mag_ramp(3), c("6 months", "12 months", "18 months")))
det <- detection %>%
  inner_join(descriptor_info, by = "descriptor") %>%
  mutate(scenario = factor(scenario_labels[scenario], levels = scenario_labels),
         label = factor(label, levels = descriptor_info$label))
det_long <- det %>%
  select(scenario, label, `Data-defined window` = old_window, `6 months` = W6, `12 months` = W12, `18 months` = W18) %>%
  pivot_longer(-c(scenario, label), names_to = "window", values_to = "rate") %>%
  mutate(window = factor(window, levels = names(window_cols)))
figS4 <- ggplot(det_long, aes(scenario, rate)) +
  geom_point(data = det, aes(y = any_significant), shape = 95, size = 6, colour = "#c9c8c2") +
  geom_point(aes(colour = window), position = position_dodge(width = 0.7), size = 1.9) +
  facet_wrap(~label, ncol = 1) +
  scale_colour_manual(values = window_cols, name = NULL) +
  scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1)) +
  labs(x = NULL, y = "Detection rate", caption = "Grey bar: share of replicates with a significant break at any date.") +
  theme_base +
  theme(legend.position = "top", axis.text.x = element_text(angle = 30, hjust = 1),
        plot.caption = element_text(colour = ink2, size = 7))
save_fig(figS4, "Fig_S3_window_sensitivity", 180, 160)

# Fig. S4: step beyond a linear and a flexible trend --------------------------------------------------------
st <- read_csv("output/14_step_beyond_spline_main.csv", show_col_types = FALSE) %>%
  inner_join(descriptor_info, by = "descriptor") %>%
  mutate(label = factor(label, levels = rev(descriptor_info$label)))
figS5 <- ggplot(st, aes(y = label)) +
  geom_segment(aes(x = R2_step_spline, xend = R2_step_linear, yend = label, colour = label), linewidth = 0.8) +
  geom_point(aes(x = R2_step_linear, colour = label), shape = 21, fill = "white", size = 2.6, stroke = 0.9) +
  geom_point(aes(x = R2_step_spline, colour = label), size = 2.6) +
  geom_text(aes(x = R2_step_linear, label = sprintf("%.3f", R2_step_linear)), nudge_y = 0.32, size = 2.4, colour = ink2) +
  geom_text(aes(x = R2_step_spline, label = sprintf("%.3f", R2_step_spline)), nudge_y = 0.32, size = 2.4, colour = ink2) +
  scale_colour_manual(values = desc_cols) +
  scale_x_continuous(limits = c(0, 0.17), expand = expansion(mult = c(0, 0.02))) +
  labs(x = expression(Delta*R^2~"due to the step"), y = NULL,
       caption = "Open circle: beyond a linear trend. Filled circle: beyond a flexible (spline) trend.") +
  theme_nolegend + theme(plot.caption = element_text(colour = ink2, size = 7))
save_fig(figS5, "Fig_S4_step_vs_trend", 140, 70)

# Fig. S9: spatial heterogeneity and dominance ---------------------------------------------------------------
het <- read_csv("output/8_mosaic_heterogeneity.csv", show_col_types = FALSE) %>%
  inner_join(descriptor_info, by = "descriptor") %>%
  select(label, spatial_pre, spatial_post, spatial_p, depth_pre, depth_post) %>%
  pivot_longer(c(spatial_pre, spatial_post, depth_pre, depth_post), names_to = c("kind", "period"),
               names_sep = "_", values_to = "d") %>%
  mutate(kind = recode(kind, spatial = "Between locations", depth = "Between depths"),
         period = factor(period, levels = c("pre", "post"), labels = c("Before", "After")),
         label = factor(label, levels = descriptor_info$label))
p_lab <- het %>% filter(kind == "Between locations", period == "After") %>% mutate(txt = sprintf("p = %.3f", spatial_p))
figS9a <- ggplot(het, aes(period, d, group = kind)) +
  geom_line(aes(linetype = kind), linewidth = 0.5, colour = ink2) +
  geom_point(aes(shape = kind, colour = period, fill = period), size = 2, stroke = 0.8) +
  geom_text(data = p_lab, aes(label = txt), hjust = -0.15, size = 2.1, colour = ink2) +
  facet_wrap(~label, nrow = 1, labeller = label_wrap_gen(14)) +
  scale_colour_manual(values = c(Before = phase_cols[["pre"]], After = phase_cols[["post"]]), guide = "none") +
  scale_fill_manual(values = c(Before = phase_cols[["pre"]], After = "white"), guide = "none") +
  scale_linetype_manual(values = c("Between locations" = "solid", "Between depths" = "22"), name = NULL) +
  scale_shape_manual(values = c("Between locations" = 21, "Between depths" = 24), name = NULL) +
  scale_x_discrete(expand = expansion(add = c(0.3, 1.1))) +
  labs(x = NULL, y = "Mean dissimilarity on the same date", title = "a  Spatial heterogeneity") +
  theme_base + theme(legend.position = "bottom", legend.text = element_text(size = 7),
                      strip.text = element_text(size = 7, colour = ink, hjust = 0))
dom <- read_csv("output/8_mosaic_dominance.csv", show_col_types = FALSE) %>%
  pivot_longer(c(BP_pre, BP_post), names_to = "period", values_to = "BP") %>%
  mutate(period = factor(period, levels = c("BP_pre", "BP_post"), labels = c("Before", "After")),
         weighting = recode(weighting, cells = "Cells", biovolume = "Biovolume"))
figS9b <- ggplot(dom, aes(period, BP, group = weighting)) +
  geom_line(aes(linetype = weighting), linewidth = 0.6, colour = ink2) +
  geom_point(aes(colour = period), size = 2) +
  geom_text(data = filter(dom, period == "After"), aes(label = sprintf("%s\np = %.3f", weighting, p)),
            hjust = -0.15, size = 2.3, lineheight = 0.9) +
  scale_colour_manual(values = c(Before = phase_cols[["pre"]], After = phase_cols[["post"]]), guide = "none") +
  scale_linetype_manual(values = weight_lty, guide = "none") +
  scale_x_discrete(expand = expansion(add = c(0.3, 0.9))) +
  scale_y_continuous(limits = c(0.4, 0.8)) +
  labs(x = NULL, y = "Berger-Parker dominance", title = "b  Dominance") +
  theme_nolegend
figS9 <- (figS9a | figS9b) + plot_layout(widths = c(3.2, 1))
save_fig(figS9, "Fig_S9_heterogeneity_dominance", 180, 85)

message("Done: figures in output/figures/")
