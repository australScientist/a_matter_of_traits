# Supporting Information tables S1 to S15 as one Word document (output/Tables_S1_S15.docx).
# Run from the repository root after scripts 01 to 07:  Rscript 09_tables_SI.R
# Every number is read from output/ or data/; nothing is typed by hand except labels and captions.

suppressPackageStartupMessages({ library(tidyverse); library(flextable); library(officer) })

fmt_date <- function(x) format(as.Date(x), "%-d %b %Y")
brk      <- function(a, b) paste(fmt_date(a), "/", fmt_date(b))
f2 <- function(x) formatC(x, format = "f", digits = 2)
f3 <- function(x) formatC(x, format = "f", digits = 3)
pct <- function(x) paste0(round(100 * x), "%")
desc_lab <- c(traits_cells = "Traits, cells", traits_biovolume = "Traits, biovolume",
              traits_cells_dom = "Traits, cells, dominant class", traits_cells_klais = "Traits, cells, Klais fraction for mixotrophy",
              traits_cells_rimet = "Traits, cells, Rimet coding of mixotrophy",
              taxa_genus = "Genera, cells", taxa_genus_biovol = "Genera, biovolume", taxa_class = "Classes, cells")
main_order <- c("traits_cells", "traits_biovolume", "taxa_genus", "taxa_genus_biovol", "taxa_class")
lab <- function(x) unname(desc_lab[x])
rd  <- function(f) read_csv(f, show_col_types = FALSE)

# Formatting shared by all tables: Times New Roman 9 pt, grey header, thin borders
style <- function(ft, widths = NULL) {
  ft <- ft %>% font(fontname = "Times New Roman", part = "all") %>% fontsize(size = 9, part = "all") %>%
    bold(part = "header") %>% bg(bg = "#E7E6E6", part = "header") %>%
    border_outer(border = fp_border(color = "#808080", width = 0.5)) %>%
    border_inner(border = fp_border(color = "#808080", width = 0.5), part = "all") %>%
    align(j = 1, align = "left", part = "all") %>% align(j = -1, align = "center", part = "all") %>%
    valign(valign = "center", part = "all") %>% padding(padding = 2, part = "all")
  if (!is.null(widths)) ft <- width(ft, width = widths) else ft <- autofit(ft)
  ft
}

doc <- read_docx()
portrait  <- prop_section(page_size = page_size(orient = "portrait", width = 8.5, height = 11),
                          page_margins = page_mar(top = 1, bottom = 1, left = 1, right = 1), type = "nextPage")
landscape <- prop_section(page_size = page_size(orient = "landscape", width = 8.5, height = 11),
                          page_margins = page_mar(top = 1, bottom = 1, left = 1, right = 1), type = "nextPage")
cap_prop <- fp_text(font.family = "Times New Roman", font.size = 11)
add_table <- function(doc, number, caption, ft, note = NULL, orient = "portrait") {
  doc <- body_add_fpar(doc, fpar(ftext(paste0("Table ", number, ". "), update(cap_prop, bold = TRUE)),
                                 ftext(caption, cap_prop)))
  doc <- body_add_flextable(doc, ft, align = "left")
  if (!is.null(note)) doc <- body_add_fpar(doc, fpar(ftext(note, update(cap_prop, font.size = 9))))
  body_end_block_section(doc, block_section(if (orient == "portrait") portrait else landscape))
}
sub_head <- function(doc, text) body_add_fpar(doc, fpar(ftext(text, update(cap_prop, italic = TRUE, font.size = 10))))

# S1. Analytical methods -------------------------------------------------------------------------------------
s1 <- rd("data/metadata/analytical_methods.csv") %>%
  transmute(Group = group, Variable = variable, Unit = replace_na(unit, ""),
            `Limit of quantification` = replace_na(as.character(loq), ""), Method = method)
doc <- add_table(doc, "S1", paste(
  "Analytical methods for environmental variables. Limits of quantification in the units given; values at the",
  "limit were recorded as the limit itself. SM, Standard Methods for the Examination of Water and Wastewater",
  "(APHA). Chemical analyses by the laboratory of the drinking water utility."),
  flextable(s1) %>% merge_v(j = "Group") %>% style(c(1.3, 1.9, 0.9, 0.8, 1.6)))

# S2. Traits -------------------------------------------------------------------------------------------------
s2 <- tribble(
  ~Trait,                          ~Unit,               ~Type,        ~Transformation, ~Rationale,
  "Cell biovolume",                "µm³",               "Continuous", "log₁₀",  "Metabolism, growth and susceptibility to grazing scale with size",
  "Surface-to-volume ratio",       "µm⁻¹",              "Continuous", "none",   "Nutrient and light uptake per unit volume",
  "Cells per colony",              "cells",             "Continuous", "none",   "Colony size; protection from grazing, sinking and buoyancy",
  "Organic carbon ratio",          "pg C µm⁻³",         "Continuous", "none",   "Carbon content per unit volume (Rimet and Druart 2018)",
  "Flagella",                      "proportion",        "Binary",     "none",   "Active vertical movement",
  "Filamentous growth",            "proportion",        "Binary",     "none",   "Colony form; grazing resistance, light harvesting",
  "Maximum growth rate",           "d⁻¹",               "Continuous", "log₁₀",  "Capacity for rapid growth after disturbance",
  "Phosphate affinity",            "L µmol⁻¹ d⁻¹",      "Continuous", "log₁₀",  "Competitive ability for phosphate",
  "Nitrogen fixation",             "proportion",        "Binary",     "none",   "Independence from dissolved inorganic nitrogen",
  "Buoyancy regulation",           "proportion",        "Binary",     "none",   "Position in the water column under stratification",
  "Silica requirement",            "proportion",        "Binary",     "none",   "Dependence on silica; heavy, sinking cells",
  "Phycocyanin",                   "proportion",        "Binary",     "none",   "Use of green and orange light",
  "Chlorophyll c",                 "proportion",        "Binary",     "none",   "Use of blue-green light",
  "Potential constitutive mixotrophy", "proportion",    "Binary",     "none",   "Phagotrophy as a source of nutrients and carbon")
src <- rd("data/traits/10_SM_trait_sources.csv") %>%
  mutate(Trait = s2$Trait[match(trait, c("cell.biovolume", "surface.volume.ratio", "colony.cell.number",
                                         "organic.Carbon.ratio", "flagella", "Filament", "maximum_growth_rate",
                                         "phosphate_affinity", "nfix", "buoyancy", "silica", "fyco", "chlc",
                                         "mixotrophic"))],
         source = str_replace_all(source, "&", "and"))
s2 <- s2 %>% left_join(src %>% select(Trait, Source = source, `Assignment to genus` = level_rule), by = "Trait")
doc <- add_table(doc, "S2", paste(
  "Trait definitions, sources and assignment rules. Binary traits enter the community-weighted means as the",
  "proportion of the community (cells or biovolume) that carries the trait. Units of growth rate and phosphate",
  "affinity as in Bruggeman (2011)."),
  flextable(s2) %>% style(c(1.4, 0.8, 0.7, 0.8, 2.1, 1.2, 1.5)), orient = "landscape")

# S3. Genus-level evidence -------------------------------------------------------------------------------------
s3 <- rd("data/traits/10_SM_genus_evidence.csv") %>%
  transmute(Genus = genus, Class = class,
            `Rimet and Druart (2018): taxa averaged` = rimet_n_taxa,
            `Klais et al. (2017): genus entry (species)` = if_else(klais_genus_entry, paste0("yes (", n_species_klais, ")"), "no"),
            `Bruggeman (2011): level` = replace_na(bruggeman_level, "missing"),
            `Evidence of phagotrophy` = if_else(mixotrophy_evidence == "none", "", mixotrophy_evidence))
doc <- add_table(doc, "S3", paste(
  "Source and taxonomic level of the trait values for each genus. Rimet and Druart (2018): number of taxa of the",
  "genus whose values were averaged. Klais et al. (2017): whether the database has a genus-level entry, with the",
  "number of species it summarizes. Bruggeman (2011): level from which growth rate and phosphate affinity were",
  "taken. Evidence of phagotrophy: source documenting phagotrophy in at least one species of the genus (potential",
  "constitutive mixotrophy coded 1). Anabaena is reported as Dolichospermum in the text."),
  flextable(s3) %>% style(c(1.1, 0.7, 1.0, 1.1, 1.4, 1.2)))

# S4. All descriptors and variants -----------------------------------------------------------------------------
base <- rd("output/1_baseline.csv"); cp <- rd("output/3_change_point.csv")
svt <- rd("output/4_step_vs_trend.csv"); spl <- rd("output/14_step_beyond_spline_main.csv")
s4 <- base %>% left_join(cp, by = "descriptor") %>% left_join(svt, by = "descriptor") %>%
  left_join(spl %>% select(descriptor, R2_step_spline, p_spline), by = "descriptor") %>%
  transmute(Descriptor = lab(descriptor),
            `Periods split in mid-June 2018: R² (p)` = paste0(f3(R2.x), " (", formatC(p, format = "f", digits = 4), ")"),
            `Dispersion: p` = f3(betadisper_p),
            `PCoA axes 1, 2 (%)` = paste0(PCoA1_pct, ", ", PCoA2_pct),
            `Break: last date before / first date after` = brk(last_pre, first_post),
            `R² at break` = f3(R2.y), `p (date permutation)` = f3(p_date_perm),
            `Step beyond linear trend, July 2018: ΔR² (p)` = paste0(f3(R2_stepJul_beyond_trend), " (", f3(p_Jul), ")"),
            `Step beyond linear trend, March 2018: ΔR² (p)` = paste0(f3(R2_stepMar_beyond_trend), " (", f3(p_Mar), ")"),
            `Step at break beyond spline: ΔR² (p)` = if_else(is.na(R2_step_spline), "",
                                                             paste0(f3(R2_step_spline), " (", f3(p_spline), ")")))
doc <- add_table(doc, "S4", paste(
  "Results for all descriptors and sensitivity variants. Baseline: PERMANOVA between the periods before and after",
  "mid-June 2018 (9999 permutations) and test of homogeneity of dispersions (betadisper, 999 permutations).",
  "Break: change point with maximum pseudo-F, p from 999 permutations of dates that include the search. Steps:",
  "increase in R² when a step in July or March 2018 is added to a model with sampling unit and a linear trend,",
  "and when a step at the detected break is added to a model with sampling unit and a natural cubic spline (3 df;",
  "main descriptors only); permutations within sampling units. Dominant class: binary traits coded by the dominant",
  "class instead of proportions; Klais fraction and Rimet coding: alternative codings of mixotrophy."),
  flextable(s4) %>% style(c(1.4, 0.9, 0.75, 0.6, 1.2, 0.5, 0.75, 0.85, 0.85, 0.85)), orient = "landscape")

# S5. Leave one genus out ---------------------------------------------------------------------------------------
full_F <- cp %>% filter(descriptor %in% c("traits_cells", "traits_biovolume")) %>%
  transmute(weighting = if_else(descriptor == "traits_cells", "cells", "biovolume"), F_full = maxF, lp_full = last_pre)
s5 <- rd("data/precomputed/7_leave_one_genus_out.csv") %>% left_join(full_F, by = "weighting") %>%
  mutate(value = paste0(fmt_date(last_pre), if_else(last_pre != lp_full, " *", ""), "; ", formatC(maxF, format = "f", digits = 1),
                        " (", sprintf("%+.1f", maxF - F_full), ")")) %>%
  select(Genus = dropped, weighting, value) %>%
  pivot_wider(names_from = weighting, values_from = value) %>% arrange(Genus) %>%
  rename(`Traits, cells: last date before break; max pseudo-F (change)` = cells,
         `Traits, biovolume: last date before break; max pseudo-F (change)` = biovolume)
doc <- add_table(doc, "S5", paste0(
  "Leave-one-genus-out analysis for trait descriptors. Each genus was removed and the change-point search repeated. ",
  "Maximum pseudo-F with all genera: ", full_F$F_full[full_F$weighting == "cells"], " (cells) and ",
  full_F$F_full[full_F$weighting == "biovolume"], " (biovolume). Asterisk, break moved."),
  flextable(s5) %>% style(c(1.4, 2.5, 2.5)))

# S6. Spatial heterogeneity, variation partitioning and dominance --------------------------------------------
het <- rd("output/8_mosaic_heterogeneity.csv")
s6 <- het %>% transmute(
  Descriptor = lab(descriptor),
  `Between locations, same date: before → after (p)` = paste0(f3(spatial_pre), " → ", f3(spatial_post), " (", f3(spatial_p), ")"),
  `Between depths, same location: before → after` = paste0(f3(depth_pre), " → ", f3(depth_post)),
  `Same unit, consecutive dates: before → after` = paste0(f3(turnover_pre), " → ", f3(turnover_post)),
  `R² date / location, before` = paste0(f3(R2_date_pre), " / ", f3(R2_location_pre)),
  `R² date / location, after` = paste0(f3(R2_date_post), " / ", f3(R2_location_post)))
dom <- rd("output/8_mosaic_dominance.csv") %>%
  transmute(Weighting = str_to_sentence(weighting), `Berger–Parker index, before` = f3(BP_pre),
            `Berger–Parker index, after` = f3(BP_post), p = f3(p))
doc <- body_add_fpar(doc, fpar(ftext("Table S6. ", update(cap_prop, bold = TRUE)), ftext(paste(
  "Spatial heterogeneity, variation partitioning and dominance before and after the transition, at the five",
  "long-term locations. (a) Mean dissimilarity between samples and R² of date and location in PERMANOVA models",
  "fitted within each period; p from permutation of period labels among dates. (b) Mean Berger–Parker dominance",
  "per sample, by cells and by biovolume."), cap_prop)))
doc <- sub_head(doc, "(a) Dissimilarity and variation partitioning")
doc <- body_add_flextable(doc, flextable(s6) %>% style(c(1.3, 1.6, 1.3, 1.3, 1.1, 1.1)), align = "left")
doc <- sub_head(doc, "(b) Dominance")
doc <- body_add_flextable(doc, flextable(dom) %>% style(c(1.2, 1.4, 1.4, 0.6)), align = "left")
doc <- body_end_block_section(doc, block_section(landscape))

# S7. Detection under thinning ----------------------------------------------------------------------------------
scen_order <- c("all data", "2/3 dates", "1/2 dates", "1/3 dates", "3 locations", "1 location", "1 location, 1/2 dates")
s7 <- rd("output/11_detection_apriori_window.csv") %>%
  mutate(scenario = factor(scenario, scen_order), descriptor = factor(descriptor, main_order)) %>%
  arrange(scenario, descriptor) %>%
  transmute(Scenario = as.character(scenario), Descriptor = lab(as.character(descriptor)), Replicates = n_reps,
            `Significant break at any date` = f2(any_significant), `Window from the full data` = f2(old_window),
            `6 months` = f2(W6), `12 months` = f2(W12), `18 months` = f2(W18))
doc <- add_table(doc, "S7", paste(
  "Detection under thinning: proportion of replicates with a significant break (p < 0.05, 999 date permutations),",
  "and with a significant break whose interval overlaps a detection window: the window derived from the full data",
  "(8 February to 13 August 2018), used in earlier versions and circular, and windows of 6, 12 and 18 months from",
  "the onset of the refilling (1 December 2017), defined a priori. All windows are evaluated on the same replicates;",
  "the 12-month window is the one reported in the main text. Spatial thinning drew whole locations from the five",
  "long-term locations."),
  flextable(s7) %>% merge_v(j = "Scenario") %>% style(c(1.2, 1.3, 0.7, 0.8, 0.8, 0.6, 0.6, 0.6)))

# S8. False alarms --------------------------------------------------------------------------------------------
stab <- rd("output/12_false_alarm_stable_stretches.csv") %>%
  mutate(stretch = factor(stretch, c("pre", "post"), c("Before the refilling", "After the shift")),
         descriptor = factor(descriptor, main_order)) %>% arrange(stretch, descriptor) %>%
  transmute(Stretch = as.character(stretch), Descriptor = lab(as.character(descriptor)), Dates = n_dates,
            `Minimum dates per side` = min_seg, `Break` = brk(last_pre, first_post), `R²` = f3(R2),
            `R² of the main break` = f3(R2_main_break), p = f3(p))
far <- rd("output/12_false_alarm_rates.csv") %>%
  mutate(scenario = factor(scenario, scen_order), descriptor = factor(descriptor, main_order)) %>%
  arrange(scenario, descriptor) %>% select(scenario, descriptor, false_alarm_rate) %>%
  mutate(descriptor = lab(as.character(descriptor))) %>%
  pivot_wider(names_from = descriptor, values_from = false_alarm_rate) %>%
  mutate(across(-scenario, f2)) %>% rename(Scenario = scenario) %>% mutate(Scenario = as.character(Scenario))
lt <- rd("output/12_post_stretch_long_term.csv") %>% mutate(descriptor = factor(descriptor, main_order)) %>%
  arrange(descriptor) %>%
  transmute(Descriptor = lab(as.character(descriptor)), Break = brk(last_pre, first_post), `R²` = f3(R2), p = f3(p))
byloc <- rd("output/12_post_stretch_by_location.csv") %>%
  transmute(Location = location, Descriptor = lab(descriptor), Break = brk(last_pre, first_post), `R²` = f3(R2), p = f3(p))
shares <- rd("data/precomputed/12_post_stretch_biovolume_shares.csv") %>%
  mutate(across(-genus, pct)) %>% rename(Genus = genus)
doc <- body_add_fpar(doc, fpar(ftext("Table S8. ", update(cap_prop, bold = TRUE)), ftext(paste(
  "False alarms in stretches without a regime shift. (a) Change-point search in the stretch before the refilling",
  "(December 2016 to November 2017) and in the stretch after the shift (from 13 August 2018), where no regime shift",
  "is expected; R² of the main break given for comparison. (b) Proportion of thinned replicates of the post-shift",
  "stretch with a significant break (scenarios as in Table S7). (c) Post-shift stretch with the five long-term",
  "locations only and (d) by location, for biovolume-weighted descriptors. (e) Mean share of biovolume of the main",
  "genera in the post-shift stretch, before and after the replacement of Closterium by Ceratium in late 2018. All",
  "p-values from 999 date permutations."), cap_prop)))
doc <- sub_head(doc, "(a) Stretches without an expected shift")
doc <- body_add_flextable(doc, flextable(stab) %>% merge_v(j = "Stretch") %>% style(c(1.2, 1.3, 0.5, 0.7, 1.7, 0.5, 0.7, 0.5)), align = "left")
doc <- sub_head(doc, "(b) False alarm rate in thinned post-shift data")
doc <- body_add_flextable(doc, flextable(far) %>% style(c(1.4, 1, 1, 1, 1, 1)), align = "left")
doc <- sub_head(doc, "(c) Post-shift stretch, five long-term locations")
doc <- body_add_flextable(doc, flextable(lt) %>% style(c(1.6, 1.8, 0.6, 0.6)), align = "left")
doc <- sub_head(doc, "(d) Post-shift stretch by location")
doc <- body_add_flextable(doc, flextable(byloc) %>% merge_v(j = "Location") %>% style(c(0.8, 1.6, 1.8, 0.6, 0.6)), align = "left")
doc <- sub_head(doc, "(e) Share of biovolume in the post-shift stretch")
doc <- body_add_flextable(doc, flextable(shares) %>% style(c(1.3, 1.3, 1.5)), align = "left")
doc <- body_end_block_section(doc, block_section(landscape))

# S9. Planted shifts --------------------------------------------------------------------------------------------
s9 <- rd("data/precomputed/14_planted_shift_summary.csv") %>%
  mutate(shape = factor(shape, c("none", "step", "ramp"), c("None", "Step", "Gradual drift")),
         descriptor = factor(descriptor, main_order)) %>% arrange(shape, alpha, descriptor) %>%
  transmute(`Planted change` = as.character(shape), `Size (fraction of the observed shift)` = f2(alpha),
            Descriptor = lab(as.character(descriptor)), Series = n,
            `Significant break` = f2(any_significant),
            `Detected at the planted date ± 1 date` = if_else(is.na(detected_near), "", f2(detected_near)),
            `Step beyond linear trend` = f2(step_beyond_trend), `Step beyond spline` = f2(step_beyond_spline),
            `Median R² of the step` = f3(median_R2_step))
doc <- add_table(doc, "S9", paste(
  "Planted shifts: proportion of series with a significant break, with a break within one sampling date of the",
  "planted date, and with a significant step at the break beyond a linear and a flexible (natural cubic spline,",
  "3 df) trend. Changes planted into the post-shift stretch after its 6th, 9th or 12th date, as steps or as gradual",
  "drifts of the same final size; None, the stretch without a planted change. Series built with the 21 combinations",
  "of five of the seven long-term sampling units; 999 date permutations."),
  flextable(s9) %>% merge_v(j = c("Planted change", "Size (fraction of the observed shift)")) %>%
    style(c(0.9, 0.9, 1.3, 0.5, 0.7, 0.9, 0.8, 0.8, 0.8)), orient = "landscape")

# S10. Break-date uncertainty -----------------------------------------------------------------------------------
boot <- rd("output/13_break_date_bootstrap.csv"); jack <- rd("output/13_break_date_jackknife.csv")
s10 <- boot %>% left_join(jack, by = "descriptor") %>% mutate(descriptor = factor(descriptor, main_order)) %>%
  arrange(descriptor) %>%
  transmute(Descriptor = lab(as.character(descriptor)), `Most frequent break` = brk(modal_last_pre, modal_first_post),
            `Share of resamples` = f3(prop_modal), `Breaks covering 95%` = n_breaks_95,
            `Range of the 95% set` = paste(fmt_date(set95_from), "to", fmt_date(set95_to)),
            `Within 12 months of the refilling` = f3(prop_in_W12),
            `Jackknife: dates removed without moving the break` = paste(n_same, "of", n_dropped),
            `Jackknife: alternative break (date removed)` = if_else(is.na(alternatives), "",
              paste0("after ", fmt_date(alternatives), " (", fmt_date(dropped_when_moved), ")")))
doc <- add_table(doc, "S10", paste(
  "Uncertainty of the break date. Bootstrap: the five long-term locations resampled with replacement (999",
  "resamples, both depths of a location kept together) and the search repeated. Jackknife: one sampling date",
  "removed at a time (27 series). Sampling units resampled instead of locations are reported in Table S15."),
  flextable(s10) %>% style(c(1.3, 1.6, 0.7, 0.7, 1.5, 0.8, 1.0, 1.4)), orient = "landscape")

# S11. Cell volumes ---------------------------------------------------------------------------------------------
vol <- rd("data/precomputed/15_volume_perturbation.csv") %>% mutate(descriptor = factor(descriptor, main_order)) %>%
  arrange(desc(error), descriptor) %>%
  transmute(`Error in cell volume` = error, Descriptor = lab(as.character(descriptor)), Replicates = n,
            `Significant break` = f2(significant), `Same break as with database volumes` = f2(same_break),
            `Median R² at break` = f3(median_R2))
col <- rd("data/precomputed/15_colony_volume_results.csv") %>%
  transmute(Descriptor = lab(descriptor), `Break with colony volume` = brk(last_pre, first_post),
            `Last date before break with cell volume` = fmt_date(original_last_pre), `R²` = f3(R2), p = f3(p),
            `Step beyond linear trend: ΔR² (p)` = paste0(f3(R2_step_linear), " (", f3(p_step_linear), ")"),
            `Step beyond spline: ΔR² (p)` = paste0(f3(R2_step_spline), " (", f3(p_step_spline), ")"))
doc <- body_add_fpar(doc, fpar(ftext("Table S11. ", update(cap_prop, bold = TRUE)), ftext(paste(
  "Sensitivity to estimated cell volumes. (a) The volume of each genus multiplied by an independent random factor",
  "drawn log-uniformly between 1/1.5 and 1.5 (moderate) or 1/2 and 2 (large), applied consistently to biovolume",
  "weights, cell volume and surface-to-volume ratio; 100 replicates, 999 date permutations each. (b) Cell volume",
  "of colonial genera replaced by the volume of the colony envelope including mucilage divided by the number of",
  "cells per colony (Rimet and Druart 2018)."), cap_prop)))
doc <- sub_head(doc, "(a) Random errors in cell volume")
doc <- body_add_flextable(doc, flextable(vol) %>% merge_v(j = 1) %>% style(c(1.2, 1.5, 0.8, 0.8, 1.2, 0.9)), align = "left")
doc <- sub_head(doc, "(b) Colony volume including mucilage")
doc <- body_add_flextable(doc, flextable(col) %>% style(c(1.4, 1.7, 1.0, 0.5, 0.5, 1.2, 1.2)), align = "left")
doc <- body_end_block_section(doc, block_section(landscape))

# S12. Distance metric ------------------------------------------------------------------------------------------
metric_lab <- c(gower = "Gower (main analysis)", eucl_z = "Euclidean, standardized", eucl_range = "Euclidean, range-scaled",
                manh_z = "Manhattan, standardized")
s12 <- rd("output/18_metric_sensitivity.csv") %>%
  transmute(Weighting = str_to_sentence(weights), Distance = unname(metric_lab[metric]),
            `R² between periods` = f3(R2_period), Break = brk(last_pre, first_post),
            `Spearman correlation with the Gower matrix` = f3(rho_vs_gower))
doc <- add_table(doc, "S12", paste(
  "Sensitivity to the distance metric in trait space: R² between the periods before and after mid-June 2018,",
  "break dates and rank correlation of each dissimilarity matrix with the Gower matrix."),
  flextable(s12) %>% merge_v(j = "Weighting") %>% style(c(1.0, 1.7, 0.9, 1.7, 1.2)))

# S13. Clusters -------------------------------------------------------------------------------------------------
s13 <- rd("output/21_cluster_agreement.csv")
doc <- add_table(doc, "S13", paste(
  "Agreement between two-group hierarchical clusters (UPGMA) and the periods before and after the transition, for",
  "traits and genera weighted by cells and by biovolume. Samples of the transition are excluded. Clusters are",
  "defined from the same dissimilarities they describe and are not a test."),
  flextable(s13) %>% style(c(1.6, 1.4, 1.6, 1.4)))

# S14. Assumptions of the change-point test ---------------------------------------------------------------------
s14 <- rd("output/22_table_S14.csv") %>% mutate(descriptor = factor(descriptor, main_order)) %>% arrange(descriptor) %>%
  transmute(Descriptor = lab(as.character(descriptor)),
            `Dates permuted one by one` = f3(`free permutation of dates`), `Blocks of 2 dates` = f3(`blocks of 2 dates`),
            `Blocks of 3 dates` = f3(`blocks of 3 dates`), `Blocks of 4 dates` = f3(`blocks of 4 dates`),
            `Welch statistic: break` = brk(welch_last_pre, welch_first_post), `Welch statistic: p` = f3(welch_p),
            `Dispersion at the break: p` = f3(betadisper_p_at_break))
doc <- add_table(doc, "S14", paste(
  "Checks of the assumptions of the change-point test. Columns 2 to 5: p-value of the change point with maximum",
  "pseudo-F (Table 1) under four null distributions, with 999 permutations each: dates permuted freely, as in the",
  "main analysis, and blocks of two, three or four consecutive dates permuted while keeping the order within each",
  "block, which preserves dependence between neighbouring dates (blocks start at a random offset in every",
  "permutation; with 27 dates, blocks of two, three and four dates span about 1.5, 3 and 4 months). The p-values",
  "with free permutation differ slightly from those of Table 1 because they come from an independent set of",
  "permutations. Columns 6 and 7: change point found when the whole search is repeated with a statistic robust to",
  "unequal dispersions, the two-group Welch statistic on distances (Alekseyenko 2016), and its p-value from 999",
  "free permutations of dates. Column 8: test of homogeneity of multivariate dispersions between the samples",
  "before and after the change point of Table 1 (betadisper, 999 permutations; Anderson 2006)."),
  flextable(s14) %>% style(c(1.4, 0.8, 0.7, 0.7, 0.7, 1.7, 0.7, 0.8)), orient = "landscape")

# S15. Sampling units instead of locations ----------------------------------------------------------------------
if (file.exists("output/23_bootstrap_units.csv")) {
  sim <- rd("output/23_unit_similarity.csv")
  type_order <- c("same location, other depth", "between inflows", "other location, same depth",
                  "other location, other depth", "reservoir and inflow", "same unit, next date")
  s15a <- sim %>% mutate(type = factor(type, type_order), descriptor = factor(descriptor, main_order)) %>%
    arrange(type, descriptor) %>%
    transmute(`Pair of samples` = str_to_sentence(as.character(type)), Descriptor = lab(as.character(descriptor)),
              Pairs = n_pairs, `Median dissimilarity` = f3(median_dissimilarity),
              `Mean coupling` = if_else(is.na(mean_coupling), "", f2(mean_coupling)))
  det_old <- rd("output/11_detection_apriori_window.csv") %>%
    filter(scenario %in% c("3 locations", "1 location", "1 location, 1/2 dates")) %>%
    transmute(scenario = recode(scenario, "3 locations" = "3", "1 location" = "1", "1 location, 1/2 dates" = "1, 1/2 dates"),
              descriptor, by_location = W12)
  det_new <- rd("output/23_detection_units.csv") %>%
    transmute(scenario = recode(scenario, "3 units" = "3", "1 unit" = "1", "1 unit, 1/2 dates" = "1, 1/2 dates"),
              descriptor, by_unit = W12)
  s15b <- det_old %>% left_join(det_new, by = c("scenario", "descriptor")) %>%
    mutate(scenario = factor(scenario, c("3", "1", "1, 1/2 dates")), descriptor = factor(descriptor, main_order)) %>%
    arrange(scenario, descriptor) %>%
    transmute(`Locations or units kept` = as.character(scenario), Descriptor = lab(as.character(descriptor)),
              `Locations (both depths together)` = f2(by_location), `Sampling units` = f2(by_unit))
  bu <- rd("output/23_bootstrap_units.csv"); bl <- rd("output/13_break_date_bootstrap.csv")
  s15c <- bl %>% select(descriptor, modal_last_pre, modal_first_post, prop_modal, n_breaks_95) %>%
    left_join(bu, by = "descriptor", suffix = c("_loc", "_unit")) %>%
    mutate(descriptor = factor(descriptor, main_order)) %>% arrange(descriptor) %>%
    transmute(Descriptor = lab(as.character(descriptor)),
              `Locations: most frequent break (share)` = paste0(brk(modal_last_pre_loc, modal_first_post_loc), " (", f3(prop_modal_loc), ")"),
              `Locations: breaks covering 95%` = n_breaks_95_loc,
              `Units: most frequent break (share)` = paste0(brk(modal_last_pre_unit, modal_first_post_unit), " (", f3(prop_modal_unit), ")"),
              `Units: breaks covering 95%` = n_breaks_95_unit)
  doc <- body_add_fpar(doc, fpar(ftext("Table S15. ", update(cap_prop, bold = TRUE)), ftext(paste0(
    "Sampling units instead of locations as resampling units. (a) Similarity and coupling of the seven long-term ",
    "sampling units: median dissimilarity between samples of the same date, by type of pair, and between ",
    "consecutive dates of the same unit; coupling is the mean cosine between the deviations of two samples from ",
    "the centroid of their date in principal coordinates space (20 axes). Without coupling its expected value is ",
    "about ", f2(first(sim$coupling_if_independent)), "; positive or less negative values mean that the two samples ",
    "deviate from the date mean in the same direction. (b) Detection within 12 months of the onset of the refilling ",
    "when locations (both depths at CLM and DLM together, as in Table S7) or single sampling units are drawn; 150 ",
    "replicates, 999 date permutations. (c) Bootstrap of the break date resampling the five locations (Table S10) or ",
    "the seven sampling units (999 resamples)."), cap_prop)))
  doc <- sub_head(doc, "(a) Similarity and coupling of sampling units")
  doc <- body_add_flextable(doc, flextable(s15a) %>% merge_v(j = 1) %>% style(c(1.8, 1.3, 0.6, 0.9, 0.8)), align = "left")
  doc <- sub_head(doc, "(b) Detection when locations or sampling units are thinned")
  doc <- body_add_flextable(doc, flextable(s15b) %>% merge_v(j = 1) %>% style(c(1.1, 1.3, 1.3, 1.0)), align = "left")
  doc <- sub_head(doc, "(c) Break-date bootstrap")
  doc <- body_add_flextable(doc, flextable(s15c) %>% style(c(1.3, 2.2, 0.8, 2.2, 0.8)), align = "left")
  doc <- body_set_default_section(doc, landscape)   # last section: no trailing blank page
}

print(doc, target = "output/Tables_S1_S15.docx")
