# Detecting a phytoplankton regime shift from traits: reproducibility package

Code and derived data for

> Marconi, G., et al. A matter of traits: an objective procedure to detect regime shifts in phytoplankton communities from short, irregular monitoring records *Limnology and Oceanography: Methods*. [[VERIFY: final title, authors, year and DOI]]

The phytoplankton counts and the environmental monitoring data of Los Molinos reservoir belong to R. Bazán (Universidad Nacional de Córdoba) and are not public. This package contains the **derived data on which the detection procedure operates**: the dissimilarity matrices between samples, the community-weighted trait means and the sample metadata. With them, every analysis of the paper that does not need the raw counts can be rerun, and its output is identical to that of the original analysis.

To apply the procedure to other data, use the R package `traitshift` github.com/australScientist/TraitShift (zenodo DOI), which includes simulated example data.

## How to run

From the root of the repository, in R 4.3 or later with the packages `tidyverse`, `vegan`, `patchwork` and `ggrepel`:

```sh
Rscript 01_main_analysis.R   # about 20 to 40 minutes
Rscript 02_robustness.R      # a few hours; needs the output of 01
Rscript 03_figures.R         # about 1 minute; needs the outputs of 01 and 02
Rscript 04_refilling.R       # a few seconds; public hydrological data only
Rscript 05_clusters.R        # a few seconds
Rscript 06_dependence_dispersion.R  # about 2 minutes
Rscript 07_units_sensitivity.R      # about 1 hour
Rscript 08_null_maxF.R              # about 15 minutes; permutation null behind the p-values of Table 1
Rscript 10_temporal_memory.R        # about 1 to 2 hours; temporal memory and calibrated change-point test
Rscript 09_tables_SI.R              # a few seconds; Tables S1 to S16 as one Word file (needs flextable, officer)
```

Run `08_null_maxF.R` before `03_figures.R` (Fig. 3b reads its output).

Results are written to `output/` and figures to `output/figures/`. The files already in `output/` are the results of the run used for the paper, so a new run can be compared with them; the scripts overwrite them.

## What reproduces what

| Manuscript item | Script | Output |
|---|---|---|
| Table 1: break dates, pseudo-F, R², p from date permutation | `01_main_analysis.R` | `3_change_point.csv` |
| Table 1: R² with and without the step, beyond linear and flexible trends | `01_main_analysis.R` | `19_step_trend_model_R2.csv` |
| Table 1: p of the step at the break, beyond linear and flexible trends (whole dates permuted) | `02_robustness.R` | `14_step_beyond_spline_main.csv` |
| Table S4: baseline PERMANOVA, dispersion and PCoA for all descriptors | `01_main_analysis.R` | `1_baseline.csv`, `3_change_point.csv`, `4_step_vs_trend.csv` |
| Table S12: sensitivity to the distance metric | `01_main_analysis.R` | `18_metric_sensitivity.csv` |
| Table S6, Fig. S10: spatial heterogeneity, variation partitioning and dominance | `01_main_analysis.R`, `03_figures.R` | `8_mosaic_heterogeneity.csv`, `8_mosaic_dominance.csv` |
| Fig. 3d, Fig. S4, Table S7: detection with thinned data and an a priori window | `02_robustness.R` | `11_detection_apriori_window.csv` |
| Table S8, Fig. S3b: false alarms in stretches without a shift | `02_robustness.R` | `12_*.csv` |
| Table S10, Fig. S3c: uncertainty of the break date | `02_robustness.R` | `13_*.csv` |
| Fig. S5: step beyond a linear and a flexible trend | `02_robustness.R` | `14_step_beyond_spline_main.csv` |
| Figs. 3, 4, S3, S4, S5 and S10 | `03_figures.R` | `output/figures/` |
| Assessment, "Environmental context": size of the 2017–18 refilling and of earlier refills | `04_refilling.R` | `refills.csv` |
| Table S13, Fig. S12: descriptive hierarchical clustering | `05_clusters.R` | `21_cluster_agreement.csv`, `output/figures/Fig_S11_clusters` |
| Fig. 3b: permutation null of the maximum pseudo-F (the permutations behind the p-values of Table 1) | `08_null_maxF.R` | `24_null_maxF.csv` |
| Table S15: sampling units instead of locations as resampling units | `07_units_sensitivity.R` | `23_unit_similarity.csv`, `23_detection_units.csv`, `23_bootstrap_units.csv` |
| Table S16: temporal memory and change-point test calibrated for it | `10_temporal_memory.R` | `25_temporal_memory.csv`, `25_calibrated_change_point.csv` |
| Tables S1 to S16 (Word file) | `09_tables_SI.R` | `Tables_S1_S16.docx` |
| Table S14: change-point test with block permutation of dates (temporal autocorrelation) and with a dispersion-robust statistic | `06_dependence_dispersion.R` | `22_table_S14.csv`, `22_block_permutation.csv`, `22_welch_statistic.csv` |

Figure numbers follow the manuscript (version 4, where the site map is Fig. S1); some file names in `output/figures/` keep the earlier numbering (for example, `Fig_S2_robustness` is Fig. S3 and `Fig_S11_clusters` is Fig. S12).

## What needs the raw counts

These analyses cannot be rerun from the derived data. Their results are provided in `data/precomputed/` so that the figures that use them can be drawn and the numbers checked.

| Analysis | Why it needs the counts | Provided result |
|---|---|---|
| Trait null model (Fig. 3c) | It shuffles traits among genera and recomputes the community-weighted means | `2_null_model.csv`, `2_null_distributions.rds` |
| Leave one genus out (Table S5, Fig. S9) | It removes one genus from the community before computing the descriptors | `7_leave_one_genus_out.csv` |
| Known-truth tests with planted shifts (Fig. 5, Table S9) | It mixes genus compositions of different dates | `14_planted_shift_summary.csv` |
| Sensitivity to cell and colony volume (Fig. S8, Table S11) | It perturbs the volume of each genus | `15_volume_perturbation.csv`, `15_colony_volume_results.csv` |
| Genus shares and densities (Results, Fig. S6) | They are the counts | `12_post_stretch_biovolume_shares.csv` (mean shares in the post-shift stretch) |

Environmental variables (Figs. 2, S2, S7 and S11, carbonate calculations) are part of the monitoring data and are not included either. The inflow series in Fig. 2 belongs to the provincial water authority (APRHI) and can be requested from it.

## Data

| File | Content |
|---|---|
| `data/samples.csv` | One row per sample: sample number, date, sampling unit and location (names as in Fig. 1 of the paper), and the Berger-Parker dominance by cells and by biovolume (share of the most abundant genus) |
| `data/unit_levels.txt` | Order of the sampling units in the original analysis (see below) |
| `data/distances/*.csv` | Dissimilarity between every pair of the 183 samples, one matrix per descriptor; row and column names are sample numbers |
| `data/cwm_cells.csv`, `data/cwm_biovolume.csv` | Community-weighted means of the traits, with cell and biovolume weights |
| `data/traits/` | Genus trait table, trait types, mixotrophy coding and trait sources (Tables S2 and S3) |
| `data/hydrology/` | Daily water level 2011–2025 compiled from the reports of the Government of Córdoba, and the 2003 level–volume curve |
| `data/precomputed/` | Results of the analyses that need the counts |

The descriptors are:

| Name | Description |
|---|---|
| `traits_cells` | Gower distance on cell-weighted community trait means (Hellinger weights) |
| `traits_biovolume` | Same, biovolume-weighted |
| `traits_cells_dom` | Cell-weighted, with binary traits coded as presence in the dominant fraction, as in the first submission |
| `traits_cells_klais`, `traits_cells_rimet` | Cell-weighted, with two alternative codings of mixotrophy |
| `taxa_genus` | Euclidean distance on Hellinger-transformed cell densities of genera |
| `taxa_genus_biovol` | Same, on biovolumes |
| `taxa_class` | Same, on cell densities aggregated to class |

Sampling units: inlets of the San Pedro (TSP) and Los Espinillos (TLE) rivers, confluence of Los Reartes and Del Medio (TLR), drinking-water intakes on the San Pedro (WISP) and Los Espinillos (WILE) arms, centre of the reservoir (CLM1 below the surface, CLM2 at the lower limit of the euphotic zone) and dam (DLM1, DLM2). WISP and WILE were sampled from February 2018 onward; the other five locations span the whole record.

## Tests that permute whole dates (revision of 27 November 2026)

The samples of one date share variation, so every permutation test moves all the samples of a date together. An earlier version tested the step versus trend, the baseline PERMANOVA and the dispersion tests by permuting single samples (within sampling units for the step), which treats samples of the same date as independent and gave p-values that were far too small. Those tests now permute whole dates (`step_test()`, `period_test()` and `dispersion_test()` in `R/functions.R`) and run with their own seeds, so every other result is unchanged. The planted-shift and colony-volume results in `data/precomputed/` were recomputed from the counts with the same test.

## Exact reproduction

Every permutation test uses the seed and the order of calculation of the original analysis, so the p-values are reproduced exactly, not only within Monte Carlo error. Three details make this possible:

- The trait null model, which needs the counts, drew random numbers before the change-point search. `01_main_analysis.R` advances the random stream by the same number of draws.
- The order of the sampling units and of the five long-term locations is that of the original analysis, because random draws of locations depend on it.
- `01_main_analysis.R` also runs a preliminary thinning analysis that is not reported in the paper, because the spatial-heterogeneity tests that follow it depend on its random draws.

The outputs of all six scripts were checked against those of the original analysis: every table is identical, and the F profiles differ only in the 14th decimal because the distance matrices are stored as text. Verified with R 4.3.3, vegan 2.6-4 and permute 0.9-7. Other versions should give the same estimates; p-values could differ within Monte Carlo error if the permutation algorithms change.

## License

The code is released under the MIT licence (`LICENSE`). The derived data in `data/` and the results in `output/` are released under CC BY 4.0 (https://creativecommons.org/licenses/by/4.0/); please cite the paper when you use them.

## Contact

Guido Marconi, Universidad Nacional de Córdoba [[VERIFY: e-mail and ORCID]]
