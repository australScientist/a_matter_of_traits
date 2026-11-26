# Semantic colour palette shared by all figures of the paper ------------------------------------------
#
# Rule: each colour means one thing across the paper.
#   - Blue and orange are RESERVED for the period: before and after the shift. Mid grey is the
#     transition. No other element uses blue or orange.
#   - Violet is the MAGNITUDE scale (similarity, detection rate): from near white to dark violet.
#   - Cyan is the REFILLING of the reservoir (December 2017 to February 2018).
#   - Descriptors are coded by family (traits black, genera magenta, classes gold) and by weighting
#     (cells solid line, biovolume dashed line).
# Categorical colours come from a palette validated for colour-vision deficiency.

ink      <- "#0b0b0b"
ink2     <- "#52514e"
grid_col <- "#e6e5e1"

# Period
phase_cols   <- c(pre = "#2a78d6", transition = "#8a8984", post = "#eb6834")
phase_light  <- c(pre = "#cde2fb", transition = "#e3e2de", post = "#fbdccd")  # light fills
phase_labels <- c(pre = "Before (to Feb 2018)", transition = "Transition", post = "After (from Aug 2018)")

# Time ramps within each period (Fig. 4d-f): hue gives the period, lightness the date (light = early)
phase_ramp <- list(pre = c("#a9cbf3", "#174a8b"), post = c("#f7c3a8", "#9c3a12"))

# Magnitude (similarity, detection rate)
mag_low  <- "#f6f4fb"
mag_high <- "#3b2a8f"
mag_ramp <- function(n) grDevices::colorRampPalette(c("#cfc8ef", mag_high))(n)

# Refilling
flood_fill <- "#d4eef2"
flood_ink  <- "#0e7c8c"

# Descriptors: family by colour, weighting by line type
desc_tbl <- tibble::tribble(
  ~descriptor,         ~label,                        ~colour,   ~linetype,
  "traits_cells",      "Traits, cell-weighted",       "#0b0b0b", "solid",
  "traits_biovolume",  "Traits, biovolume-weighted",  "#0b0b0b", "22",
  "taxa_genus",        "Genera, cell-weighted",       "#c2407f", "solid",
  "taxa_genus_biovol", "Genera, biovolume-weighted",  "#c2407f", "22",
  "taxa_class",        "Classes, cell-weighted",      "#b07d00", "solid")
desc_cols  <- setNames(desc_tbl$colour,   desc_tbl$label)
desc_lty   <- setNames(desc_tbl$linetype, desc_tbl$label)
weight_lty <- c(Cells = "solid", Biovolume = "22")

# Genera (no blue or orange). "Other" in light grey.
genus_palette <- c("#7a5c00", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#6250d6", "#e34948", "#c9c8c2")
