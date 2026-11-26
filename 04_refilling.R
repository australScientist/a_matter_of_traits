# Size of the refilling of December 2017 to February 2018 compared with the refills of 2012, 2014 and 2021
# (Assessment, "Environmental context"). Uses public data only.
# Run from the repository root:  Rscript 04_refilling.R   (a few seconds)
#
# Daily water levels (gauge reading; spillway crest at 53.00 m) were compiled from the reports of the
# Government of Córdoba (https://www.cba.gov.ar/nivel-de-diques-y-embalses/). On every date the
# median of all reports was taken; range_between_sources_m gives the spread between reports.
# Volumes come from the 2003 level-volume curve of the reservoir (Catalini and Caamaño Nelli,
# Instituto Nacional del Agua), interpolated with a monotone spline. Later sedimentation would make the
# percentages slightly larger without changing the comparison between years.
# The inflow series used in Fig. 2 belongs to the provincial water authority (APRHI) and is not included.

suppressPackageStartupMessages(library(tidyverse))
dir.create("output", showWarnings = FALSE)

level <- read_csv("data/hydrology/los_molinos_level_daily.csv", show_col_types = FALSE)
curve <- read_csv("data/hydrology/level_volume_curve_2003.csv", show_col_types = FALSE)
volume <- splinefun(curve$level_m, curve$volume_hm3, method = "monoH.FC")

# Lowest level in a window, then the highest level after it
refill <- function(from, to) {
  x  <- level %>% filter(date >= as.Date(from), date <= as.Date(to))
  lo <- x %>% slice_min(level_m, n = 1, with_ties = FALSE)
  hi <- x %>% filter(date > lo$date) %>% slice_max(level_m, n = 1, with_ties = FALSE)
  tibble(from = lo$date, level_from = lo$level_m, to = hi$date, level_to = hi$level_m,
         rise_m = hi$level_m - lo$level_m, volume_from = volume(lo$level_m), volume_to = volume(hi$level_m),
         volume_increase_pct = 100 * (volume(hi$level_m) / volume(lo$level_m) - 1))
}
refills <- bind_rows(`2011-12` = refill("2011-09-01", "2012-06-30"), `2013-14` = refill("2013-09-01", "2014-06-30"),
                     `2017-18` = refill("2017-09-01", "2018-06-30"), `2020-21` = refill("2020-09-01", "2021-06-30"),
                     .id = "refill") %>% mutate(across(where(is.double), ~ round(.x, 2)))
write_csv(refills, "output/refills.csv")
print(refills)
