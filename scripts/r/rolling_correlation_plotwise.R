################################################################################
# Rolling-Window Vegetation-Index (VI) vs. Tree-Ring Index (TRI) Analysis
# ─────────────────────────────────────────────────────────────────────────────
# This script
#   1. extracts Landsat surface-reflectance pixels (10 m buffer) for each plot
#   2. computes daily vegetation indices (VIs)
#   3. builds one detrended tree-ring chronology (TRI) per plot
#   4. rolls 1- to 24-day cumulative sums of every VI
#   5. finds, for every plot × VI, the (window, DOY) pair with the strongest
#      absolute Pearson correlation to TRI
#   6. draws per-plot time-series grids, global summary tables, and several
#      distribution / heat-map figures
################################################################################

# ─────────────────────────────────────────────────────────────────────────────
# 0. LIBRARIES & HELPER FUNCTIONS
# ─────────────────────────────────────────────────────────────────────────────
library(dplyr)      # data wrangling
library(tidyr)      # pivot_* helpers
library(lubridate)  # date-time helpers (year(), yday(), …)
library(zoo)        # rollapply() for rolling sums
library(purrr)      # map_dfc(), walk()
library(ggplot2)    # all plotting
library(dplR)
library(terra)
library(sf)
library(zoo)        # for rolling means
library(signal)     # for Savitzky-Golay smoothing
library(readxl)
library(tibble)
library(RColorBrewer) 
library(data.table)
library(fuzzyjoin)

# custom extractors (paths relative to your repo)
source("scripts/r/functions/extract_TRW.R")          # tree-ring widths
source("scripts/r/functions/extract_SDC_buffered.R") # Landsat reflectance
source("scripts/r/functions/extract_spi_by_plot.R")
source("scripts/r/functions/extract_predictors.R")

# ─────────────────────────────────────────────────────────────────────────────
# 1. SATELLITE DATA  – daily VIs per plot
# ─────────────────────────────────────────────────────────────────────────────
## 1A  raw reflectance pixels --------------------------------------------------
SDC_df <- read.csv("data/analysis_ready_data/SDC_extracted.csv")
names(SDC_df)[1] <- "plot"                          # first col → plot ID

## 1B  add calendar fields -----------------------------------------------------
SDC_df <- SDC_df |>
  mutate(year = year(date),           # calendar year
         doy  = yday(date),
         date = as.Date(date))          # ② same for the SDC table)           # day-of-year (1–365/366)

par_tbl <- read.csv("data/analysis_ready_data/daily_PAR_mean_by_plot.csv") |>
  mutate(date = as.Date(date))          # ① guarantee Date class


setDT(SDC_df)[ , date := as.IDate(date)]     # IDate = fast Date
setDT(par_tbl)[, date := as.IDate(date)]

par_tbl[ , par_mean :=
           frollmean(par_mean,             # data.table’s fast roll
                     n = 16,               # 8-day window
                     align = "right",
                     na.rm = TRUE),
         by = plot]

setkey(par_tbl, plot, date)                  # key on both cols
SDC_raw <- par_tbl[SDC_df, roll = "nearest"]   # keeps SDC order


## 1C  daily mean vegetation indices ------------------------------------------

library(forecast)

# make sure every table uses Date (not IDate or POSIXct)
SDC_raw  <- SDC_raw  %>% mutate(date = as.Date(date))
band_cols <- c("blue","green","red","nir","swir1","swir2")

SDC_clean <- SDC_raw %>% 
  arrange(plot, date) %>%          # keep each plot's dates in order
  group_by(plot) %>% 
  mutate(across(
    all_of(band_cols),
    ~ {
      tsx <- ts(.x, frequency = 46)          # ≈ 8-day composites
      as.numeric(
        tsclean(tsx,
                iterate         = 2,
                replace.missing = FALSE)
      )
    }
  )) %>% 
  ungroup()



# ── 2. compute indices for each row (no grouping needed) ------------------
VI_raw <- SDC_clean %>%
  mutate(
    NDVI    = (nir - red) / (nir + red),
    GNDVI   = (nir - green) / (nir + red),
    EVI     = 2.5 * (nir - red) / (nir + 6 * red - 7.5 * blue + 1),
    GRVI    = (green - red) / (green + red),
    #OSAVI   = (nir - red) / (nir + red + 0.16),
    NIRv    = NDVI * nir,
    NIRvP = NIRv * par_mean,
    kNDVI   = (1 - exp(-(nir - red)^2 / (2 * 0.2^2))) / (1 + exp(-(nir - red)^2 / (2 * 0.2^2))),
    #VMI     = (nir - swir1) * (nir - swir2) / (nir + swir1 + swir2),
    SIPI    = ((nir - blue) / (nir - red)) * -1,
    NMDI    = (nir - (swir1 -swir2))/(nir+(swir1+swir2)),
    GVMI    = ((nir + 0.10) - (swir2 + 0.02)) /((nir + 0.10) + (swir2 + 0.02)),
    CIG     = (nir / green) - 1,
    #MSI     = nir / swir1,
    NDWI    = (nir - swir1) / (nir + swir1),
    year    = year(date),
    doy     = yday(date)
  ) %>%
  select( -par_mean, - blue, -green, -nir, -red, -swir1, -swir2)          # drop helper column if you like

# make sure every table uses Date (not IDate or POSIXct)
VI_raw  <- VI_raw  %>% mutate(date = as.Date(date))
index_cols <- names(VI_raw)[seq(6,17,1)]

VI_clean <- VI_raw %>% 
  arrange(plot, date) %>%          # keep each plot's dates in order
  group_by(plot) %>% 
  mutate(across(
    all_of(index_cols),
    ~ {
      tsx <- ts(.x, frequency = 46)          # ≈ 8-day composites
      as.numeric(
        tsclean(tsx,
                iterate         = 2,
                replace.missing = FALSE)
      )
    }
  )) %>% 
  ungroup()

# ─────────────────────────────────────────────────────────────────────────────
# 2. TREE-RING CHRONOLOGIES  – detrended TRI per plot
# ─────────────────────────────────────────────────────────────────────────────
chronologies_by_plot <- read.csv("data/analysis_ready_data/TRI_chronologies.csv") %>%
  dplyr::filter(year < 2023)

# ─────────────────────────────────────────────────────────────────────────────
# 3. ROLLING-SUM VIs  (1–24-day windows, daily resolution)
# ─────────────────────────────────────────────────────────────────────────────
accum_windows <- 1:91                          # candidate window lengths

## 3A  long table: one VI per row ------------------------------------------------
SDC_long <- VI_clean |>
  pivot_longer(
    -c(plot, date, year, doy, area),                  # keep date + calendar cols
    names_to  = "VI", values_to = "value"
  )

## 3B  rolling sums per window ---------------------------------------------------
SDC_cum_continuous <- SDC_long |>
  group_by(plot, VI) |>
  group_modify(~{
    bind_cols(
      .x,
      map_dfc(accum_windows, \(w)
              tibble(!!paste0("cum", w) :=
                       zoo::rollapply(.x$value, w, sum, align = "right", fill = NA)))
    )
  }) |>
  ungroup()

## 3C  tidy: one row = one window length -----------------------------------------
SDC_cum_long <- SDC_cum_continuous |>
  pivot_longer(starts_with("cum"),
               names_to  = "window", names_prefix = "cum",
               values_to = "cum_value") |>
  mutate(window = as.integer(window))

# ─────────────────────────────────────────────────────────────────────────────
# 4. JOIN TRI  +  tag forest type (coniferous / deciduous)
# ─────────────────────────────────────────────────────────────────────────────
cum_TRI_filtered <- SDC_cum_long |>
  left_join(chronologies_by_plot, by = c("plot","year")) |>
  dplyr::filter(year >= 2000) |>
  mutate(forest_type = case_when(
    plot %in% c("PF01","SF03","SF06","SF07","SF10", "SF11",
                "SF13","SF15","SF16", "SF21") ~ "coniferous",
    plot %in% c("PF02","PF03","SF01", "SF02","SF04","SF05","SF08",
                "SF09","SF12", "SF14","SF17", "SF18","SF19","SF20") ~ "deciduous",
    TRUE                                                  ~ NA_character_
  ))
# ─────────────────────────────────────────────────────────────────────────────
# 5. CORRELATION HEAT-TABLE  &  best (window, DOY)
# ─────────────────────────────────────────────────────────────────────────────
heatmap_df <- cum_TRI_filtered |>
  group_by(plot, forest_type, VI, window, doy) |>
  summarise(
    correlation = if (sum(!is.na(cum_value) & !is.na(TRI)) > 1)
      cor(cum_value, TRI, use = "complete.obs")
    else NA_real_,
    .groups = "drop"
  )

heat_curr <- heatmap_df %>% 
  mutate(lag = 0)

# ——— 2. build the lag‐1 correlations ——————————————————————
# Start from your cum‐values + TRI joined table:
heat_prev <- cum_TRI_filtered %>%
  # pull out only what we need, and rename:
  select(
    plot, forest_type, VI, window, doy,
    cum_prev = cum_value,
    year_prev = year
  ) %>%
  # shift the year so that cum_prev comes from Y−1
  mutate(year = year_prev + 1) %>%
  # bring in *only* the current‐year TRI as TRI_current
  inner_join(
    chronologies_by_plot %>% 
      select(plot, year, TRI_current = TRI),
    by = c("plot", "year")
  ) %>%
  # now correlate the lagged cum_prev vs. TRI_current
  group_by(plot, forest_type, VI, window, doy) %>%
  summarise(
    correlation = if (sum(!is.na(cum_prev) & !is.na(TRI_current)) > 1)
      cor(cum_prev, TRI_current, use = "complete.obs")
    else NA_real_,
    .groups = "drop"
  ) %>%
  mutate(lag = 1)

# 3. Bind them and plot as before:
heat_all <- bind_rows(heat_curr, heat_prev)

saveRDS(heat_all, "corr_table_vi_tri.rds")

best_combo <- heat_all |>
  dplyr::filter(!is.na(correlation)) |>           # keep only positive r
  group_by(plot, VI) |>
  slice_max(correlation,               # largest positive r
            n = 1, with_ties = FALSE) |>
  ungroup()


# ─────────────────────────────────────────────────────────────────────────────
# 6. TIME-SERIES TABLE for plotting  (only best windows)
# ─────────────────────────────────────────────────────────────────────────────
plot_tbl <- SDC_cum_long |>
  inner_join(best_combo |>
               select(plot, VI, win_days = window, doy),
             by = c("plot","VI","window"="win_days","doy")) |>
  left_join(select(chronologies_by_plot, plot, year, TRI),
            by = c("plot","year")) |>
  rename(VI_val = cum_value) |>
  select(-value) |>
  pivot_longer(c(VI_val, TRI), names_to = "series", values_to = "value") |>
  group_by(plot, VI, series) |>
  mutate(value = scale(value)[,1]) |>
  ungroup()

plot_meta <- read_xlsx("data/vector_data/plot_metadata.xlsx")
plot_tbl <- left_join(plot_tbl, plot_meta, by = "plot")       # ▶︎ B

# ─────────────────────────────────────────────────────────────────────────────
# 7. PER-PLOT TIME-SERIES GRID  (one PNG per plot)
# ─────────────────────────────────────────────────────────────────────────────
dir.create("figures/VI_TRI_by_plot", showWarnings = FALSE)

walk(unique(plot_tbl$plot), \(pl){
  
  plt   <- dplyr::filter(plot_tbl, plot == pl)
  stats <- dplyr::filter(best_combo,   plot == pl)
  overall <- mean(abs(stats$correlation), na.rm = TRUE)
  
  ## ▲  meta for this plot
  species <- unique(plt$species)
  area    <- unique(plt$area)
  
  ## facet labels (unchanged)
  labels <- transmute(stats,
                      VI,
                      facet_lab = sprintf("%s\nr = %.2f | win = %dd | DOY = %d",
                                          VI, correlation, window, doy))
  plt <- left_join(plt, labels, by = "VI")
  
  ## plot
  p <- ggplot(plt, aes(year, value, colour = series)) +
    geom_line(linewidth = .4, na.rm = TRUE) +
    facet_wrap(~facet_lab, ncol = 4, scales = "free_y",
               strip.position = "top") +
    scale_colour_manual(values = c("VI_val" = "steelblue",
                                   "TRI"    = "black"),
                        breaks = c("VI_val","TRI"),
                        labels = c("VI","TRI")) +
    labs(title = sprintf("Best window & DOY — plot %s  (%s, %s)  |mean |r| = %.2f",
                         pl, species, area, overall),
         x = NULL, y = "z-score", colour = NULL) +
    theme(legend.position = "top",
          axis.text.x     = element_text(angle = 90,
                                         vjust = .5, hjust = 1),
          panel.spacing   = unit(3, "mm"),
          strip.placement = "outside")
  
  ## ▲  filename starts with species
  ggsave(file.path("figures/VI_TRI_by_plot",
                   sprintf("%s_%s_VI_TRI_best_windows.png", species, pl)),
         plot = p, width = 11, height = 8, dpi = 300, bg = "white")
})

dir.create("figures/VI_TRI_by_plot", showWarnings = FALSE)

walk(unique(plot_tbl$plot), \(pl){
  
  ## ── meta & helpers ─────────────────────────────────────────────
  stats_pl  <- dplyr::filter(best_combo, plot == pl)
  
  ## **1. pick the single “best” VI (|r| max) for this plot**
  best_stat <- slice_max(stats_pl, abs(correlation), n = 1, with_ties = FALSE)
  best_vi   <- best_stat$VI
  best_r    <- best_stat$correlation
  
  species <- unique(plot_tbl$species[plot_tbl$plot == pl])
  area    <- unique(plot_tbl$area   [plot_tbl$plot == pl])
  
  ## ── data for plotting ──────────────────────────────────────────
  ## **2. keep TRI + the best VI only**
  plt <- dplyr::filter(plot_tbl, plot == pl,
                series == "TRI" | VI == best_vi)
  
  ## single-facet label
  labels <- best_stat |> 
    transmute(VI,
              facet_lab = sprintf("%s\nr = %.2f | win = %dd | DOY = %d",
                                  VI, correlation, window, doy))
  plt <- left_join(plt, labels, by = "VI")
  
  ## ── plot ───────────────────────────────────────────────────────
  p <- ggplot(plt, aes(year, value, colour = series)) +
    geom_line(linewidth = .4, na.rm = TRUE) +
    scale_colour_manual(values = c("VI_val" = "steelblue",
                                   "TRI"    = "black"),
                        breaks = c("VI_val","TRI"),
                        labels = c(best_vi,"TRI")) +
    labs(title = sprintf("Best window %s & DOY %s — plot %s  (%s, %s)  | r = %.2f",
                         best_stat$window, best_stat$doy, pl, species, area, best_r),
         x = NULL, y = "z-score", colour = NULL) +
    theme(legend.position = "top",
          axis.text.x     = element_text(angle = 90, vjust = .5, hjust = 1),
          panel.spacing   = unit(3, "mm"),
          strip.placement = "outside")
  
  ## ── output ─────────────────────────────────────────────────────
  ## **3. include the best VI name in the file name**
  ggsave(file.path("figures/VI_TRI_by_plot",
                   sprintf("%s_%s_%s_TRI_best_window.png",
                           species, pl, best_vi)),
         plot = p, width = 11, height = 6, dpi = 300, bg = "white")
})


# ─────────────────────────────────────────────────────────────────────────────
# 8. SUMMARY TABLE  (mean |r|, window, DOY per VI × forest type)
# ─────────────────────────────────────────────────────────────────────────────
summary_tbl <- best_combo |>
  group_by(VI) |>
  summarise(
    plots_used  = n(),
    mean_abs_r  = mean(abs(correlation)),
    mean_window = mean(window), sd_window = sd(window),
    mean_doy    = mean(doy),    sd_doy    = sd(doy),
    .groups     = "drop") |>
  arrange(desc(mean_abs_r))
print(summary_tbl)

## --- 10B  output folder ------------------------------------------------------
dir.create("figures/VI_TRI_heatmaps_by_plot", showWarnings = FALSE)
heatmap_df <- dplyr::left_join(heatmap_df, plot_meta, by = "plot")

## --- 10C  loop over plots ----------------------------------------------------
for (pl in unique(heatmap_df$plot)) {
  
  df <- heatmap_df %>% 
    dplyr::filter(plot == pl, !is.na(correlation)) %>% 
    mutate(
      month_date   = as.Date(doy - 1, origin = "2000-01-01"),
      window_month = window * 8 / 30.437
    )
  if (nrow(df) == 0) next
  
  rng    <- range(df$correlation, na.rm = TRUE)          # data limits
  breaks <- seq(floor(rng[1] / 0.1) * 0.1,   # round down to nearest 0.1
                ceiling(rng[2] / 0.1) * 0.1, # round up   to nearest 0.1
                by = 0.1)                    # ← constant 0.1 gap
  
  # q_breaks <- quantile(df$correlation, probs = seq(0, 1, length.out = 9), na.rm = TRUE)
  
  # how many bands (= how many colours do we need)?
  n_bands <- length(breaks) - 1         
  
  # build a long Spectral palette
  spectral_long <- colorRampPalette(
    rev(RColorBrewer::brewer.pal(11, "Spectral"))
  )(n_bands)
  
  p <- ggplot(df, aes(month_date, window_month, z = correlation)) +
    geom_contour_filled(
      colour = "grey", size = 0.2, na.rm = TRUE, breaks = breaks
    ) +
    facet_wrap(~VI, ncol = 4) +
    scale_x_date(date_breaks = "1 month", date_labels = "%b", expand = c(0, 0)) +
    scale_y_continuous(
      breaks = seq(0, ceiling(max(df$window_month, na.rm = TRUE)), 1),
      expand = c(0, 0)
    ) +
    labs(
      title = sprintf("VI–TRI rolling correlation — plot %s  (%s, %s)",
                      pl, unique(df$species), unique(df$area)),
      x = "Month", y = "Window length (months)"
    ) +
    theme(strip.text = element_text(face = "bold")) +
    scale_fill_manual(values = spectral_long)   # discrete scale
  
  ggsave(
    filename = file.path("figures/VI_TRI_heatmaps_by_plot",
                         sprintf("%s_%s_VI_TRI_contour_%s.png",
                                 unique(df$species), unique(df$area), pl)),
    plot   = p,
    width  = 14, height = 12, dpi = 300, bg = "white"
  )
}

vi_order <- c("NDVI", "EVI",  "kNDVI", "NIRv", "NIRvP",
              "NMDI", "NDWI", "GVMI", "SIPI", "GRVI", "CIG", "GNDVI")

library(ggrepel)
dir.create("figures/VI_TRI_heatmaps_by_plot", showWarnings = FALSE)
heat_all <- dplyr::left_join(heat_all, plot_meta, by = "plot")
for (pl in unique(heat_all$plot)) {
  
  df <- heat_all %>%
    dplyr::filter(plot == pl, !is.na(correlation)) %>%
    mutate(
      VI          = factor(VI, levels = vi_order),   #
      ext_doy      = if_else(lag == 1, doy, doy + 365),
      month_date   = as.Date(ext_doy - 1, origin = "2000-01-01"),
      window_month = window * 8 / 30.437
    )
  if (nrow(df) == 0) next
  
  ## ── strongest positive r per VI ───────────────────────────────
  best_pts <- df %>%                               # ← NEW
    dplyr::filter(correlation > 0) %>%                    # keep positives only
    group_by(VI) %>%                               # one row per facet
    slice_max(correlation, n = 1, with_ties = FALSE) %>% 
    ungroup()                                      # ← NEW
  
  ## colour-band setup (unchanged) --------------------------------
  rng     <- range(df$correlation, na.rm = TRUE)
  breaks  <- seq(floor(rng[1]/0.1)*0.1,
                 ceiling(rng[2]/0.1)*0.1,
                 by = 0.1)
  n_bands <- length(breaks) - 1
  spectral_long <- colorRampPalette(
    rev(RColorBrewer::brewer.pal(11, "RdYlBu"))
  )(n_bands)
  
  ## ── plot ──────────────────────────────────────────────────────
  p <- ggplot(df, aes(month_date, window_month, z = correlation)) +
    geom_contour_filled(colour = "gray50", size = 0.2,
                        na.rm = TRUE, breaks = breaks) +
    facet_wrap(~VI, ncol = 4) +
    scale_x_date(date_breaks = "2 month",
                 date_labels = "%b",
                 expand = c(0, 0)) +
    scale_y_continuous(breaks = seq(0, ceiling(max(df$window_month, na.rm = TRUE)), 2),
                       expand = c(0, 0)) +
    
    ## — highlight the max-r point — -------------------------------
  geom_point(                                    # ← NEW
    data   = best_pts,
    aes(month_date, window_month),
    size   = 2,
    shape  = 21,
    stroke = 0.8,
    fill   = "yellow",
    colour = "black"
    ) +
    geom_label_repel(
      data   = best_pts,
      aes(month_date, window_month,
          label = sprintf("r = %.2f", correlation)),
      size        = 3,
      fill        = alpha("white", 0.3),   # 60 % opacity white
      colour      = "black",               # text colour
      label.size  = 0,                     # no border line
      label.r     = unit(0.1, "lines"),    # corner radius
      box.padding   = 0.3,
      point.padding = 0.2,
      segment.color = NA                   # no leader line
    )+
    
    labs(title = sprintf("VI–TRI rolling correlation — plot %s  (%s, %s)",
                         pl, unique(df$species), unique(df$area)),fill = expression("Pearson "~italic(r)),
         x = "Month (prev year → current year)",
         y = "Window length (months)") +
    theme(strip.text = element_text(face = "bold")) +
    scale_fill_manual(values = spectral_long)
  
  ggsave(
    filename = file.path(
      "figures/VI_TRI_heatmaps_by_plot",
      sprintf("%s_%s_VI_TRI_contour_%s.png",
              unique(df$species), unique(df$area), pl)),
    plot   = p,
    width  = 16, height = 8, dpi = 300, bg = "white"
  )
}


for (sp in unique(heat_all$species)) {
  
  ## ── prepare data ──────────────────────────────────────────────
  df <- heat_all %>%                                             # start with full table
    dplyr::filter(species == sp, !is.na(correlation)) %>%               # keep one species
    mutate(
      ext_doy      = if_else(lag == 1, doy, doy + 365),          # prev-year days →  366-730
      month_date   = as.Date(ext_doy - 1, origin = "2000-01-01"),# dummy calendar on year 2000
      window_month = window * 8 / 30.437                         # your old 8-day step in months
    ) %>%
    ## average the correlations across every plot that feeds this species
    group_by(VI, window_month, month_date) %>%
    summarise(correlation = mean(correlation, na.rm = TRUE), .groups = "drop")
  
  if (nrow(df) == 0) next                                         # safety
  
  ## ── strongest positive r per VI (after averaging) ─────────────
  best_pts <- df %>%
    dplyr::filter(correlation > 0) %>%                                   # keep positives only
    group_by(VI) %>%
    slice_max(correlation, n = 1, with_ties = FALSE) %>%
    ungroup()
  
  ## ── colour-band setup (unchanged) ─────────────────────────────
  rng     <- range(df$correlation, na.rm = TRUE)
  breaks  <- seq(floor(rng[1] / 0.1) * 0.1,
                 ceiling(rng[2] / 0.1) * 0.1,
                 by = 0.1)
  n_bands <- length(breaks) - 1
  spectral_long <- colorRampPalette(
    rev(brewer.pal(11, "Spectral"))
  )(n_bands)
  
  ## ── plot ──────────────────────────────────────────────────────
  p <- ggplot(df, aes(month_date, window_month, z = correlation)) +
    geom_contour_filled(colour = "grey", size = 0.2,
                        na.rm = TRUE, breaks = breaks) +
    facet_wrap(~VI, ncol = 4) +
    scale_x_date(date_breaks = "2 month",
                 date_labels = "%b",
                 expand = c(0, 0)) +
    scale_y_continuous(breaks = seq(0, ceiling(max(df$window_month, na.rm = TRUE)), 2),
                       expand = c(0, 0)) +
    geom_point(                                # highlight best point per VI
      data   = best_pts,
      aes(month_date, window_month),
      size   = 2,
      shape  = 21,
      stroke = 0.8,
      fill   = "yellow",
      colour = "black"
    ) +
    geom_text(
      data  = best_pts,
      aes(month_date, window_month,
          label = sprintf("r = %.2f", correlation)),
      vjust = -0.5,
      size  = 3
    ) +
    labs(title = sprintf("VI–TRI rolling correlation — species %s",
                         sp),
         subtitle = "Averaged across all plots where the species occurs",
         x = "Month (prev year → current year)",
         y = "Window length (months)") +
    theme(strip.text = element_text(face = "bold")) +
    scale_fill_manual(values = spectral_long)
  
  ggsave(
    filename = file.path(
      "figures/VI_TRI_heatmaps_by_species",
      sprintf("%s_VI_TRI_contour.png", sp)),
    plot   = p,
    width  = 16, height = 9, dpi = 300, bg = "white"
  )
}

list.files("figures/VI_TRI_heatmaps_by_plot/")

