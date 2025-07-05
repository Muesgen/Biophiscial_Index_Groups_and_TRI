# ─────────────────────────────────────────────────────────────────────────────
# 0. Libraries (only what we actually use)
# ─────────────────────────────────────────────────────────────────────────────
library(dplyr)
library(tidyr)
library(lubridate)
library(data.table)
library(readr)
library(forecast)
library(ggplot2)
library(purrr)     # for walk()

# ─────────────────────────────────────────────────────────────────────────────
# 1. Load reflectance + PAR  →  SDC_raw
# ─────────────────────────────────────────────────────────────────────────────
band_cols <- c("blue","green","red","nir","swir1","swir2")

SDC_raw <- read_csv("data/analysis_ready_data/SDC_extracted.csv") |>
  rename(plot = 1) |>
  mutate(date = as.Date(date),
         year = year(date),
         doy  = yday(date))

par_tbl <- read_csv("data/analysis_ready_data/daily_PAR_mean_by_plot.csv") |>
  mutate(date = as.Date(date))

setDT(par_tbl)[ , par_mean :=
                  frollmean(par_mean, n = 16, align = "right",
                            na.rm = TRUE), by = plot]
setkey(setDT(par_tbl), plot, date)
setDT(SDC_raw)[ , date := as.IDate(date)]

SDC_raw <- par_tbl[SDC_raw, roll = "nearest"] |>
  as.data.frame() |>
  mutate(date = as.Date(date))

# ─────────────────────────────────────────────────────────────────────────────
# 2. Clean the six bands  →  SDC_band
# ─────────────────────────────────────────────────────────────────────────────
SDC_band <- SDC_raw |>
  arrange(plot, date) |>
  group_by(plot) |>
  mutate(across(all_of(band_cols), ~ {
    tsx <- ts(.x, frequency = 46)
    as.numeric(tsclean(tsx, iterate = 2,
                       replace.missing = FALSE,
                       lambda = NULL))
  })) |>
  ungroup()

# ─────────────────────────────────────────────────────────────────────────────
# 3. Compute VIs from cleaned bands  →  SDC_vi
# ─────────────────────────────────────────────────────────────────────────────
SDC_vi <- SDC_band |>
  mutate(
    NDVI  = (nir - red) / (nir + red),
    GNDVI = (nir - green) / (nir + red),
    EVI   = 2.5 * (nir - red) / (nir + 6*red - 7.5*blue + 1),
    GRVI  = (green - red) / (green + red),
    NIRv  = NDVI * nir,
    NIRvP = NIRv * par_mean,
    kNDVI = (1 - exp(-(nir - red)^2 / (2 * 0.2^2))) /
      (1 + exp(-(nir - red)^2 / (2 * 0.2^2))),
    SIPI  = ((nir - blue) / (nir - red)) * -1,
    NMDI  = (nir - (swir1 - swir2)) / (nir + (swir1 + swir2)),
    GVMI  = ((nir + 0.10) - (swir2 + 0.02)) /
      ((nir + 0.10) + (swir2 + 0.02)),
    CIG   = (nir / green) - 1,
    NDWI  = (nir - swir1) / (nir + swir1)
  )

# ─────────────────────────────────────────────────────────────────────────────
# 4. Clean each VI again  →  SDC_vi_clean
# ─────────────────────────────────────────────────────────────────────────────
vi_cols <- c("NDVI","GNDVI","EVI","GRVI",
             "NIRv","NIRvP","kNDVI","SIPI",
             "NMDI","GVMI","CIG","NDWI")

SDC_vi_clean <- SDC_vi |>
  arrange(plot, date) |>
  group_by(plot) |>
  mutate(across(all_of(vi_cols), ~ {
    tsx <- ts(.x, frequency = 46)
    as.numeric(tsclean(tsx, iterate = 2,
                       replace.missing = FALSE,
                       lambda = NULL))
  })) |>
  ungroup()

# ─────────────────────────────────────────────────────────────────────────────
# 5A.  Reflectance: raw vs clean audit + overlay
# ─────────────────────────────────────────────────────────────────────────────
compare_band <- full_join(
  pivot_longer(SDC_raw,  band_cols, names_to="band", values_to="raw"),
  pivot_longer(SDC_band, band_cols, names_to="band", values_to="clean"),
  by = c("plot","date","band")
) |>
  mutate(changed = raw != clean & !is.na(raw))

band_stats <- compare_band |>
  count(plot, band, changed) |>
  pivot_wider(names_from=changed, values_from=n, values_fill = 0) |>
  mutate(pct_changed = TRUE*100/(TRUE+FALSE))

print(band_stats)

dir.create("figures/compare_raw_clean_reflectance", showWarnings = FALSE, recursive = TRUE)

walk(unique(compare_band$plot), function(pl) {
  ggplot(dplyr::filter(compare_band, plot == pl), aes(date)) +
    geom_line(aes(y=raw),   colour="#B2B2B2", linewidth=.3) +
    geom_line(aes(y=clean), colour="#D95F02", linewidth=.3) +
    facet_wrap(~band, ncol=1, scales="free_y") +
    labs(title=sprintf("Reflectance raw vs clean — plot %s", pl),
         subtitle="Grey = original  ·  Orange = after tsclean()") +
    theme_minimal(9) +
    theme(strip.text = element_text(face="bold")) ->
    p
  ggsave(sprintf("figures/compare_raw_clean_reflectance/%s_raw_vs_clean.png", pl),
         p, width=11, height=8, units="in", bg="white")
})

# ─────────────────────────────────────────────────────────────────────────────
# 5B.  VIs: pre- vs post-clean audit + overlay
# ─────────────────────────────────────────────────────────────────────────────
compare_vi <- full_join(
  pivot_longer(SDC_vi,       vi_cols, names_to="VI", values_to="pre"),
  pivot_longer(SDC_vi_clean, vi_cols, names_to="VI", values_to="post"),
  by = c("plot","date","VI")
) |>
  mutate(changed = pre != post & !is.na(pre))

vi_stats <- compare_vi |>
  count(plot, VI, changed) |>
  pivot_wider(names_from=changed, values_from=n, values_fill = 0) |>
  mutate(pct_changed = TRUE*100/(TRUE+FALSE))

print(vi_stats)

dir.create("figures/compare_raw_clean_VI", showWarnings = FALSE, recursive = TRUE)

walk(unique(compare_vi$plot), function(pl) {
  ggplot(dplyr::filter(compare_vi, plot == pl), aes(date)) +
    geom_line(aes(y=pre),  colour="#B2B2B2", linewidth=.3) +
    geom_line(aes(y=post), colour="#D95F02", linewidth=.3) +
    facet_wrap(~VI, ncol=1, scales="free_y") +
    labs(title=sprintf("Vegetation indices raw vs clean — plot %s", pl),
         subtitle="Grey = from cleaned bands · Orange = extra VI-level clean") +
    theme_minimal(9) +
    theme(strip.text = element_text(face="bold")) ->
    p
  ggsave(sprintf("figures/compare_raw_clean_VI/%s_VI_raw_vs_clean.png", pl),
         p, width=11, height=18, units="in", bg="white")
})
