# ─────────────────────────────────────────────────────────────
#  POINT-BASED CLIMATE & ELEVATION SUMMARY → AGGREGATED BY PLOT
#  Study area: Five regions in Germany (DE)
#
#  This script:
#    ✓ extracts climate variables at tree sampling points
#    ✓ aggregates monthly climate data to plot level
#    ✓ outputs one row per plot per month with mean climate values
#
#  Input data:
#    - Points: tree sampling points (with plot membership)
#    - Climate rasters: temperature, precipitation, evapotranspiration, soil moisture
#    - Areas shapefile: to add area names to points
#
#  Output:
#    - CSV: mean monthly climate data per plot
#
#  Author:  Konstantin Engelmayer
#  Date:  2025-06-12
# ─────────────────────────────────────────────────────────────

# 1 ── libraries ──────────────────────────────────────────────
library(sf)         # vector data handling
library(terra)      # raster data handling
library(dplyr)      # data manipulation
library(tidyr)      # reshaping data
library(lubridate)  # date handling

# 2 ── paths (edit if needed) ─────────────────────────────────
shp_file        <- "data/vector_data/all_merged.shp"               # polygons (area names)
points_file     <- "data/vector_data/trees_all_plots.gpkg"         # sampling points
raster_dir      <- "data/climate/dwd_4_variables/"

temp_file       <- file.path(raster_dir, "air_temp_mean_1991_2024_merged.tif")   # °C × 10
prec_file       <- file.path(raster_dir, "precipitation_1991_2024_merged.tif")   # mm
evapo_file      <- file.path(raster_dir, "evapo_p_1991_2024_merged.tif")         # mm
soilmoist_file  <- file.path(raster_dir, "soi_moist_1991_2024_merged.tif")       # %

# 3 ── vectors & rasters ──────────────────────────────────────

## 3.1 Read area polygons and assign area names
areas <- st_read(shp_file)
areas$area <- areas$name <- c("Eifel", "Kellerwald", "Koenigsforst",
                              "Calderner Wald", "Lindenberger Wald")

## 3.2 Read sampling points (keep only `plot` column)
points <- st_read(points_file) %>% 
  select(plot)

## 3.3 Load climate rasters
Tmean   <- rast(temp_file)      / 10    # Convert from 0.1 °C to °C
Pmm     <- rast(prec_file)             # Precipitation (mm)
ETp     <- rast(evapo_file)            # Potential evapotranspiration (mm)
SoilMs  <- rast(soilmoist_file)        # Soil moisture (%, or index)

## 3.4 Label raster layers with "YYYY-MM"
start_date  <- ymd("1991-01-01")
layer_names <- format(seq(start_date, by = "1 month", length.out = nlyr(Tmean)), "%Y-%m")

names(Tmean)  <- layer_names
names(Pmm)    <- layer_names
names(ETp)    <- layer_names
names(SoilMs) <- layer_names

# 4 ── harmonise CRS & attach area names to points ────────────
crs_target <- crs(Tmean)

points <- st_transform(points, crs_target) |>
  st_join(st_transform(areas, crs_target), left = FALSE)

# 5 ── extract raster values at points ────────────────────────
pts_v <- vect(points)

temp_df  <- terra::extract(Tmean,   pts_v, ID = FALSE) |> as_tibble()
prec_df  <- terra::extract(Pmm,     pts_v, ID = FALSE) |> as_tibble()
evapo_df <- terra::extract(ETp,     pts_v, ID = FALSE) |> as_tibble()
soil_df  <- terra::extract(SoilMs,  pts_v, ID = FALSE) |> as_tibble()

# 6 ── tidy to long format (one row per point ⋅ month) ────────

## 6.1 Create point metadata table
pt_meta <- points |> 
  st_drop_geometry() |> 
  mutate(pt_id = row_number()) |> 
  select(pt_id, area, plot)

## 6.2 Reshape extracted data to long format

temp_long <- temp_df |> 
  mutate(pt_id = row_number()) |>
  pivot_longer(-pt_id, names_to = "date", values_to = "temp")

prec_long <- prec_df |> 
  mutate(pt_id = row_number()) |>
  pivot_longer(-pt_id, names_to = "date", values_to = "prec")

evapo_long <- evapo_df |> 
  mutate(pt_id = row_number()) |>
  pivot_longer(-pt_id, names_to = "date", values_to = "evapo")

soil_long <- soil_df |> 
  mutate(pt_id = row_number()) |>
  pivot_longer(-pt_id, names_to = "date", values_to = "soil_moist")

# 7 ── merge all variables into a single climate table ────────
clim_long <- temp_long |> 
  left_join(prec_long,  by = c("pt_id", "date")) |>
  left_join(evapo_long, by = c("pt_id", "date")) |>
  left_join(soil_long,  by = c("pt_id", "date")) |>
  mutate(date  = ymd(paste0(date, "-01")),
         month = month(date)) |>
  left_join(pt_meta, by = "pt_id")

# 8 ── aggregate to mean per plot and month ───────────────────
clim_plot_month <- clim_long |> 
  group_by(plot, area, date, month) |> 
  summarise(temp       = mean(temp, na.rm = TRUE),
            prec       = mean(prec, na.rm = TRUE),
            evapo      = mean(evapo, na.rm = TRUE),
            soil_moist = mean(soil_moist, na.rm = TRUE),
            .groups    = "drop")

# 9 ── export final table ─────────────────────────────────────
write.csv(clim_plot_month, "data/analysis_ready_data/climate_df.csv", row.names = FALSE)

# ─────────────────────────────────────────────────────────────
# END OF SCRIPT
# ─────────────────────────────────────────────────────────────

# ---- libraries -----------------------------------------------------------
library(kableExtra)
library(scales)   # for rescale()

# ---- read phenology ------------------------------------------------------
phenology <- read.csv("data/analysis_ready_data/phenology_df.csv")

# ---- (0) helper: map plots to areas, median SOS/EOS by area --------------
plot_area <- clim_plot_month %>%
  distinct(plot, area)                                # 1 row per plot

sos_eos_area <- phenology %>%
  left_join(plot_area, by = "plot") %>%               # add area
  dplyr::filter(!is.na(SOS_doy), !is.na(EOS_doy)) %>%
  group_by(area) %>%
  summarise(
    median_SOS = median(SOS_doy, na.rm = TRUE),
    median_EOS = median(EOS_doy, na.rm = TRUE),
    .groups    = "drop"
  )

# ---- (1) add year + mid-month DOY to climate data ------------------------
clim_tbl <- clim_plot_month %>%
  mutate(
    year    = year(date),
    doy_mid = yday(date) + days_in_month(date) / 2
  )

# ---- (2) keep rows inside the vegetation period -------------------------
veg_clim <- clim_tbl %>%
  inner_join(phenology, by = c("plot", "year")) %>%   # SOS/EOS added
  dplyr::filter(!is.na(SOS_doy), !is.na(EOS_doy)) %>%
  dplyr::filter(doy_mid >= SOS_doy & doy_mid <= EOS_doy)

# ---- (3) plot-year aggregates -------------------------------------------
plot_year <- veg_clim %>%
  group_by(area, plot, year) %>%
  summarise(
    temp_mean_plot_year = mean(temp, na.rm = TRUE),   # °C
    prec_sum_plot_year  = sum(prec,  na.rm = TRUE),   # mm
    .groups = "drop"
  )

# ---- (4) area-year averages of plot values ------------------------------
area_year <- plot_year %>%
  group_by(area, year) %>%
  summarise(
    mean_temp_yr     = mean(temp_mean_plot_year, na.rm = TRUE),
    mean_prec_sum_yr = mean(prec_sum_plot_year,  na.rm = TRUE),
    .groups = "drop"
  )

# ---- (5) long-term means per area ---------------------------------------
area_means <- area_year %>%
  group_by(area) %>%
  summarise(
    mean_temp_veg     = mean(mean_temp_yr,     na.rm = TRUE), # °C
    mean_sum_prec_veg = mean(mean_prec_sum_yr, na.rm = TRUE), # mm
    .groups = "drop"
  ) %>%
  left_join(sos_eos_area, by = "area")                # add medians

# ---- (6) one-colour shading helpers -------------------------------------
shade_vec <- function(x, low = 10, high = 60) {
  rng <- range(x, na.rm = TRUE)
  rescale(x, to = c(low, high), from = rng)
}

area_means_colour <- area_means %>%
  mutate(
    temp_col = shade_vec(mean_temp_veg),
    prec_col = shade_vec(mean_sum_prec_veg),
    sos_col  = shade_vec(median_SOS),   # ← NEW
    eos_col  = shade_vec(median_EOS)    # ← NEW
  )

# ---- (7) build the LaTeX table ------------------------------------------
veg_table_tex <- kbl(
  area_means_colour %>%
    transmute(
      Area = area,
      `Median\\newline SOS (DOY)` =
        cell_spec(sprintf("%.0f", median_SOS),
                  "latex",
                  color      = "black",
                  background = sprintf("yellow!%0.f", sos_col)),   # ← shaded
      `Median\\newline EOS (DOY)` =
        cell_spec(sprintf("%.0f", median_EOS),
                  "latex",
                  color      = "black",
                  background = sprintf("yellow!%0.f", eos_col)),   # ← shaded
      `Vegetation\\textendash period\\newline mean $T$ (°C)` =
        cell_spec(sprintf("%.1f", mean_temp_veg),
                  "latex",
                  color      = "black",
                  background = sprintf("yellow!%0.f", temp_col)),
      `Vegetation\\textendash period\\newline mean $\\Sigma P$ (mm)` =
        cell_spec(sprintf("%.0f", mean_sum_prec_veg),
                  "latex",
                  color      = "black",
                  background = sprintf("yellow!%0.f", prec_col))
    ),
  format   = "latex",
  booktabs = TRUE,
  escape   = FALSE,
  caption  = "Vegetation-period climate and phenology by study area."
) %>%
  kable_styling(position = "center",
                latex_options = "hold_position")

# ---- (8) print LaTeX code to console ------------------------------------
cat(veg_table_tex)

