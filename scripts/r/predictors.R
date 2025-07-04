# 1 ── libraries ──────────────────────────────────────────────
library(sf)         # vector data handling
library(terra)      # raster data handling
library(dplyr)      # data manipulation
library(tidyr)      # reshaping data
library(lubridate)  # date handling

# 2 ── paths (edit if needed) ─────────────────────────────────
shp_file        <- "data/vector_data/all_merged.shp"               # polygons (area names)
points_file     <- "data/vector_data/trees_all_plots.gpkg"         # sampling points

# 3 ── vectors & rasters ──────────────────────────────────────

## 3.1 Read area polygons and assign area names
areas <- st_read(shp_file)
areas$area <- areas$name <- c("Eifel", "Kellerwald", "Koenigsforst",
                              "Calderner Wald", "Lindenberger Wald")

## 3.2 Read sampling points (keep only `plot` column)
points <- st_read(points_file) %>% 
  select(plot)

# 5 ── extract raster values at points ────────────────────────
pts_v <- vect(points)

## 6.1 Create point metadata table
pt_meta <- points |> 
  st_drop_geometry() |> 
  mutate(pt_id = row_number()) |> 
  select(pt_id, plot)


# ─────────────────────────────────────────────────────────────
# 5 bis ─  Topographic predictors at points  ──────────────────

## 5.1  load the 4-band raster we wrote earlier --------------
predictors  <- rast("data/satellite_data/predictors/predictors.tif")
names(predictors) <- c("slope", "aspect", "dem", "ai")  # rename for clarity

## 5.2  extract at the same sampling points ------------------
topo_df <- terra::extract(predictors, pts_v, ID = FALSE) |>
  as_tibble() |>
  mutate(pt_id = row_number())        # same id key as climate tables

# ─────────────────────────────────────────────────────────────
# 6 bis ─  summarise to AREA level  ──────────────────────────

library(dplyr)

topo_area <- topo_df |>
  left_join(pt_meta, by = "pt_id") |>
  mutate(
    aspect_rad = aspect * pi / 180      # degrees → radians
  ) |>
  group_by(plot) |>
  summarise(
    dem_mean   = mean(dem,   na.rm = TRUE),      # m a.s.l.
    slope_mean = mean(slope, na.rm = TRUE),      # degrees
    # circular mean of aspect (convert back to 0-360°)
    aspect_mean = {
      y <- mean(sin(aspect_rad), na.rm = TRUE)
      x <- mean(cos(aspect_rad), na.rm = TRUE)
      (atan2(y, x) * 180 / pi) %% 360
    },
    ai_mean    = mean(ai,    na.rm = TRUE),      # if you need aridity index
    .groups = "drop"
  )

# ─────────────────────────────────────────────────────────────
# 7 bis ─  merge with the existing area-level climate table ──

# optional: save
write.csv(topo_area,
          "data/analysis_ready_data/topo_ai_df.csv",
          row.names = FALSE)
