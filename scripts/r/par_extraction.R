library(terra)      # fast raster IO
library(sf)         # vector data
library(dplyr)
library(tidyr)
library(stringr)

## ------------------------------------------------------------------
## 1  Load & stack the PAR rasters
## ------------------------------------------------------------------
par_paths <- list.files("data/satellite_data/daily_par",
                        pattern = "\\.tif$", full.names = TRUE) |>
  sort()                     # make sure years are in order

# read each yearly file and concatenate the layers
par_stack <- do.call(c, lapply(par_paths, rast))

# clean layer names: "YYYY_MM_DD"
names(par_stack) <- str_remove(names(par_stack), "_PAR_daily")

## ------------------------------------------------------------------
## 2  Load tree points & re-project to the raster CRS
## ------------------------------------------------------------------
points <- read_sf("data/vector_data/trees_all_plots.gpkg") |>
  st_transform(crs(par_stack))

points_tv <- vect(points)   # terra’s SpatVector

## ------------------------------------------------------------------
## 3  Extract daily PAR at every tree
##     → matrix: n_trees × n_days
## ------------------------------------------------------------------
par_vals <- terra::extract(par_stack, points_tv, ID = FALSE)

## ------------------------------------------------------------------
## 4  Bind tree attributes, reshape to long,
##     and compute the plot-level daily mean
## ------------------------------------------------------------------
df_long <- cbind(st_drop_geometry(points), par_vals) |>
  pivot_longer(
    cols = starts_with("20"),                # all PAR layers
    names_to   = "date",
    values_to  = "par"
  ) |>
  mutate(date = as.Date(date, format = "%Y_%m_%d"))

daily_plot_mean <- df_long |>
  group_by(plot, date) |>
  summarise(par_mean = mean(par, na.rm = TRUE), .groups = "drop") |>
  na.omit()

## ------------------------------------------------------------------
## 5  Save or inspect
## ------------------------------------------------------------------
write.csv(daily_plot_mean,
          "data/analysis_ready_data/daily_PAR_mean_by_plot.csv",
          row.names = FALSE)
