# ── packages ───────────────────────────────────────────────────────────────
library(terra)      # raster I/O & math
library(sf)         # polygons
library(pbapply)    # progress bar without parallel
library(stringr)

# ── 0. paths & polygons ────────────────────────────────────────────────────
in_root  <- "data/satellite_data/SDC"
out_root <- "data/satellite_data/kNDVI"
dir.create(out_root, recursive = TRUE, showWarnings = FALSE)

shp_all  <- read_sf("data/vector_data/all_merged.shp") |>
  st_transform(crs = 32632)

lin_eif_koe <- shp_all[c(1,3,5), ]
lah_kel     <- shp_all[c(2,4), ]

# two first-level sub-dirs (assumed area folders)
subdirs <- list.dirs(in_root, full.names = TRUE, recursive = FALSE)
lah_kel_paths     <- list.files(subdirs[1], pattern = "\\.tif$", recursive = TRUE, full.names = TRUE)
lin_eif_koe_paths <- list.files(subdirs[2], pattern = "\\.tif$", recursive = TRUE, full.names = TRUE)

terraOptions(progress = 0, threads = 1)   # single-thread, no terra bar

# ── 1. helper: clip, compute kNDVI, save ───────────────────────────────────
process_scene <- function(path, clip_shape) {
  
  # read only Red (band 4) and NIR (band 5) and clip immediately
  pair <- rast(path, lyrs = c(4, 5)) |>
    crop(clip_shape)  * 1e-4
    
  names(pair) <- c("red", "nir")
  
  red <- ifel(is.finite(pair$red), pair$red, NA)
  nir <- ifel(is.finite(pair$nir), pair$nir, NA)
  
  sigma <- global(abs(nir - red), fun = "mean", na.rm = TRUE)[[1]]
  
  kndvi <- (1 - exp(-(nir - red)^2 / (2 * sigma^2))) /
    (1 + exp(-(nir - red)^2 / (2 * sigma^2)))
  kndvi <- ifel(is.finite(kndvi), kndvi, NA)
  names(kndvi) <- "kNDVI"
  
  # mirrored output directory
  rel_dir <- dirname(sub(in_root, "", path, fixed = TRUE))
  out_dir <- file.path(out_root, rel_dir)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  out_path <- file.path(out_dir, paste0("kNDVI_", basename(path)))
  writeRaster(kndvi, out_path, overwrite = TRUE)
}

# ── 2. run sequentially with pbapply progress bar ──────────────────────────
pboptions(type = "timer")

cat("\n▸ Processing LAH-KEL area …\n")
pblapply(lah_kel_paths, function(p) process_scene(p, lah_kel))

cat("\n▸ Processing LIN-EIF-KOE area …\n")
pblapply(lin_eif_koe_paths, function(p) process_scene(p, lin_eif_koe))
