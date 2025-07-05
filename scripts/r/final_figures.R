library(dplyr)
library(ggplot2)
library(readxl)

heat_all <- readRDS("corr_table_vi_tri.rds")
plot_meta <- read_xlsx("data/vector_data/plot_metadata.xlsx")
heat_all <- dplyr::left_join(heat_all, plot_meta, by = "plot")

head(heat_all)

# ── 1. keep the 10 strongest positive r per plot/VI ───────────────────────
best_combo <- heat_all %>% 
  dplyr::filter(!is.na(correlation)) %>% 
  group_by(plot, VI) %>% 
  slice_max(correlation, n = 10, with_ties = FALSE) %>% 
  ungroup()

# ── 2. median r by VI  &  order factor levels low → high ──────────────────
med_tbl <- best_combo %>% 
  group_by(VI) %>% 
  summarise(median_r = median(correlation), .groups = "drop") %>% 
  arrange(median_r)                       # ascending (low → high)

best_combo <- best_combo %>% 
  mutate(VI = factor(VI, levels = med_tbl$VI))   # apply ordering

# ── 3. box-plot with species-coloured jitter points ───────────────────────
p1 <- ggplot(best_combo, aes(VI, correlation)) +
        geom_boxplot(outlier.shape = NA) +                      # hide default outliers
        geom_text(data = med_tbl,
                  aes(VI, median_r, label = sprintf("%.2f", median_r)),
                  vjust = -0.6, size = 3) +
        scale_colour_brewer(palette = "Dark2", name = "Species") +
        labs(title    = "Top-10 VI–TRI Correlations per Plot, by VI",
             x        = "Vegetation index (ordered by median r)",
             y        = "Pearson correlation (r)") +
        theme_minimal(base_size = 11) +
        theme(plot.title        = element_text(face = "bold"),
          panel.grid.major.x = element_line(colour = "grey85"),
              panel.grid.major.y = element_blank(),
              legend.position = "right")

ggsave("figures/final_figures/VI_top_10.png", p1, height = 4, width = 11, bg = "white")

library(scales)   

# ── 1. species-level median r (one row per species × VI) ──────────────────────
species_level <- best_combo %>%               
  group_by(species, VI) %>% 
  summarise(species_median = median(correlation), .groups = "drop")

# ── 2. counts of plots per species ────────────────────────────────────────────
species_counts <- best_combo %>% 
  count(species, name = "n") %>%                 # raw plot count
  mutate(
    n = n/120,     
  )

# ── 3. named vector for legend labels ─────────────────────────────────────────
species_labs <- setNames(
  paste0(species_counts$species, " (n = ", species_counts$n, ")"),
  species_counts$species
)
# e.g. "FASY" → "FASY (n = 9)"

# ── 4. build the figure ───────────────────────────────────────────────────────
p_species <- ggplot() +
  # background boxplots
  geom_boxplot(
    data    = best_combo,
    aes(x = VI, y = correlation),
    outlier.shape = NA,
    width   = .6,
    colour  = "grey40",
    fill    = "grey90"
  ) +
  # species-level median lines
  geom_line(
    data    = species_level,
    aes(x = VI, y = species_median, group = species, colour = species),
    linewidth = .7, alpha = .2
  ) +
  # species-level median points
  geom_point(
    data    = species_level,
    aes(x = VI, y = species_median, colour = species),
    size = 2
  ) +
  # overall median-r labels on box tops
  geom_text(
    data = med_tbl,
    aes(VI, median_r, label = sprintf("%.2f", median_r)),
    vjust = -0.6, size = 3
  ) +
  scale_colour_brewer(
    palette = "Dark2",
    name    = "Species",
    labels  = species_labs          # ← counts appear here
  ) +
  labs(
    title    = "Top-10 VI–TRI correlations (species medians)",
    subtitle = "Lines connect each species’ median correlation across vegetation indices",
    x        = "Vegetation index (ordered by overall median r)",
    y        = "Pearson correlation (r)"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    plot.title         = element_text(face = "bold"),
    panel.grid.major.x = element_line(colour = "grey85"),
    panel.grid.major.y = element_blank(),
    legend.position    = "right"
  )

# ── 5. display & save ─────────────────────────────────────────────────────────
print(p_species)

ggsave(
  filename = "figures/final_figures/VI_top10_species_medians.png",
  plot     = p_species,
  width    = 11,
  height   = 5,
  dpi      = 300,
  bg       = "white"
)



library(kableExtra)
library(scales)  # make sure it's loaded

# ------------------------------------------------------------------
# 1. Load data (unchanged)
heat_all  <- readRDS("corr_table_vi_tri.rds")
plot_meta <- readxl::read_xlsx("data/vector_data/plot_metadata.xlsx")
heat_all  <- dplyr::left_join(heat_all, plot_meta, by = "plot")

# ------------------------------------------------------------------
# 2. Keep 10 strongest VI–TRI correlations per plot and get medians
top10 <- heat_all %>%
  dplyr::filter(!is.na(correlation)) %>%
  group_by(plot, VI) %>%
  slice_max(correlation, n = 10, with_ties = FALSE) %>%
  ungroup()

medians <- top10 %>%
  group_by(plot, VI, species) %>%
  summarise(median_r = median(correlation), .groups = "drop")

# ------------------------------------------------------------------
# 3. Re-order & pivot for the kable
vi_order <- c("NDVI","EVI","kNDVI","NIRv","NIRvP","NMDI",
              "NDWI","GVMI","SIPI","GRVI","CIG","GNDVI")

table_data <- medians %>%
  tidyr::pivot_wider(names_from = VI, values_from = median_r) %>%
  dplyr::select(plot, species, dplyr::all_of(vi_order)) %>%
  arrange(species, plot)

# ------------------------------------------------------------------
# 4. Build a *ggplot default red* gradient for the cells
gp_red <- "#F8766D"                                        # ggplot2 default red
all_vi_values <- unlist(table_data[vi_order], use.names = FALSE)

shade_fun <- scales::col_numeric(
  palette = c("#FFFFFF", gp_red),                          # white → red
  domain  = range(all_vi_values, na.rm = TRUE)
)

for (vi in vi_order) {
  table_data[[vi]] <- mapply(function(val) {
    bg_col <- shade_fun(val)
    cell_spec(sprintf("%.2f", val),
              format      = "latex",
              background  = bg_col,
              color       = "black")
  }, table_data[[vi]])
}

# ------------------------------------------------------------------
# 5. Produce the LaTeX table (no grey striping)
kbl(
  table_data,
  format    = "latex",
  booktabs  = TRUE,
  escape    = FALSE,
  caption   = paste(
    "Median of the ten strongest VI–TRI correlations per plot.",
    "Cell shading — deeper", gp_red, "= stronger correlation."
  ),
  col.names = c("Plot", "Species", vi_order)
) %>%
  kable_styling(latex_options = "hold_position")   # <- no "striped"

# ── 4. print the medians for reference ────────────────────────────────────
med_tbl

best_combo <- heat_all %>% 
  dplyr::filter(!is.na(correlation)) %>% 
  group_by(plot, VI) %>% 
  slice_max(correlation, n = 100, with_ties = FALSE) %>% 
  ungroup()

# ── 3. median r, DOY and window by VI and lag ─────────────────────────────
med_details <- best_combo %>%
  group_by(VI, lag) %>%
  summarise(
    median_r      = median(correlation, na.rm = TRUE),
    median_doy    = median(doy, na.rm = TRUE),
    median_window = median(window, na.rm = TRUE),
    count         = n(),  # number of records (i.e., top correlations in this lag)
    .groups       = "drop"
  ) %>%
  arrange(lag, median_r)


library(lubridate)

# keep only the VIs of interest, lag-0
vi_keep <- c("NDVI","kNDVI", "GVMI","NMDI", "SIPI", "GRVI")
best_combo <- readRDS("corr_table_vi_tri.rds") %>%
  dplyr::filter(!is.na(correlation), lag == 0, VI %in% vi_keep) %>%
  group_by(plot, VI) %>% slice_max(correlation, n = 100) %>% ungroup() %>%
  mutate(
    VI  = factor(VI, levels = c("NDVI","EVI","kNDVI","NIRv","NIRvP",
                                "NMDI","NDWI","GVMI","SIPI","GRVI","CIG", "GNDVI")),
    month = month.abb[month(as.Date(doy, origin = "2000-01-01"))],
    win_month = window * 8 / 30.437            # 8-day steps → calendar months
  )

# ── month ticks (one per month, middle of each) ───────────────────────────
x_breaks  <- yday(as.Date(paste0("2000-", 1:12, "-15")))  # 12 numeric DOYs

# build a labels vector of length 12, blanks for even months
x_labels <- month.abb           # "Jan" … "Dec"
x_labels[seq(2, 12, by = 2)] <- ""   # blank out Feb, Apr, …

# ── y-breaks unchanged ────────────────────────────────────────────────────
y_breaks <- seq(0, 24, by = 4)

# ── ggplot ────────────────────────────────────────────────────────────────
p2 <- ggplot(best_combo, aes(doy, win_month)) +
  geom_density_2d_filled(contour_var = "ndensity", alpha = 0.85) +
  geom_vline(xintercept = x_breaks, colour = "grey60",
             linewidth = 0.3, alpha = 0.3) +
  geom_hline(yintercept = y_breaks, colour = "grey60",
             linewidth = 0.3, alpha = 0.3) +
  facet_wrap(forest_type ~ VI, nrow = 2) +
  scale_x_continuous(breaks = x_breaks, labels = x_labels, expand = c(0, 0)) +
  scale_y_continuous(breaks = y_breaks,  expand = c(0, 0)) +
  labs(
    title = "Density of Top 100 VI–TRI Correlations per Plot, by VI and Forest Type",
    x     = "Month",
    y     = "Accumulation window (months)",
    fill  = "Density"
  ) +
  theme(
    strip.text        = element_text(face = "bold"),
    legend.position   = "bottom",
    legend.direction  = "horizontal"
  )+
  scale_fill_viridis_d(
    name   = "Density",
    labels = function(z) sub("\\(([^,]+),([^]]+)\\]", "\\1 – \\2", z)  # "(0,0.1]" → "0 – 0.1"
  )

p2
ggsave("figures/final_figures/VI_top_100_timing.png", p2, height = 6, width = 13, bg = "white")





library(readxl)
library(dplyr)
library(ggplot2)
library(lubridate)
library(RColorBrewer)

dir.create("figures/VI_TRI_heatmaps_by_VI", showWarnings = FALSE)

plot_meta <- read_xlsx("data/vector_data/plot_metadata.xlsx")
heat_all  <- left_join(heat_all, plot_meta, by = "plot")

# plots to omit ------------------------------------------------------------
omit_plots <- c("SF04","SF08","PF03","SF06","SF10",
                "SF17","SF12","SF09")
unique(heat_all$plot)

# custom VI order ----------------------------------------------------------
vi_order <- c("NDVI","EVI","OSAVI","kNDVI","NIRv","NIRvP",
              "NMDI","NDWI","GVMI","SIPI","GRVI","CIG")

# -------------------------------------------------------------------------
## 2 -─ loop over every VI -------------------------------------------------
for (vi in vi_order) {
  
  df <- heat_all %>% 
    filter(VI == vi,
           !plot %in% omit_plots,
           !is.na(correlation)) %>% 
    mutate(
      ext_doy      = if_else(lag == 1, doy, doy + 365),
      month_date   = as.Date(ext_doy - 1, origin = "2000-01-01"),
      window_month = window * 8 / 30.437
    )
  length(unique(df$plot))
  if (nrow(df) == 0) next                       # nothing to plot
  
  ## 2A – facet order & labels  -------------------------------------------
  plot_levels <- df %>%                         # order = species → plot
    distinct(plot, species) %>% 
    arrange(species, plot) %>% 
    pull(plot)
  
  df  <- df %>% mutate(plot = factor(plot, levels = plot_levels))
  
  facet_lab <- df %>% 
    distinct(plot, species, area) %>% 
    mutate(label = sprintf("%s – %s – %s", species, area, plot)) %>% 
    select(plot, label) %>%          # ① keep only two columns
    deframe()                        # ② convert to named vector
  # named vector (names = plot)
  
  ## 2B – strongest point per plot  ---------------------------------------
  best_pts <- df %>% 
    filter(correlation > 0) %>% 
    group_by(plot) %>% 
    slice_max(correlation, n = 1, with_ties = FALSE) %>% 
    ungroup()
  
  ## 2C – colour bands  ----------------------------------------------------
  rng     <- range(df$correlation, na.rm = TRUE)
  breaks  <- seq(floor(rng[1]/0.1)*0.1,
                 ceiling(rng[2]/0.1)*0.1,
                 by = 0.1)
  spectral_long <- colorRampPalette(
    rev(brewer.pal(11, "Spectral"))
  )(length(breaks) - 1)
  
  ## 2D – plot  ------------------------------------------------------------
  p <- ggplot(df, aes(month_date, window_month, z = correlation)) +
    geom_contour_filled(colour = "grey40", size = 0.2,
                        breaks = breaks, na.rm = TRUE) +
    facet_wrap(~ plot, ncol = 4,
               labeller = labeller(plot = facet_lab)) +
    scale_x_date(date_breaks = "2 month", date_labels = "%b",
                 expand = c(0, 0)) +
    scale_y_continuous(breaks = seq(0,
                                    ceiling(max(df$window_month, na.rm = TRUE)),
                                    by = 4),
                       expand = c(0, 0)) +
    
    geom_point(data = best_pts,
               aes(month_date, window_month),
               size = 2, shape = 21, stroke = 0.8,
               fill = "yellow", colour = "black") +
    geom_text(data = best_pts,
              aes(month_date, window_month,
                  label = sprintf("r = %.2f", correlation)),
              vjust = -0.6, size = 3) +
    
    labs(title    = sprintf("VI–TRI rolling correlations for %s (2000 – 2022)", vi),
         subtitle = "Facets ordered by species — strip shows species – area",
         x        = "Month (previous year → current year)",
         y        = "Window length (months)") +
    theme_minimal(base_size = 11) +
    theme(strip.text = element_text(face = "bold")) +
    scale_fill_manual(values = spectral_long, name = "r")
  
  ggsave(file.path("figures/VI_TRI_heatmaps_by_VI",
                   sprintf("VI_TRI_contour_%s.png", vi)),
         p, width = 14, height = 8, dpi = 300, bg = "white")
}

# ───────────────────────── 0. PACKAGES ───────────────────────────────────────
# install.packages(c("dplyr", "tidyr", "dendextend", "pals"))   # first time only
library(dplyr)
library(tidyr)
library(dendextend)
library(pals)
library(readxl)
# ───────────────────────── 1. LOAD DATA ──────────────────────────────────────
heat_all <- readRDS("corr_table_vi_tri.rds")     # adjust path if required
plot_meta <- read_xlsx("data/vector_data/plot_metadata.xlsx")
heat_all  <- left_join(heat_all, plot_meta, by = "plot")

# ───────────────────────── 2. HELPERS ────────────────────────────────────────
## make a coloured dendrogram from a correlation matrix -----------------------
make_dend <- function(corr_mat, h_cut, palette_fun = glasbey) {
  d   <- as.dist(1 - corr_mat)            # dissimilarity = 1 − r
  hc  <- hclust(d, method = "average")
  dend <- as.dendrogram(hc)
  
  k    <- length(unique(cutree(hc, h = h_cut)))
  cols <- palette_fun(k)
  
  dend <- color_branches(dend, h = h_cut, col = cols)
  labels_colors(dend) <- get_leaves_branches_col(dend)
  labels_cex(dend)    <- 0.9
  hang.dendrogram(dend, hang = 0.1)
}

## find right-margin lines needed for labels ----------------------------------
right_lines <- function(dend) {
  w_in   <- max(strwidth(labels(dend), units = "inches"))
  char_w <- par("cin")[1]                 # width of one margin line (inches)
  ceiling(w_in / char_w) + 1
}

# ───────────────────────── 3. DENDROGRAM 1  (plots × plots | kNDVI) ──────────
kndvi_mat <- heat_all                                                  %>% 
  dplyr::filter(tolower(VI) == "evi")                                       %>% 
  mutate(time_id = paste0("w", window, "_d", doy, "_lag", lag))        %>% 
  group_by(time_id, plot)                                              %>% 
  summarise(corr = mean(correlation, na.rm = TRUE), .groups = "drop")  %>% 
  pivot_wider(names_from = plot, values_from = corr, values_fill = NA) %>% 
  select(-time_id)                                                     %>% 
  cor(use = "pairwise.complete.obs")

dend_plot <- make_dend(kndvi_mat, h_cut = 0.20)

## --------- ► neue Labels: "Plot – Species"  ------------------
labs         <- labels(dend_plot)                             # aktuelle Plots
sp_vec       <- plot_meta$species[match(labs, plot_meta$plot)]
new_labels   <- paste(labs, sp_vec, sep = " – ")

labels(dend_plot) <- new_labels                               # ersetzen
labels_colors(dend_plot) <- get_leaves_branches_col(dend_plot) # Farben beibehalten

# ───────────────────────── 4. DENDROGRAM 2  (VI × VI | all plots) ────────────
vi_mat <- heat_all                                                     %>% 
  mutate(rec_id = paste(plot, "w", window, "d", doy, "lag", lag, sep = "_")) %>% 
  group_by(rec_id, VI)                                                 %>% 
  summarise(corr = mean(correlation, na.rm = TRUE), .groups = "drop")  %>% 
  pivot_wider(names_from = VI, values_from = corr, values_fill = NA)   %>% 
  select(-rec_id)                                                      %>% 
  cor(use = "pairwise.complete.obs")

dend_vi <- make_dend(vi_mat, h_cut = 0.20)

# ───────────────────────── 5. PLOT GRID (base graphics) ──────────────────────
png("figures/final_figures/dendrogram_grid.png", width = 4000, height = 1750, res = 300)  # adjust size/res as needed

op <- par(no.readonly = TRUE)                 # save current settings
par(mfrow = c(1, 2))                          # 1 row × 2 columns

## — left panel: VI similarity tree ------------------------------------------
par(mar = c(4, 2, 2, right_lines(dend_vi)))
plot(dend_vi, horiz = TRUE,
     main = "VIs Clustered by TRI Correlation Pattern",
     xlab  = "1 − Pearson r")
abline(v = 0.20, lty = 2)

## — right panel: plot similarity tree -----------------------------------------
par(mar = c(4, 2, 2, right_lines(dend_plot)))
plot(dend_plot, horiz = TRUE,
     main = "Plots Clustered by kNDVI ↔ TRI Correlation Pattern",
     xlab = "1 − Pearson r")
abline(v = 0.20, lty = 2)


par(op)       
dev.off()

# ───────────────────────── 0. PACKAGES ───────────────────────────────────────
library(dplyr)
library(tidyr)
library(ggcorrplot)
library(patchwork)   # simple plot-grid with |

# ───────────────────────── 1. LOAD DATA ──────────────────────────────────────
heat_all <- readRDS("corr_table_vi_tri.rds")   # adjust path as needed

# ───────────────────────── 2. PLOT × PLOT CORR  (kNDVI only) ────────────────
corr_plot <- heat_all %>% 
  dplyr::filter(tolower(VI) == "kndvi") %>% 
  mutate(id = paste0("w", window, "_d", doy, "_lag", lag)) %>% 
  group_by(id, plot) %>% 
  summarise(corr = mean(correlation, na.rm = TRUE), .groups = "drop") %>% 
  pivot_wider(names_from = plot, values_from = corr, values_fill = NA) %>% 
  select(-id) %>% 
  cor(use = "pairwise.complete.obs")

p1 <- ggcorrplot(corr_plot, type = "upper", lab = TRUE, 
                 title = "Plot × Plot (kNDVI)")

# ───────────────────────── 3. VI × VI CORR  (all plots) ──────────────────────
corr_vi <- heat_all %>% 
  mutate(id = paste(plot, "w", window, "d", doy, "lag", lag, sep = "_")) %>% 
  group_by(id, VI) %>% 
  summarise(corr = mean(correlation, na.rm = TRUE), .groups = "drop") %>% 
  pivot_wider(names_from = VI, values_from = corr, values_fill = NA) %>% 
  select(-id) %>% 
  cor(use = "pairwise.complete.obs")

p2 <- ggcorrplot(corr_vi, type = "upper", lab = TRUE, 
                 title = "VI × VI (all plots)")
ggsave("figures/final_figures/plot_plot_kndvi_cor.png", p1, bg = "white", height = 11 , width = 12)
ggsave("figures/final_figures/Vi_VI_cor.png", p2, bg = "white", height = 8 , width = 9)


# ggsave("two_coreplots.png", p1 | p2, width = 18, height = 9, units = "cm", dpi = 300)

