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
  filename = "figures/final_figures/VI_top_10.png",
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





# ───────────────────────── 0. PACKAGES ───────────────────────────────────────
library(dplyr)
library(ggplot2)
library(ggrepel)
library(readxl)
library(RColorBrewer)

# ───────────────────────── 1. LOAD DATA (adjust paths) ───────────────────────
heat_all  <- readRDS("corr_table_vi_tri.rds")
plot_meta <- read_xlsx("data/vector_data/plot_metadata.xlsx")
heat_all  <- left_join(heat_all, plot_meta, by = "plot")

# ───────────────────────── 2. PREPARE DATA — EVI for SF02 & SF07 ─────────────
df <- heat_all %>% 
  filter(plot %in% c("SF02", "SF07"),
         tolower(VI) == "evi",
         !is.na(correlation)) %>% 
  mutate(
    ext_doy      = if_else(lag == 1, doy, doy + 365),
    month_date   = as.Date(ext_doy - 1, origin = "2000-01-01"),
    window_month = window * 8 / 30.437,
    
    ## ─ NEW ─ build facet label: 'SF02 – Species – Area'
    plot_lab = sprintf("%s – %s – %s", plot, species, area)
  )

# ── strongest positive r per plot (now uses plot_lab) ────────────────────────
best_pts <- df %>% 
  filter(correlation > 0) %>% 
  group_by(plot_lab) %>% 
  slice_max(correlation, n = 1, with_ties = FALSE) %>% 
  ungroup()

# ── colour-band setup (unchanged) ────────────────────────────────────────────
rng     <- range(df$correlation, na.rm = TRUE)
breaks  <- seq(floor(rng[1] / 0.1) * 0.1,
               ceiling(rng[2] / 0.1) * 0.1,
               by = 0.1)
n_bands <- length(breaks) - 1
spectral_long <- colorRampPalette(
  rev(RColorBrewer::brewer.pal(11, "RdYlBu"))
)(n_bands)

x_grid <- seq(min(df$month_date), max(df$month_date), by = "2 month")
y_grid <- seq(0, ceiling(max(df$window_month, na.rm = TRUE)), by = 2)

# ───────────────────────── 3. BUILD THE GRID PLOT ────────────────────────────
p <- ggplot(df, aes(month_date, window_month, z = correlation)) +
  geom_contour_filled(colour = "gray50", size = 0.2,
                      na.rm = TRUE, breaks = breaks) +
  facet_wrap(~plot_lab, ncol = 2) +             # ← uses new label
  scale_x_date(date_breaks = "2 month",
               date_labels = "%b",
               expand = c(0, 0)) +
  scale_y_continuous(breaks = seq(0,
                                  ceiling(max(df$window_month, na.rm = TRUE)), 2),
                     expand = c(0, 0)) +
  
  ## highlight the max-r points ------------------------------------------------
geom_point(
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
    fill        = alpha("white", 0.3),
    colour      = "black",
    label.size  = 0,
    label.r     = unit(0.1, "lines"),
    box.padding   = 0.3,
    point.padding = 0.2,
    segment.color = NA
  ) +
  
  labs(
    title = "EVI–TRI rolling correlation — plots SF02 & SF07",
    fill  = expression("Pearson "~italic(r)),
    x     = "Month (prev year → current year)",
    y     = "Window length (months)"
  ) +
  scale_fill_manual(values = spectral_long) +
  theme(title  = element_text(face = "bold"),
        plot.title  = element_text(hjust = 0.5),
        strip.text  = element_text(face = "bold"))+
  geom_vline(xintercept = x_grid,
             colour = "grey50", linetype = "dotted", size = 0.3) +
  geom_hline(yintercept = y_grid,
             colour = "grey50", linetype = "dotted", size = 0.3)

print(p)
# ───────────────────────── 4. SAVE THE FIGURE ────────────────────────────────

ggsave(
  filename = "figures/final_figures/EVI_TRI_SF02_SF07_grid.png",
  plot     = p,
  width    = 12, height = 5, dpi = 300, bg = "white"
)


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

# ───────────────────────── 5. EXPORT SEPARATE DENDROGRAMS ───────────────────

# Common device specs
fig_w  <- 2000   # pixels   (≈ 170 mm at 300 dpi)
fig_h  <- 1750   # pixels   (≈ 150 mm at 300 dpi)
fig_res <- 300   # dpi

# Keep a copy of the current par settings
op <- par(no.readonly = TRUE)

## ---------- 5a. VI-similarity tree (dend_vi) -------------------------------
png("figures/final_figures/dendrogram_vi.png",
    width = fig_w, height = fig_h, res = fig_res)

par(mfrow = c(1, 1),                      # single panel
    mar = c(4, 2, 2, right_lines(dend_vi)))  # dynamic right margin

plot(dend_vi, horiz = TRUE,
     main = "VIs Clustered by TRI Correlation Pattern",
     xlab = "1 − Pearson r")
abline(v = 0.20, lty = 2)

dev.off()                                 # close first device

## ---------- 5b. Plot-similarity tree (dend_plot) ---------------------------
png("figures/final_figures/dendrogram_plot.png",
    width = fig_w, height = fig_h, res = fig_res)

par(mfrow = c(1, 1),
    mar = c(4, 2, 2, right_lines(dend_plot)))

plot(dend_plot, horiz = TRUE,
     main = "Plots Clustered by EVI ↔ TRI Correlation Pattern",
     xlab = "1 − Pearson r")
abline(v = 0.20, lty = 2)

dev.off()                                 # close second device

# Restore original graphics settings
par(op)


# ───────────────────────── 0. PACKAGES ───────────────────────────────────────
library(dplyr)
library(tidyr)
library(ggcorrplot)
library(patchwork)   # simple plot-grid with |

# ───────────────────────── 1. LOAD DATA ──────────────────────────────────────
heat_all <- readRDS("corr_table_vi_tri.rds")   # adjust path as needed

# ───────────────────────── 2. PLOT × PLOT CORR  (kNDVI only) ────────────────
corr_plot <- heat_all %>% 
  dplyr::filter(tolower(VI) == "evi") %>% 
  mutate(id = paste0("w", window, "_d", doy, "_lag", lag)) %>% 
  group_by(id, plot) %>% 
  summarise(corr = mean(correlation, na.rm = TRUE), .groups = "drop") %>% 
  pivot_wider(names_from = plot, values_from = corr, values_fill = NA) %>% 
  select(-id) %>% 
  cor(use = "pairwise.complete.obs")

p1 <- ggcorrplot(corr_plot, type = "upper", lab = TRUE, 
                 title = "Plot × Plot (EVI)")

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

