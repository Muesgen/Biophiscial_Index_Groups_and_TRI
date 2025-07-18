# ── 0. Libraries ───────────────────────────────────────────────────────
library(tidyverse)   # dplyr, ggplot2, purrr, …
library(readxl)      # read_xlsx()
library(scales)      # pretty_breaks(), etc.

# ── 1. Data ────────────────────────────────────────────────────────────
heat_all  <- readRDS("corr_table_vi_tri.rds")
plot_meta <- read_xlsx("data/vector_data/plot_metadata.xlsx")

heat_all <- left_join(heat_all, plot_meta, by = "plot")

# ── 2. Top-10 |r| per VI × plot  (lag is ignored) ──────────────────────
best_combo <- heat_all %>% 
  filter(!is.na(correlation)) %>% 
  group_by(plot, VI) %>%                     # lag NOT in the group
  slice_max(abs(correlation), n = 10,        # strongest |r|
            with_ties = FALSE) %>% 
  ungroup()

# ── 3. Overall median |r| per VI  (for ordering & labels) ──────────────
med_tbl <- best_combo %>% 
  group_by(VI) %>% 
  summarise(median_r = median(abs(correlation)), .groups = "drop")

# put VIs in descending-median order
vi_order <- med_tbl %>% arrange((median_r)) %>% pull(VI)
best_combo   <- mutate(best_combo,   VI = factor(VI, levels = vi_order))
med_tbl      <- mutate(med_tbl,      VI = factor(VI, levels = vi_order))

# ── 4. Species-level summaries ─────────────────────────────────────────
species_level <- best_combo %>% 
  group_by(species, VI) %>% 
  summarise(species_median = median(abs(correlation)), .groups = "drop")

species_counts <- best_combo %>%       # unique plots per species
  distinct(plot, species) %>% 
  count(species, name = "n") %>% 
  mutate(n = n )                  # << adjust divisor if 120≠total-plots

species_labs <- setNames(
  paste0(species_counts$species, " (n = ", species_counts$n, ")"),
  species_counts$species
)

# ── 5. Plot ────────────────────────────────────────────────────────────
p_species <- ggplot() +
  geom_boxplot(
    data    = best_combo,
    aes(VI, abs(correlation)),
    outlier.shape = NA, width = .6,
    colour = "grey40",  fill = "grey90"
  ) +
  geom_line(
    data = species_level,
    aes(VI, species_median, group = species, colour = species),
    linewidth = .7, alpha = .2
  ) +
  geom_point(
    data = species_level,
    aes(VI, species_median, colour = species),
    size = 2
  ) +
  geom_text(                               # median |r| on box tops
    data = med_tbl,
    aes(VI, median_r, label = sprintf("%.2f", median_r)),
    vjust = -0.6, size = 3
  ) +
  scale_colour_brewer(
    palette = "Dark2", name = "Species", labels = species_labs
  ) +
  labs(
    title    = "Top-10 VI–TRI correlations by VI",
    subtitle = "Lines connect each species’ median |r| across vegetation indices",
    x        = "Vegetation index (ordered by overall median |r|)",
    y        = "Pearson correlation (|r|)"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    plot.title         = element_text(face = "bold"),
    panel.grid.major.x = element_line(colour = "grey85"),
    panel.grid.major.y = element_blank(),
    legend.position    = "right"
  )

print(p_species)

# ── 6. Save ────────────────────────────────────────────────────────────
ggsave(
  filename = "figures/final_figures/VI_top10_species_medians.png",
  plot     = p_species,
  width    = 11, height = 5, dpi = 300, bg = "white"
)




# ───────────────────────── 0.  PACKAGES ─────────────────────────
library(dplyr)
library(tidyr)
library(readxl)
library(kableExtra)
library(scales)

# ───────────────────────── 1.  LOAD DATA ───────────────────────
heat_all  <- readRDS("corr_table_vi_tri.rds")
plot_meta <- readxl::read_xlsx("data/vector_data/plot_metadata.xlsx")
heat_all  <- left_join(heat_all, plot_meta, by = "plot")

# ───────────────────────── 2.  TOP-10 MEDIANS ──────────────────
top10 <- heat_all %>%
  filter(!is.na(correlation)) %>%
  group_by(plot, VI) %>%
  slice_max(abs(correlation),               # ← use magnitude, not sign
            n = 10, with_ties = FALSE) %>%
  ungroup()

medians <- top10 %>%
  group_by(plot, VI, species) %>%
  summarise(median_r = median(abs(correlation)), .groups = "drop")

# ───────────────────────── 3.  PIVOT WIDER  ─────────────────────
vi_order <- c("NDVI","GNDVI","CIG","EVI","kNDVI","NIRv","NIRvP",
              "NMDI","NDWI","GVMI","SIPI","GRVI")

table_data <- medians %>%
  pivot_wider(names_from = VI, values_from = median_r) %>%
  select(plot, species, all_of(vi_order))

# ───── 3b.  ORDER & COLOUR ROWS LIKE THE DENDROGRAM ────────────
leaf_labels <- labels(dend_plot)
leaf_cols   <- labels_colors(dend_plot)            # already "#RRGGBB"
plot_ids    <- sub("^([^ ]+).*", "\\1", leaf_labels)
plot_col_vec <- setNames(leaf_cols, plot_ids)

table_data <- table_data %>%
  mutate(plot = factor(plot, levels = rev(plot_ids))) %>%
  arrange(plot)

# ───── 3c.  SIGNIFICANT-COUNT MATRIX ───────────────────────────
sig_counts <- heat_all %>%
  filter(p_value < 0.05) %>%
  group_by(plot, VI) %>%
  summarise(n_sig = n(), .groups = "drop") %>%
  pivot_wider(names_from = VI, values_from = n_sig) %>%
  select(plot, all_of(vi_order))

# ───────── 4.  RED GRADIENT (while numerics intact) ────────────
gp_red        <- "#F8766D"
all_vi_values <- unlist(table_data[vi_order], use.names = FALSE)
shade_fun <- col_numeric(c("#FFFFFF", gp_red),
                         domain = range(all_vi_values, na.rm = TRUE))

# ───────── 5.  SHADE VI CELLS  ─────────────────────────────────
for (vi in vi_order) {
  cnt <- sig_counts[[vi]][match(table_data$plot, sig_counts$plot)]
  cnt[is.na(cnt)] <- 0
  
  table_data[[vi]] <- mapply(function(val, n_sig){
    if (n_sig < 10){
      cell_spec(sprintf("%.2f", val), format = "latex",
                background = "#D9D9D9", color = "black")   ## ← changed
    } else {
      cell_spec(sprintf("%.2f", val), format = "latex",
                background = shade_fun(val), color = "black") ## ← keep "#"
    }
  }, table_data[[vi]], cnt, SIMPLIFY = FALSE)
}

# ───────── 6.  COLOUR PLOT & SPECIES LABELS  ───────────────────
table_data <- table_data %>%
  mutate(row_hex = plot_col_vec[as.character(plot)],             # with "#"
         plot    = cell_spec(as.character(plot),
                             format = "latex", color = row_hex),
         species = cell_spec(as.character(species),
                             format = "latex", color = row_hex)) %>%
  select(-row_hex)

# ───────── 7.  RENDER LaTeX TABLE  ─────────────────────────────
old_fmt <- getOption("knitr.table.format")
options(knitr.table.format = "latex")

kbl(table_data,
    format   = "latex",
    booktabs = TRUE,
    escape   = FALSE,
    caption  = paste(
      "Median of the ten strongest VI–TRI correlations per plot.",
      "Cells with fewer than 10 significant correlations are greyed out;",
      "deeper", gp_red, "= stronger correlation."
    ),
    col.names = c("Plot","Species", vi_order)) %>%
  kable_styling(latex_options = "hold_position")

options(knitr.table.format = old_fmt)




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
library(lubridate)
library(scales)

## -------------------------------------------------------------
## 0. Settings
sig_thr   <- 0.05   # p-value threshold
## -------------------------------------------------------------
## 1. DATA — two plots, one VI
df <- heat_all %>% 
  dplyr::filter(plot %in% c("SF12", "SF21"),
         tolower(VI) == "sipi",
         !is.na(correlation)) %>% 
  mutate(
    sig          = p_value < sig_thr,
    ext_doy      = if_else(lag == 1, doy, doy + 365),
    month_date   = as.Date(ext_doy - 1, origin = "2000-01-01"),
    window_month = window * 8 / 30.437,
    plot_lab     = sprintf("%s – %s – %s", plot, species, area)
  )

df_nsig <- df  %>% dplyr::filter(!sig)          # non-significant rows
df_sig  <- df  %>% dplyr::filter( sig)          # significant rows

## strongest |r| among sig. values, one per panel
best_pts <- df_sig %>% 
  group_by(plot_lab) %>% 
  slice_max(abs(correlation), n = 1, with_ties = FALSE) %>% 
  ungroup()

## -------------------------------------------------------------
## 2. Colour scale (same bin setup as before)
rng     <- range(df$correlation, na.rm = TRUE)
breaks  <- seq(floor(rng[1]/0.1)*0.1, ceiling(rng[2]/0.1)*0.1, 0.1)
nbins   <- length(breaks) - 1
pal     <- colorRampPalette(rev(RColorBrewer::brewer.pal(11,"RdYlBu")))(nbins)

mid_vals <- (head(breaks,-1) + tail(breaks,-1)) / 2
lab_vals <- sprintf("%.2f", mid_vals)
lab_vals[lab_vals == "-0.00"] <- "0.00"

## grid for dotted guides
x_grid <- seq(min(df$month_date), max(df$month_date), by = "2 month")
y_grid <- seq(0, ceiling(max(df$window_month, na.rm = TRUE)), by = 2)

## -------------------------------------------------------------
## 3. PLOT — fade first, opaque second
p <- ggplot() +
  ## significant layer (opaque)
  geom_contour_filled(
    data   = df,
    aes(month_date, window_month,
        z    = correlation,
        fill = after_stat(as.numeric(level))),
    breaks = breaks,
    colour = NA, linewidth = 0,
    na.rm  = TRUE
  ) +
  ## non-significant layer (faded)
  geom_contour_filled(
    data   = df_nsig,
    aes(month_date, window_month,
        z    = correlation,
        fill = after_stat(as.numeric(level))),
    breaks = breaks,
    fill = "white",
    alpha  = 0.5,
    colour = NA, linewidth = 0,
    na.rm  = TRUE
  ) +
  
  ## thin grey isolines for context
  geom_contour(
    data   = df,
    aes(month_date, window_month, z = correlation),
    colour = "grey70", size = 0.15,
    breaks = breaks
  ) +
  ## highlight best significant |r|
  { if (nrow(best_pts) > 0)
    list(
      geom_point(
        data = best_pts,
        aes(month_date, window_month),
        size = 2, shape = 21, stroke = .8,
        fill = "yellow", colour = "black"
      ),
      geom_label_repel(
        data   = best_pts,
        aes(month_date, window_month,
            label = sprintf("r = %.2f", correlation)),
        size = 3,
        fill = scales::alpha("white", 0.3),
        colour = "black",
        label.size = 0,
        label.r    = unit(0.1, "lines"),
        box.padding = 0.3,
        point.padding = 0.2,
        segment.color = NA
      )
    )
  } +
  ## axes, facets, palette
  facet_wrap(~ plot_lab, ncol = 2) +
  scale_x_date(
    date_breaks = "2 month",
    labels      = date_format("%b", locale = "en"),   # English month abbrev.
    expand      = c(0, 0)
  ) +
  scale_y_continuous(
    breaks = y_grid, expand = c(0,0)
  ) +
  geom_vline(xintercept = x_grid, colour = "grey50",
             linetype = "dotted", size = 0.3) +
  geom_hline(yintercept = y_grid, colour = "grey50",
             linetype = "dotted", size = 0.3) +
  scale_fill_stepsn(
    colours = pal,
    limits  = c(1, nbins),
    breaks  = seq_len(nbins),
    labels  = lab_vals,
    name    = expression("Pearson "~italic(r)),
    guide   = guide_colourbar(
      barheight = unit(5, "cm"),
      barwidth  = unit(0.5, "cm"),
      ticks.colour = "black"
    )
  ) +
  labs(
    title = "SIPI–TRI rolling correlation — plots SF12 & SF21",
    x     = "Month (prev year → current year)",
    y     = "Window length (months)"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    plot.title  = element_text(face = "bold", hjust = 0.5),
    strip.text  = element_text(face = "bold")
  )

print(p)
# ───────────────────────── 4. SAVE THE FIGURE ────────────────────────────────

ggsave(
  filename = "figures/final_figures/SIPI_TRI_SF12_SF21_grid.png",
  plot     = p,
  width    = 12, height = 5, dpi = 300, bg = "white"
)

library(dplyr)
library(ggplot2)
library(ggrepel)
library(readxl)
library(RColorBrewer)
library(lubridate)
library(scales)

## -------------------------------------------------------------
## 0. Settings
sig_thr   <- 0.05   # p-value threshold
## -------------------------------------------------------------
## 1. DATA — two plots, one VI
df <- heat_all %>% 
  dplyr::filter(plot %in% c("SF04", "SF08"),
                tolower(VI) == "gndvi",
                !is.na(correlation)) %>% 
  mutate(
    sig          = p_value < sig_thr,
    ext_doy      = if_else(lag == 1, doy, doy + 365),
    month_date   = as.Date(ext_doy - 1, origin = "2000-01-01"),
    window_month = window * 8 / 30.437,
    plot_lab     = sprintf("%s – %s – %s", plot, species, area)
  )

df_nsig <- df  %>% dplyr::filter(!sig)          # non-significant rows
df_sig  <- df  %>% dplyr::filter( sig)          # significant rows

## strongest |r| among sig. values, one per panel
best_pts <- df_sig %>% 
  group_by(plot_lab) %>% 
  slice_max(abs(correlation), n = 1, with_ties = FALSE) %>% 
  ungroup()

## -------------------------------------------------------------
## 2. Colour scale (same bin setup as before)
rng     <- range(df$correlation, na.rm = TRUE)
breaks  <- seq(floor(rng[1]/0.1)*0.1, ceiling(rng[2]/0.1)*0.1, 0.1)
nbins   <- length(breaks) - 1
pal     <- colorRampPalette(rev(RColorBrewer::brewer.pal(11,"RdYlBu")))(nbins)

mid_vals <- (head(breaks,-1) + tail(breaks,-1)) / 2
lab_vals <- sprintf("%.2f", mid_vals)
lab_vals[lab_vals == "-0.00"] <- "0.00"

## grid for dotted guides
x_grid <- seq(min(df$month_date), max(df$month_date), by = "2 month")
y_grid <- seq(0, ceiling(max(df$window_month, na.rm = TRUE)), by = 2)

## -------------------------------------------------------------
## 3. PLOT — fade first, opaque second
p <- ggplot() +
  ## significant layer (opaque)
  geom_contour_filled(
    data   = df,
    aes(month_date, window_month,
        z    = correlation,
        fill = after_stat(as.numeric(level))),
    breaks = breaks,
    colour = NA, linewidth = 0,
    na.rm  = TRUE
  ) +
  ## non-significant layer (faded)
  geom_contour_filled(
    data   = df_nsig,
    aes(month_date, window_month,
        z    = correlation,
        fill = after_stat(as.numeric(level))),
    breaks = breaks,
    fill = "white",
    alpha  = 0.5,
    colour = NA, linewidth = 0,
    na.rm  = TRUE
  ) +
  
  ## thin grey isolines for context
  geom_contour(
    data   = df,
    aes(month_date, window_month, z = correlation),
    colour = "grey70", size = 0.15,
    breaks = breaks
  ) +
  ## highlight best significant |r|
  { if (nrow(best_pts) > 0)
    list(
      geom_point(
        data = best_pts,
        aes(month_date, window_month),
        size = 2, shape = 21, stroke = .8,
        fill = "yellow", colour = "black"
      ),
      geom_label_repel(
        data   = best_pts,
        aes(month_date, window_month,
            label = sprintf("r = %.2f", correlation)),
        size = 3,
        fill = scales::alpha("white", 0.3),
        colour = "black",
        label.size = 0,
        label.r    = unit(0.1, "lines"),
        box.padding = 0.3,
        point.padding = 0.2,
        segment.color = NA
      )
    )
  } +
  ## axes, facets, palette
  facet_wrap(~ plot_lab, ncol = 2) +
  scale_x_date(
    date_breaks = "2 month",
    labels      = date_format("%b", locale = "en"),   # English month abbrev.
    expand      = c(0, 0)
  ) +
  scale_y_continuous(
    breaks = y_grid, expand = c(0,0)
  ) +
  geom_vline(xintercept = x_grid, colour = "grey50",
             linetype = "dotted", size = 0.3) +
  geom_hline(yintercept = y_grid, colour = "grey50",
             linetype = "dotted", size = 0.3) +
  scale_fill_stepsn(
    colours = pal,
    limits  = c(1, nbins),
    breaks  = seq_len(nbins),
    labels  = lab_vals,
    name    = expression("Pearson "~italic(r)),
    guide   = guide_colourbar(
      barheight = unit(5, "cm"),
      barwidth  = unit(0.5, "cm"),
      ticks.colour = "black"
    )
  ) +
  labs(
    title = "GNDVI–TRI rolling correlation — plots SF04 & SF08",
    x     = "Month (prev year → current year)",
    y     = "Window length (months)"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    plot.title  = element_text(face = "bold", hjust = 0.5),
    strip.text  = element_text(face = "bold")
  )

print(p)
# ───────────────────────── 4. SAVE THE FIGURE ────────────────────────────────

ggsave(
  filename = "figures/final_figures/GNDVI_TRI_SF04_SF08_grid.png",
  plot     = p,
  width    = 12, height = 5, dpi = 300, bg = "white"
)
# ───────────────────────── 0. PACKAGES ───────────────────────────────────────
library(dplyr)
library(ggplot2)
library(ggrepel)
library(readxl)
library(RColorBrewer)

sig_thr <- 0.05                       # p-value threshold
vi_set  <- c("kNDVI")         # two indices
plot_id <- "SF01"                     # one plot

# ───────────────────────── 1. DATA — one plot, two VIs ──────────────────────
df <- heat_all %>% 
  filter(plot == plot_id,
         VI %in% vi_set,
         !is.na(correlation)) %>% 
  mutate(
    sig          = p_value < sig_thr,
    ext_doy      = if_else(lag == 1, doy, doy + 365),
    month_date   = as.Date(ext_doy - 1, origin = "2000-01-01"),
    window_month = window * 8 / 30.437,
    VI           = factor(VI, levels = vi_set), 
    panel_lab    = paste(VI),
  )

df_nsig <- filter(df, !sig)
df_sig  <- filter(df,  sig)

best_pts <- df_sig %>%                      # one “best” per panel
  group_by(panel_lab) %>% 
  slice_max(abs(correlation), n = 1, with_ties = FALSE) %>% 
  ungroup()

# ───────────────────────── 2. COLOUR SCALE (same bins) ──────────────────────
rng    <- range(df$correlation, na.rm = TRUE)
breaks <- seq(floor(rng[1]/0.1)*0.1, ceiling(rng[2]/0.1)*0.1, 0.1)
nbins  <- length(breaks) - 1
pal    <- colorRampPalette(rev(brewer.pal(11, "RdYlBu")))(nbins)

lab_vals <- sprintf("%.2f",
                    (head(breaks, -1) + tail(breaks, -1)) / 2)
lab_vals[lab_vals == "-0.00"] <- "0.00"

x_grid <- seq(min(df$month_date), max(df$month_date), by = "2 month")
y_grid <- seq(0, ceiling(max(df$window_month)), by = 2)

# ───────────────────────── 3. PLOT ───────────────────────────────────────────
p <- ggplot() +
  geom_contour_filled(
    data   = df,
    aes(month_date, window_month,
        z    = correlation,
        fill = after_stat(as.numeric(level))),
    breaks = breaks, colour = NA
  ) +
  geom_contour_filled(
    data   = df_nsig,
    aes(month_date, window_month,
        z    = correlation,
        fill = after_stat(as.numeric(level))),
    breaks = breaks, fill = "white", alpha = 0.5, colour = NA
  ) +
  geom_contour(
    data   = df,
    aes(month_date, window_month, z = correlation),
    colour = "grey70", size = 0.15,
    breaks = breaks
  ) +
  geom_point(
    data = best_pts,
    aes(month_date, window_month),
    size = 2, shape = 21, stroke = .8,
    fill = "yellow", colour = "black"
  ) +
  geom_label_repel(
    data = best_pts,
    aes(month_date, window_month,
        label = sprintf("r = %.2f", correlation)),
    size = 3,
    fill = scales::alpha("white", 0.3),
    colour = "black",
    label.size = 0,
    label.r    = unit(0.1, "lines"),
    segment.color = NA,
    box.padding   = 0.3,
    point.padding = 0.2
  ) +
  #(~ VI, ncol = 2) +                # one column, two rows
  scale_x_date(
    date_breaks = "2 month",
    labels      = date_format("%b", locale = "en"),   # English month abbrev.
    expand      = c(0, 0)
  ) +
  scale_y_continuous(breaks = y_grid, expand = c(0,0)) +
  geom_vline(xintercept = x_grid, colour = "grey50", linetype = "dotted") +
  geom_hline(yintercept = y_grid, colour = "grey50", linetype = "dotted") +
  scale_fill_stepsn(
    colours = pal,
    limits  = c(1, nbins),
    breaks  = seq_len(nbins),
    labels  = lab_vals,
    name    = expression("Pearson "~italic(r)),
    guide   = guide_colourbar(barheight = unit(5,"cm"),
                              barwidth  = unit(0.5,"cm"),
                              ticks.colour = "black")
  ) +
  labs(
    title = "kNDVI - TRI Rolling Correlation — Plot SF01",
    x     = "Month (prev year → current year)",
    y     = "Window length (months)"
  ) +
  theme_minimal(base_size = 11) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5),
        strip.text  = element_text(face = "bold"))

print(p)

ggsave(
  filename = "figures/final_figures/kNDVI_SF01.png",
  plot     = p,
  width    = 6, height = 4, dpi = 300, bg = "white"
)
#######################################
library(dplyr)
library(ggplot2)
library(ggrepel)
library(readxl)
library(RColorBrewer)

sig_thr <- 0.05                       # p-value threshold
vi_set  <- c("kNDVI", "GVMI")         # two indices
plot_id <- "SF18"                     # one plot

# ───────────────────────── 1. DATA — one plot, two VIs ──────────────────────
df <- heat_all %>% 
  filter(plot == plot_id,
         VI %in% vi_set,
         !is.na(correlation)) %>% 
  mutate(
    sig          = p_value < sig_thr,
    ext_doy      = if_else(lag == 1, doy, doy + 365),
    month_date   = as.Date(ext_doy - 1, origin = "2000-01-01"),
    window_month = window * 8 / 30.437,
    VI           = factor(VI, levels = vi_set), 
    panel_lab    = paste(VI),
  )

df_nsig <- filter(df, !sig)
df_sig  <- filter(df,  sig)

best_pts <- df_sig %>%                      # one “best” per panel
  group_by(panel_lab) %>% 
  slice_max(abs(correlation), n = 1, with_ties = FALSE) %>% 
  ungroup()

# ───────────────────────── 2. COLOUR SCALE (same bins) ──────────────────────
rng    <- range(df$correlation, na.rm = TRUE)
breaks <- seq(floor(rng[1]/0.1)*0.1, ceiling(rng[2]/0.1)*0.1, 0.1)
nbins  <- length(breaks) - 1
pal    <- colorRampPalette(rev(brewer.pal(11, "RdYlBu")))(nbins)

lab_vals <- sprintf("%.2f",
                    (head(breaks, -1) + tail(breaks, -1)) / 2)
lab_vals[lab_vals == "-0.00"] <- "0.00"

x_grid <- seq(min(df$month_date), max(df$month_date), by = "2 month")
y_grid <- seq(0, ceiling(max(df$window_month)), by = 2)

# ───────────────────────── 3. PLOT ───────────────────────────────────────────
p <- ggplot() +
  geom_contour_filled(
    data   = df,
    aes(month_date, window_month,
        z    = correlation,
        fill = after_stat(as.numeric(level))),
    breaks = breaks, colour = NA
  ) +
  geom_contour_filled(
    data   = df_nsig,
    aes(month_date, window_month,
        z    = correlation,
        fill = after_stat(as.numeric(level))),
    breaks = breaks, fill = "white", alpha = 0.5, colour = NA
  ) +
  geom_contour(
    data   = df,
    aes(month_date, window_month, z = correlation),
    colour = "grey70", size = 0.15,
    breaks = breaks
  ) +
  geom_point(
    data = best_pts,
    aes(month_date, window_month),
    size = 2, shape = 21, stroke = .8,
    fill = "yellow", colour = "black"
  ) +
  geom_label_repel(
    data = best_pts,
    aes(month_date, window_month,
        label = sprintf("r = %.2f", correlation)),
    size = 3,
    fill = scales::alpha("white", 0.3),
    colour = "black",
    label.size = 0,
    label.r    = unit(0.1, "lines"),
    segment.color = NA,
    box.padding   = 0.3,
    point.padding = 0.2
  ) +
  facet_grid(~ VI) +                # one column, two rows
  scale_x_date(
    date_breaks = "2 month",
    labels      = date_format("%b", locale = "en"),   # English month abbrev.
    expand      = c(0, 0)
  ) +
  scale_y_continuous(breaks = y_grid, expand = c(0,0)) +
  geom_vline(xintercept = x_grid, colour = "grey50", linetype = "dotted") +
  geom_hline(yintercept = y_grid, colour = "grey50", linetype = "dotted") +
  scale_fill_stepsn(
    colours = pal,
    limits  = c(1, nbins),
    breaks  = seq_len(nbins),
    labels  = lab_vals,
    name    = expression("Pearson "~italic(r)),
    guide   = guide_colourbar(barheight = unit(5,"cm"),
                              barwidth  = unit(0.5,"cm"),
                              ticks.colour = "black")
  ) +
  labs(
    title = "kNDVI & GVMI Rolling Correlation with TRI — Plot SF18 (FASY-Kellerwald)",
    x     = "Month (prev year → current year)",
    y     = "Window length (months)"
  ) +
  theme_minimal(base_size = 11) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5),
        strip.text  = element_text(face = "bold"))

print(p)

ggsave(
  filename = "figures/final_figures/SF18_kNDVI_GVMI_grid.png",
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
  dplyr::filter(tolower(VI) == "kndvi")                                       %>% 
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

library(ggplot2)
# 1. filggplot2# 1. filter to only significant correlations
heat_sig <- heat_all %>% 
  filter(p_value < 0.05)%>%
  mutate(window_days = window * 8)

# 2. plot a histogram of window‐lengths, faceted by VI
ggplot(heat_sig, aes(x = window_days)) +
  geom_histogram(binwidth = 5, fill = "gray80", color = "black") +
  facet_wrap(~ VI, scales = "free_y", ncol = 3) +
  labs(
    x = "Accumulation window (days)",
    y = "Count of significant correlations",
    title = "Distribution of Significant VI–TRI Correlation Windows",
    subtitle = "Each panel shows the window‐lengths for one vegetation index"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    strip.text = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )


library(ggplot2)

ggplot(heat_sig, aes(x = window_days, y = correlation)) +
  # 1) bin in 2D, mapping fill to the raw count
  stat_bin2d(
    aes(fill = after_stat(count)),
    binwidth = c(40, 0.1),      # e.g. 40-day × 0.1-corr bins
    color    = "grey70"         # outline colour for each tile
  ) +
  # 2) label each tile with its count
  geom_text(
    stat    = "bin2d",
    aes(label = after_stat(count)),
    size    = 2,
    color   = "black"
  ) +
  # 3) facet by VI
  facet_wrap(~ VI, scales = "free") +
  # 4) a continuous fill scale over count
  scale_fill_viridis_c(
    option = "magma",
    name   = "N observations"
  ) +
  labs(
    x        = "Accumulation window (days)",
    y        = "Correlation",
    title    = "Counts of Correlations per Window×Strength Bin",
    subtitle = "Each tile’s colour & label = number of points in that bin"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    strip.text       = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )
ggplot(heat_sig, aes(x = correlation)) +
  geom_histogram(binwidth = 0.1, fill = "grey80", color = "black") +
  geom_text(
    stat    = "bin",
    aes(label = after_stat(count), y = after_stat(count)),
    vjust   = -0.5,
    size    = 3
  ) +
  facet_wrap(~ VI, scales = "free_y") +
  labs(
    x = "Correlation strength",
    y = "Count",
    title = "Histogram of Significant Correlations by Strength Bin"
  ) +
  theme_minimal()

# ───────────────────── 0.  PACKAGES ───────────────────────────
library(dplyr)
library(tidyr)
library(ggplot2)
library(RColorBrewer)
library(readxl)
library(scales)

# ───────────────────── 1.  SETTINGS & DATA ────────────────────
vi_target <- "kNDVI"                                   # change if desired
plot_set  <- c("SF19","SF18","SF05","SF20",
               "SF14","SF02","PF03","SF15","SF21")

heat_all  <- readRDS("corr_table_vi_tri.rds")
plot_meta <- read_xlsx("data/vector_data/plot_metadata.xlsx")
heat_all  <- left_join(heat_all, plot_meta, by = "plot")

# ───────────────────── 2.  MEDIAN ACROSS 9 PLOTS ─────────────
df_med <- heat_all %>%
  filter(plot %in% plot_set,
         tolower(VI) == tolower(vi_target),
         !is.na(correlation)) %>%
  group_by(window, doy, lag) %>%                           # one row per combo
  summarise(med_r = median(correlation), .groups = "drop") %>%
  mutate(
    ext_doy      = if_else(lag == 1, doy, doy + 365),
    month_date   = as.Date(ext_doy - 1, origin = "2000-01-01"),
    window_month = window * 8 / 30.437
  )

# ───────────────────── 3.  COLOUR BINS (same rules) ──────────
rng    <- range(df_med$med_r, na.rm = TRUE)
breaks <- seq(floor(rng[1]/0.1)*0.1, ceiling(rng[2]/0.1)*0.1, 0.1)
nbins  <- length(breaks) - 1
pal    <- colorRampPalette(rev(brewer.pal(11,"RdYlBu")))(nbins)

mid_vals <- (head(breaks,-1) + tail(breaks,-1)) / 2
lab_vals <- sprintf("%.2f", mid_vals)
lab_vals[lab_vals == "-0.00"] <- "0.00"

x_grid <- seq(min(df_med$month_date), max(df_med$month_date), by = "2 month")
y_grid <- seq(0, ceiling(max(df_med$window_month)), by = 2)

# ───────────────────── 4.  HEAT-MAP ──────────────────────────
p1 <- ggplot(df_med,
       aes(month_date, window_month, z = med_r,
           fill = after_stat(as.numeric(level)))) +
  geom_contour_filled(breaks = breaks, colour = NA) +
  geom_contour(breaks = breaks, colour = "grey70", size = 0.15) +
  scale_x_date(
    date_breaks = "2 month",
    labels      = date_format("%b", locale = "en"),   # English month abbrev.
    expand      = c(0, 0)
  ) +
  scale_y_continuous(breaks = y_grid, expand = c(0,0)) +
  geom_vline(xintercept = x_grid, colour = "grey50", linetype = "dotted") +
  geom_hline(yintercept = y_grid, colour = "grey50", linetype = "dotted") +
  scale_fill_stepsn(
    colours = pal,
    limits  = c(1, nbins),
    breaks  = seq_len(nbins),
    labels  = lab_vals,
    name    = expression("Median "~italic(r)),
    guide   = guide_colourbar(barheight = unit(5,"cm"),
                              barwidth  = unit(0.5,"cm"),
                              ticks.colour = "black")
  ) +
  labs(title = paste(vi_target, "– TRI Median Rolling Correlation of the FASY Cluster"),
       x = "Month (previous → current year)",
       y = "Window length (months)") +
  theme_minimal(base_size = 11) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5))
p1

ggsave(
  filename = "figures/final_figures/kNDVI_FASY_Cluster.png",
  plot     = p1,
  width    = 6, height = 4, dpi = 300, bg = "white"
)

