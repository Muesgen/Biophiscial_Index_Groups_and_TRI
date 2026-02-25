################################################################################
# Script: Site Mean Curves
# Date: 25.02.206
#
# This script processes tree ring width (TRW) data:
# 1. Loads required libraries and functions.
# 2. Loads detrended site mean vures from extract_TRW.R function.
# 3. Detrends TRW using several methods and computes ring-width indices (RWI).
# 4. Creates species Specific plots for Site Mean Curves 
################################################################################

library(dplyr)
library(ggplot2)
library(scales)
library(viridis)
library(tidyr)
library(ggrepel)

SMC <- read.csv2("~/data/analysis_ready_data/TRI_chronologies.csv", sep=",")

SMC <- SMC %>%
  mutate(Species = case_when(
    plot %in% c("PF01") ~ "PCAB",
    plot %in% c("PF02","SF12","SF17") ~ "QUPE",
    plot %in% c("PF03","SF02","SF04","SF05","SF08","SF14","SF18","SF19","SF20") ~ "FASY",
    plot %in% c("SF01","SF09") ~ "QURO",
    plot %in% c("SF03","SF06","SF10","SF16") ~ "PISY",
    plot %in% c("SF07") ~ "PINI",
    plot %in% c("SF15","SF21","SF11") ~ "PSME",
    plot %in% c("SF13") ~ "LADE",
    TRUE ~ NA_character_
  ))
head(SMC)

# Clean and Filter
SMC_00_22 <- SMC %>%
  mutate(
    year = as.integer(year),
    TRI  = as.numeric(as.character(TRI)),
    Species = factor(Species),
    plot    = factor(plot)
  ) %>%
  filter(year >= 2000, year <= 2022) %>%
  drop_na(year, TRI, Species, plot) %>%
  arrange(Species, plot, year)

# Label options
end_pts <- SMC_00_22 %>%
  group_by(Species, plot) %>%
  slice_max(order_by = year, n = 1, with_ties = FALSE) %>%
  ungroup()

x_max <- max(SMC_00_22$year, na.rm = TRUE)

# Create Plot
p <- ggplot(SMC_00_22, aes(x = year, y = TRI, group = plot, color = plot)) +
  
  # White label area on the right
  annotate("rect",
           xmin = x_max + 0.2, xmax = x_max + 2,
           ymin = -Inf, ymax = Inf,
           fill = "white", color = NA) +
  geom_line(linewidth = 0.55, alpha = 0.9) +
  
  geom_text_repel(
    data = end_pts,
    aes(label = plot, color = plot),
    nudge_x = 0.8,
    direction = "y",
    hjust = 0,          # linke Textkante = Ankerpunkt
    vjust = 0.5,        # vertikal zentriert
    size = 3,
    min.segment.length = 0,
    segment.size = 0.25,
    segment.color = "black",
    box.padding = 0.15,   # kleiner → Linie trifft Mitte
    point.padding = 0.15,
    max.overlaps = Inf,
    show.legend = FALSE
  )+
  
  facet_wrap(~ Species, ncol = 2, scales = "fixed") +
  
  # Scientific color palette
  scale_color_viridis_d(end = 0.95, guide = "none") +
  
  scale_x_continuous(
    limits = c(2000, x_max + 2),
    breaks = seq(2000, 2022, by = 5),
    minor_breaks = seq(2000, 2022, by = 1),
    expand = expansion(mult = c(0.01, 0.01))
  ) +
  scale_y_continuous(
    breaks = pretty_breaks(n = 5),
    expand = expansion(mult = c(0.05, 0.05))
  ) +
  
  labs(
    x = "Year",
    y = "Site Mean Curves - (TRI)"
  ) +
  
  coord_cartesian(clip = "off") +
  
  theme_classic(base_size = 11) +
  theme(
    strip.background = element_rect(fill = "grey95", color = "grey70", linewidth = 0.3),
    strip.text       = element_text(face = "bold"),
    panel.grid.major = element_line(color = "grey90", linewidth = 0.25),
    panel.grid.minor = element_line(color = "grey95", linewidth = 0.2),
    axis.title       = element_text(face = "bold", size = 12),
    axis.line        = element_line(color = "black", linewidth = 0.4),
    axis.ticks       = element_line(color = "black", linewidth = 0.3),
    plot.margin      = margin(t = 6, r = 40, b = 6, l = 6)
  )

print(p)

# Save Plot

ggsave(
  filename = "SMC_TRI_2000_2022_bySpecies.tiff",
  plot = p,
  width = 360, height = 320, units = "mm",
  dpi = 800, compression = "lzw"
)
ggsave(
  filename = "SMC_TRI_2000_2022_bySpecies.png",
  plot = p,
  width = 360, height = 320, units = "mm",
  dpi = 800, bg = "white"
)

