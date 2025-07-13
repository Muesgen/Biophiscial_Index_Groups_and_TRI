################################################################################
## Dendro-climate workflow:                                                   ##
##  • read phenology, monthly climate, and tree-ring (TRI) data               ##
##  • derive growing-season drought & climate indices                         ##
##  • join with TRI and compute plot-level correlations (ρ + p-values)        ##
################################################################################

## ────────────────────────────────────────────────────────────────────────────
## 0.  Libraries
## ────────────────────────────────────────────────────────────────────────────
library(dplyr)       # data wrangling
library(lubridate)   # dates (year(), yday(), ...)
library(tidyr)       # unnest_wider()
library(purrr)       # map-style helpers (implicit in across())
library(broom)       # tidying cor.test()
library(SPEI)        # Standardised Precipitation-Evapotranspiration Index
# install.packages("SPEI")   # run once if missing
library(ggplot2)     # <-- optional plotting at the end

## ────────────────────────────────────────────────────────────────────────────
## 1.  Read data
## ────────────────────────────────────────────────────────────────────────────
pheno    <- read.csv("data/analysis_ready_data/phenology_df.csv")

climate  <- read.csv("data/analysis_ready_data/climate_df.csv")

TRI      <- read.csv("data/analysis_ready_data/TRI_chronologies.csv")

heat_all <- readRDS("corr_table_vi_tri.rds")

## ────────────────────────────────────────────────────────────────────────────
## 2.  Add monthly water-balance & SPEI
## ────────────────────────────────────────────────────────────────────────────
climate_wb <- climate %>%
  mutate(
    wb   = prec - evapo,                   # monthly climatic water balance
    year = year(date),
    doy  = yday(date) + 14                 # mid-month day-of-year
  ) %>%
  group_by(plot) %>%
  arrange(date) %>%
  mutate(                                   # rolling SPEI inside each plot
    spei3 = as.numeric(spei(wb,  3)$fitted),
    spei6 = as.numeric(spei(wb,  6)$fitted)
  ) %>%
  ungroup()

## ────────────────────────────────────────────────────────────────────────────
## 3.  Restrict to growing-season months (SOS–EOS per year & plot)
## ────────────────────────────────────────────────────────────────────────────
gs_climate <- climate_wb %>%
  inner_join(pheno, by = c("plot", "year")) %>%   # add SOS/EOS
  filter(doy >= SOS_doy, doy <= EOS_doy)          # keep only GS months

## ────────────────────────────────────────────────────────────────────────────
## 4.  Annual growing-season aggregates & drought indices
## ────────────────────────────────────────────────────────────────────────────
clim_annual <- gs_climate %>%
  group_by(plot, year) %>%
  summarise(
    gs_temp   = mean(temp,  na.rm = TRUE),         # °C
    gs_prec   = sum(prec,  na.rm = TRUE),          # mm
    gs_pet    = sum(evapo, na.rm = TRUE),          # mm
    gs_cwb    = gs_prec - gs_pet,                  # mm  (P − PET)
    gs_ppet   = gs_prec / gs_pet,                  # P : PET ratio
    gs_cmi    = (gs_prec - gs_pet) / gs_pet,       # Climatic Moisture Index
    gs_demart = gs_prec / (gs_temp + 10),          # De Martonne aridity
    gs_soil   = mean(soil_moist, na.rm = TRUE),    # %
    gs_spei3  = mean(spei3, na.rm = TRUE),         # dimensionless
    gs_spei6  = mean(spei6, na.rm = TRUE),
    .groups   = "drop"
  ) %>%
  group_by(plot) %>%
  mutate(gs_cwb_z = (gs_cwb - mean(gs_cwb)) / sd(gs_cwb)) %>%  # z-score
  ungroup()

## ────────────────────────────────────────────────────────────────────────────
## 5.  Join with tree-ring data (TRI)                                         ##
##      TRI must contain at least: plot, tree, year, TRI                      ##
## ────────────────────────────────────────────────────────────────────────────
tri_clim <- TRI %>%
  inner_join(clim_annual, by = c("plot", "year")) %>%
  na.omit()

## ────────────────────────────────────────────────────────────────────────────
## 6.  Plot-level correlations ρ & p using reframe() (no list columns)
## ────────────────────────────────────────────────────────────────────────────

vars <- c("gs_temp", "gs_prec", "gs_pet", "gs_cwb", "gs_cwb_z",
          "gs_ppet", "gs_cmi", "gs_demart",
          "gs_spei3", "gs_spei6", "gs_soil")

corr_tbl <- tri_clim %>%
  group_by(plot) %>%
  reframe(
    ## correlation estimates (ρ)
    across(
      all_of(vars),
      ~ cor(.x, TRI, use = "pairwise.complete.obs"),
      .names = "r_{.col}"
    ),
    ## p-values
    across(
      all_of(vars),
      ~ cor.test(.x, TRI, use = "pairwise.complete.obs")$p.value,
      .names = "p_{.col}"
    )
  )

print(corr_tbl)

## ────────────────────────────────────────────────────────────────────────────
## 9.  Median of the top-10 VI–TRI correlations  (per plot & VI)
## ────────────────────────────────────────────────────────────────────────────
best_vi <- heat_all %>%                                 # has: plot, VI, correlation, …
  filter(!is.na(correlation)) %>%
  group_by(plot, VI) %>%
  slice_max(abs(correlation), n = 10, with_ties = FALSE) %>%
  summarise(median_vi_r = median(abs(correlation)), .groups = "drop")

## ────────────────────────────────────────────────────────────────────────────
## 10.  Climate–TRI correlations in long form (absolute value)
## ────────────────────────────────────────────────────────────────────────────
clim_long <- corr_tbl %>%                               # has r_<index> columns
  select(plot, starts_with("r_")) %>%
  pivot_longer(-plot,
               names_to  = "clim_index",
               values_to = "clim_r") %>%
  mutate(
    clim_index = sub("^r_", "", clim_index),
    clim_r     = abs(clim_r)                            # strength only
  )

## ────────────────────────────────────────────────────────────────────────────
## 11.  Test every VI × climate-index combination
## ────────────────────────────────────────────────────────────────────────────
combo_strength <- best_vi %>%
  inner_join(clim_long, by = "plot") %>%                # add climate r per plot
  group_by(VI, clim_index) %>%
  summarise(
    rho          = cor(median_vi_r, clim_r,
                       use = "pairwise.complete.obs"),  # Pearson across plots
    n_plots      = n_distinct(plot),
    .groups      = "drop"
  ) %>%
  mutate(abs_rho = abs(rho)) %>%                        # correlation strength
  arrange(desc(abs_rho))

print(combo_strength)

## ────────────────────────────────────────────────────────────────────────────
## 12.  Which combo is best?
## ────────────────────────────────────────────────────────────────────────────
best_pair <- combo_strength %>% slice_max(abs_rho, n = 1)

cat("\nTop combination:\n")
print(best_pair)


## ────────────────────────────────────────────────────────────────────────────
## 0-bis.  Annual Aridity Index (AI = P / PET)
## ────────────────────────────────────────────────────────────────────────────

ai_tbl <- climate %>%                               # uses full-year climate
  mutate(year = year(date)) %>%                     
  group_by(plot, year) %>%
  summarise(
    P_ann   = sum(prec,  na.rm = TRUE),             # mm
    PET_ann = sum(evapo, na.rm = TRUE),             # mm
    AI      = P_ann / PET_ann,
    .groups = "drop"
  ) %>%                                             # long-term mean AI per plot
  group_by(plot) %>%
  summarise(
    ai_mean = mean(AI, na.rm = TRUE),               # our dryness proxy
    ai_sd   = sd(AI,   na.rm = TRUE),
    .groups = "drop"
  )

## ────────────────────────────────────────────────────────────────────────────
## 13.  Does dryness strengthen VI–TRI coupling?
## ────────────────────────────────────────────────────────────────────────────

vi_ai <- best_vi %>%              # from step 9
  inner_join(ai_tbl, by = "plot") # add long-term AI

## A) one overall test across all VIs
overall <- cor.test(
  vi_ai$ai_mean, vi_ai$median_vi_r, use = "pairwise.complete.obs"
)

## B) or per VI, if you want to see whether the pattern holds for each index
vi_level <- vi_ai %>%
  group_by(VI) %>%
  summarise(
    r_ai   = cor(ai_mean, median_vi_r, use = "pairwise.complete.obs"),
    p_ai   = cor.test(ai_mean, median_vi_r)$p.value,
    n_plots = n(),
    .groups = "drop"
  ) %>%
  arrange(r_ai)                   # drier (= lower AI) → stronger if r_ai is NEGATIVE

print(overall)
print(vi_level)

ggplot(vi_ai, aes(ai_mean, median_vi_r, colour = VI)) +
  geom_point(size = 3) +
  geom_smooth(method = "lm", se = FALSE) +
  scale_x_reverse() +                                 # left = dry, right = wet
  labs(
    x = "Long-term Aridity Index (P / PET)",
    y = "Median |VI – TRI| correlation (ρ)",
    title = "Are drier plots more tightly coupled?"
  ) +
  theme_classic(base_size = 13)

library(ggplot2)

ggplot(vi_ai, aes(ai_mean, median_vi_r)) +
  geom_point(size = 3, colour = "steelblue") +         # same colour everywhere
  geom_smooth(method = "lm", se = FALSE, colour = "black") +
  scale_x_reverse() +                                  # left = dry, right = wet
  facet_wrap(~ VI, ncol = 2, scales = "free_y") +      # one panel per VI
  labs(
    x = "Long-term Aridity Index (P / PET)",
    y = "Median |VI – TRI| correlation (ρ)",
    title = "Is VI ↔ TRI coupling stronger in drier plots?",
    subtitle = "One panel per vegetation index"
  ) +
  theme_classic(base_size = 13)
