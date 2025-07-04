library(tidyverse)
library(lme4)      # mixed models
library(mgcv)      # non-linear smoothers if you need them
library(randomForest)  # ML comparison
pheno <- read_csv("data/analysis_ready_data/phenology_df.csv")

# Average phenology across years per plot
pheno_mean <- pheno %>%
  group_by(plot) %>%
  summarise(
    SOS_doy = mean(SOS_doy, na.rm = TRUE),
    EOS_doy = mean(EOS_doy, na.rm = TRUE)
  )

# Join into heatmap_df
heatmap_df_test <- heatmap_df %>%
  left_join(pheno_mean, by = "plot")

##  a)  Fisher z-transform so residuals are ~Gaussian
heatmap_df_test <- heatmap_df_test %>% 
  mutate(z_r  = atanh(correlation),          # Fisher-z
         absr = abs(correlation),            # handy for classification
         n    = 23,                          # years in each correlation
         weight = n - 3)                     # var(z) ≈ 1/(n-3)

##  b)  Window attributes
heatmap_df_test <- heatmap_df_test %>% 
  mutate(win_len   = window,
         win_mid   = doy - (window/2),       # centre DOY of the window
         phase_rel = case_when(              # phenological phase
           win_mid < SOS_doy        ~ "pre-green",
           win_mid <= EOS_doy       ~ "in-season",
           TRUE                     ~ "post-senescence")
  )

##  c)  Static topo-site predictors
topo   <- read_csv("data/analysis_ready_data/topo_ai_df.csv")  # DEM, slope, aspect
clim <- read_csv("data/analysis_ready_data/climate_df.csv")

heatmap_df_test <- heatmap_df_test %>% 
  left_join(topo,   by = "plot")

## 1A ─ add year + mid-month DOY  ----------------------------------------
clim_month <- clim %>%                                   # 9 744 rows
  mutate(
    date_mid = date + days(14),            # 1st → ≈15th (optional)
    year     = year(date_mid),
    doy_mid  = yday(date_mid),
    doy      = ((doy_mid - 1) %/% 8) * 8 + 1   # 1,9,17,25,… ≤ DOY_mid
  ) %>% 
  select(plot, year, doy, temp, prec, evapo, soil_moist) %>% 
  distinct()                                   # one row per month

## 1B ─ build an 8-day grid & fill  (no warnings) ------------------------
## ── 2. make 8-day grid & forward-fill inside each plot-year ───────────
clim_8d <- clim_month %>% 
  group_by(plot, year) %>% 
  group_modify(~{
    end <- if (leap_year(.y$year)) 366 else 365
    grid <- tibble(doy = seq(1, end, by = 8))
    
    grid %>% 
      left_join(.x, by = "doy") %>%            # now 12 rows line up
      arrange(doy) %>% 
      fill(temp, prec, evapo, soil_moist,      # value persists *until next month*
           .direction = "downup")
  }) %>% 
  ungroup()


library(zoo)
library(purrr)

accum_windows <- 1:91              # 1 × 8 d  …  91 × 8 d

## 2A ─ long→wide helper  -----------------------------------------------
roll_one_window <- function(w) {
  clim_8d %>%
    group_by(plot, year) %>%
    arrange(doy) %>%
    mutate(
      temp  = rollapply(temp,  w, mean, align = "right", fill = NA, na.rm = TRUE),
      prec  = rollapply(prec,  w, mean, align = "right", fill = NA, na.rm = TRUE),
      evapo = rollapply(evapo, w, mean, align = "right", fill = NA, na.rm = TRUE),
      soil_moist = rollapply(soil_moist, w, mean, align = "right", fill = NA, na.rm = TRUE)
    ) %>%
    ungroup() %>%
    mutate(window = w)
}

## 2B ─ loop over all windows and bind rows ------------------------------
clim_roll <- map_dfr(accum_windows, roll_one_window)

clim_means <- clim_roll %>% 
  group_by(plot, window, doy) %>%      # *no* year here
  summarise(
    mean_temp       = mean(temp,       na.rm = TRUE),
    mean_prec       = mean(prec,       na.rm = TRUE),
    mean_evapo      = mean(evapo,      na.rm = TRUE),
    mean_soil_moist = mean(soil_moist, na.rm = TRUE),
    .groups = "drop"
  )

heatmap_df_test <- heatmap_df_test %>%
  left_join(as_tibble(clim_means), by = c("plot", "window", "doy"))


ggplot(heatmap_df_test,
       aes(win_mid, win_len, fill = z_r)) +
  geom_tile() +
  scale_fill_gradient2(limits = c(-2,2), name = "Fisher z") +
  facet_wrap(~ VI, nrow = 2) +
  labs(x = "Window mid-DOY", y = "Window length (×8 d)") +
  theme_classic()

ggplot(heatmap_df_test,
       aes(phase_rel, absr, colour = phase_rel)) +
  geom_boxplot() + facet_wrap(~ VI)
cor.test(heatmap_df_test$absr, heatmap_df_test$mean_prec)


pred_vars <- c("win_len", "win_mid",
               "mean_temp", "mean_prec", "mean_evapo", "mean_soil_moist",
               "dem_mean", "slope_mean", "aspect_mean", "ai_mean")

heat_df <- heatmap_df_test %>%
  mutate(across(all_of(pred_vars), scale),                 # z-score
         VI          = factor(VI),
         forest_type = factor(forest_type))

###############################################################################
# 3.  Fit the maximal fixed-effects model  -------------------------------------
#    • Random intercept for plot (23 yrs collapsed ⇒ very few higher-level units)
#    • Weights = (n-3) as in your earlier pipeline
###############################################################################
m_clean <- lmer(
  absr ~ win_len + win_mid +
    mean_evapo +                     # keep the strongest water-stress term
    dem_mean + slope_mean + ai_mean +# topo/climate indices, low VIF
    (1 + win_len | plot),            # random intercept & slope
  data    = heat_df,
  weights = weight,
  control = lmerControl(optimizer = "bobyqa",
                        optCtrl   = list(maxfun = 1e6)),
  REML    = FALSE
)

summary(m_full)

###############################################################################
# 4.  Basic diagnostics  -------------------------------------------------------
# 4a. multicollinearity check
library(performance)# R², ICC, VIF helpers

check_collinearity(m_full)          # VIFs; > 5 means trouble

# 4b. marginal vs conditional R² (fixed-only vs fixed+random)
r2_full <- r2(m_full)
print(r2_full)

# 4c. heteroscedasticity plot
plot(m_full, sqrt(abs(resid(.))) ~ fitted(.), type = c("p","smooth"))

###############################################################################
# 5.  Optional random-slope upgrade  ------------------------------------------
#    Add a random slope for window length, which is often plot-specific
###############################################################################
m_rs <- update(m_full, . ~ . + (win_len | plot))
anova(m_full, m_rs)     # if p < 0.05 keep the extra complexity

###############################################################################
# 6.  Optional term trimming  --------------------------------------------------
#    Drop non-significant climate vars en bloc and compare AIC
###############################################################################
m_trim <- update(m_full, . ~ . - mean_evapo - mean_soil_moist)
anova(m_full, m_trim)   # ΔAIC, χ², p-value

sjPlot::plot_model(m_clean, type="pred", terms=c("win_len","mean_evapo"))
