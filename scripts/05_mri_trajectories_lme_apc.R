source(here::here("config.R"))

####Calculating rate of MRI longitudianl changes

library(lme4)
library(lmerTest)
library(tidyverse)

# ICV-normalize all bilateral volumes
lme_data <- analytical_cohort %>%
  mutate(
    across(c(feature_map$left, feature_map$right, "ST149SV"),
           ~ . / ST10CV * 100,
           .names = "{.col}_norm")
  ) %>%
  # Total Cortical Vol = CorticalGM + CorticalWM, then normalize
  mutate(TotalCortVol_norm = (ST149SV + ST152SV) / ST10CV * 100)

cat("Normalized columns created:", sum(grepl("_norm$|TotalCort", names(lme_data))), "\n")
cat("Rows:", nrow(lme_data), "| Participants:", n_distinct(lme_data$RID), "\n")

#####
within_12m <- analytical_cohort %>%
  filter(time_years <= 1.0)

cat("Data within first 12 months:\n")
cat("  Total observations:", nrow(within_12m), "\n")
cat("  Participants:", n_distinct(within_12m$RID), "\n\n")

# Visits per participant within 12 months
visits_12m <- within_12m %>%
  count(RID, name = "n_visits_12m")

cat("Visits per participant within 12 months:\n")
print(table(visits_12m$n_visits_12m))

cat("\nParticipants with < 2 visits in 12 months:", 
    sum(visits_12m$n_visits_12m < 2), "\n")
cat("Participants with >= 2 visits in 12 months:", 
    sum(visits_12m$n_visits_12m >= 2), "\n")
cat("Participants with >= 3 visits in 12 months:", 
    sum(visits_12m$n_visits_12m >= 3), "\n")

for (window in c(1.0, 1.5, 2.0)) {
  within_window <- analytical_cohort %>%
    filter(time_years <= window) %>%
    count(RID, name = "n_visits")
  
  cat("Window:", window * 12, "months\n")
  cat("  Participants with data:", nrow(within_window), "\n")
  cat("  >= 2 visits:", sum(within_window$n_visits >= 2), 
      "(", round(sum(within_window$n_visits >= 2)/479*100, 1), "%)\n")
  cat("  >= 3 visits:", sum(within_window$n_visits >= 3),
      "(", round(sum(within_window$n_visits >= 3)/479*100, 1), "%)\n")
  cat("  Visit distribution:\n")
  print(table(within_window$n_visits))
  cat("\n")
}

####### check no of converters before 12, 18 and 24 months
early_18m <- conversion_timeline %>%
  filter(trajectory_group %in% c("CN→MCI", "CN→Dementia", "MCI→Dementia")) %>%
  mutate(
    baseline_group = if_else(str_starts(trajectory_group, "CN"), "CN", "MCI"),
    converted_before_18m = conversion_time_years < 1.5
  )

cat("Conversions before 18 months:\n\n")
early_18m %>%
  group_by(baseline_group) %>%
  summarise(
    total_converters = n(),
    before_18m = sum(converted_before_18m),
    after_18m  = sum(!converted_before_18m),
    .groups = "drop"
  ) %>%
  print()

cat("\nDetailed: conversion time distribution for all converters:\n")
early_18m %>%
  group_by(baseline_group) %>%
  summarise(
    n = n(),
    before_12m = sum(conversion_time_years < 1),
    m12_to_18m = sum(conversion_time_years >= 1 & conversion_time_years < 1.5),
    after_18m  = sum(conversion_time_years >= 1.5),
    .groups = "drop"
  ) %>%
  print()


#####lets remove those with single visits within 12 months
# Filter to first 12 months only, exclude participants with <2 visits
lme_12m <- lme_data %>%
  filter(time_years <= 1.0) %>%
  group_by(RID) %>%
  filter(n() >= 2) %>%
  ungroup()

cat("12-month analytical dataset:\n")
cat("  Participants:", n_distinct(lme_12m$RID), "\n")
cat("  Observations:", nrow(lme_12m), "\n")
cat("  Visits per participant:\n")
print(table(lme_12m %>% count(RID) %>% pull(n)))

cat("\nFitting 32 LME models on 12-month window...\n\n")

### LME model 
slope_list_12m <- list()
model_summary_12m <- list()

for (i in seq_along(norm_cols)) {
  region_name <- names(norm_cols)[i]
  col_name    <- norm_cols[i]
  
  formula_str <- paste0(col_name, " ~ time_years + (1 + time_years | RID)")
  
  fit <- tryCatch(
    lmer(as.formula(formula_str), data = lme_12m, REML = FALSE,
         control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 20000))),
    error = function(e) NULL
  )
  
  if (is.null(fit)) {
    cat("  ❌", region_name, "— failed\n")
    next
  }
  
  re <- ranef(fit)$RID
  fe <- fixef(fit)
  
  slopes <- tibble(
    RID = as.numeric(rownames(re)),
    !!paste0("slope_", region_name) := re$time_years + fe["time_years"]
  )
  
  slope_list_12m[[region_name]] <- slopes
  
  model_summary_12m[[region_name]] <- tibble(
    region = region_name,
    fixed_slope = round(fe["time_years"], 4)
  )
  
  cat("  ✅", region_name, "| Fixed slope:", round(fe["time_years"], 4), "%ICV/year\n")
}

cat("\nModels completed:", length(slope_list_12m), "/", length(norm_cols), "\n")

individual_slopes_12m <- slope_list_12m %>%
  reduce(full_join, by = "RID")

cat("Individual slopes dataset (12-month window):\n")
cat("  Participants:", nrow(individual_slopes_12m), "\n")
cat("  Slope columns:", ncol(individual_slopes_12m) - 1, "\n")

########
# Combine all slopes
individual_slopes_12m <- slope_list_12m %>%
  reduce(full_join, by = "RID")

cat("Individual slopes (12-month window):\n")
cat("  Participants:", nrow(individual_slopes_12m), "\n")
cat("  Slope columns:", ncol(individual_slopes_12m) - 1, "\n\n")

# Check which models had singular fits (unreliable random slope estimation)
singular_models <- c("L_Hippocampus", "L_Amygdala", "L_Thalamus", 
                     "L_NucAccumbens", "R_Thalamus")

cat("Models with singular fit (low between-person slope variation):\n")
print(singular_models)

# Summary statistics for all slopes
slope_stats <- individual_slopes_12m %>%
  select(starts_with("slope_")) %>%
  summarise(across(everything(), list(
    mean = ~ round(mean(.x, na.rm = TRUE), 5),
    sd   = ~ round(sd(.x, na.rm = TRUE), 5),
    min  = ~ round(min(.x, na.rm = TRUE), 5),
    max  = ~ round(max(.x, na.rm = TRUE), 5)
  ))) %>%
  pivot_longer(everything(),
               names_to = c("region", "stat"),
               names_sep = "_(?=[^_]+$)") %>%
  pivot_wider(names_from = stat, values_from = value) %>%
  mutate(
    region = str_remove(region, "^slope_"),
    range = max - min,
    singular = region %in% singular_models
  ) %>%
  arrange(desc(sd))

cat("\nSlope quality summary (sorted by between-person variability):\n\n")
print(slope_stats, n = Inf)

#####comparing the feature trajectories incase we have taken 18 months as our window interval
# Converters before 18 months (includes the 23 already excluded at <12m)
converters_before_18m <- conversion_timeline %>%
  filter(trajectory_group %in% c("CN→MCI", "CN→Dementia", "MCI→Dementia"),
         conversion_time_years < 1.5) %>%
  pull(RID)

cat("Converters before 18 months:", length(converters_before_18m), "\n")
cat("  (includes", length(early_converter_rids), "already excluded at <12m)\n")
cat("  Additional exclusions (12-18m):", length(setdiff(converters_before_18m, early_converter_rids)), "\n\n")

# Build 18-month analytical cohort from trajectory_cohort (which already excludes baseline Dementia + MCI→CN)
analytical_18m <- trajectory_cohort %>%
  filter(!RID %in% converters_before_18m)

# Filter to 18-month window, require ≥2 visits
lme_18m_data <- analytical_18m %>%
  filter(time_years <= 1.5) %>%
  group_by(RID) %>%
  filter(n() >= 2) %>%
  ungroup()

# ICV-normalize
lme_18m_data <- lme_18m_data %>%
  mutate(
    across(c(feature_map$left, feature_map$right, "ST149SV"),
           ~ . / ST10CV * 100,
           .names = "{.col}_norm"),
    TotalCortVol_norm = (ST149SV + ST152SV) / ST10CV * 100
  )

cat("18-month cohort:\n")
cat("  Participants:", n_distinct(lme_18m_data$RID), "\n")
cat("  Observations:", nrow(lme_18m_data), "\n")
cat("  Visits per participant:\n")
print(table(lme_18m_data %>% count(RID) %>% pull(n)))

##fitting LME model on the 18 months cohort
cat("\nFitting 32 LME models on 18-month window...\n\n")

slope_list_18m <- list()
model_summary_18m <- list()

for (i in seq_along(norm_cols)) {
  region_name <- names(norm_cols)[i]
  col_name    <- norm_cols[i]
  
  formula_str <- paste0(col_name, " ~ time_years + (1 + time_years | RID)")
  
  fit <- tryCatch(
    lmer(as.formula(formula_str), data = lme_18m_data, REML = FALSE,
         control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 20000))),
    error = function(e) NULL
  )
  
  if (is.null(fit)) {
    cat("  ❌", region_name, "— failed\n")
    next
  }
  
  re <- ranef(fit)$RID
  fe <- fixef(fit)
  
  slopes <- tibble(
    RID = as.numeric(rownames(re)),
    !!paste0("slope_", region_name) := re$time_years + fe["time_years"]
  )
  
  slope_list_18m[[region_name]] <- slopes
  
  model_summary_18m[[region_name]] <- tibble(
    region = region_name,
    fixed_slope = round(fe["time_years"], 4),
    singular = isSingular(fit)
  )
  
  sing <- if (isSingular(fit)) " [SINGULAR]" else ""
  cat("  ✅", region_name, "| Fixed slope:", round(fe["time_years"], 4), sing, "\n")
}

cat("\nModels completed:", length(slope_list_18m), "/", length(norm_cols), "\n")


# Combine 18m slopes
individual_slopes_18m <- slope_list_18m %>% reduce(full_join, by = "RID")

# Build comparison table
compare_12_18 <- bind_rows(
  slope_stats %>% mutate(window = "12m"),
  individual_slopes_18m %>%
    select(starts_with("slope_")) %>%
    summarise(across(everything(), list(
      mean = ~ round(mean(.x, na.rm = TRUE), 5),
      sd   = ~ round(sd(.x, na.rm = TRUE), 5)
    ))) %>%
    pivot_longer(everything(),
                 names_to = c("region", "stat"),
                 names_sep = "_(?=[^_]+$)") %>%
    pivot_wider(names_from = stat, values_from = value) %>%
    mutate(region = str_remove(region, "^slope_"),
           window = "18m")
) %>%
  select(region, window, mean, sd) %>%
  pivot_wider(names_from = window, values_from = c(mean, sd),
              names_glue = "{.value}_{window}") %>%
  mutate(sd_improvement = round((sd_18m - sd_12m) / sd_12m * 100, 1)) %>%
  arrange(desc(sd_18m))

cat("12m vs 18m Slope Comparison (sorted by 18m SD):\n\n")
print(compare_12_18, n = Inf)

###########
# Rebuild trajectory group lookup from conversion_status
traj_lookup <- conversion_status %>%
  mutate(
    traj_group = case_when(
      trajectory_group == "CN-Stable" ~ "CN-Stable",
      trajectory_group %in% c("CN→MCI", "CN→Dementia") ~ "CN-Progressor",
      trajectory_group == "MCI-Stable" ~ "MCI-Stable",
      trajectory_group == "MCI→Dementia" ~ "MCI-Progressor",
      TRUE ~ NA_character_
    )
  ) %>%
  select(RID, traj_group)

# 12-month group breakdown
cat("12-MONTH GROUP BREAKDOWN (N=460):\n")
lme_12m %>%
  distinct(RID) %>%
  left_join(traj_lookup, by = "RID") %>%
  count(traj_group) %>%
  print()

# 18-month group breakdown
cat("\n18-MONTH GROUP BREAKDOWN:\n")
lme_18m_data %>%
  distinct(RID) %>%
  left_join(traj_lookup, by = "RID") %>%
  count(traj_group) %>%
  print()

# Check third ventricle column exists
cat("\nThird Ventricle (ST127SV) in master_visit:", "ST127SV" %in% names(master_visit), "\n")

##########Including Third ventricle as a feature##########
# Normalize third ventricle in both datasets
lme_12m <- lme_12m %>% mutate(ST127SV_norm = ST127SV / ST10CV * 100)
lme_18m_data <- lme_18m_data %>% mutate(ST127SV_norm = ST127SV / ST10CV * 100)

# Fit 12-month LME
fit_3v_12m <- lmer(ST127SV_norm ~ time_years + (1 + time_years | RID),
                   data = lme_12m, REML = FALSE,
                   control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 20000)))

fe_12 <- fixef(fit_3v_12m)
re_12 <- ranef(fit_3v_12m)$RID
sing_12 <- isSingular(fit_3v_12m)

cat("Third Ventricle — 12m: Fixed slope:", round(fe_12["time_years"], 4),
    if (sing_12) "[SINGULAR]" else "", "\n")

# Fit 18-month LME
fit_3v_18m <- lmer(ST127SV_norm ~ time_years + (1 + time_years | RID),
                   data = lme_18m_data, REML = FALSE,
                   control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 20000)))

fe_18 <- fixef(fit_3v_18m)
re_18 <- ranef(fit_3v_18m)$RID
sing_18 <- isSingular(fit_3v_18m)

cat("Third Ventricle — 18m: Fixed slope:", round(fe_18["time_years"], 4),
    if (sing_18) "[SINGULAR]" else "", "\n")

# Add to 12m slope list
slopes_3v_12m <- tibble(
  RID = as.numeric(rownames(re_12)),
  slope_ThirdVentricle = re_12$time_years + fe_12["time_years"]
)
individual_slopes_12m <- individual_slopes_12m %>% left_join(slopes_3v_12m, by = "RID")

# Add to 18m slope list
slopes_3v_18m <- tibble(
  RID = as.numeric(rownames(re_18)),
  slope_ThirdVentricle = re_18$time_years + fe_18["time_years"]
)
individual_slopes_18m <- individual_slopes_18m %>% left_join(slopes_3v_18m, by = "RID")

cat("\nUpdated slope datasets:\n")
cat("  12m: ", nrow(individual_slopes_12m), "participants,", ncol(individual_slopes_12m)-1, "slopes\n")
cat("  18m: ", nrow(individual_slopes_18m), "participants,", ncol(individual_slopes_18m)-1, "slopes\n")


########Slope correlation matrix for 12months#######
slope_matrix <- individual_slopes_12m %>%
  select(starts_with("slope_")) %>%
  rename_with(~ str_remove(., "slope_"))

cor_slopes <- cor(slope_matrix, use = "pairwise.complete.obs")

# Plot
library(corrplot)
png("slope_correlation_12m.png", width = 14, height = 14, units = "in", res = 300)
corrplot(cor_slopes, method = "color", type = "lower",
         tl.col = "black", tl.cex = 0.55, tl.srt = 45,
         addCoef.col = "black", number.cex = 0.3,
         col = colorRampPalette(c("#2166AC", "white", "#B2182B"))(200),
         title = "Slope Correlations — 12-Month Window (N=460)",
         mar = c(0, 0, 2, 0))
dev.off()
cat("Plot saved: slope_correlation_12m.png\n\n")

# Highly correlated pairs (|r| > 0.5)
cat("Slope pairs with |r| > 0.5:\n")
high_cor <- which(abs(cor_slopes) > 0.5 & upper.tri(cor_slopes), arr.ind = TRUE)
if (nrow(high_cor) > 0) {
  for (k in 1:nrow(high_cor)) {
    r1 <- rownames(cor_slopes)[high_cor[k, 1]]
    r2 <- colnames(cor_slopes)[high_cor[k, 2]]
    rval <- round(cor_slopes[high_cor[k, 1], high_cor[k, 2]], 3)
    cat("  ", r1, "↔", r2, ": r =", rval, "\n")
  }
} else {
  cat("  None\n")
}

# L-R slope correlations
cat("\nLeft-Right slope correlations:\n")
lr_cors <- feature_map %>%
  mutate(cor = map2_dbl(
    paste0("slope_L_", region), paste0("slope_R_", region),
    ~ cor(individual_slopes_12m[[.x]], individual_slopes_12m[[.y]], use = "pairwise.complete.obs")
  )) %>%
  select(region, cor) %>%
  arrange(cor)
print(lr_cors, n = Inf)

############run correlation for the 18 months#######
# Add third ventricle to 18m (already done above)
# Get slope matrix
slope_matrix_18m <- individual_slopes_18m %>%
  select(starts_with("slope_")) %>%
  rename_with(~ str_remove(., "slope_"))

cor_slopes_18m <- cor(slope_matrix_18m, use = "pairwise.complete.obs")

# Highly correlated pairs (|r| > 0.5)
cat("18-MONTH: Slope pairs with |r| > 0.5:\n\n")
high_cor_18 <- which(abs(cor_slopes_18m) > 0.5 & upper.tri(cor_slopes_18m), arr.ind = TRUE)
if (nrow(high_cor_18) > 0) {
  for (k in 1:nrow(high_cor_18)) {
    r1 <- rownames(cor_slopes_18m)[high_cor_18[k, 1]]
    r2 <- colnames(cor_slopes_18m)[high_cor_18[k, 2]]
    rval <- round(cor_slopes_18m[high_cor_18[k, 1], high_cor_18[k, 2]], 3)
    cat("  ", r1, "↔", r2, ": r =", rval, "\n")
  }
}

# L-R slope correlations at 18m
cat("\n\n18-MONTH: Left-Right slope correlations:\n\n")
lr_cors_18m <- feature_map %>%
  mutate(cor_18m = map2_dbl(
    paste0("slope_L_", region), paste0("slope_R_", region),
    ~ cor(individual_slopes_18m[[.x]], individual_slopes_18m[[.y]], use = "pairwise.complete.obs")
  )) %>%
  select(region, cor_18m) %>%
  left_join(lr_cors %>% rename(cor_12m = cor), by = "region") %>%
  mutate(change = round(cor_18m - cor_12m, 3)) %>%
  arrange(cor_12m)

cat("Comparison of L-R slope correlations: 12m vs 18m\n\n")
print(lr_cors_18m, n = Inf)

png("slope_correlation_18m.png", width = 14, height = 14, units = "in", res = 300)
corrplot(cor_slopes_18m, method = "color", type = "lower",
         tl.col = "black", tl.cex = 0.55, tl.srt = 45,
         addCoef.col = "black", number.cex = 0.3,
         col = colorRampPalette(c("#2166AC", "white", "#B2182B"))(200),
         title = "Slope Correlations — 18-Month Window (N=459)",
         mar = c(0, 0, 2, 0))

dev.off()
cat("Plot saved: slope_correlation_18m.png\n\n")


##########Final feature selection steps#######
# 12-MONTH: Combined quality assessment

# Slope quality stats
quality_12m <- individual_slopes_12m %>%
  select(starts_with("slope_")) %>%
  summarise(across(everything(), list(
    mean = ~ round(mean(.x, na.rm = TRUE), 5),
    sd   = ~ round(sd(.x, na.rm = TRUE), 5)
  ))) %>%
  pivot_longer(everything(),
               names_to = c("region", "stat"),
               names_sep = "_(?=[^_]+$)") %>%
  pivot_wider(names_from = stat, values_from = value) %>%
  mutate(region = str_remove(region, "^slope_"))

# Singular fit flags (from 12m models)
singular_12m <- c("L_Hippocampus", "L_Amygdala", "L_Thalamus", 
                  "L_NucAccumbens", "R_Thalamus")

# L-R correlations (from earlier)
lr_12m <- lr_cors %>%
  mutate(
    L_region = paste0("L_", region),
    R_region = paste0("R_", region)
  ) %>%
  rename(lr_cor_12m = cor)

# Combine into one summary
summary_12m <- quality_12m %>%
  mutate(
    singular = region %in% singular_12m,
    tier = case_when(
      sd > 0.005  ~ "Tier 1 (High)",
      sd > 0.001  ~ "Tier 2 (Moderate)",
      TRUE        ~ "Tier 3 (Low)"
    )
  ) %>%
  arrange(desc(sd))

cat(strrep("=", 60), "\n")
cat("12-MONTH WINDOW: Slope Quality Summary (N=460)\n")
cat(strrep("=", 60), "\n\n")
print(summary_12m, n = Inf)

cat("\n\nL-R Slope Correlations (12m):\n")
print(lr_cors, n = Inf)

cat("\n\nTier counts:\n")
summary_12m %>% count(tier) %>% print()


#########Build complete summary with anatomical categories#######

# Build complete summary with anatomical categories
full_summary_12m <- summary_12m %>%
  mutate(
    anatomy = case_when(
      region %in% c("L_LatVentricle", "R_LatVentricle", "L_InfLatVentricle",
                    "R_InfLatVentricle", "ThirdVentricle") ~ "Ventricular",
      region %in% c("L_Hippocampus", "R_Hippocampus", "L_Amygdala", "R_Amygdala",
                    "L_Thalamus", "R_Thalamus", "L_Putamen", "R_Putamen",
                    "L_Caudate", "R_Caudate", "L_NucAccumbens", "R_NucAccumbens",
                    "L_Pallidum", "R_Pallidum") ~ "Subcortical",
      region %in% c("L_Entorhinal", "R_Entorhinal", "L_Fusiform", "R_Fusiform",
                    "L_MiddleTemporal", "R_MiddleTemporal", "L_Precuneus", "R_Precuneus",
                    "L_Insula", "R_Insula", "L_LateralOccipital", "R_LateralOccipital") ~ "Cortical",
      region %in% c("CorticalGM", "TotalCortVol") ~ "Global",
      TRUE ~ "Other"
    ),
    # Add L-R correlation where applicable
    base_region = str_remove(region, "^[LR]_"),
    hemisphere = case_when(
      str_starts(region, "L_") ~ "L",
      str_starts(region, "R_") ~ "R",
      TRUE ~ "Bilateral"
    )
  ) %>%
  left_join(
    lr_cors %>% rename(base_region = region, lr_cor = cor),
    by = "base_region"
  ) %>%
  select(anatomy, region, hemisphere, mean, sd, singular, tier, lr_cor) %>%
  arrange(anatomy, desc(sd))

cat("12-MONTH WINDOW: Complete Feature Quality by Anatomical Category\n\n")

for (cat_name in c("Ventricular", "Cortical", "Subcortical", "Global")) {
  cat("\n", strrep("-", 50), "\n")
  cat(cat_name, "Regions:\n")
  cat(strrep("-", 50), "\n")
  full_summary_12m %>%
    filter(anatomy == cat_name) %>%
    print(n = Inf, width = Inf)
}

# Count by tier and anatomy
cat("\n\nFeature counts by Tier and Anatomy:\n")
full_summary_12m %>%
  count(anatomy, tier) %>%
  pivot_wider(names_from = tier, values_from = n, values_fill = 0) %>%
  print()

##########
# 18-month slope quality stats
quality_18m <- individual_slopes_18m %>%
  select(starts_with("slope_")) %>%
  summarise(across(everything(), list(
    mean = ~ round(mean(.x, na.rm = TRUE), 5),
    sd   = ~ round(sd(.x, na.rm = TRUE), 5)
  ))) %>%
  pivot_longer(everything(),
               names_to = c("region", "stat"),
               names_sep = "_(?=[^_]+$)") %>%
  pivot_wider(names_from = stat, values_from = value) %>%
  mutate(region = str_remove(region, "^slope_"))

# 18m singular fits (from model_summary_18m)
singular_18m <- model_summary_18m %>%
  bind_rows() %>%
  filter(singular == TRUE) %>%
  pull(region)

cat("18m singular models:", paste(singular_18m, collapse = ", "), "\n\n")

# 18m L-R correlations
lr_cors_18m_full <- feature_map %>%
  mutate(lr_cor = map2_dbl(
    paste0("slope_L_", region), paste0("slope_R_", region),
    ~ cor(individual_slopes_18m[[.x]], individual_slopes_18m[[.y]], use = "pairwise.complete.obs")
  )) %>%
  select(region, lr_cor)

# Build combined summary
summary_18m <- quality_18m %>%
  mutate(
    singular = region %in% singular_18m,
    tier = case_when(
      sd > 0.005  ~ "Tier 1 (High)",
      sd > 0.001  ~ "Tier 2 (Moderate)",
      TRUE        ~ "Tier 3 (Low)"
    ),
    anatomy = case_when(
      region %in% c("L_LatVentricle", "R_LatVentricle", "L_InfLatVentricle",
                    "R_InfLatVentricle", "ThirdVentricle") ~ "Ventricular",
      region %in% c("L_Hippocampus", "R_Hippocampus", "L_Amygdala", "R_Amygdala",
                    "L_Thalamus", "R_Thalamus", "L_Putamen", "R_Putamen",
                    "L_Caudate", "R_Caudate", "L_NucAccumbens", "R_NucAccumbens",
                    "L_Pallidum", "R_Pallidum") ~ "Subcortical",
      region %in% c("L_Entorhinal", "R_Entorhinal", "L_Fusiform", "R_Fusiform",
                    "L_MiddleTemporal", "R_MiddleTemporal", "L_Precuneus", "R_Precuneus",
                    "L_Insula", "R_Insula", "L_LateralOccipital", "R_LateralOccipital") ~ "Cortical",
      region %in% c("CorticalGM", "TotalCortVol") ~ "Global",
      TRUE ~ "Other"
    ),
    base_region = str_remove(region, "^[LR]_")
  ) %>%
  left_join(lr_cors_18m_full %>% rename(base_region = region, lr_cor = lr_cor), by = "base_region") %>%
  select(anatomy, region, mean, sd, singular, tier, lr_cor) %>%
  arrange(anatomy, desc(sd))

cat("18-MONTH WINDOW: Complete Feature Quality by Anatomical Category (N=459)\n\n")

for (cat_name in c("Ventricular", "Cortical", "Subcortical", "Global")) {
  cat("\n", strrep("-", 50), "\n")
  cat(cat_name, "Regions:\n")
  cat(strrep("-", 50), "\n")
  summary_18m %>%
    filter(anatomy == cat_name) %>%
    print(n = Inf, width = Inf)
}

cat("\n\n18m Tier counts by Anatomy:\n")
summary_18m %>%
  count(anatomy, tier) %>%
  pivot_wider(names_from = tier, values_from = n, values_fill = 0) %>%
  print()

# Compare tiers between windows
tier_compare <- full_summary_12m %>%
  select(region, anatomy, tier_12m = tier, sd_12m = sd, singular_12m = singular) %>%
  left_join(
    summary_18m %>% select(region, tier_18m = tier, sd_18m = sd, singular_18m = singular),
    by = "region"
  ) %>%
  mutate(
    tier_changed = tier_12m != tier_18m,
    sd_change_pct = round((sd_18m - sd_12m) / sd_12m * 100, 1)
  ) %>%
  arrange(anatomy, region)

cat("\n\n12m vs 18m: Features that CHANGED tier:\n\n")
tier_compare %>%
  filter(tier_changed) %>%
  select(anatomy, region, tier_12m, tier_18m, sd_12m, sd_18m, sd_change_pct) %>%
  print(n = Inf, width = Inf)

cat("\n\nFull comparison:\n\n")
tier_compare %>%
  select(anatomy, region, tier_12m, tier_18m, sd_change_pct, singular_12m, singular_18m) %>%
  print(n = Inf, width = Inf)


#################Removing tier 3 features###############
# For 12-month window: features to keep (Tier 1 or Tier 2 at 12m)
keep_12m <- full_summary_12m %>%
  filter(tier %in% c("Tier 1 (High)", "Tier 2 (Moderate)")) %>%
  pull(region)

cat("Features passing Tier 1/2 at 12 months:", length(keep_12m), "\n")
print(sort(keep_12m))

# For 18-month window
keep_18m <- summary_18m %>%
  filter(tier %in% c("Tier 1 (High)", "Tier 2 (Moderate)")) %>%
  pull(region)

cat("\nFeatures passing Tier 1/2 at 18 months:", length(keep_18m), "\n")
print(sort(keep_18m))

# Features that are Tier 1/2 in BOTH windows (most robust)
keep_both <- intersect(keep_12m, keep_18m)
cat("\nFeatures in Tier 1/2 at BOTH windows:", length(keep_both), "\n")
print(sort(keep_both))

# Features that improved from Tier 3 → Tier 2 at 18m (window-dependent)
improved_at_18m <- setdiff(keep_18m, keep_12m)
cat("\nFeatures that improved to Tier 2 only at 18m:\n")
print(improved_at_18m)

# Features that dropped from Tier 2 → Tier 3 at 18m
dropped_at_18m <- setdiff(keep_12m, keep_18m)
cat("\nFeatures that dropped to Tier 3 at 18m:\n")
print(dropped_at_18m)

cat("Tier 1+2 features WITH singular fit at 12 months:\n\n")
full_summary_12m %>%
  filter(tier %in% c("Tier 1 (High)", "Tier 2 (Moderate)"),
         singular == TRUE) %>%
  select(anatomy, region, mean, sd, tier) %>%
  print(n = Inf)

cat("\nTier 1+2 features WITH singular fit at 18 months:\n\n")
summary_18m %>%
  filter(tier %in% c("Tier 1 (High)", "Tier 2 (Moderate)"),
         singular == TRUE) %>%
  select(anatomy, region, mean, sd, tier) %>%
  print(n = Inf)


########correlation among each features in 12 months#######
vent_12m <- full_summary_12m %>%
  filter(anatomy == "Ventricular", tier %in% c("Tier 1 (High)", "Tier 2 (Moderate)")) %>%
  pull(region)

cat("VENTRICULAR (12m):", paste(vent_12m, collapse = ", "), "\n\n")

individual_slopes_12m %>%
  select(all_of(paste0("slope_", vent_12m))) %>%
  rename_with(~ str_remove(., "slope_")) %>%
  cor(use = "pairwise.complete.obs") %>%
  round(3) %>%
  print()

cort_12m <- full_summary_12m %>%
  filter(anatomy == "Cortical", tier %in% c("Tier 1 (High)", "Tier 2 (Moderate)")) %>%
  pull(region)

cat("\nCORTICAL (12m):", length(cort_12m), "features\n")
cat(paste(cort_12m, collapse = ", "), "\n\n")

individual_slopes_12m %>%
  select(all_of(paste0("slope_", cort_12m))) %>%
  rename_with(~ str_remove(., "slope_")) %>%
  cor(use = "pairwise.complete.obs") %>%
  round(3) %>%
  print()

sub_12m <- full_summary_12m %>%
  filter(anatomy == "Subcortical", tier %in% c("Tier 1 (High)", "Tier 2 (Moderate)")) %>%
  pull(region)

cat("\nSUBCORTICAL (12m):", length(sub_12m), "features\n")
cat(paste(sub_12m, collapse = ", "), "\n\n")

individual_slopes_12m %>%
  select(all_of(paste0("slope_", sub_12m))) %>%
  rename_with(~ str_remove(., "slope_")) %>%
  cor(use = "pairwise.complete.obs") %>%
  round(3) %>%
  print()

global_12m <- full_summary_12m %>%
  filter(anatomy == "Global") %>%
  pull(region)

all_kept_12m <- c(vent_12m, cort_12m, sub_12m, global_12m)

cat("\nGLOBAL (12m): CorticalGM and TotalCortVol\n\n")

all_cors_12m <- individual_slopes_12m %>%
  select(all_of(paste0("slope_", all_kept_12m))) %>%
  rename_with(~ str_remove(., "slope_")) %>%
  cor(use = "pairwise.complete.obs")

cat("CorticalGM correlations with other features:\n")
sort(round(all_cors_12m["CorticalGM", ], 3), decreasing = TRUE) %>% print()

cat("\nTotalCortVol correlations with other features:\n")
sort(round(all_cors_12m["TotalCortVol", ], 3), decreasing = TRUE) %>% print()


###########CORRELATION AMONG FEATURES IN 18M#############
# 18-MONTH VENTRICULAR
vent_18m <- summary_18m %>%
  filter(anatomy == "Ventricular", tier %in% c("Tier 1 (High)", "Tier 2 (Moderate)")) %>%
  pull(region)

cat("VENTRICULAR (18m):", paste(vent_18m, collapse = ", "), "\n\n")
individual_slopes_18m %>%
  select(all_of(paste0("slope_", vent_18m))) %>%
  rename_with(~ str_remove(., "slope_")) %>%
  cor(use = "pairwise.complete.obs") %>% round(3) %>% print()

# 18-MONTH CORTICAL
cort_18m <- summary_18m %>%
  filter(anatomy == "Cortical", tier %in% c("Tier 1 (High)", "Tier 2 (Moderate)")) %>%
  pull(region)

cat("\nCORTICAL (18m):", length(cort_18m), "features\n")
cat(paste(cort_18m, collapse = ", "), "\n\n")
individual_slopes_18m %>%
  select(all_of(paste0("slope_", cort_18m))) %>%
  rename_with(~ str_remove(., "slope_")) %>%
  cor(use = "pairwise.complete.obs") %>% round(3) %>% print()

# 18-MONTH SUBCORTICAL
sub_18m <- summary_18m %>%
  filter(anatomy == "Subcortical", tier %in% c("Tier 1 (High)", "Tier 2 (Moderate)")) %>%
  pull(region)

cat("\nSUBCORTICAL (18m):", length(sub_18m), "features\n")
cat(paste(sub_18m, collapse = ", "), "\n\n")
individual_slopes_18m %>%
  select(all_of(paste0("slope_", sub_18m))) %>%
  rename_with(~ str_remove(., "slope_")) %>%
  cor(use = "pairwise.complete.obs") %>% round(3) %>% print()

# 18-MONTH GLOBAL
cat("\nGLOBAL (18m): CorticalGM ↔ TotalCortVol\n")
cor(individual_slopes_18m$slope_CorticalGM,
    individual_slopes_18m$slope_TotalCortVol,
    use = "pairwise.complete.obs") %>% round(3) %>% cat("r =", ., "\n")

################Better way of doing this feature extraction analyisis########
# All 33 normalised features, including Third Ventricle
norm_cols_all <- c(norm_cols, ThirdVentricle = "ST127SV_norm")

metrics_12m <- map_dfr(names(norm_cols_all), function(region_name) {
  
  col_name <- norm_cols_all[[region_name]]
  
  fit <- lmer(as.formula(paste0(col_name, " ~ time_years + (1 + time_years | RID)")),
              data = lme_12m, REML = TRUE,
              control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 20000)))
  
  vc        <- as.data.frame(VarCorr(fit))
  tau_slope <- vc$sdcor[vc$grp == "RID" & vc$var1 == "time_years" & is.na(vc$var2)]
  sigma_res <- sigma(fit)
  
  # Baseline mean volume (first visit per participant)
  bl_mean <- lme_12m %>%
    group_by(RID) %>% slice_min(time_years, n = 1, with_ties = FALSE) %>% ungroup() %>%
    pull(all_of(col_name)) %>% mean(na.rm = TRUE)
  
  # Per-participant slope reliability
  rel <- lme_12m %>%
    group_by(RID) %>%
    summarise(ss_t = sum((time_years - mean(time_years))^2), .groups = "drop") %>%
    mutate(rel_i = tau_slope^2 / (tau_slope^2 + sigma_res^2 / ss_t))
  
  tibble(
    region          = region_name,
    APC_pct_yr      = round(fixef(fit)["time_years"] / bl_mean * 100, 3),
    between_SD_pct  = round(tau_slope / bl_mean * 100, 3),
    reliability     = round(mean(rel$rel_i), 3),
    singular        = isSingular(fit)
  )
})

metrics_12m %>% arrange(desc(reliability)) %>% print(n = Inf)


#for 18 months window
metrics_18m <- map_dfr(names(norm_cols_all), function(region_name) {
  
  col_name <- norm_cols_all[[region_name]]
  
  fit <- lmer(as.formula(paste0(col_name, " ~ time_years + (1 + time_years | RID)")),
              data = lme_18m_data, REML = TRUE,
              control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 20000)))
  
  vc        <- as.data.frame(VarCorr(fit))
  tau_slope <- vc$sdcor[vc$grp == "RID" & vc$var1 == "time_years" & is.na(vc$var2)]
  sigma_res <- sigma(fit)
  
  # Baseline mean volume (first visit per participant)
  bl_mean <- lme_18m_data %>%
    group_by(RID) %>% slice_min(time_years, n = 1, with_ties = FALSE) %>% ungroup() %>%
    pull(all_of(col_name)) %>% mean(na.rm = TRUE)
  
  # Per-participant slope reliability
  rel <- lme_18m_data %>%
    group_by(RID) %>%
    summarise(ss_t = sum((time_years - mean(time_years))^2), .groups = "drop") %>%
    mutate(rel_i = tau_slope^2 / (tau_slope^2 + sigma_res^2 / ss_t))
  
  tibble(
    region          = region_name,
    APC_pct_yr      = round(fixef(fit)["time_years"] / bl_mean * 100, 3),
    between_SD_pct  = round(tau_slope / bl_mean * 100, 3),
    reliability     = round(mean(rel$rel_i), 3),
    singular        = isSingular(fit)
  )
})

metrics_18m %>% arrange(desc(reliability)) %>% print(n = Inf)

##########
# Add bilateral composite volumes (L + R) and a total ventricular volume, ICV-normalised
add_composites <- function(df) {
  for (i in seq_len(nrow(feature_map))) {
    df[[paste0(feature_map$region[i], "_bi")]] <-
      (df[[feature_map$left[i]]] + df[[feature_map$right[i]]]) / df$ST10CV * 100
  }
  df$Ventricles_all <- (df$ST37SV + df$ST96SV + df$ST30SV + df$ST89SV + df$ST127SV) / df$ST10CV * 100
  df
}

lme_12m_comp <- add_composites(lme_12m)
lme_18m_comp <- add_composites(lme_18m_data)

composite_cols <- c(paste0(feature_map$region, "_bi"), "Ventricles_all",
                    "ST127SV_norm", "ST149SV_norm", "TotalCortVol_norm")

# Same metrics as before: APC, scaled between-person SD, slope reliability
lme_metrics <- function(df, col) {
  fit <- lmer(as.formula(paste0(col, " ~ time_years + (1 + time_years | RID)")),
              data = df, REML = TRUE,
              control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 20000)))
  vc        <- as.data.frame(VarCorr(fit))
  tau_slope <- vc$sdcor[vc$grp == "RID" & vc$var1 == "time_years" & is.na(vc$var2)]
  bl_mean   <- df %>% group_by(RID) %>% slice_min(time_years, n = 1, with_ties = FALSE) %>%
    ungroup() %>% pull(all_of(col)) %>% mean(na.rm = TRUE)
  ss_t <- df %>% group_by(RID) %>%
    summarise(ss = sum((time_years - mean(time_years))^2), .groups = "drop") %>% pull(ss)
  tibble(region         = col,
         APC_pct_yr     = round(fixef(fit)["time_years"] / bl_mean * 100, 3),
         between_SD_pct = round(tau_slope / bl_mean * 100, 3),
         reliability    = round(mean(tau_slope^2 / (tau_slope^2 + sigma(fit)^2 / ss_t)), 3),
         singular       = isSingular(fit))
}

comp_12m <- map_dfr(composite_cols, ~ lme_metrics(lme_12m_comp, .x))
comp_18m <- map_dfr(composite_cols, ~ lme_metrics(lme_18m_comp, .x))

# Side-by-side comparison of composite reliability across windows
comp_compare <- comp_12m %>%
  rename_with(~ paste0(.x, "_12m"), -region) %>%
  left_join(comp_18m %>% rename_with(~ paste0(.x, "_18m"), -region), by = "region") %>%
  select(region, APC_pct_yr_12m, APC_pct_yr_18m,
         reliability_12m, reliability_18m, singular_12m, singular_18m) %>%
  arrange(desc(reliability_18m))

print(comp_compare, n = Inf, width = Inf)

########redoing analyis based onbasleine diagnosis grouping##########
# Refit composites with time × baseline diagnosis; also compute reliability within CN and MCI separately
lme_metrics_dx <- function(df, col) {
  df  <- df %>% mutate(BL_DIAGNOSIS = factor(BL_DIAGNOSIS, levels = c("CN", "MCI")))
  fit <- lmer(as.formula(paste0(col, " ~ time_years * BL_DIAGNOSIS + (1 + time_years | RID)")),
              data = df, REML = TRUE,
              control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 20000)))
  vc        <- as.data.frame(VarCorr(fit))
  tau_slope <- vc$sdcor[vc$grp == "RID" & vc$var1 == "time_years" & is.na(vc$var2)]
  sig       <- sigma(fit)
  bl        <- df %>% group_by(RID) %>% slice_min(time_years, n = 1, with_ties = FALSE) %>% ungroup()
  bl_mean   <- mean(bl[[col]], na.rm = TRUE)
  
  # Per-person reliability, then averaged overall and within each baseline group
  rel <- df %>% group_by(RID, BL_DIAGNOSIS) %>%
    summarise(ss_t = sum((time_years - mean(time_years))^2), .groups = "drop") %>%
    mutate(rel_i = tau_slope^2 / (tau_slope^2 + sig^2 / ss_t))
  
  fe <- fixef(fit)
  tibble(region          = col,
         APC_CN          = round(fe["time_years"] / bl_mean * 100, 3),
         APC_MCI         = round((fe["time_years"] + fe["time_years:BL_DIAGNOSISMCI"]) / bl_mean * 100, 3),
         reliability_all = round(mean(rel$rel_i), 3),
         reliability_CN  = round(mean(rel$rel_i[rel$BL_DIAGNOSIS == "CN"]), 3),
         reliability_MCI = round(mean(rel$rel_i[rel$BL_DIAGNOSIS == "MCI"]), 3),
         singular        = isSingular(fit))
}

comp_dx_12m <- map_dfr(composite_cols, ~ lme_metrics_dx(lme_12m_comp, .x))
comp_dx_18m <- map_dfr(composite_cols, ~ lme_metrics_dx(lme_18m_comp, .x))

# Readable labels for the three coded features
label_fix <- c(ST127SV_norm = "ThirdVentricle", ST149SV_norm = "CorticalGM",
               TotalCortVol_norm = "TotalCortVol")
tidy_labels <- function(x) x %>% mutate(region = recode(region, !!!label_fix))

comp_dx_12m %>% tidy_labels() %>% arrange(desc(reliability_all)) %>% print(n = Inf, width = Inf)
comp_dx_18m %>% tidy_labels() %>% arrange(desc(reliability_all)) %>% print(n = Inf, width = Inf)


##########
apply_rule <- function(metrics, window) {
  metrics %>%
    filter(region != "Ventricles_all") %>%                     # redundant sum of other ventricle measures
    mutate(is_ventricle = str_detect(region, "Ventric|ST127SV"),
           sign_ok      = if_else(is_ventricle, APC_CN > 0 & APC_MCI > 0, APC_CN < 0 & APC_MCI < 0),
           mci_faster   = abs(APC_MCI) >= abs(APC_CN),         # descriptive only, not a filter
           keep         = reliability_all >= 0.20 & !singular & sign_ok,
           window       = window)
}

rule_12m <- apply_rule(comp_dx_12m, "12m")
rule_18m <- apply_rule(comp_dx_18m, "18m")

# Side-by-side view of reliability, checks and keep decision in each window
rule_compare <- rule_12m %>%
  select(region, rel_12m = reliability_all, sign_ok_12m = sign_ok,
         mci_faster_12m = mci_faster, keep_12m = keep) %>%
  left_join(rule_18m %>% select(region, rel_18m = reliability_all, sign_ok_18m = sign_ok,
                                mci_faster_18m = mci_faster, keep_18m = keep), by = "region") %>%
  tidy_labels() %>%
  arrange(desc(keep_18m), desc(rel_18m))

print(rule_compare, n = Inf, width = Inf)

# Each window's own selected list
selected_12m <- rule_12m %>% filter(keep) %>% pull(region)
selected_18m <- rule_18m %>% filter(keep) %>% pull(region)

cat("\n12m selected (", length(selected_12m), "):", paste(recode(selected_12m, !!!label_fix), collapse = ", "), "\n")
cat("18m selected (", length(selected_18m), "):", paste(recode(selected_18m, !!!label_fix), collapse = ", "), "\n")
cat("In both:", paste(recode(intersect(selected_12m, selected_18m), !!!label_fix), collapse = ", "), "\n")

##Extract individual slopes (group slope + personal deviation) for each window's selected features, then check row counts, missing values and redundancy
extract_slopes_dx <- function(df, col) {
  df  <- df %>% mutate(BL_DIAGNOSIS = factor(BL_DIAGNOSIS, levels = c("CN", "MCI")))
  fit <- lmer(as.formula(paste0(col, " ~ time_years * BL_DIAGNOSIS + (1 + time_years | RID)")),
              data = df, REML = TRUE,
              control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 20000)))
  fe <- fixef(fit)
  
  # Personal deviation from the group slope (random slope / BLUP)
  re <- ranef(fit)$RID %>%
    rownames_to_column("RID") %>%
    transmute(RID = as.numeric(RID), re_slope = time_years)
  
  df %>% distinct(RID, BL_DIAGNOSIS) %>%
    left_join(re, by = "RID") %>%
    transmute(RID,
              !!paste0("slope_", str_remove(recode(col, !!!label_fix), "_bi$")) :=
                fe["time_years"] +                                                     # CN average slope
                if_else(BL_DIAGNOSIS == "MCI", fe["time_years:BL_DIAGNOSISMCI"], 0) +  # extra for MCI
                re_slope)                                                              # personal deviation
}

final_slopes_12m <- map(selected_12m, ~ extract_slopes_dx(lme_12m_comp, .x)) %>% reduce(left_join, by = "RID")
final_slopes_18m <- map(selected_18m, ~ extract_slopes_dx(lme_18m_comp, .x)) %>% reduce(left_join, by = "RID")

# Integrity: one row per participant, no missing slopes (expect 460 and 459)
cat("12m:", nrow(final_slopes_12m), "rows,", n_distinct(final_slopes_12m$RID), "RIDs, NA =", sum(is.na(final_slopes_12m)), "\n")
cat("18m:", nrow(final_slopes_18m), "rows,", n_distinct(final_slopes_18m$RID), "RIDs, NA =", sum(is.na(final_slopes_18m)), "\n")

# Redundancy among the final slopes in each window
cat("\n12m correlations:\n")
final_slopes_12m %>% select(-RID) %>% rename_with(~ str_remove(., "slope_")) %>% cor() %>% round(2) %>% print()
cat("\n18m correlations:\n")
final_slopes_18m %>% select(-RID) %>% rename_with(~ str_remove(., "slope_")) %>% cor() %>% round(2) %>% print()

#####Build the final metrics table from the shared-spread model: APC for each diagnosis, between-person SD, reliability and the keep decision for both windows

table_metrics <- function(df, col, window) {
  df  <- df %>% mutate(BL_DIAGNOSIS = factor(BL_DIAGNOSIS, levels = c("CN", "MCI")))
  fit <- lmer(as.formula(paste0(col, " ~ time_years * BL_DIAGNOSIS + (1 + time_years | RID)")),
              data = df, REML = TRUE,
              control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 20000)))
  vc  <- as.data.frame(VarCorr(fit))
  tau <- vc$sdcor[vc$grp == "RID" & vc$var1 == "time_years" & is.na(vc$var2)]
  fe  <- fixef(fit)
  
  # Baseline mean volume per diagnosis, and overall
  bl  <- df %>% group_by(RID) %>% slice_min(time_years, n = 1, with_ties = FALSE) %>% ungroup()
  bm  <- bl %>% group_by(BL_DIAGNOSIS) %>% summarise(m = mean(.data[[col]], na.rm = TRUE)) %>% deframe()
  ss  <- df %>% group_by(RID) %>% summarise(ss = sum((time_years - mean(time_years))^2)) %>% pull(ss)
  
  tibble(window     = window,
         region     = col,
         APC_CN     = round(fe["time_years"] / bm["CN"] * 100, 3),
         APC_MCI    = round((fe["time_years"] + fe["time_years:BL_DIAGNOSISMCI"]) / bm["MCI"] * 100, 3),
         between_SD = round(tau / mean(bl[[col]], na.rm = TRUE) * 100, 3),   # shared across CN and MCI
         reliability = round(mean(tau^2 / (tau^2 + sigma(fit)^2 / ss)), 3),
         singular   = isSingular(fit))
}

composites_eval <- setdiff(composite_cols, "Ventricles_all")

feature_table <- bind_rows(
  map_dfr(composites_eval, ~ table_metrics(lme_12m_comp, .x, "12m")),
  map_dfr(composites_eval, ~ table_metrics(lme_18m_comp, .x, "18m"))
) %>%
  tidy_labels() %>%
  left_join(bind_rows(rule_12m, rule_18m) %>% tidy_labels() %>% select(window, region, keep),
            by = c("window", "region")) %>%
  mutate(region = str_remove(region, "_bi$")) %>%
  pivot_wider(names_from = window,
              values_from = c(APC_CN, APC_MCI, between_SD, reliability, singular, keep),
              names_glue = "{.value}_{window}") %>%
  select(region, ends_with("_12m"), ends_with("_18m")) %>%
  arrange(desc(keep_18m), desc(reliability_18m))

print(feature_table, n = Inf, width = Inf)
write_csv(feature_table, "MRI_feature_metrics_final.csv")

######creating the table of reliable features in each window for the results section######
# Readable region names for the manuscript
region_labels <- c(InfLatVentricle = "Inferior lateral ventricle", LatVentricle = "Lateral ventricle",
                   ThirdVentricle = "Third ventricle", Caudate = "Caudate", Hippocampus = "Hippocampus",
                   TotalCortVol = "Total cortical volume", Fusiform = "Fusiform", Amygdala = "Amygdala",
                   MiddleTemporal = "Middle temporal", Thalamus = "Thalamus", Precuneus = "Precuneus",
                   LateralOccipital = "Lateral occipital")

# One table per window: keep only that window's reliable features
make_results_table <- function(win) {
  feature_table %>%
    filter(.data[[paste0("keep_", win)]]) %>%
    transmute(Region                     = recode(region, !!!region_labels),
              `APC CN (%/yr)`            = round(.data[[paste0("APC_CN_", win)]], 2),
              `APC MCI (%/yr)`           = round(.data[[paste0("APC_MCI_", win)]], 2),
              `Between-person SD (%/yr)` = round(.data[[paste0("between_SD_", win)]], 2),
              Reliability                = round(.data[[paste0("reliability_", win)]], 3)) %>%
    arrange(desc(Reliability))
}

results_12m <- make_results_table("12m")
results_18m <- make_results_table("18m")

print(results_12m); print(results_18m, n = Inf)

write_csv(results_12m, "Table_MRI_features_12m.csv")
write_csv(results_18m, "Table_MRI_features_18m.csv")

#########correlation between fetures#######
# Find feature pairs with |r| >= 0.80 and mark the less reliable one for dropping
flag_redundant <- function(slopes, rule_tbl, cutoff = 0.80) {
  rel <- rule_tbl %>% tidy_labels() %>%
    transmute(feature = str_remove(region, "_bi$"), reliability = reliability_all)
  
  slopes %>% select(-RID) %>% rename_with(~ str_remove(., "slope_")) %>%
    cor() %>% as.data.frame() %>% rownames_to_column("f1") %>%
    pivot_longer(-f1, names_to = "f2", values_to = "r") %>%
    filter(f1 < f2, abs(r) >= cutoff) %>%                        # each pair once
    left_join(rel, by = c("f1" = "feature")) %>% rename(rel_f1 = reliability) %>%
    left_join(rel, by = c("f2" = "feature")) %>% rename(rel_f2 = reliability) %>%
    mutate(drop = if_else(rel_f1 >= rel_f2, f2, f1), r = round(r, 4))
}

redund_12m <- flag_redundant(final_slopes_12m, rule_12m)
redund_18m <- flag_redundant(final_slopes_18m, rule_18m)
print(redund_12m); print(redund_18m)

# Two candidate feature lists per window: all reliable features vs redundancy-pruned
mri_all_12m    <- str_remove(names(final_slopes_12m)[-1], "slope_")
mri_all_18m    <- str_remove(names(final_slopes_18m)[-1], "slope_")
mri_pruned_12m <- setdiff(mri_all_12m, redund_12m$drop)
mri_pruned_18m <- setdiff(mri_all_18m, redund_18m$drop)

cat("\n12m all (", length(mri_all_12m), "):", paste(mri_all_12m, collapse = ", "),
    "\n12m pruned (", length(mri_pruned_12m), "):", paste(mri_pruned_12m, collapse = ", "),
    "\n18m all (", length(mri_all_18m), "):", paste(mri_all_18m, collapse = ", "),
    "\n18m pruned (", length(mri_pruned_18m), "):", paste(mri_pruned_18m, collapse = ", "), "\n")

# Save the slope tables for the W-score step
write_csv(final_slopes_12m, "MRI_slopes_12m.csv")
write_csv(final_slopes_18m, "MRI_slopes_18m.csv")


#########
# Lower-triangle correlation heatmap with readable names and r values
theme_pub <- theme_classic(base_size = 12) +
  theme(panel.grid       = element_blank(),
        panel.background = element_rect(fill = "white", colour = NA),
        plot.background  = element_rect(fill = "white", colour = NA),
        axis.line        = element_line(colour = "black"),
        axis.ticks       = element_line(colour = "black"),
        axis.text        = element_text(colour = "black"),
        plot.title       = element_text(face = "bold", hjust = 0.5),
        legend.title     = element_text(face = "bold"))

class(theme_pub)   # should show "theme" "gg"

plot_slope_cor <- function(slopes, title) {
  cm <- slopes %>% select(-RID) %>%
    rename_with(~ recode(str_remove(., "slope_"), !!!region_labels)) %>%
    cor()
  
  ord <- colnames(cm)
  cm[upper.tri(cm)] <- NA                                     # keep lower triangle only
  
  cm %>% as.data.frame() %>% rownames_to_column("f1") %>%
    pivot_longer(-f1, names_to = "f2", values_to = "r") %>%
    filter(!is.na(r)) %>%
    mutate(f1 = factor(f1, levels = rev(ord)), f2 = factor(f2, levels = ord)) %>%
    ggplot(aes(f2, f1, fill = r)) +
    geom_tile(colour = "white") +
    geom_text(aes(label = sprintf("%.2f", r)), size = 3.2) +
    scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                         midpoint = 0, limits = c(-1, 1), name = "Pearson r") +
    labs(title = title, x = NULL, y = NULL) +
    theme_pub +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          axis.line = element_blank(), axis.ticks = element_blank())
}

p_cor_12m <- plot_slope_cor(final_slopes_12m, "Correlation of individual MRI slopes: 12-month window")
p_cor_18m <- plot_slope_cor(final_slopes_18m, "Correlation of individual MRI slopes: 18-month window")

p_cor_12m
p_cor_18m

ggsave("Fig_MRI_slope_correlation_12m.png", p_cor_12m, width = 6,  height = 5, dpi = 300)
ggsave("Fig_MRI_slope_correlation_18m.png", p_cor_18m, width = 10, height = 8, dpi = 300)

#######checking whether to select the temporooccipital composite or only middle temporal#########
add_to <- function(df) {
  df %>% mutate(TemporoOccipital_bi = (ST40CV + ST99CV + ST26CV + ST85CV + ST35CV + ST94CV) / ST10CV * 100)
}
lme_12m_comp <- add_to(lme_12m_comp)
lme_18m_comp <- add_to(lme_18m_comp)

# Readable names: new primary composite, temporoparietal kept for sensitivity
region_labels <- c(region_labels[!names(region_labels) %in% c("Temporoparietal", "TemporoOccipital")],
                   TemporoOccipital = "Temporo-occipital composite",
                   Temporoparietal  = "Temporoparietal composite")

# Same metrics and same rule as every other feature, for both composites
composite_metrics <- function(col, name) {
  bind_rows(table_metrics(lme_12m_comp, col, "12m"),
            table_metrics(lme_18m_comp, col, "18m")) %>%
    mutate(region = name,
           keep   = reliability >= 0.20 & !singular & APC_CN < 0 & APC_MCI < 0) %>%
    pivot_wider(names_from = window,
                values_from = c(APC_CN, APC_MCI, between_SD, reliability, singular, keep),
                names_glue = "{.value}_{window}")
}

to_metrics <- composite_metrics("TemporoOccipital_bi", "TemporoOccipital")
tp_metrics <- composite_metrics("Temporoparietal_bi",  "Temporoparietal")
print(bind_rows(to_metrics, tp_metrics), width = Inf)

lme_12m_comp <- lme_12m_comp %>% select(-any_of("Temporoparietal_bi"))
lme_18m_comp <- lme_18m_comp %>% select(-any_of("Temporoparietal_bi"))
region_labels <- region_labels[!names(region_labels) %in% c("Temporoparietal", "TemporoOccipital")]
rm(list = intersect(ls(), c("tp_slopes_18m", "tp_metrics", "feature_table_tp", "final_slopes_18m_tp")))

# Check: both should print FALSE
"Temporoparietal_bi" %in% names(lme_18m_comp); "Temporoparietal" %in% names(region_labels)


add_to <- function(df) {
  df %>% mutate(TemporoOccipital_bi = (ST40CV + ST99CV + ST26CV + ST85CV + ST35CV + ST94CV) / ST10CV * 100)
}
lme_12m_comp <- add_to(lme_12m_comp)
lme_18m_comp <- add_to(lme_18m_comp)

region_labels <- c(region_labels, TemporoOccipital = "Temporo-occipital composite")

# Same rule: reliability >= 0.20, non-singular fit, tissue shrinking in both groups
to_metrics <- bind_rows(table_metrics(lme_12m_comp, "TemporoOccipital_bi", "12m"),
                        table_metrics(lme_18m_comp, "TemporoOccipital_bi", "18m")) %>%
  mutate(region = "TemporoOccipital",
         keep   = reliability >= 0.20 & !singular & APC_CN < 0 & APC_MCI < 0) %>%
  pivot_wider(names_from = window,
              values_from = c(APC_CN, APC_MCI, between_SD, reliability, singular, keep),
              names_glue = "{.value}_{window}")
print(to_metrics, width = Inf)


feature_table_comp <- bind_rows(feature_table, to_metrics)

#########
make_results_table <- function(tbl, win) {
  tbl %>%
    filter(.data[[paste0("keep_", win)]]) %>%
    transmute(Region                     = recode(region, !!!region_labels),
              `APC CN (%/yr)`            = round(.data[[paste0("APC_CN_", win)]], 2),
              `APC MCI (%/yr)`           = round(.data[[paste0("APC_MCI_", win)]], 2),
              `Between-person SD (%/yr)` = round(.data[[paste0("between_SD_", win)]], 2),
              Reliability                = round(.data[[paste0("reliability_", win)]], 3)) %>%
    arrange(desc(Reliability))
}

results_12m <- make_results_table(feature_table_comp, "12m")
results_18m <- make_results_table(feature_table_comp, "18m")
print(results_12m); print(results_18m, n = Inf)

write_csv(feature_table_comp, "MRI_feature_metrics_final.csv")
write_csv(results_12m, "Table_MRI_features_12m.csv")
write_csv(results_18m, "Table_MRI_features_18m.csv")

############
to_slopes_12m <- extract_slopes_dx(lme_12m_comp, "TemporoOccipital_bi")
to_slopes_18m <- extract_slopes_dx(lme_18m_comp, "TemporoOccipital_bi")

# The composite goes into a window's heatmap only if it passed the rule in that window
slopes_12m_plot <- if (to_metrics$keep_12m) left_join(final_slopes_12m, to_slopes_12m, by = "RID") else final_slopes_12m
slopes_18m_plot <- if (to_metrics$keep_18m) left_join(final_slopes_18m, to_slopes_18m, by = "RID") else final_slopes_18m

# Integrity check: expect 460 and 459 rows, 0 NA
cat("12m:", nrow(slopes_12m_plot), "rows, NA =", sum(is.na(slopes_12m_plot)),
    "| 18m:", nrow(slopes_18m_plot), "rows, NA =", sum(is.na(slopes_18m_plot)), "\n")

p_cor_12m <- plot_slope_cor(slopes_12m_plot, "Correlation of individual MRI slopes: 12-month window")
p_cor_18m <- plot_slope_cor(slopes_18m_plot, "Correlation of individual MRI slopes: 18-month window")
p_cor_12m; p_cor_18m
ggsave("Fig_MRI_slope_correlation_12m.png", p_cor_12m, width = 6,  height = 5, dpi = 300)
ggsave("Fig_MRI_slope_correlation_18m.png", p_cor_18m, width = 11, height = 9, dpi = 300)

wscore_mri_12m <- final_slopes_12m %>% select(RID, slope_InfLatVentricle, slope_Hippocampus)

wscore_mri_18m <- final_slopes_18m %>%
  select(RID, slope_InfLatVentricle, slope_Hippocampus, slope_TotalCortVol) %>%
  left_join(to_slopes_18m, by = "RID") %>%
  relocate(slope_TemporoOccipital, .after = slope_Hippocampus)

cat("12m:", nrow(wscore_mri_12m), "rows, NA =", sum(is.na(wscore_mri_12m)), "\n")
cat("18m:", nrow(wscore_mri_18m), "rows, NA =", sum(is.na(wscore_mri_18m)), "\n")

# Overlap among the four 18-month W-score features (flag any |r| >= 0.80)
wscore_mri_18m %>% select(-RID) %>%
  rename_with(~ recode(str_remove(., "slope_"), !!!region_labels)) %>%
  cor() %>% round(2) %>% print()

write_csv(wscore_mri_12m, "Wscore_MRI_slopes_12m.csv")
write_csv(wscore_mri_18m, "Wscore_MRI_slopes_18m.csv")




