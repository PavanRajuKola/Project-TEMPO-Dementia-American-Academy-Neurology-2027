source(here::here("config.R"))

# Check what diagnosis variables exist in master_visit
cat("Diagnosis variables available:\n")
cat("  DIAGNOSIS (observed):", "DIAGNOSIS" %in% names(master_visit), "\n")
cat("  DIAGNOSIS_FINAL (observed + imputed):", "DIAGNOSIS_FINAL" %in% names(master_visit), "\n")
cat("  DIAGNOSIS_SOURCE:", "DIAGNOSIS_SOURCE" %in% names(master_visit), "\n\n")

# Use baseline_dx (already created: earliest visit per participant)
# Check it still exists and has the expected structure
cat("baseline_dx exists:", exists("baseline_dx"), "\n")
cat("Rows:", nrow(baseline_dx), "| Unique RIDs:", n_distinct(baseline_dx$RID), "\n\n")

# Baseline diagnosis distribution
cat("Baseline DIAGNOSIS_FINAL distribution:\n")
baseline_dx %>%
  count(DIAGNOSIS_FINAL, sort = TRUE) %>%
  mutate(pct = round(n / sum(n) * 100, 1)) %>%
  print()

# Create a baseline diagnosis lookup (one per participant)
bl_dx_lookup <- baseline_dx %>%
  select(RID, BL_DIAGNOSIS = DIAGNOSIS_FINAL)

# Attach to master_visit
master_visit <- master_visit %>%
  left_join(bl_dx_lookup, by = "RID")

# Compute time from baseline (years) for each participant
master_visit <- master_visit %>%
  group_by(RID) %>%
  mutate(
    baseline_date = min(EXAMDATE, na.rm = TRUE),
    time_years    = as.numeric(EXAMDATE - baseline_date) / 365.25
  ) %>%
  ungroup()

# Integrity checks
cat("1. Baseline diagnosis coverage:\n")
cat("   Participants with BL_DIAGNOSIS:", sum(!is.na(master_visit$BL_DIAGNOSIS)), "/", nrow(master_visit), "\n\n")

cat("2. Duplicate RID × VISCODE2 check:\n")
dup_check <- master_visit %>% count(RID, VISCODE2) %>% filter(n > 1)
cat("   Duplicates:", nrow(dup_check), "\n\n")

cat("3. Negative follow-up times:\n")
neg_time <- sum(master_visit$time_years < 0, na.rm = TRUE)
cat("   Observations with time < 0:", neg_time, "\n\n")

cat("4. Observations per participant:\n")
obs_per_pt <- master_visit %>%
  count(RID, name = "n_obs")
print(summary(obs_per_pt$n_obs))

cat("\n5. Follow-up duration (years) per participant:\n")
fu_per_pt <- master_visit %>%
  group_by(RID) %>%
  summarise(fu_years = max(time_years), .groups = "drop")
print(summary(fu_per_pt$fu_years))

cat("\n6. Group sizes and follow-up by baseline diagnosis:\n")
master_visit %>%
  group_by(BL_DIAGNOSIS) %>%
  summarise(
    n_participants = n_distinct(RID),
    n_observations = n(),
    median_visits  = median(table(RID)),
    median_fu_yrs  = round(median(time_years[time_years == max(time_years)]), 1),
    .groups = "drop"
  ) %>%
  print()

###
# Compute bilateral hippocampal volume normalized by ICV
plot_data <- master_visit %>%
  mutate(
    Hippocampus_bilateral = (ST29SV + ST88SV),
    Hippocampus_norm = Hippocampus_bilateral / ST10CV * 100,
    BL_DIAGNOSIS = factor(BL_DIAGNOSIS, levels = c("CN", "MCI", "Dementia"))
  ) %>%
  filter(!is.na(BL_DIAGNOSIS))

# Spaghetti plot with group-level smooth
p_hippo <- ggplot(plot_data, aes(x = time_years, y = Hippocampus_norm,
                                 group = RID, color = BL_DIAGNOSIS)) +
  geom_line(alpha = 0.1, linewidth = 0.3) +
  geom_smooth(aes(group = BL_DIAGNOSIS), method = "loess", se = TRUE,
              linewidth = 1.5, alpha = 0.2) +
  scale_color_manual(values = c("CN" = "#4CAF50", "MCI" = "#FF9800", "Dementia" = "#F44336")) +
  labs(
    title = "Longitudinal Hippocampal Volume Trajectories by Baseline Diagnosis",
    subtitle = paste0("N = ", n_distinct(plot_data$RID),
                      " participants | Individual trajectories with group-level LOESS smooth"),
    x = "Time from Baseline (Years)",
    y = "Bilateral Hippocampal Volume (% of ICV)",
    color = "Baseline\nDiagnosis"
  ) +
  theme_pub +
  theme(legend.position = "top")

ggsave("trajectory_hippocampus.png", p_hippo,
       width = 8, height = 6, dpi = 300, bg = "white")
p_hippo

########
# Prepare data with all four key regions normalized by ICV
traj_data <- master_visit %>%
  mutate(
    BL_DIAGNOSIS = factor(BL_DIAGNOSIS, levels = c("CN", "MCI", "Dementia")),
    Hippocampus  = (ST29SV + ST88SV) / ST10CV * 100,
    Amygdala     = (ST12SV + ST71SV) / ST10CV * 100,
    Entorhinal   = (ST24CV + ST83CV) / ST10CV * 100,
    LatVentricle = (ST37SV + ST96SV) / ST10CV * 100
  ) %>%
  filter(!is.na(BL_DIAGNOSIS)) %>%
  select(RID, time_years, BL_DIAGNOSIS,
         Hippocampus, Amygdala, Entorhinal, LatVentricle) %>%
  pivot_longer(cols = c(Hippocampus, Amygdala, Entorhinal, LatVentricle),
               names_to = "Region", values_to = "Volume_pct_ICV") %>%
  mutate(Region = factor(Region, levels = c("Hippocampus", "Amygdala",
                                            "Entorhinal", "LatVentricle"),
                         labels = c("Hippocampus", "Amygdala",
                                    "Entorhinal Cortex", "Lateral Ventricles")))

p_multi <- ggplot(traj_data, aes(x = time_years, y = Volume_pct_ICV,
                                 group = RID, color = BL_DIAGNOSIS)) +
  geom_line(alpha = 0.08, linewidth = 0.3) +
  geom_smooth(aes(group = BL_DIAGNOSIS), method = "loess", se = TRUE,
              linewidth = 1.3, alpha = 0.2) +
  facet_wrap(~ Region, scales = "free_y", ncol = 2) +
  scale_color_manual(values = c("CN" = "#4CAF50", "MCI" = "#FF9800",
                                "Dementia" = "#F44336")) +
  labs(
    title = "Longitudinal MRI Volume Trajectories by Baseline Diagnosis",
    subtitle = paste0("N = ", n_distinct(master_visit$RID),
                      " participants | ICV-normalized | Individual + LOESS smooth"),
    x = "Time from Baseline (Years)",
    y = "Volume (% of ICV)",
    color = "Baseline Diagnosis"
  ) +
  theme_pub +
  theme(
    legend.position = "top",
    strip.text = element_text(face = "bold", size = 12)
  )

ggsave("trajectory_4regions.png", p_multi,
       width = 11, height = 9, dpi = 300, bg = "white")
p_multi

#########
# For each participant, get their baseline and last known diagnosis
conversion_status <- master_visit %>%
  filter(!is.na(DIAGNOSIS_FINAL)) %>%
  arrange(RID, EXAMDATE) %>%
  group_by(RID) %>%
  summarise(
    BL_DX   = first(DIAGNOSIS_FINAL),
    LAST_DX = last(DIAGNOSIS_FINAL),
    .groups = "drop"
  ) %>%
  left_join(bl_dx_lookup, by = "RID") %>%
  mutate(
    trajectory_group = case_when(
      BL_DIAGNOSIS == "CN" & LAST_DX == "CN"       ~ "CN-Stable",
      BL_DIAGNOSIS == "CN" & LAST_DX == "MCI"      ~ "CN→MCI",
      BL_DIAGNOSIS == "CN" & LAST_DX == "Dementia"  ~ "CN→Dementia",
      BL_DIAGNOSIS == "MCI" & LAST_DX == "MCI"      ~ "MCI-Stable",
      BL_DIAGNOSIS == "MCI" & LAST_DX == "Dementia"  ~ "MCI→Dementia",
      BL_DIAGNOSIS == "MCI" & LAST_DX == "CN"       ~ "MCI→CN",
      BL_DIAGNOSIS == "Dementia"                     ~ "Dementia",
      TRUE                                           ~ "Other"
    )
  )

cat("Trajectory group distribution:\n\n")
conversion_status %>%
  count(trajectory_group, sort = TRUE) %>%
  mutate(pct = round(n / sum(n) * 100, 1)) %>%
  print(n = Inf)

######
# First join trajectory_group from conversion_status
master_visit <- master_visit %>%
  left_join(conversion_status %>% select(RID, trajectory_group), by = "RID")

# Verify it attached
cat("trajectory_group in master_visit:", "trajectory_group" %in% names(master_visit), "\n")
master_visit %>% count(trajectory_group, sort = TRUE) %>% print()


# Define the 4 trajectory groups for the focused study
trajectory_cohort <- master_visit %>%
  filter(trajectory_group %in% c("CN-Stable", "CN→MCI", "CN→Dementia",
                                 "MCI-Stable", "MCI→Dementia")) %>%
  mutate(
    traj_group = case_when(
      trajectory_group == "CN-Stable"    ~ "CN-Stable",
      trajectory_group %in% c("CN→MCI", "CN→Dementia") ~ "CN-Progressor",
      trajectory_group == "MCI-Stable"   ~ "MCI-Stable",
      trajectory_group == "MCI→Dementia" ~ "MCI-Progressor"
    ),
    traj_group = factor(traj_group, levels = c("CN-Stable", "CN-Progressor",
                                               "MCI-Stable", "MCI-Progressor"))
  )

# Verify group sizes
cat("Trajectory cohort:\n")
cat("  Total participants:", n_distinct(trajectory_cohort$RID), "\n")
cat("  Total observations:", nrow(trajectory_cohort), "\n\n")

trajectory_cohort %>%
  group_by(traj_group) %>%
  summarise(
    n_participants = n_distinct(RID),
    n_observations = n(),
    median_visits  = median(table(RID)),
    median_fu_yrs  = round(median(time_years[time_years == max(time_years)]), 1),
    .groups = "drop"
  ) %>%
  print()

############
plot_data_4g <- trajectory_cohort %>%
  mutate(Hippocampus_norm = (ST29SV + ST88SV) / ST10CV * 100)

p_4g_hippo <- ggplot(plot_data_4g, aes(x = time_years, y = Hippocampus_norm,
                                       group = RID, color = traj_group)) +
  geom_line(alpha = 0.1, linewidth = 0.3) +
  geom_smooth(aes(group = traj_group), method = "loess", se = TRUE,
              linewidth = 1.3, alpha = 0.15) +
  scale_color_manual(values = c("CN-Stable"      = "#4CAF50",
                                "CN-Progressor"  = "#1B5E20",
                                "MCI-Stable"     = "#FF9800",
                                "MCI-Progressor" = "#E65100")) +
  labs(
    title = "Hippocampal Volume Trajectories: Stable vs Progressor",
    subtitle = paste0("N = ", n_distinct(plot_data_4g$RID),
                      " participants | ICV-normalized"),
    x = "Time from Baseline (Years)",
    y = "Bilateral Hippocampal Volume (% of ICV)",
    color = "Trajectory Group"
  ) +
  theme_pub +
  theme(legend.position = "top")

ggsave("trajectory_hippo_4groups.png", p_4g_hippo,
       width = 8, height = 6, dpi = 300, bg = "white")
p_4g_hippo


########Multi panel plot
plot_multi_4g <- trajectory_cohort %>%
  mutate(
    Hippocampus  = (ST29SV + ST88SV) / ST10CV * 100,
    Amygdala     = (ST12SV + ST71SV) / ST10CV * 100,
    Entorhinal   = (ST24CV + ST83CV) / ST10CV * 100,
    LatVentricle = (ST37SV + ST96SV) / ST10CV * 100
  ) %>%
  pivot_longer(cols = c(Hippocampus, Amygdala, Entorhinal, LatVentricle),
               names_to = "Region", values_to = "Volume_pct") %>%
  mutate(Region = factor(Region, levels = c("Hippocampus", "Amygdala",
                                            "Entorhinal", "LatVentricle"),
                         labels = c("Hippocampus", "Amygdala",
                                    "Entorhinal Cortex", "Lateral Ventricles")))

p_multi_4g <- ggplot(plot_multi_4g, aes(x = time_years, y = Volume_pct,
                                        group = RID, color = traj_group)) +
  geom_line(alpha = 0.07, linewidth = 0.3) +
  geom_smooth(aes(group = traj_group), method = "loess", se = TRUE,
              linewidth = 1.2, alpha = 0.15) +
  facet_wrap(~ Region, scales = "free_y", ncol = 2) +
  scale_color_manual(values = c("CN-Stable"      = "#4CAF50",
                                "CN-Progressor"  = "#1B5E20",
                                "MCI-Stable"     = "#FF9800",
                                "MCI-Progressor" = "#E65100")) +
  labs(
    title = "Longitudinal MRI Trajectories: Stable vs Progressor Groups",
    subtitle = paste0("N = ", n_distinct(trajectory_cohort$RID),
                      " participants | ICV-normalized"),
    x = "Time from Baseline (Years)",
    y = "Volume (% of ICV)",
    color = "Trajectory Group"
  ) +
  theme_pub +
  theme(
    legend.position = "top",
    strip.text = element_text(face = "bold", size = 12)
  )

ggsave("trajectory_4regions_4groups.png", p_multi_4g,
       width = 11, height = 9, dpi = 300, bg = "white")
p_multi_4g

############
# For each progressor, identify the visit where diagnosis first changed
conversion_timeline <- master_visit %>%
  filter(!is.na(DIAGNOSIS_FINAL), !is.na(EXAMDATE)) %>%
  arrange(RID, EXAMDATE) %>%
  group_by(RID) %>%
  mutate(
    first_dx = first(DIAGNOSIS_FINAL),
    dx_changed = DIAGNOSIS_FINAL != first_dx
  ) %>%
  filter(dx_changed) %>%
  slice(1) %>%
  ungroup() %>%
  select(RID, conversion_date = EXAMDATE, conversion_dx = DIAGNOSIS_FINAL,
         conversion_time_years = time_years)

# Merge with trajectory group info
conversion_timeline <- conversion_timeline %>%
  left_join(conversion_status %>% select(RID, trajectory_group), by = "RID")

# Only CN and MCI progressors
progressors <- conversion_timeline %>%
  filter(trajectory_group %in% c("CN→MCI", "CN→Dementia", "MCI→Dementia"))

cat("Conversion timeline for progressors:\n\n")
cat("N progressors:", nrow(progressors), "\n\n")

progressors %>%
  group_by(trajectory_group) %>%
  summarise(
    n = n(),
    median_conversion_yrs = round(median(conversion_time_years), 1),
    min_conversion_yrs    = round(min(conversion_time_years), 1),
    max_conversion_yrs    = round(max(conversion_time_years), 1),
    q25_conversion_yrs    = round(quantile(conversion_time_years, 0.25), 1),
    q75_conversion_yrs    = round(quantile(conversion_time_years, 0.75), 1),
    .groups = "drop"
  ) %>%
  print()

######
# For all 502 participants in trajectory_cohort: follow-up and visit distribution
fu_summary <- trajectory_cohort %>%
  group_by(RID) %>%
  summarise(
    traj_group = first(traj_group),
    n_visits = n(),
    fu_years = max(time_years),
    .groups = "drop"
  )

# How many have < 1 year follow-up?
cat("Participants with < 1 year total follow-up:\n")
fu_summary %>%
  filter(fu_years < 1) %>%
  count(traj_group) %>%
  print()

cat("\nParticipants with < 1.5 years total follow-up:\n")
fu_summary %>%
  filter(fu_years < 1.5) %>%
  count(traj_group) %>%
  print()

cat("\nFull distribution of follow-up by group:\n")
fu_summary %>%
  group_by(traj_group) %>%
  summarise(
    n = n(),
    min_fu  = round(min(fu_years), 1),
    q25_fu  = round(quantile(fu_years, 0.25), 1),
    median_fu = round(median(fu_years), 1),
    q75_fu  = round(quantile(fu_years, 0.75), 1),
    max_fu  = round(max(fu_years), 1),
    under_1yr = sum(fu_years < 1),
    .groups = "drop"
  ) %>%
  print()

########
early_converters <- conversion_timeline %>%
  filter(trajectory_group %in% c("CN→MCI", "CN→Dementia", "MCI→Dementia")) %>%
  mutate(
    baseline_group = if_else(str_starts(trajectory_group, "CN"), "CN", "MCI"),
    converted_before_1yr = conversion_time_years < 1
  )

cat("Conversions before 1 year:\n\n")
early_converters %>%
  group_by(baseline_group) %>%
  summarise(
    total_converters = n(),
    before_1yr = sum(converted_before_1yr),
    after_1yr  = sum(!converted_before_1yr),
    .groups = "drop"
  ) %>%
  print()

cat("\nDetailed: participants who converted before 1 year:\n")
early_converters %>%
  filter(converted_before_1yr) %>%
  select(RID, baseline_group, trajectory_group, 
         conversion_time_years, conversion_dx) %>%
  arrange(conversion_time_years) %>%
  print(n = Inf)


#############
# Flag early converters in trajectory_cohort
early_converter_rids <- early_converters %>%
  filter(converted_before_1yr) %>%
  pull(RID)

# Get baseline data for all 502 participants
char_data <- trajectory_cohort %>%
  arrange(RID, EXAMDATE) %>%
  distinct(RID, .keep_all = TRUE) %>%
  mutate(excluded = RID %in% early_converter_rids)

cat("Excluded (early converters):", sum(char_data$excluded), "\n")
cat("Retained:", sum(!char_data$excluded), "\n\n")

# Characteristics comparison
char_table <- char_data %>%
  group_by(excluded) %>%
  summarise(
    N = n(),
    Age_mean = round(mean(as.numeric(format(EXAMDATE, "%Y")) - PTDOBYY, na.rm = TRUE), 1),
    Age_sd   = round(sd(as.numeric(format(EXAMDATE, "%Y")) - PTDOBYY, na.rm = TRUE), 1),
    Female_n = sum(PTGENDER == "Female", na.rm = TRUE),
    Female_pct = round(sum(PTGENDER == "Female", na.rm = TRUE) / n() * 100, 1),
    Educ_mean = round(mean(PTEDUCAT, na.rm = TRUE), 1),
    Educ_sd   = round(sd(PTEDUCAT, na.rm = TRUE), 1),
    MMSE_mean = round(mean(MMSE_MMSCORE, na.rm = TRUE), 1),
    MMSE_sd   = round(sd(MMSE_MMSCORE, na.rm = TRUE), 1),
    APOE4_n   = sum(APOE4_carrier == TRUE, na.rm = TRUE),
    APOE4_pct = round(sum(APOE4_carrier == TRUE, na.rm = TRUE) / n() * 100, 1),
    CN_n      = sum(BL_DIAGNOSIS == "CN"),
    MCI_n     = sum(BL_DIAGNOSIS == "MCI"),
    .groups = "drop"
  ) %>%
  mutate(Group = if_else(excluded, "Excluded (early converters)", "Retained"))

cat("Characteristics: Excluded vs Retained\n\n")
char_table %>%
  select(Group, N, Age_mean, Age_sd, Female_n, Female_pct,
         Educ_mean, Educ_sd, MMSE_mean, MMSE_sd,
         APOE4_n, APOE4_pct, CN_n, MCI_n) %>%
  print(width = Inf)

View(char_table)
