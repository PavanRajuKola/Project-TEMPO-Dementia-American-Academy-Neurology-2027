source(here::here("config.R"))

####Working on calculating cognitive score decline#######
tab_dir <- TABLE_DIR
fig_dir <- FIGURE_DIR

cog_candidates <- names(master_visit) %>%
  str_subset(regex("^(ADAS|MMSE|CDR|FAQ|NB)_", ignore_case = TRUE)) %>%
  str_subset(regex("TOT|SCORE|SB|AVTOT|AVDEL|LDEL|TRAB|TRAA|CATANIM|DIGIT|FORGET",
                   ignore_case = TRUE))

master_visit %>%
  summarise(across(all_of(cog_candidates), ~ sum(!is.na(.)))) %>%
  pivot_longer(everything(), names_to = "column", values_to = "non_missing") %>%
  arrange(desc(non_missing)) %>%
  print(n = Inf)


##########
# Cohort RIDs and baseline diagnosis from the MRI LME data
cohort_12m <- lme_12m      %>% distinct(RID, BL_DIAGNOSIS)
cohort_18m <- lme_18m_data %>% distinct(RID, BL_DIAGNOSIS)

# ADAS-13 for these RIDs, timed from each person's first MRI date (same time origin as the MRI slopes)
get_adas <- function(cohort, max_years) {
  master_visit %>%
    filter(RID %in% cohort$RID) %>%
    group_by(RID) %>%
    mutate(time_years = as.numeric(as.Date(EXAMDATE) - min(as.Date(EXAMDATE))) / 365.25) %>%
    ungroup() %>%
    filter(time_years <= max_years, !is.na(ADAS_TOTAL13)) %>%
    select(RID, VISCODE2, time_years, ADAS13 = ADAS_TOTAL13) %>%
    left_join(cohort, by = "RID")
}

adas_12m <- get_adas(cohort_12m, 1.0)
adas_18m <- get_adas(cohort_18m, 1.5)

# Integrity: one row per RID per visit
cat("Duplicate RID-visit rows: 12m =", sum(duplicated(adas_12m[c("RID", "VISCODE2")])),
    "| 18m =", sum(duplicated(adas_18m[c("RID", "VISCODE2")])), "\n")

# Coverage: number of ADAS-13 scores per person within the window
adas_coverage <- function(adas, cohort, label) {
  per <- adas %>% count(RID, name = "n_scores")
  cat("\n", label, ": cohort =", nrow(cohort),
      "| any ADAS =", nrow(per),
      "| >= 2 scores =", sum(per$n_scores >= 2),
      "| no ADAS =", nrow(cohort) - nrow(per), "\n")
  per %>% count(n_scores, name = "participants") %>% print()
  adas %>% count(VISCODE2) %>% print()                         # should be bl / m06 / m12 (/ m18)
  adas %>% group_by(RID) %>% filter(n() >= 2) %>% ungroup() %>%
    distinct(RID, BL_DIAGNOSIS) %>% count(BL_DIAGNOSIS) %>% print()
}
adas_coverage(adas_12m, cohort_12m, "12-month")
adas_coverage(adas_18m, cohort_18m, "18-month")

##########
# Time of each MRI visit code from the first MRI (does m12 fall just after 1.0 year?)
master_visit %>%
  filter(RID %in% cohort_18m$RID) %>%
  group_by(RID) %>%
  mutate(time_years = as.numeric(as.Date(EXAMDATE) - min(as.Date(EXAMDATE))) / 365.25) %>%
  ungroup() %>%
  filter(VISCODE2 %in% c("scmri", "bl", "m03", "m06", "m12", "m18", "m24")) %>%
  group_by(VISCODE2) %>%
  summarise(n = n(), median_yr = round(median(time_years), 3),
            pct_over_1yr = round(100 * mean(time_years > 1.0), 1), .groups = "drop") %>%
  print()

# Visit composition of the 12-month MRI window actually used for the slopes
lme_12m %>% count(VISCODE2) %>% print()

# Raw ADAS for the cohort: which visit codes exist vs which reached master_visit
adas_raw <- read_csv(file.path(data_dir, "ADAS.csv"), show_col_types = FALSE)
adas_raw %>% filter(RID %in% cohort_18m$RID, !is.na(TOTAL13)) %>%
  count(VISCODE2, name = "in_ADAS_csv") %>%
  left_join(master_visit %>% filter(RID %in% cohort_18m$RID, !is.na(ADAS_TOTAL13)) %>%
              count(VISCODE2, name = "in_master_visit"), by = "VISCODE2") %>%
  filter(VISCODE2 %in% c("sc", "bl", "m06", "m12", "m18", "m24")) %>%
  print()

######
# Gap between first ADAS test and first MRI per person (expect mostly 0 to 60 days)
adas_12m %>% group_by(RID) %>% slice_min(time_years, n = 1, with_ties = FALSE) %>% ungroup() %>%
  summarise(median_days = median(time_years * 365.25),
            q05 = quantile(time_years * 365.25, 0.05), q95 = quantile(time_years * 365.25, 0.95),
            before_MRI = sum(time_years < 0), more_than_90d = sum(abs(time_years * 365.25) > 90)) %>%
  print()

###############date wise calculation of cognitive traj rates"###########
# ADAS-13 with its own test date; same dedup rule as before
adas_dated <- adas_raw %>%
  filter(!is.na(TOTAL13), !is.na(VISDATE)) %>%
  arrange(RID, VISCODE2, VISCODE) %>%
  distinct(RID, VISCODE2, .keep_all = TRUE) %>%
  transmute(RID, ADAS_VISCODE2 = VISCODE2, ADAS_DATE = as.Date(VISDATE), ADAS13 = TOTAL13)

# First MRI date and baseline diagnosis per participant
first_mri <- function(lme_df) {
  lme_df %>% group_by(RID) %>%
    summarise(MRI_T0 = min(as.Date(EXAMDATE)), BL_DIAGNOSIS = first(BL_DIAGNOSIS), .groups = "drop")
}

# ADAS timed from first MRI, cut at the window's MRI cut-off (one-to-many: one person, several tests)
# Same model and reliability formula as for MRI; individual slope = group slope + personal deviation
adas_lme <- function(adas, window) {
  d <- adas %>% group_by(RID) %>% filter(n() >= 2) %>% ungroup() %>%
    mutate(BL_DIAGNOSIS = factor(BL_DIAGNOSIS, levels = c("CN", "MCI")))
  
  fit <- lmer(ADAS13 ~ time_years * BL_DIAGNOSIS + (1 + time_years | RID), data = d, REML = TRUE,
              control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 20000)))
  vc  <- as.data.frame(VarCorr(fit))
  tau <- vc$sdcor[vc$grp == "RID" & vc$var1 == "time_years" & is.na(vc$var2)]
  fe  <- fixef(fit)
  ss  <- d %>% group_by(RID) %>% summarise(ss = sum((time_years - mean(time_years))^2)) %>% pull(ss)
  bl_sd <- d %>% group_by(RID) %>% slice_min(time_years, n = 1, with_ties = FALSE) %>%
    ungroup() %>% pull(ADAS13) %>% sd()
  
  metrics <- tibble(window = window, N = n_distinct(d$RID), obs = nrow(d),
                    slope_CN_pts_yr   = round(fe["time_years"], 3),
                    slope_MCI_pts_yr  = round(fe["time_years"] + fe["time_years:BL_DIAGNOSISMCI"], 3),
                    between_SD_pts_yr = round(tau, 3),
                    std_between_SD    = round(tau / bl_sd, 3),     # SD units per year
                    reliability       = round(mean(tau^2 / (tau^2 + sigma(fit)^2 / ss)), 3),
                    singular          = isSingular(fit))
  
  re <- ranef(fit)$RID %>% rownames_to_column("RID") %>%
    transmute(RID = as.numeric(RID), re_slope = time_years)
  slopes <- d %>% distinct(RID, BL_DIAGNOSIS) %>% left_join(re, by = "RID") %>%
    transmute(RID, slope_ADAS13 = fe["time_years"] +
                if_else(BL_DIAGNOSIS == "MCI", fe["time_years:BL_DIAGNOSISMCI"], 0) + re_slope)
  
  list(metrics = metrics, slopes = slopes)
}

adas_res_6m  <- adas_lme(adas_6m,  "6m")
adas_res_12m <- adas_lme(adas_12m, "12m")

# Do the early-m12 ADAS people match the early-m12 MRI people?
early_adas <- adas_6m %>% filter(ADAS_VISCODE2 == "m12") %>% distinct(RID) %>% pull()
early_mri  <- lme_12m %>% filter(VISCODE2 == "m12") %>% distinct(RID) %>% pull()
cat("Early m12 - ADAS:", length(early_adas), "| MRI:", length(early_mri),
    "| in both:", length(intersect(early_adas, early_mri)), "\n")

# Sensitivity: 6-month ADAS model without the early-m12 participants
adas_res_6m_sens <- adas_lme(adas_6m %>% filter(!RID %in% early_adas), "6m_excl_early_m12")

bind_rows(adas_res_6m$metrics, adas_res_6m_sens$metrics, adas_res_12m$metrics) %>% print(width = Inf)

cat("6m slopes:",  nrow(adas_res_6m$slopes),  "| NA =", sum(is.na(adas_res_6m$slopes)), "\n")
cat("12m slopes:", nrow(adas_res_12m$slopes), "| NA =", sum(is.na(adas_res_12m$slopes)), "\n")

# Reliability and singular columns that were cut off
bind_rows(adas_res_6m$metrics, adas_res_12m$metrics) %>%
  select(window, N, slope_CN_pts_yr, slope_MCI_pts_yr, reliability, singular) %>% print()

# Raw ADAS-13 mean by visit and diagnosis (12-month window): does MCI really improve?
adas_12m %>%
  filter(ADAS_VISCODE2 %in% c("bl", "m06", "m12")) %>%
  group_by(BL_DIAGNOSIS, ADAS_VISCODE2) %>%
  summarise(n = n(), mean_ADAS13 = round(mean(ADAS13), 2), sd = round(sd(ADAS13), 2), .groups = "drop") %>%
  print()

##########summary of MRI and ADAS timing by starting-visit group, with readable labels######
# Readable names for the starting-visit groups and for visit codes
start_group <- function(code) case_when(
  code == "scmri" ~ "Screening MRI",
  code == "sc"    ~ "Screening visit",
  code == "bl"    ~ "Baseline visit",
  code == "m03"   ~ "Month 3",
  TRUE            ~ "Rollover from earlier ADNI phase")

visit_label <- function(x) case_when(
  x == "scmri" ~ "Screening MRI", x == "sc" ~ "Screening", x == "bl" ~ "Baseline",
  str_detect(x, "^m\\d+") ~ paste("Month", as.integer(str_remove(x, "m"))), TRUE ~ x)

ordinal <- function(i) c("1st", "2nd", "3rd", "4th", "5th")[i]

# One row per test, timed from each person's own first MRI (12-month window cohort)
timing_long <- bind_rows(
  lme_18m_data %>% group_by(RID) %>%
    mutate(days = as.numeric(as.Date(EXAMDATE) - min(as.Date(EXAMDATE)))) %>% ungroup() %>%
    transmute(RID, test = "MRI scan", code = VISCODE2, days),
  adas_12m %>% transmute(RID, test = "ADAS-Cog 13", code = ADAS_VISCODE2, days = time_years * 365.25)
) %>%
  left_join(first_code, by = "RID") %>%
  mutate(group = start_group(first_MRI_code)) %>%
  group_by(RID, test) %>% mutate(order = rank(days, ties.method = "first")) %>% ungroup() %>%
  mutate(visit = if_else(group == "Rollover from earlier ADNI phase",
                         paste(ordinal(order), "visit"), visit_label(code)))

# Group sizes for the axis labels (e.g. "Screening MRI (n = 415)")
group_n <- first_code %>% mutate(group = start_group(first_MRI_code)) %>% count(group)

timing_summary <- timing_long %>%
  group_by(group, test, visit) %>%
  summarise(n = n(), median_days = round(median(days)),
            q25 = round(quantile(days, 0.25)), q75 = round(quantile(days, 0.75)), .groups = "drop") %>%
  filter(n >= 2 | group != "Screening MRI") %>%             # drop stray single tests in the main group
  left_join(group_n %>% rename(group_n = n), by = "group") %>%
  mutate(group_lab = paste0(group, " (n = ", group_n, ")"))

print(timing_summary, n = Inf)

########
group_levels <- c("Rollover from earlier ADNI phase", "Month 3", "Baseline visit",
                  "Screening visit", "Screening MRI")
lab_levels <- timing_summary %>% distinct(group, group_lab) %>%
  arrange(match(group, group_levels)) %>% pull(group_lab)

plot_df <- timing_summary %>%
  mutate(y = match(group_lab, lab_levels) + if_else(test == "MRI scan", 0.15, -0.15))

p_timing <- ggplot(plot_df, aes(median_days, y, colour = test, shape = test)) +
  annotate("rect", xmin = -91, xmax = 548, ymin = -Inf, ymax = Inf, fill = "grey92") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_errorbarh(aes(xmin = q25, xmax = q75), height = 0.1) +
  geom_point(size = 3) +
  geom_text(aes(label = visit), vjust = if_else(plot_df$test == "MRI scan", -0.9, 1.8),
            size = 2.8, show.legend = FALSE) +
  scale_y_continuous(breaks = seq_along(lab_levels), labels = lab_levels,
                     limits = c(0.5, length(lab_levels) + 0.5)) +
  scale_x_continuous(breaks = seq(-90, 540, 90)) +
  scale_colour_manual(values = c("MRI scan" = "#2166AC", "ADAS-Cog 13" = "#B2182B"), name = NULL) +
  scale_shape_manual(values = c("MRI scan" = 16, "ADAS-Cog 13" = 17), name = NULL) +
  labs(title = "Timing of MRI scans and ADAS-Cog 13 tests from each participant's first MRI",
       subtitle = "12-month window (N = 459). Shaded: assessments included (−3 months to 18 months). Dots: median; bars: 25th–75th percentile.",
       x = "Days from participant's first MRI", y = "First MRI visit") +
  theme_pub +
  theme(legend.position = "top")

p_timing
ggsave("Fig_visit_timing_by_start_visit.png", p_timing, width = 10, height = 6, dpi = 300)

#######
timing_table <- timing_summary %>%
  mutate(timing = paste0(median_days, " (", q25, " to ", q75, ")")) %>%
  select(`First MRI visit` = group_lab, Visit = visit, test, n, timing) %>%
  pivot_wider(names_from = test, values_from = c(n, timing), names_glue = "{test}: {.value}") %>%
  arrange(match(str_remove(`First MRI visit`, " \\(n.*"), rev(group_levels)), `Visit`)

print(timing_table, n = Inf, width = Inf)
write_csv(timing_table, "Table_visit_timing_by_start_visit.csv")


##########
timing_long <- timing_long %>%
  mutate(window_member = case_when(days <= 365 ~ "Both windows",
                                   days <= 548 ~ "18-month window only",
                                   TRUE        ~ "Outside"))

window_counts <- timing_long %>%
  count(group, test, visit, window_member) %>%
  pivot_wider(names_from = window_member, values_from = n, values_fill = 0)
print(window_counts, n = Inf)
write_csv(window_counts, file.path(tab_dir, "Table_assessments_by_window_membership.csv"))

# Points and labels: one label per visit, placed above the MRI point
plot_df <- timing_summary %>%
  filter(!(group == "Screening MRI" & n < 10)) %>%
  mutate(y = match(group_lab, lab_levels) + if_else(test == "MRI scan", 0.12, -0.12),
         months = median_days / 30.44, q25_m = q25 / 30.44, q75_m = q75 / 30.44)

visit_labels <- plot_df %>%
  group_by(group_lab, visit) %>%
  summarise(months = median(months), y = max(y) + 0.22, .groups = "drop")

# Arrow positions above the top group
n_groups <- length(lab_levels)
y_arrow_12 <- n_groups + 0.75
y_arrow_18 <- n_groups + 1.15

# Arrows start at -3 months (cognitive tests up to 3 months before the first MRI are included)
p_timing <- ggplot(plot_df, aes(months, y, colour = test, shape = test)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_vline(xintercept = 12, linetype = "dashed", colour = "grey50") +
  geom_errorbarh(aes(xmin = q25_m, xmax = q75_m), height = 0.08) +
  geom_point(size = 3) +
  geom_text(data = visit_labels, aes(months, y, label = visit), inherit.aes = FALSE,
            size = 2.8, colour = "black") +
  annotate("segment", x = -3, xend = 12, y = y_arrow_12, yend = y_arrow_12,
           arrow = arrow(length = unit(0.18, "cm"), ends = "both"), linewidth = 0.6) +
  annotate("label", x = 4.5, y = y_arrow_12, label = "12-month observation window",
           size = 3.1, label.size = 0, fill = "white") +
  annotate("segment", x = -3, xend = 18, y = y_arrow_18, yend = y_arrow_18,
           arrow = arrow(length = unit(0.18, "cm"), ends = "both"), linewidth = 0.6) +
  annotate("label", x = 7.5, y = y_arrow_18, label = "18-month observation window",
           size = 3.1, label.size = 0, fill = "white") +
  scale_y_continuous(breaks = seq_along(lab_levels), labels = lab_levels,
                     limits = c(0.6, y_arrow_18 + 0.3)) +
  scale_x_continuous(breaks = seq(-3, 18, 3), limits = c(-3.5, 18.5)) +
  scale_colour_manual(values = c("MRI scan" = "#2166AC", "ADAS-Cog 13" = "#B2182B"), name = NULL) +
  scale_shape_manual(values = c("MRI scan" = 16, "ADAS-Cog 13" = 17), name = NULL) +
  labs(title = "Timing of MRI and ADAS-Cog 13 assessments within the observation windows",
       x = "Months from first MRI", y = "Visit of first MRI (participants)") +
  theme_pub +
  theme(legend.position = "top")

p_timing
ggsave(file.path(fig_dir, "Fig_visit_timing_by_start_visit.png"), p_timing, width = 10, height = 6, dpi = 300)
###############

# Name each assessment by its timing, not by its ADNI code
pooled <- timing_long %>%
  mutate(months = days / 30.44,
         visit_t = case_when(months < 1.5 ~ "Baseline",
                             months < 4.5 ~ "Month 3",
                             months < 9   ~ "Month 6",
                             TRUE         ~ "Month 12"),
         test = factor(test, levels = c("ADAS-Cog 13", "MRI scan")))

# Median, IQR and number of assessments per visit and test
pooled_sum <- pooled %>%
  group_by(test, visit_t) %>%
  summarise(n = n(), med = median(months), q25 = quantile(months, 0.25),
            q75 = quantile(months, 0.75), .groups = "drop") %>%
  mutate(y = as.numeric(test), lab = paste0(visit_t, "\n(n = ", n, ")"))

y_a12 <- 2.75; y_a18 <- 3.15

p_pooled <- ggplot(pooled, aes(months, as.numeric(test), colour = test)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_jitter(height = 0.12, width = 0, alpha = 0.25, size = 1.2) +
  geom_errorbarh(data = pooled_sum, aes(xmin = q25, xmax = q75, y = y, x = med),
                 inherit.aes = FALSE, height = 0.2, colour = "black") +
  geom_point(data = pooled_sum, aes(med, y), inherit.aes = FALSE, size = 2.5, colour = "black") +
  geom_text(data = pooled_sum, aes(med, y + 0.33, label = lab), inherit.aes = FALSE, size = 2.8) +
  annotate("segment", x = -3, xend = 12, y = y_a12, yend = y_a12,
           arrow = arrow(length = unit(0.18, "cm"), ends = "both"), linewidth = 0.6) +
  annotate("label", x = 4.5, y = y_a12, label = "12-month observation window",
           size = 3.1, label.size = 0, fill = "white") +
  annotate("segment", x = -3, xend = 18, y = y_a18, yend = y_a18,
           arrow = arrow(length = unit(0.18, "cm"), ends = "both"), linewidth = 0.6) +
  annotate("label", x = 7.5, y = y_a18, label = "18-month observation window",
           size = 3.1, label.size = 0, fill = "white") +
  scale_y_continuous(breaks = 1:2, labels = levels(pooled$test), limits = c(0.6, 3.4)) +
  scale_x_continuous(breaks = seq(-3, 18, 3), limits = c(-3.5, 18.5)) +
  scale_colour_manual(values = c("MRI scan" = "#2166AC", "ADAS-Cog 13" = "#B2182B"), guide = "none") +
  labs(title = paste0("Timing of MRI and ADAS-Cog 13 assessments (N = ", n_distinct(pooled$RID), ")"),
       x = "Months from first MRI", y = "Assessment") +
  theme_pub

p_pooled
ggsave(file.path(fig_dir, "Fig_visit_timing_pooled.png"), p_pooled, width = 10, height = 5, dpi = 300)

###########Density timeline of MRI and ADAS-Cog 13 assessments, with visits named by order and window arrows underneath
# Nearest scheduled visit: MRI 0/3/6/12 months, ADAS 0/6/12 months
nearest_visit <- function(m, targets) targets[apply(abs(outer(m, targets, "-")), 1, which.min)]

pooled <- timing_long %>%
  mutate(months = days / 30.44) %>%
  mutate(target = if_else(test == "MRI scan",
                          nearest_visit(months, c(0, 3, 6, 12)),
                          nearest_visit(months, c(0, 6, 12))),
         visit = if_else(target == 0, "Baseline", paste("Month", target)))

# Integrity: nobody should have two assessments with the same visit label
dups <- pooled %>% count(RID, test, visit) %>% filter(n > 1)
cat("Duplicate person-visit labels:", nrow(dups), "\n")

visit_sum <- pooled %>%
  group_by(test, visit, target) %>%
  summarise(n = n(), med = median(months), .groups = "drop") %>%
  arrange(test, target)
print(visit_sum)

#######
# Density curves scaled to height 1; ADAS row at y = 0, MRI row at y = 1.4
dens_df <- pooled %>%
  group_by(test) %>%
  group_modify(~ { d <- density(.x$months, adjust = 0.35, from = -3.5, to = 18.5)
  tibble(x = d$x, h = d$y / max(d$y)) }) %>%
  ungroup() %>%
  mutate(base = if_else(test == "MRI scan", 1.4, 0), top = base + h)

row_base <- c("MRI scan" = 1.4, "ADAS-Cog 13" = 0)
test_cols <- c("MRI scan" = "#2166AC", "ADAS-Cog 13" = "#B2182B")

lab_df <- visit_sum %>%
  mutate(base = row_base[test], lab = paste0(visit, "\n(n = ", n, ")"))

arrow_df <- tibble(window = c("12-month observation window", "18-month observation window"),
                   xend = c(12, 18), y = c(2.75, 3.15),
                   fill = c("#FDE0C5", "#D9EAF7"))

# Visit labels without counts
lab_df <- visit_sum %>%
  mutate(base = row_base[test], lab = visit)

p_timeline <- ggplot() +
  geom_ribbon(data = dens_df, aes(x, ymin = base, ymax = top, fill = test), alpha = 0.35, colour = NA) +
  geom_line(data = dens_df, aes(x, top, colour = test), linewidth = 0.5) +
  geom_segment(data = lab_df, aes(x = med, xend = med, y = base, yend = base + 1.05, colour = test),
               linetype = "dashed", linewidth = 0.4) +
  geom_text(data = lab_df, aes(med, base + 1.15, label = lab), size = 3) +
  geom_vline(xintercept = 0, colour = "grey60", linewidth = 0.3) +
  geom_vline(xintercept = 12, colour = "grey60", linewidth = 0.3) +
  geom_segment(data = arrow_df, aes(x = -3, xend = xend, y = y, yend = y),
               arrow = arrow(length = unit(0.2, "cm"), ends = "both"), linewidth = 0.7) +
  geom_label(data = arrow_df, aes(x = (xend - 3) / 2, y = y, label = window, fill = I(fill)),
             size = 3.8, fontface = "bold", label.size = 0.4, label.padding = unit(0.25, "lines")) +
  annotate("text", x = -3.4, y = row_base + 0.35, label = names(row_base),
           hjust = 1, fontface = "bold", size = 3.6) +
  scale_fill_manual(values = test_cols, guide = "none") +
  scale_colour_manual(values = test_cols, guide = "none") +
  scale_x_continuous(breaks = seq(-3, 18, 3)) +
  scale_y_continuous(limits = c(0, 3.4), expand = c(0, 0)) +
  coord_cartesian(xlim = c(-3.5, 18.5), clip = "off") +
  labs(title = "Timing of MRI and ADAS-Cog 13 assessments within the observation windows",
       x = "Months from first MRI", y = NULL) +
  theme_pub +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.line.y = element_blank(),
        plot.margin = margin(10, 10, 10, 75))

p_timeline
ggsave(file.path(fig_dir, "Fig_visit_timing_density.png"), p_timeline, width = 11, height = 6, dpi = 300)

##########LME models check for cognitive rates#########
# Retest indicator: 0 at a person's first ADAS test in the window, 1 at later tests
if (exists("adas_6m") && exists("adas_12m")) {
  adas_win12 <- adas_6m          # 12-month observation window (N = 460 cohort)
  adas_win18 <- adas_12m         # 18-month observation window (N = 459 cohort)
} else {
  adas_win12 <- get_adas_dated(lme_12m,      1.0)
  adas_win18 <- get_adas_dated(lme_18m_data, 1.5)
}

# Integrity check: expect about 460/459 participants and 0 duplicate RID-visit rows
for (nm in c("adas_win12", "adas_win18")) {
  d <- get(nm)
  cat(nm, ": rows =", nrow(d), "| participants =", n_distinct(d$RID),
      "| duplicates =", sum(duplicated(d[c("RID", "ADAS_VISCODE2")])), "\n")
}




add_retest <- function(adas) {
  adas %>% group_by(RID) %>% filter(n() >= 2) %>%
    mutate(retest = as.numeric(time_years > min(time_years))) %>% ungroup() %>%
    mutate(BL_DIAGNOSIS = factor(BL_DIAGNOSIS, levels = c("CN", "MCI")))
}

cog_models <- function(adas, window) {
  d <- add_retest(adas)
  ctrl <- lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 20000))
  forms <- list(
    `Model 1: standard`        = ADAS13 ~ time_years * BL_DIAGNOSIS + (1 + time_years | RID),
    `Model 2: practice effect` = ADAS13 ~ time_years * BL_DIAGNOSIS + retest * BL_DIAGNOSIS + (1 + time_years | RID))
  
  ss <- d %>% group_by(RID) %>% summarise(ss = sum((time_years - mean(time_years))^2)) %>% pull(ss)
  
  map_dfr(names(forms), function(m) {
    fit_ml   <- lmer(forms[[m]], data = d, REML = FALSE, control = ctrl)   # for AIC/BIC
    fit_reml <- lmer(forms[[m]], data = d, REML = TRUE,  control = ctrl)   # for variance components
    fe  <- fixef(fit_reml)
    vc  <- as.data.frame(VarCorr(fit_reml))
    tau <- vc$sdcor[vc$grp == "RID" & vc$var1 == "time_years" & is.na(vc$var2)]
    tibble(window = window, model = m, N = n_distinct(d$RID), obs = nrow(d),
           slope_CN  = round(fe["time_years"], 3),
           slope_MCI = round(fe["time_years"] + fe["time_years:BL_DIAGNOSISMCI"], 3),
           practice_CN  = if ("retest" %in% names(fe)) round(fe["retest"], 3) else NA_real_,
           practice_MCI = if ("retest" %in% names(fe)) round(fe["retest"] + fe["BL_DIAGNOSISMCI:retest"], 3) else NA_real_,
           reliability = round(mean(tau^2 / (tau^2 + sigma(fit_reml)^2 / ss)), 3),
           singular = isSingular(fit_reml),
           AIC = round(AIC(fit_ml), 1), BIC = round(BIC(fit_ml), 1))
  })
}

cog_compare <- bind_rows(cog_models(adas_win12, "12-month window"),
                         cog_models(adas_win18, "18-month window"))
print(cog_compare, width = Inf)
write_csv(cog_compare, file.path(tab_dir, "Table_ADAS13_standard_vs_practice_models.csv"))

# Can the practice effect be told apart from time?
for (nm in c("adas_win12", "adas_win18")) {
  d <- add_retest(get(nm))
  cat(nm, ": r(retest, time) =", round(cor(d$retest, d$time_years), 3), "\n")
}

####compares all four practice-effect forms for ADAS-13 using the 24-month feasibility data.######
# Practice-effect terms: first-visit indicator, number of prior tests, square root of prior tests
add_pe_terms <- function(adas) {
  adas %>% group_by(RID) %>% filter(n() >= 2) %>% arrange(time_years, .by_group = TRUE) %>%
    mutate(n_prior = row_number() - 1,
           retest  = as.numeric(n_prior > 0),
           sqrt_prior = sqrt(n_prior)) %>% ungroup() %>%
    mutate(BL_DIAGNOSIS = factor(BL_DIAGNOSIS, levels = c("CN", "MCI")))
}

pe_compare <- function(adas, score, window) {
  d <- add_pe_terms(adas) %>% rename(y = all_of(score))
  ctrl <- lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 20000))
  base <- "y ~ time_years * BL_DIAGNOSIS"
  forms <- c(`None`           = base,
             `Jump (first visit)` = paste(base, "+ retest * BL_DIAGNOSIS"),
             `Hop (n prior)`      = paste(base, "+ n_prior * BL_DIAGNOSIS"),
             `Skip (sqrt prior)`  = paste(base, "+ sqrt_prior * BL_DIAGNOSIS"))
  ss <- d %>% group_by(RID) %>% summarise(ss = sum((time_years - mean(time_years))^2)) %>% pull(ss)
  
  map_dfr(names(forms), function(m) {
    f  <- as.formula(paste(forms[[m]], "+ (1 + time_years | RID)"))
    ml <- lmer(f, data = d, REML = FALSE, control = ctrl)
    re <- lmer(f, data = d, REML = TRUE,  control = ctrl)
    fe <- fixef(re); vc <- as.data.frame(VarCorr(re))
    tau <- vc$sdcor[vc$grp == "RID" & vc$var1 == "time_years" & is.na(vc$var2)]
    pe_term <- intersect(c("retest", "n_prior", "sqrt_prior"), names(fe))
    tibble(window = window, score = score, practice_form = m, N = n_distinct(d$RID),
           slope_CN = round(fe["time_years"], 3),
           slope_MCI = round(fe["time_years"] + fe["time_years:BL_DIAGNOSISMCI"], 3),
           practice_CN = if (length(pe_term)) round(fe[pe_term], 3) else NA_real_,
           r_practice_time = if (length(pe_term)) round(cor(d[[pe_term]], d$time_years), 3) else NA_real_,
           reliability = round(mean(tau^2 / (tau^2 + sigma(re)^2 / ss)), 3),
           singular = isSingular(re), AIC = round(AIC(ml), 1), BIC = round(BIC(ml), 1))
  })
}

pe_24 <- pe_compare(adas_win24_check, "ADAS13", "24-month (feasibility)")
print(pe_24, width = Inf)
write_csv(pe_24, file.path(tab_dir, "Table_ADAS13_practice_effect_forms_24m.csv"))













































