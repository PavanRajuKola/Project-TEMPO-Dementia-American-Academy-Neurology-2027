source(here::here("config.R"))

###Cognitive domain wise & Behaviour/functional domain wise trajectory
rda_dir <- ADNIMERGE2_DIR
load_rda <- function(name) { e <- new.env(); load(file.path(rda_dir, paste0(name, ".rda")), envir = e); get(name, envir = e) }

uw    <- load_rda("UWNPSYCHSUM")
npiq  <- load_rda("NPIQ")
ecogs <- load_rda("ECOGSP")

# Confirm the columns we expect exist
list(uw    = c("RID", "VISCODE2", "EXAMDATE", "ADNI_MEM", "ADNI_EF", "ADNI_LAN", "ADNI_VS") %in% names(uw),
     npiq  = c("RID", "VISCODE2", "VISDATE", "NPISCORE") %in% names(npiq),
     ecogs = c("RID", "VISCODE2", "VISDATE", "EcogSPTotal") %in% names(ecogs)) %>% print()

# Duplicate key check (RID + VISCODE2) and coverage in the 18-month cohort
check_tbl <- function(df, name) {
  tibble(file = name,
         rows = nrow(df),
         dup_RID_visit = sum(duplicated(df[c("RID", "VISCODE2")])),
         cohort_RIDs_present = sum(unique(lme_18m_data$RID) %in% df$RID))
}
bind_rows(check_tbl(uw, "UWNPSYCHSUM"), check_tbl(npiq, "NPIQ"), check_tbl(ecogs, "ECOGSP")) %>% print()

#####
dup_key <- ecogs %>% count(RID, VISCODE2) %>% filter(n > 1)

ecogs %>% semi_join(dup_key, by = c("RID", "VISCODE2")) %>%
  select(RID, VISCODE, VISCODE2, VISDATE, EcogSPTotal) %>%
  print()

# Is this participant in our cohort?
dup_key$RID %in% lme_18m_data$RID

#######
faq <- load_rda("FAQ")
c("RID", "VISCODE2", "VISDATE", "FAQTOTAL") %in% names(faq)   # confirm before using

get_dated <- function(df, date_col, score_col, label) {
  df %>% filter(!is.na(.data[[score_col]]), !is.na(.data[[date_col]])) %>%
    arrange(RID, VISCODE2, VISCODE) %>% distinct(RID, VISCODE2, .keep_all = TRUE) %>%
    transmute(RID, VISCODE2, test_date = as.Date(.data[[date_col]]),
              score = as.numeric(.data[[score_col]]), measure = label)
}

domain_long <- bind_rows(
  get_dated(adas_raw, "VISDATE",  "TOTAL13",     "ADAS-Cog 13"),
  get_dated(uw,       "EXAMDATE", "ADNI_MEM",    "Memory (ADNI-MEM)"),
  get_dated(uw,       "EXAMDATE", "ADNI_EF",     "Executive function (ADNI-EF)"),
  get_dated(uw,       "EXAMDATE", "ADNI_LAN",    "Language (ADNI-LAN)"),
  get_dated(uw,       "EXAMDATE", "ADNI_VS",     "Visuospatial (ADNI-VS)"),
  get_dated(faq,      "VISDATE",  "FAQTOTAL",    "Daily function (FAQ)"),
  get_dated(npiq,     "VISDATE",  "NPISCORE",    "Behaviour (NPI-Q)"),
  get_dated(ecogs,    "VISDATE",  "EcogSPTotal", "Everyday cognition (ECog-SP)")
)

# Integrity: no duplicate person-visit rows within a measure
domain_long %>% count(measure, RID, VISCODE2) %>% filter(n > 1) %>% nrow() %>% cat("Duplicates:", ., "\n")
domain_long %>% count(measure, name = "assessments") %>% print()

########Restrict to each window, count assessments per person, measure first-to-second change, and fit the LME for every measure
worse_if_higher <- c("ADAS-Cog 13", "Daily function (FAQ)", "Behaviour (NPI-Q)", "Everyday cognition (ECog-SP)")
ctrl <- lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 20000))

window_data <- function(lme_df, max_years) {
  first_mri(lme_df) %>%
    inner_join(domain_long, by = "RID", relationship = "one-to-many") %>%
    mutate(time_years = as.numeric(test_date - MRI_T0) / 365.25,
           BL_DIAGNOSIS = factor(BL_DIAGNOSIS, levels = c("CN", "MCI"))) %>%
    filter(time_years >= -0.25, time_years <= max_years)
}

domain_metrics <- function(d, window) {
  d %>% group_split(measure) %>% map_dfr(function(x) {
    m <- x$measure[1]
    x <- x %>% mutate(y = if (m %in% worse_if_higher) -score else score) %>%   # negative = getting worse
      group_by(RID) %>% filter(n() >= 2) %>% arrange(time_years, .by_group = TRUE) %>% ungroup()
    bl    <- x %>% group_by(RID) %>% slice(1) %>% ungroup()
    bl_sd <- sd(bl$y)
    per   <- x %>% count(RID)
    
    # First-to-second assessment change, in baseline SD units
    first_change <- x %>% group_by(RID, BL_DIAGNOSIS) %>%
      summarise(chg = (y[2] - y[1]) / bl_sd, .groups = "drop")
    
    fit <- tryCatch(lmer(y ~ time_years * BL_DIAGNOSIS + (1 + time_years | RID),
                         data = x, REML = TRUE, control = ctrl), error = function(e) NULL)
    if (is.null(fit)) return(tibble(window = window, measure = m, N = nrow(per), note = "model failed"))
    
    fe  <- fixef(fit); vc <- as.data.frame(VarCorr(fit))
    tau <- vc$sdcor[vc$grp == "RID" & vc$var1 == "time_years" & is.na(vc$var2)]
    ss  <- x %>% group_by(RID) %>% summarise(ss = sum((time_years - mean(time_years))^2)) %>% pull(ss)
    s_cn <- fe["time_years"] / bl_sd
    s_mci <- (fe["time_years"] + fe["time_years:BL_DIAGNOSISMCI"]) / bl_sd
    
    tibble(window = window, measure = m,
           N = nrow(per), with_2 = sum(per$n == 2), with_3plus = sum(per$n >= 3),
           pct_at_floor_bl = round(100 * mean(bl$score == min(x$score)), 1),
           std_slope_CN  = round(s_cn, 3), std_slope_MCI = round(s_mci, 3),
           MCI_minus_CN  = round(s_mci - s_cn, 3),
           first_change_CN  = round(mean(first_change$chg[first_change$BL_DIAGNOSIS == "CN"]), 3),
           first_change_MCI = round(mean(first_change$chg[first_change$BL_DIAGNOSIS == "MCI"]), 3),
           pct_improved_first_CN = round(100 * mean(first_change$chg[first_change$BL_DIAGNOSIS == "CN"] > 0), 1),
           reliability = round(mean(tau^2 / (tau^2 + sigma(fit)^2 / ss)), 3),
           singular = isSingular(fit))
  })
}

domain_results <- bind_rows(
  domain_metrics(window_data(lme_12m,      1.0), "12-month window"),
  domain_metrics(window_data(lme_18m_data, 1.5), "18-month window")
) %>% arrange(window, desc(reliability))

print(domain_results, n = Inf, width = Inf)
write_csv(domain_results, file.path(tab_dir, "Table_domain_slopes_reliability_12m_18m.csv"))


#######24 months check
domain_24 <- domain_metrics(window_data(lme_18m_data, 2.2), "24-month (feasibility)") %>%
  arrange(desc(reliability))

print(domain_24, n = Inf, width = Inf)
write_csv(domain_24, file.path(tab_dir, "Table_domain_slopes_reliability_24m_feasibility.csv"))

#########Verify the conversion time origin, then build the third-window cohort with the same exclusion rule
forward_groups <- c("CN→MCI", "CN→Dementia", "MCI→Dementia")

# Converters before a cut-off, counting forward conversions only
converters_before <- function(cutoff_years) {
  conversion_timeline %>%
    filter(trajectory_group %in% forward_groups, conversion_time_years < cutoff_years) %>%
    pull(RID) %>% unique()
}

# Check 1: does this rule reproduce the existing 18-month exclusion list?
check_18 <- converters_before(1.5)
cat("Rebuilt 18m list:", length(check_18), "| existing:", length(converters_before_18m),
    "| not matching:", length(setdiff(union(check_18, converters_before_18m),
                                      intersect(check_18, converters_before_18m))), "\n")

# Check 2: trajectory_cohort has every column used by the 18-month MRI data
missing_cols <- setdiff(names(lme_18m_data), names(trajectory_cohort))
cat("Columns in lme_18m_data but not in trajectory_cohort:", length(missing_cols), "\n")
print(head(missing_cols, 20))

# Build the third-window cohort: exclude converters before 2.2 years, keep scans <= 2.2 years, require >= 2 scans
converters_before_24m <- converters_before(2.2)

lme_24m_data <- trajectory_cohort %>%
  filter(!RID %in% converters_before_24m, time_years <= 2.2) %>%
  group_by(RID) %>% filter(n() >= 2) %>% ungroup()

# Integrity and flowchart numbers
cat("\nConverters before 2.2 years:", length(converters_before_24m), "\n")
conversion_timeline %>% filter(RID %in% converters_before_24m) %>%
  mutate(when = cut(conversion_time_years, c(0, 1, 1.5, 2.2), right = FALSE,
                    labels = c("< 12 m", "12–18 m", "18–26 m"))) %>%
  count(trajectory_group, when) %>% print()

cat("\nThird-window cohort: participants =", n_distinct(lme_24m_data$RID),
    "| scans =", nrow(lme_24m_data),
    "| duplicate RID-visit rows =", sum(duplicated(lme_24m_data[c("RID", "VISCODE2")])), "\n")

lme_24m_data %>% count(RID, name = "scans") %>% count(scans, name = "participants") %>% print()
lme_24m_data %>% distinct(RID, BL_DIAGNOSIS) %>% count(BL_DIAGNOSIS) %>% print()
lme_24m_data %>% count(VISCODE2) %>% print()

# Median follow-up (for the window name and Methods)
lme_24m_data %>% group_by(RID) %>%
  summarise(fu_m = as.numeric(max(as.Date(EXAMDATE)) - min(as.Date(EXAMDATE))) / 30.44) %>%
  summarise(median = round(median(fu_m), 1), q25 = round(quantile(fu_m, 0.25), 1),
            q75 = round(quantile(fu_m, 0.75), 1), min = round(min(fu_m), 1), max = round(max(fu_m), 1)) %>%
  print()

#####
print(missing_cols)

# Proceed only if every missing column is an ICV-normalised volume
stopifnot(all(str_detect(missing_cols, "_norm$")))

# Rebuild: TotalCortVol = cortical GM + cortical WM; every other column = its raw FreeSurfer volume
add_norm <- function(df) {
  df <- df %>% mutate(TotalCortVol_raw = ST149SV + ST152SV)
  for (col in missing_cols) {
    raw <- if (col == "TotalCortVol_norm") "TotalCortVol_raw" else str_remove(col, "_norm$")
    df[[col]] <- df[[raw]] / df$ST10CV * 100
  }
  df %>% select(-TotalCortVol_raw)
}

# Verify the formula: recompute in lme_18m_data and compare with the stored values
recheck <- add_norm(lme_18m_data %>% select(-all_of(missing_cols)))
max_diff <- map_dbl(missing_cols, ~ max(abs(recheck[[.x]] - lme_18m_data[[.x]]), na.rm = TRUE))
cat("Largest difference from stored values:", signif(max(max_diff), 3), "\n")   # expect ~0
stopifnot(max(max_diff) < 1e-8)

# Add the columns to the 24-month data and confirm nothing is missing now
lme_24m_data <- add_norm(lme_24m_data)
cat("Still missing:", length(setdiff(names(lme_18m_data), names(lme_24m_data))),
    "| rows:", nrow(lme_24m_data), "| participants:", n_distinct(lme_24m_data$RID), "\n")

#######
# Trajectory groups remaining in the 26-month cohort
lme_24m_data %>% distinct(RID) %>%
  left_join(conversion_timeline %>% select(RID, trajectory_group), by = "RID") %>%
  left_join(lme_24m_data %>% distinct(RID, BL_DIAGNOSIS), by = "RID") %>%
  mutate(group = case_when(is.na(trajectory_group) ~ paste0(BL_DIAGNOSIS, "-Stable"),
                           TRUE ~ trajectory_group)) %>%
  count(group) %>% print()

# The single 'nv' scan
lme_24m_data %>% filter(VISCODE2 == "nv") %>%
  select(RID, VISCODE2, EXAMDATE, time_years) %>% print()


#############
# Same composites and same rule as the other windows
lme_24m_comp <- add_to(add_composites(lme_24m_data))

comp_dx_24m <- map_dfr(c(composites_eval, "TemporoOccipital_bi"),
                       ~ lme_metrics_dx(lme_24m_comp, .x)) %>% tidy_labels()

rule_24m <- comp_dx_24m %>%
  mutate(is_ventricle = str_detect(region, "Ventric"),
         sign_ok = if_else(is_ventricle, APC_CN > 0 & APC_MCI > 0, APC_CN < 0 & APC_MCI < 0),
         keep = reliability_all >= 0.20 & !singular & sign_ok)

# Reliability across all three windows, side by side
rel_3win <- rule_24m %>%
  select(region, rel_26m = reliability_all, keep_26m = keep, APC_CN_26m = APC_CN, APC_MCI_26m = APC_MCI) %>%
  left_join(rule_compare %>% select(region, rel_12m, rel_18m), by = "region") %>%
  relocate(rel_12m, rel_18m, .after = region) %>%
  arrange(desc(rel_26m))

print(rel_3win, n = Inf, width = Inf)
write_csv(rel_3win, file.path(tab_dir, "Table_MRI_reliability_three_windows.csv"))

#########
# Rename the 24-month reliability columns and save
rel_3win <- rel_3win %>% rename(rel_24m = rel_26m, keep_24m = keep_26m,
                                APC_CN_24m = APC_CN_26m, APC_MCI_24m = APC_MCI_26m)
write_csv(rel_3win, file.path(tab_dir, "Table_MRI_reliability_three_windows.csv"))

# Cognitive and functional measures in the 24-month cohort (converters before 2.2 years excluded)
domain_24m <- domain_metrics(window_data(lme_24m_data, 2.2), "24-month window") %>%
  filter(measure %in% c("Memory (ADNI-MEM)", "ADAS-Cog 13", "Language (ADNI-LAN)",
                        "Executive function (ADNI-EF)", "Visuospatial (ADNI-VS)", "Daily function (FAQ)")) %>%
  arrange(desc(reliability))

# Compare with the earlier 24-month check, which still included converters
domain_24m %>%
  select(measure, N, with_3plus, std_slope_CN, std_slope_MCI, MCI_minus_CN, reliability, singular) %>%
  left_join(domain_24 %>% select(measure, reliability_with_converters = reliability), by = "measure") %>%
  print(width = Inf)

write_csv(domain_24m, file.path(tab_dir, "Table_domain_slopes_reliability_24m.csv"))


