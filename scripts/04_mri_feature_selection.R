source(here::here("config.R"))

#Feature selection

st_cols <- names(master_visit)[grepl("^ST", names(master_visit))]

cat("Total ST columns in master_visit:", length(st_cols), "\n\n")

# Categorize by suffix
st_types <- tibble(column = st_cols) %>%
  mutate(suffix = str_extract(column, "[A-Z]+$")) %>%
  count(suffix, sort = TRUE)

cat("Column types by suffix:\n")
print(st_types)

# Print ALL ST columns grouped by suffix
cat("\n\nAll SV (Subcortical Volume) columns:\n")
print(sort(st_cols[grepl("SV$", st_cols)]))

cat("\n\nAll CV (Cortical Volume) columns:\n")
print(sort(st_cols[grepl("CV$", st_cols)]))

cat("\n\nAll TA (Thickness Average) columns:\n")
print(sort(st_cols[grepl("TA$", st_cols)]))

cat("\n\nAll HS (Hippocampal Subfield) columns:\n")
print(sort(st_cols[grepl("HS$", st_cols)]))


# Bilateral regional features
feature_map <- tribble(
  ~region,              ~left,       ~right,
  "Hippocampus",        "ST29SV",    "ST88SV",
  "Amygdala",           "ST12SV",    "ST71SV",
  "Thalamus",           "ST61SV",    "ST120SV",
  "Putamen",            "ST53SV",    "ST112SV",
  "Caudate",            "ST16SV",    "ST75SV",
  "NucAccumbens",       "ST11SV",    "ST70SV",
  "Pallidum",           "ST42SV",    "ST101SV",
  "LatVentricle",       "ST37SV",    "ST96SV",
  "InfLatVentricle",    "ST30SV",    "ST89SV",
  "Entorhinal",         "ST24CV",    "ST83CV",
  "Fusiform",           "ST26CV",    "ST85CV",
  "MiddleTemporal",     "ST40CV",    "ST99CV",
  "Precuneus",          "ST52CV",    "ST111CV",
  "Insula",             "ST129CV",   "ST130CV"
)

# Global columns needed
global_cols <- c("ST149SV", "ST152SV", "ST10CV")

# Verify all columns exist
all_cols <- c(feature_map$left, feature_map$right, global_cols)
cat("Columns found:", sum(all_cols %in% names(master_visit)), "/", length(all_cols), "\n")
missing <- all_cols[!all_cols %in% names(master_visit)]
if (length(missing) > 0) { print(missing) } else { cat("All columns present ✓\n") }

# Compute Total Cortical Volume = Cortical GM + Cortical WM
master_visit <- master_visit %>%
  mutate(TotalCorticalVol = ST149SV + ST152SV)

cat("\nTotalCorticalVol computed (ST149SV + ST152SV)\n")
cat("  Mean:", round(mean(master_visit$TotalCorticalVol, na.rm = TRUE), 0), "mm³\n")
cat("  Cortical GM mean:", round(mean(master_visit$ST149SV, na.rm = TRUE), 0), "mm³\n")
cat("  Cortical WM mean:", round(mean(master_visit$ST152SV, na.rm = TRUE), 0), "mm³\n")

# All feature columns including computed TotalCorticalVol
all_feature_cols <- c(feature_map$left, feature_map$right,
                      "ST149SV", "TotalCorticalVol", "ST10CV")

name_lookup <- c(
  setNames(paste0("L_", feature_map$region), feature_map$left),
  setNames(paste0("R_", feature_map$region), feature_map$right),
  ST149SV          = "CorticalGM",
  TotalCorticalVol = "TotalCorticalVol",
  ST10CV           = "ICV"
)

missingness <- tibble(
  code   = all_feature_cols,
  region = name_lookup[all_feature_cols],
  n_total   = nrow(master_visit),
  n_missing = map_int(all_feature_cols, ~ sum(is.na(master_visit[[.x]]))),
  pct_missing = round(n_missing / n_total * 100, 2)
) %>%
  arrange(desc(pct_missing))

print(missingness, n = Inf)

cat("\nSummary:\n")
cat("  Total feature columns:", length(all_feature_cols), "\n")
cat("  Features with 0% missing:", sum(missingness$pct_missing == 0), "\n")
cat("  Features with >0% missing:", sum(missingness$pct_missing > 0), "\n")
cat("  Max missingness:", max(missingness$pct_missing), "%\n")

#correlation test between our 31 features
library(corrplot)

bilateral_cols <- c(feature_map$left, feature_map$right)

# Build renaming vector: old name → new name
rename_vec <- c(
  setNames(paste0("L.", feature_map$region), feature_map$left),
  setNames(paste0("R.", feature_map$region), feature_map$right)
)

cor_data <- baseline_mri %>%
  select(all_of(bilateral_cols))

# Rename columns for readability
names(cor_data) <- rename_vec[names(cor_data)]

cor_matrix <- cor(cor_data, use = "pairwise.complete.obs")

png("mri_feature_correlation.png", width = 12, height = 12, units = "in", res = 300)
corrplot(cor_matrix, method = "color", type = "lower",
         tl.col = "black", tl.cex = 0.65, tl.srt = 45,
         addCoef.col = "black", number.cex = 0.4,
         col = colorRampPalette(c("#2166AC", "white", "#B2182B"))(200),
         title = "Correlation Matrix — Baseline MRI Volumes (602 participants)",
         mar = c(0, 0, 2, 0))
dev.off()

cat("Plot saved: mri_feature_correlation.png\n")

#checking the correlation coefficient values
lr_cors <- feature_map %>%
  mutate(cor = map2_dbl(left, right, ~ cor(
    baseline_mri[[.x]], baseline_mri[[.y]], use = "pairwise.complete.obs"
  ))) %>%
  select(region, cor) %>%
  arrange(cor)

cat("Left-Right Hemisphere Correlations:\n\n")
print(lr_cors, n = Inf)

cat("\nAll L-R correlations > 0.90:", all(lr_cors$cor > 0.90), "\n")
cat("Min L-R correlation:", round(min(lr_cors$cor), 3),
    "(", lr_cors$region[which.min(lr_cors$cor)], ")\n")



########### Verify occipital columns############
cat("Left Lateral Occipital (ST35CV):", "ST35CV" %in% names(master_visit), "\n")
cat("Right Lateral Occipital (ST94CV):", "ST94CV" %in% names(master_visit), "\n\n")

# Complete feature map with occipital added
feature_map <- tribble(
  ~region,              ~left,       ~right,
  "Hippocampus",        "ST29SV",    "ST88SV",
  "Amygdala",           "ST12SV",    "ST71SV",
  "Thalamus",           "ST61SV",    "ST120SV",
  "Putamen",            "ST53SV",    "ST112SV",
  "Caudate",            "ST16SV",    "ST75SV",
  "NucAccumbens",       "ST11SV",    "ST70SV",
  "Pallidum",           "ST42SV",    "ST101SV",
  "LatVentricle",       "ST37SV",    "ST96SV",
  "InfLatVentricle",    "ST30SV",    "ST89SV",
  "Entorhinal",         "ST24CV",    "ST83CV",
  "Fusiform",           "ST26CV",    "ST85CV",
  "MiddleTemporal",     "ST40CV",    "ST99CV",
  "Precuneus",          "ST52CV",    "ST111CV",
  "Insula",             "ST129CV",   "ST130CV",
  "LateralOccipital",   "ST35CV",    "ST94CV"
)

global_cols <- c(CorticalGM = "ST149SV", CorticalWM = "ST152SV", ICV = "ST10CV")

# Verify all exist
all_cols <- c(feature_map$left, feature_map$right, global_cols)
cat("Feature columns found:", sum(all_cols %in% names(master_visit)), "/", length(all_cols), "\n")

# Confirm analytical cohort: exclude baseline Dementia, MCI→CN, and early converters
analytical_cohort <- trajectory_cohort %>%
  filter(!RID %in% early_converter_rids)

cat("\nAnalytical cohort:\n")
cat("  Participants:", n_distinct(analytical_cohort$RID), "\n")
cat("  Observations:", nrow(analytical_cohort), "\n\n")

analytical_cohort %>%
  group_by(traj_group) %>%
  summarise(n_participants = n_distinct(RID), n_obs = n(), .groups = "drop") %>%
  print()




