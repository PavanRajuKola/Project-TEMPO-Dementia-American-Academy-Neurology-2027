source(here::here("config.R"))

############
 ============================================================
# ADNI FreeSurfer DATA AUDIT
 ============================================================

# Path to your ADNIMERGE output folder
data_dir <- ADNIMERGE_DIR

# Check that the folder exists
dir.exists(data_dir)

# List all files
files <- list.files(
  path = data_dir,
  recursive = TRUE,
  full.names = TRUE
)

# Display files
files

 ============================================================
# 1. LOAD ALL FREESURFER TABLES
============================================================

data_dir <- ADNIMERGE_DIR

fs_files <- c(
  "UCSFFSX.csv",
  "UCSFFSX51.csv",
  "UCSFFSL.csv",
  "UCSFFSL51.csv"
)

fs_list <- lapply(fs_files, function(f) {
  read.csv(
    file.path(data_dir, f),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
})

names(fs_list) <- fs_files

# ============================================================
# 2. BASIC STRUCTURE OF EACH FREESURFER DATASET
# ============================================================

for (i in seq_along(fs_list)) {
  
  x <- fs_list[[i]]

  cat("\n\n============================================\n")
  cat("FILE:", names(fs_list)[i], "\n")
  cat("============================================\n")
  
  cat("Rows:", nrow(x), "\n")
  cat("Columns:", ncol(x), "\n")
  
  cat("\nFirst 20 columns:\n")
  print(names(x)[1:min(20, ncol(x))])
  
  cat("\nID-related columns:\n")
  print(names(x)[grepl("RID|PTID|SUBJ|ID", names(x), ignore.case = TRUE)])
  
  cat("\nVisit/date-related columns:\n")
  print(names(x)[grepl("VIS|DATE|EXAM", names(x), ignore.case = TRUE)])
  
  cat("\nVersion-related columns:\n")
  print(names(x)[grepl("VERSION|FSVER|FREESURFER", names(x), ignore.case = TRUE)])
}

 ============================================================
# 3. NUMBER OF UNIQUE SUBJECTS
 ============================================================

for (i in seq_along(fs_list)) {
  
  x <- fs_list[[i]]
  
  rid_col <- names(x)[
    grepl("^RID$", names(x), ignore.case = TRUE)
  ]
  
  cat(
    "\n",
    names(fs_list)[i],
    ": ",
    ifelse(
      length(rid_col) > 0,
      length(unique(x[[rid_col[1]]])),
      "RID column not found"
    ),
    "unique RIDs\n",
    sep = ""
  )
}

# ============================================================
# AUDIT FREESURFER VERSION + FIELD STRENGTH
# ============================================================

for (nm in c("UCSFFSX.csv", "UCSFFSX51.csv",
             "UCSFFSL.csv", "UCSFFSL51.csv")) {
  
  x <- fs_list[[nm]]
  
  cat("\n\n============================================\n")
  cat(nm, "\n")
  cat("============================================\n")
  
  cat("\nVERSION:\n")
  print(table(x$VERSION, useNA = "ifany"))
  
  if ("FLDSTRENG" %in% names(x)) {
    cat("\nFIELD STRENGTH:\n")
    print(table(x$FLDSTRENG, useNA = "ifany"))
  } else {
    cat("\nNo FLDSTRENG column in this table\n")
  }
  
  cat("\nSTATUS:\n")
  print(table(x$STATUS, useNA = "ifany"))
  
  if ("OVERALLQC" %in% names(x)) {
    cat("\nOVERALL QC:\n")
    print(table(x$OVERALLQC, useNA = "ifany"))
  }
}

 ============================================================
# LONGITUDINAL AVAILABILITY IN CROSS-SECTIONAL FREESURFER
============================================================

for (nm in c("UCSFFSX.csv", "UCSFFSX51.csv")) {
  
  x <- fs_list[[nm]]
  
  # Make sure date is Date
  x$EXAMDATE <- as.Date(x$EXAMDATE)
  
  # Remove rows without RID/date
  x <- x[!is.na(x$RID) & !is.na(x$EXAMDATE), ]
  
  # Number of observations per subject
  n_visits <- aggregate(
    EXAMDATE ~ RID,
    data = x,
    FUN = function(z) length(unique(z))
  )
  
  names(n_visits)[2] <- "n_FS_visits"
  
  # Follow-up in months
  followup <- aggregate(
    EXAMDATE ~ RID,
    data = x,
    FUN = function(z) {
      z <- unique(z)
      
      if (length(z) >= 2) {
        as.numeric(max(z) - min(z)) / 30.44
      } else {
        NA
      }
    }
  )
  
  names(followup)[2] <- "followup_months"
  
  audit <- merge(
    n_visits,
    followup,
    by = "RID"
  )
  
  # Eligibility matching your current study
  eligible <- audit[
    audit$n_FS_visits >= 3 &
      audit$followup_months >= 12,
  ]
  
  cat("\n============================================\n")
  cat(nm, "\n")
  cat("============================================\n")
  
  cat("Unique subjects:", nrow(audit), "\n")
  
  cat(
    "Subjects with >=3 FS visits:",
    sum(audit$n_FS_visits >= 3),
    "\n"
  )
  
  cat(
    "Subjects with >=3 FS visits AND >=1 year:",
    nrow(eligible),
    "\n"
  )
  
  cat("\nDistribution of number of FS visits:\n")
  print(table(audit$n_FS_visits))
  
  cat("\nFollow-up summary (months):\n")
  print(summary(audit$followup_months))
}

fs51 <- fs_list[["UCSFFSX51.csv"]]

# Examine protocol variables
cat("ORIGPROT:\n")
print(table(fs51$ORIGPROT, useNA = "ifany"))

cat("\nCOLPROT:\n")
print(table(fs51$COLPROT, useNA = "ifany"))

cat("\nIMAGETYPE:\n")
print(table(fs51$IMAGETYPE, useNA = "ifany"))

cat("\nVISCODE:\n")
print(table(fs51$VISCODE, useNA = "ifany"))

cat("\nVISCODE2:\n")
print(table(fs51$VISCODE2, useNA = "ifany"))

# Look at a sample of the relevant metadata
fs51[, c(
  "RID",
  "VISCODE",
  "VISCODE2",
  "EXAMDATE",
  "ORIGPROT",
  "COLPROT",
  "IMAGETYPE",
  "LONISID",
  "LONIUID",
  "IMAGEUID",
  "STATUS",
  "OVERALLQC"
)] |> head(30)

library(dplyr)

# Count FS51 records per participant + visit
dup_visits <- fs51 %>%
  group_by(RID, VISCODE2, EXAMDATE) %>%
  summarise(
    n_records = n(),
    image_types = paste(unique(IMAGETYPE), collapse = " | "),
    colprots = paste(unique(COLPROT), collapse = " | "),
    statuses = paste(unique(STATUS), collapse = " | "),
    qc = paste(unique(OVERALLQC), collapse = " | "),
    .groups = "drop"
  ) %>%
  filter(n_records > 1)

# How many participant-visits have duplicates?
nrow(dup_visits)

# Show the duplicated visits
head(dup_visits, 30)

# Distribution of number of records per RID + visit
fs51 %>%
  count(RID, VISCODE2, EXAMDATE, name = "n_records") %>%
  count(n_records)

# Examine all duplicated participant-visits in detail
dup_detail <- fs51 %>%
  semi_join(
    dup_visits,
    by = c("RID", "VISCODE2", "EXAMDATE")
  ) %>%
  select(
    RID, VISCODE2, EXAMDATE,
    ORIGPROT, COLPROT, IMAGETYPE,
    LONISID, LONIUID, IMAGEUID,
    STATUS, OVERALLQC
  ) %>%
  arrange(RID, EXAMDATE)

print(dup_detail, n = 100)

dup_detail |> 
  as.data.frame() |> 
  head(100)

dup_detail %>%
  count(IMAGETYPE, sort = TRUE)

dup_detail %>%
  count(COLPROT, IMAGETYPE, sort = TRUE)

dup_detail %>%
  count(STATUS, OVERALLQC, sort = TRUE)

dup_detail %>%
  group_by(RID, VISCODE2, EXAMDATE) %>%
  summarise(
    n_images = n_distinct(IMAGEUID),
    n_loniuid = n_distinct(LONIUID),
    image_types = paste(unique(IMAGETYPE), collapse = " | "),
    .groups = "drop"
  ) %>%
  count(n_images, image_types, sort = TRUE)

grep(
  "hippo|entorh|amyg|thalam|ventr|caud|cort|volume|thick",
  names(fs51),
  ignore.case = TRUE,
  value = TRUE
)

data.frame(
  index = seq_along(names(fs51)),
  column = names(fs51)
)

str(fs51, max.level = 1)
#######################################

# Load the packages needed for reading, cleaning, summarising, and visualising the longitudinal FreeSurfer data
library(tidyverse)

# Define the path to the UCSFFSL51 longitudinal FreeSurfer dataset
fs_long_path <- file.path(ADNIMERGE_DIR, "UCSFFSL51.csv")

# Read the longitudinal FreeSurfer dataset into R
fs_long <- read_csv(fs_long_path, show_col_types = FALSE)

# Check the basic structure and dimensions of the longitudinal FreeSurfer dataset
dim(fs_long)
glimpse(fs_long)

# Count the number of unique participants, observations, and participant-visit combinations
fs_long %>%
  summarise(
    total_rows = n(),
    unique_RIDs = n_distinct(RID),
    unique_visits = n_distinct(RID, VISCODE2),
    unique_exams = n_distinct(RID, EXAMDATE)
  )

# Examine the main identification, visit, date, protocol, and quality-control variables
fs_long %>%
  select(RID, VISCODE, VISCODE2, EXAMDATE, ORIGPROT, COLPROT, 
         LONISID, LONIUID, IMAGEUID, RUNDATE, STATUS, OVERALLQC) %>%
  glimpse()

# Check whether any participant-visit combinations occur more than once
duplicate_visits <- fs_long %>%
  count(RID, VISCODE2, EXAMDATE, name = "n_records") %>%
  filter(n_records > 1) %>%
  arrange(desc(n_records))

duplicate_visits

# Summarise how many participant-visits have one or multiple FreeSurfer records
fs_long %>%
  count(RID, VISCODE2, EXAMDATE, name = "n_records") %>%
  count(n_records, name = "number_of_participant_visits") %>%
  arrange(n_records)

# Calculate the number of longitudinal FreeSurfer visits available for each participant
subject_visit_counts <- fs_long %>%
  distinct(RID, VISCODE2, EXAMDATE) %>%
  group_by(RID) %>%
  summarise(
    n_visits = n(),
    first_exam_date = min(EXAMDATE, na.rm = TRUE),
    last_exam_date = max(EXAMDATE, na.rm = TRUE),
    .groups = "drop"
  )

# Examine the distribution of the number of longitudinal FreeSurfer visits per participant
subject_visit_counts %>%
  count(n_visits, name = "number_of_subjects") %>%
  arrange(n_visits)

# Identify participants who have at least three longitudinal FreeSurfer visits
subjects_3plus_visits <- subject_visit_counts %>%
  filter(n_visits >= 3)

# Count participants with at least three longitudinal FreeSurfer visits
subjects_3plus_visits %>%
  summarise(
    subjects_with_3plus_visits = n()
  )

# Calculate the actual follow-up duration for each participant in months
subject_followup <- fs_long %>%
  distinct(RID, VISCODE2, EXAMDATE) %>%
  filter(!is.na(EXAMDATE)) %>%
  group_by(RID) %>%
  summarise(
    n_visits = n(),
    first_exam_date = min(EXAMDATE),
    last_exam_date = max(EXAMDATE),
    followup_months = as.numeric(
      difftime(last_exam_date, first_exam_date, units = "days")
    ) / 30.4375,
    .groups = "drop"
  )

# Summarise the distribution of longitudinal follow-up duration across participants
subject_followup %>%
  summarise(
    minimum_months = min(followup_months, na.rm = TRUE),
    q1_months = quantile(followup_months, 0.25, na.rm = TRUE),
    median_months = median(followup_months, na.rm = TRUE),
    mean_months = mean(followup_months, na.rm = TRUE),
    q3_months = quantile(followup_months, 0.75, na.rm = TRUE),
    maximum_months = max(followup_months, na.rm = TRUE)
  )

# Identify participants who satisfy both the minimum three visits and minimum one-year follow-up criteria
eligible_longitudinal_subjects <- subject_followup %>%
  filter(
    n_visits >= 3,
    followup_months >= 12
  )

# Count participants who satisfy the three-visit and one-year follow-up criteria
eligible_longitudinal_subjects %>%
  summarise(
    eligible_subjects = n()
  )

# Examine the distribution of follow-up duration among participants eligible for longitudinal analysis
eligible_longitudinal_subjects %>%
  summarise(
    minimum_months = min(followup_months),
    q1_months = quantile(followup_months, 0.25),
    median_months = median(followup_months),
    mean_months = mean(followup_months),
    q3_months = quantile(followup_months, 0.75),
    maximum_months = max(followup_months)
  )

# Examine the FreeSurfer processing status of all longitudinal records
fs_long %>%
  count(STATUS, sort = TRUE)

# Examine the overall FreeSurfer quality-control status of all longitudinal records
fs_long %>%
  count(OVERALLQC, sort = TRUE)

# Examine the combination of processing status and overall quality-control status
fs_long %>%
  count(STATUS, OVERALLQC, sort = TRUE)

# Calculate the quality-control profile specifically among participants eligible for longitudinal analysis
eligible_rids <- eligible_longitudinal_subjects %>%
  select(RID)

fs_long %>%
  semi_join(eligible_rids, by = "RID") %>%
  count(STATUS, OVERALLQC, sort = TRUE)

# Examine how many eligible participants have at least one failed or incomplete FreeSurfer record
eligible_qc_summary <- fs_long %>%
  semi_join(eligible_rids, by = "RID") %>%
  group_by(RID) %>%
  summarise(
    n_records = n(),
    n_pass = sum(OVERALLQC == "Pass", na.rm = TRUE),
    n_partial = sum(OVERALLQC == "Partial", na.rm = TRUE),
    n_fail = sum(OVERALLQC == "Fail", na.rm = TRUE),
    n_hippocampus_only = sum(OVERALLQC == "Hippocampus Only", na.rm = TRUE),
    .groups = "drop"
  )

# Summarise the participant-level quality-control profile among eligible subjects
eligible_qc_summary %>%
  summarise(
    eligible_subjects = n(),
    subjects_with_any_pass = sum(n_pass > 0),
    subjects_with_any_partial = sum(n_partial > 0),
    subjects_with_any_fail = sum(n_fail > 0),
    subjects_with_hippocampus_only = sum(n_hippocampus_only > 0)
  )

# Examine the ADNI protocols represented in the longitudinal FreeSurfer dataset
fs_long %>%
  count(ORIGPROT, COLPROT, sort = TRUE)

# Examine the longitudinal FreeSurfer visit codes represented in the dataset
fs_long %>%
  count(VISCODE2, sort = TRUE)

# Examine the number of records and unique participants contributed by each protocol
fs_long %>%
  group_by(ORIGPROT, COLPROT) %>%
  summarise(
    records = n(),
    unique_RIDs = n_distinct(RID),
    .groups = "drop"
  ) %>%
  arrange(desc(unique_RIDs))

# Check the availability of examination dates and identify records without dates
fs_long %>%
  summarise(
    total_records = n(),
    records_with_examdate = sum(!is.na(EXAMDATE)),
    records_without_examdate = sum(is.na(EXAMDATE)),
    percentage_with_examdate = mean(!is.na(EXAMDATE)) * 100
  )

# Check whether participant identifiers or visit identifiers are missing
fs_long %>%
  summarise(
    missing_RID = sum(is.na(RID)),
    missing_VISCODE2 = sum(is.na(VISCODE2)),
    missing_EXAMDATE = sum(is.na(EXAMDATE)),
    missing_IMAGEUID = sum(is.na(IMAGEUID))
  )

# Identify the FreeSurfer variables available after the metadata and quality-control columns
fs_long %>%
  select(RID, VISCODE2, EXAMDATE, everything()) %>%
  select(RID, VISCODE2, EXAMDATE, starts_with("ST")) %>%
  glimpse()

# Count the number of FreeSurfer measurement variables available in the longitudinal dataset
fs_long %>%
  select(starts_with("ST")) %>%
  summarise(
    number_of_FreeSurfer_variables = ncol(.)
  )

# Examine missingness across the FreeSurfer measurement variables
fs_long %>%
  summarise(
    across(
      starts_with("ST"),
      ~ mean(is.na(.)) * 100
    )
  ) %>%
  pivot_longer(
    cols = everything(),
    names_to = "variable",
    values_to = "missing_percentage"
  ) %>%
  arrange(desc(missing_percentage))

# Create a participant-level table containing the longitudinal coverage needed for later trajectory analysis
longitudinal_coverage <- subject_followup %>%
  mutate(
    has_3plus_visits = n_visits >= 3,
    has_1year_followup = followup_months >= 12,
    eligible_for_trajectory = has_3plus_visits & has_1year_followup
  )

# Display the final longitudinal coverage summary
longitudinal_coverage %>%
  summarise(
    total_subjects = n(),
    subjects_3plus_visits = sum(has_3plus_visits),
    subjects_1year_followup = sum(has_1year_followup),
    subjects_3plus_visits_1year = sum(eligible_for_trajectory)
  )

View(longitudinal_coverage)

# Save the participant-level longitudinal eligibility table for use in subsequent analysis
write_csv(
  longitudinal_coverage,
  file.path(ADNIMERGE_DIR, "UCSFFSL51_longitudinal_coverage.csv")
)


###### QC check #####

# Keep only FreeSurfer records from the 644 eligible participants
fs_long_644 <- fs_long %>%
  semi_join(
    eligible_longitudinal_subjects,
    by = "RID"
  )

# Examine the QC distribution before applying the MRI quality-control filter
fs_long_644 %>%
  count(OVERALLQC, sort = TRUE)

# Examine the combination of processing status and overall FreeSurfer quality control before filtering
fs_long_644 %>%
  count(STATUS, OVERALLQC, sort = TRUE)

# Examine the FreeSurfer processing status
fs_long_644 %>%
  count(STATUS, sort = TRUE)

# Examine processing status together with overall quality control
fs_long_644 %>%
  count(STATUS, OVERALLQC, sort = TRUE)

# Keep scans with acceptable overall QC
fs_qc <- fs_long_644 %>%
  filter(OVERALLQC %in% c("Pass", "Partial"))

# Check QC distribution and number of subjects/records remaining
fs_qc %>%
  count(OVERALLQC, sort = TRUE)

fs_qc %>%
  summarise(
    subjects = n_distinct(RID),
    records = n()
  )

# Reassess longitudinal eligibility after QC
subject_visit_counts <- fs_qc %>%
  count(RID, name = "n_visits")

subject_visit_counts %>%
  count(n_visits, sort = TRUE)

# Retain subjects with ≥3 QC-eligible visits
eligible_qc_subjects <- subject_visit_counts %>%
  filter(n_visits >= 3)

fs_qc_longitudinal <- fs_qc %>%
  semi_join(eligible_qc_subjects, by = "RID")

# Final cohort size
fs_qc_longitudinal %>%
  summarise(
    subjects = n_distinct(RID),
    records = n()
  )

# Publication-quality visualization of FreeSurfer QC and longitudinal visit distribution

library(ggplot2)

# Distribution of QC-eligible visits per participant

subject_visit_counts <- fs_qc %>%
  count(RID, name = "n_visits") %>%
  count(n_visits, name = "n_subjects")

ggplot(subject_visit_counts, aes(x = n_visits, y = n_subjects)) +
  geom_col(width = 0.7, fill = "grey35") +
  geom_text(
    aes(label = n_subjects),
    vjust = -0.4,
    size = 4
  ) +
  scale_x_continuous(
    breaks = seq(
      min(subject_visit_counts$n_visits),
      max(subject_visit_counts$n_visits)
    )
  ) +
  labs(
    x = "Number of QC-eligible MRI visits",
    y = "Number of participants"
  ) +
  theme_classic(base_size = 13) +
  theme(
    axis.title = element_text(size = 13),
    axis.text = element_text(size = 11),
    plot.margin = margin(10, 15, 10, 10)
  ) +
  expand_limits(y = max(subject_visit_counts$n_subjects) * 1.08)

###linking diagnosis data#####

# Load the DXSUM clinical diagnosis data
library(tidyverse)

dxsum <- read_csv(
  file.path(ADNIMERGE_DIR, "DXSUM.csv")
)

str(dxsum)
colnames(dxsum)

# Prepare clinical diagnosis history
dx_history <- dxsum %>%
  select(RID, VISCODE2, EXAMDATE, DIAGNOSIS) %>%
  filter(!is.na(EXAMDATE), !is.na(DIAGNOSIS))

# Generate all MRI–clinical pairs within ±90 days
dx_candidates <- fs_qc_longitudinal %>%
  select(RID, MRI_VISCODE2 = VISCODE2, MRI_EXAMDATE = EXAMDATE) %>%
  distinct() %>%
  left_join(dx_history, by = "RID", relationship = "many-to-many") %>%
  mutate(date_gap_days = abs(as.numeric(MRI_EXAMDATE - EXAMDATE))) %>%
  filter(date_gap_days <= 90)

# Identify the nearest clinical assessment(s) for each MRI visit
dx_nearest <- dx_candidates %>%
  group_by(RID, MRI_VISCODE2, MRI_EXAMDATE) %>%
  filter(date_gap_days == min(date_gap_days)) %>%
  mutate(
    n_nearest = n(),
    n_diagnoses_nearest = n_distinct(DIAGNOSIS),
    diagnosis_conflict = n_diagnoses_nearest > 1
  ) %>%
  ungroup()

# Summarise linkage quality before resolving any discrepancies
dx_nearest %>%
  summarise(
    MRI_visits = n_distinct(paste(RID, MRI_VISCODE2, MRI_EXAMDATE)),
    linked_MRI_visits = n_distinct(
      paste(RID, MRI_VISCODE2, MRI_EXAMDATE)
    ),
    tied_nearest_visits = n_distinct(
      paste(RID, MRI_VISCODE2, MRI_EXAMDATE)[n_nearest > 1]
    ),
    conflicting_nearest_visits = n_distinct(
      paste(RID, MRI_VISCODE2, MRI_EXAMDATE)[diagnosis_conflict]
    )
  )

# Inspect MRI visits with conflicting equally-nearest diagnoses
diagnosis_discrepancies <- dx_nearest %>%
  filter(diagnosis_conflict) %>%
  arrange(RID, MRI_EXAMDATE, date_gap_days)

# Create the unambiguous one-to-one MRI–diagnosis linkage
fs_clinical_90day_linked <- dx_nearest %>%
  filter(n_nearest == 1 | !diagnosis_conflict) %>%
  group_by(RID, MRI_VISCODE2, MRI_EXAMDATE) %>%
  slice_head(n = 1) %>%
  ungroup() %>%
  rename(
    DX_VISCODE2 = VISCODE2,
    DX_EXAMDATE = EXAMDATE
  )

# Assess diagnosis coverage after linkage
diagnosis_coverage_90day <- fs_qc_longitudinal %>%
  select(RID, VISCODE2, EXAMDATE) %>%
  distinct() %>%
  left_join(
    fs_clinical_90day_linked %>%
      select(RID, MRI_VISCODE2, MRI_EXAMDATE, DIAGNOSIS),
    by = c(
      "RID",
      "VISCODE2" = "MRI_VISCODE2",
      "EXAMDATE" = "MRI_EXAMDATE"
    )
  ) %>%
  group_by(RID) %>%
  summarise(
    n_MRI_visits = n(),
    n_with_diagnosis = sum(!is.na(DIAGNOSIS)),
    n_without_diagnosis = sum(is.na(DIAGNOSIS)),
    .groups = "drop"
  )

# Participants with complete diagnosis linkage for every MRI visit
complete_diagnosis_subjects <- diagnosis_coverage_90day %>%
  filter(n_without_diagnosis == 0) %>%
  select(RID)

# Restrict to participants who also satisfy the longitudinal criteria
complete_longitudinal <- complete_diagnosis_subjects %>%
  semi_join(
    eligible_longitudinal_subjects,
    by = "RID"
  )

# Baseline diagnosis distribution
baseline_diagnosis <- fs_clinical_90day_linked %>%
  semi_join(complete_longitudinal, by = "RID") %>%
  group_by(RID) %>%
  slice_min(MRI_EXAMDATE, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  count(DIAGNOSIS, name = "n_participants")

baseline_diagnosis

dx_nearest %>%
  count(n_nearest, diagnosis_conflict, sort = TRUE) 
#CN                   150
#Dementia              68
#MCI                  269

diagnosis_discrepancies #0

# Assess MRI-to-clinical diagnosis coverage
diagnosis_coverage <- fs_qc_longitudinal %>%
  select(RID, VISCODE2, EXAMDATE) %>%
  distinct() %>%
  left_join(
    fs_clinical_90day_linked %>%
      select(RID, MRI_VISCODE2, MRI_EXAMDATE, DIAGNOSIS),
    by = c(
      "RID",
      "VISCODE2" = "MRI_VISCODE2",
      "EXAMDATE" = "MRI_EXAMDATE"
    )
  ) %>%
  summarise(
    MRI_visits = n(),
    linked_visits = sum(!is.na(DIAGNOSIS)),
    unmatched_visits = sum(is.na(DIAGNOSIS)),
    linkage_rate = 100 * mean(!is.na(DIAGNOSIS))
  )

diagnosis_coverage

#plots f
ggplot(
  tibble(
    status = c("Linked", "Unmatched"),
    n = c(
      diagnosis_coverage$linked_visits,
      diagnosis_coverage$unmatched_visits
    )
  ),
  aes(x = status, y = n)
) +
  geom_col(width = 0.65, fill = "grey35") +
  geom_text(
    aes(label = n),
    vjust = -0.4,
    size = 4
  ) +
  labs(
    x = "MRI–clinical linkage status",
    y = "Number of MRI visits"
  ) +
  theme_classic(base_size = 13) +
  theme(
    axis.title = element_text(size = 13),
    axis.text = element_text(size = 11)
  )

# Identify MRI visits without a clinical diagnosis within ±90 days
unmatched_mri <- fs_qc_longitudinal %>%
  select(RID, VISCODE2, EXAMDATE) %>%
  distinct() %>%
  anti_join(
    fs_clinical_90day_linked %>%
      select(RID, MRI_VISCODE2, MRI_EXAMDATE),
    by = c(
      "RID",
      "VISCODE2" = "MRI_VISCODE2",
      "EXAMDATE" = "MRI_EXAMDATE"
    )
  )

# Summarise unmatched visits by participant
unmatched_summary <- unmatched_mri %>%
  count(RID, name = "unmatched_visits") %>%
  arrange(desc(unmatched_visits)) 

unmatched_summary %>%
  summarise(
    participants_with_unmatched_visits = n(),
    unmatched_visits = sum(unmatched_visits),
    max_unmatched_per_subject = max(unmatched_visits)
  )

# Distribution of the number of unmatched visits per participant
ggplot(unmatched_summary, aes(x = unmatched_visits)) +
  geom_histogram(
    binwidth = 1,
    boundary = 0.5,
    fill = "grey35",
    color = "white"
  ) +
  geom_text(
    stat = "bin",
    aes(label = after_stat(count)),
    vjust = -0.4,
    size = 4
  ) +
  labs(
    x = "Number of unmatched MRI visits",
    y = "Number of participants"
  ) +
  theme_classic(base_size = 13) +
  theme(
    axis.title = element_text(size = 13),
    axis.text = element_text(size = 11)
  )























































###linking overall data
# Before dedup
cat("Before:", nrow(fs_qc_longitudinal), "rows\n")

fs_qc_longitudinal <- fs_qc_longitudinal %>%
  arrange(RID, VISCODE2, EXAMDATE) %>%
  distinct(RID, VISCODE2, .keep_all = TRUE)

# After dedup
cat("After:", nrow(fs_qc_longitudinal), "rows\n")
cat("Unique RID+VISCODE2:", n_distinct(fs_qc_longitudinal$RID, fs_qc_longitudinal$VISCODE2), "\n")
cat("Unique RIDs:", n_distinct(fs_qc_longitudinal$RID), "\n")

# Deduplicate fs_qc_longitudinal (RID 5040 at m03 has 2 scans)
cat("Before dedup:", nrow(fs_qc_longitudinal), "rows\n")

fs_qc_longitudinal <- fs_qc_longitudinal %>%
  arrange(RID, VISCODE2, EXAMDATE) %>%
  distinct(RID, VISCODE2, .keep_all = TRUE)

cat("After dedup:", nrow(fs_qc_longitudinal), "rows\n")
cat("Unique RID+VISCODE2:", n_distinct(fs_qc_longitudinal$RID, fs_qc_longitudinal$VISCODE2), "\n")
cat("Unique RIDs:", n_distinct(fs_qc_longitudinal$RID), "\n\n")

# Join demographics (subject-level → one-to-many is expected)
master_visit <- fs_qc_longitudinal %>%
  left_join(demog, by = "RID")

# Validate the demographics merge
cat("After demographics merge:\n")
cat("  Rows:", nrow(master_visit), "(should be", nrow(fs_qc_longitudinal), ")\n")
cat("  Unique RIDs:", n_distinct(master_visit$RID), "\n")
cat("  Row count changed:", nrow(master_visit) != nrow(fs_qc_longitudinal), "\n")
cat("  Participants with demog:", sum(!is.na(master_visit$PTGENDER)), "/", nrow(master_visit), "\n")
cat("  Participants missing demog:", sum(is.na(master_visit$PTGENDER)), "\n\n")

# Prepare MMSE: deduplicate RID 830, keep only the score we need
mmse <- fread(file.path(data_dir, "MMSE.csv"), showProgress = FALSE) %>%
  as_tibble() %>%
  arrange(RID, VISCODE2, VISCODE) %>%
  distinct(RID, VISCODE2, .keep_all = TRUE) %>%
  select(RID, VISCODE2, MMSCORE)

cat("MMSE prepared:\n")
cat("  Rows:", nrow(mmse), "\n")
cat("  Unique RID+VISCODE2:", n_distinct(mmse$RID, mmse$VISCODE2), "\n")
cat("  Score range:", range(mmse$MMSCORE, na.rm = TRUE), "\n")
cat("  Missing scores:", sum(is.na(mmse$MMSCORE)), "\n\n")

# Join MMSE to master_visit
master_visit <- master_visit %>%
  left_join(mmse, by = c("RID", "VISCODE2"))

cat("After MMSE merge:\n")
cat("  Rows:", nrow(master_visit), "(should still be", nrow(fs_qc_longitudinal), ")\n")
cat("  MRI visits with MMSE:", sum(!is.na(master_visit$MMSCORE)), "/", nrow(master_visit), "\n")
cat("  MRI visits without MMSE:", sum(is.na(master_visit$MMSCORE)), "\n")

# CDR — primary score is CDRSB (CDR Sum of Boxes) and CDGLOBAL (global CDR)
cdr <- fread(file.path(data_dir, "CDR.csv"), showProgress = FALSE) %>%
  as_tibble() %>%
  arrange(RID, VISCODE2, VISCODE) %>%
  distinct(RID, VISCODE2, .keep_all = TRUE) %>%
  select(RID, VISCODE2, CDRSB, CDGLOBAL)

# ADAS — primary score is typically TOTSCORE (ADAS-Cog total) or similar
adas_raw <- fread(file.path(data_dir, "ADAS.csv"), showProgress = FALSE) %>% as_tibble()
cat("ADAS score columns:\n")
print(names(adas_raw)[!names(adas_raw) %in% c("Phase", "ID", "SITEID", "ORIGPROT",
                                              "COLPROT", "VISCODE", "VISCODE2", "PTID", "RID", "update_stamp")])

# FAQ — primary score is FAQTOTAL
faq <- fread(file.path(data_dir, "FAQ.csv"), showProgress = FALSE) %>%
  as_tibble() %>%
  arrange(RID, VISCODE2, VISCODE) %>%
  distinct(RID, VISCODE2, .keep_all = TRUE) %>%
  select(RID, VISCODE2, FAQTOTAL)

# NEUROBAT — inspect score columns before selecting
neuro_raw <- fread(file.path(data_dir, "NEUROBAT.csv"), showProgress = FALSE) %>% as_tibble()
cat("\nNEUROBAT score columns:\n")
print(names(neuro_raw)[!names(neuro_raw) %in% c("Phase", "ID", "SITEID", "ORIGPROT",
                                                "COLPROT", "VISCODE", "VISCODE2", "PTID", "RID", "update_stamp",
                                                "USERDATE", "USERDATE2", "EXAMDATE", "VISDATE")])

# Columns to drop from every cognitive file (if they exist)
drop_cols <- c("Phase", "ID", "SITEID", "ORIGPROT", "COLPROT", "PTID",
               "VISCODE", "update_stamp", "EXAMDATE", "VISDATE",
               "SOURCE", "LMSTORY")

# Helper function: deduplicate and drop admin columns
clean_cog <- function(file) {
  fread(file.path(data_dir, file), showProgress = FALSE) %>%
    as_tibble() %>%
    arrange(RID, VISCODE2, VISCODE) %>%
    distinct(RID, VISCODE2, .keep_all = TRUE) %>%
    select(-any_of(drop_cols)) %>%
    select(-matches("USERDATE|DD_CRF|LANGUAGE_CODE|HAS_QC"))
}

mmse_clean    <- clean_cog("MMSE.csv")
cdr_clean     <- clean_cog("CDR.csv")
adas_clean    <- clean_cog("ADAS.csv")
faq_clean     <- clean_cog("FAQ.csv")
neurobat_clean <- clean_cog("NEUROBAT.csv")

# Check for column name conflicts across datasets
all_names <- c(names(mmse_clean), names(cdr_clean), names(adas_clean),
               names(faq_clean), names(neurobat_clean))
non_key <- all_names[!all_names %in% c("RID", "VISCODE2")]
dup_names <- non_key[duplicated(non_key)]

cat("Column name conflicts:\n")
if (length(dup_names) > 0) {
  print(unique(dup_names))
} else {
  cat("None — safe to merge\n")
}

cat("\nColumns per dataset:\n")
cat("  MMSE:", ncol(mmse_clean), "\n")
cat("  CDR:", ncol(cdr_clean), "\n")
cat("  ADAS:", ncol(adas_clean), "\n")
cat("  FAQ:", ncol(faq_clean), "\n")
cat("  NEUROBAT:", ncol(neurobat_clean), "\n")


# SPID is an admin ID, not a score — drop it from all datasets
mmse_clean    <- mmse_clean %>% select(-any_of("SPID"))
cdr_clean     <- cdr_clean %>% select(-any_of("SPID"))
adas_clean    <- adas_clean %>% select(-any_of("SPID"))
faq_clean     <- faq_clean %>% select(-any_of("SPID"))
neurobat_clean <- neurobat_clean %>% select(-any_of("SPID"))

# Prefix cognitive columns to avoid future ambiguity (except join keys)
# This makes it clear which test each column came from
mmse_clean <- mmse_clean %>% rename_with(~ paste0("MMSE_", .), -c(RID, VISCODE2))
cdr_clean  <- cdr_clean %>% rename_with(~ paste0("CDR_", .), -c(RID, VISCODE2))
adas_clean <- adas_clean %>% rename_with(~ paste0("ADAS_", .), -c(RID, VISCODE2))
faq_clean  <- faq_clean %>% rename_with(~ paste0("FAQ_", .), -c(RID, VISCODE2))
neurobat_clean <- neurobat_clean %>% rename_with(~ paste0("NB_", .), -c(RID, VISCODE2))

# Rebuild master_visit from clean anchor
master_visit <- fs_qc_longitudinal %>%
  left_join(demog, by = "RID") %>%
  left_join(mmse_clean, by = c("RID", "VISCODE2")) %>%
  left_join(cdr_clean, by = c("RID", "VISCODE2")) %>%
  left_join(adas_clean, by = c("RID", "VISCODE2")) %>%
  left_join(faq_clean, by = c("RID", "VISCODE2")) %>%
  left_join(neurobat_clean, by = c("RID", "VISCODE2"))

# Validate
cat("master_visit dimensions:", nrow(master_visit), "rows x", ncol(master_visit), "cols\n")
cat("Expected rows: 2954\n")
cat("Row count correct:", nrow(master_visit) == 2954, "\n")
cat("Unique RIDs:", n_distinct(master_visit$RID), "\n")
cat("Unique RID+VISCODE2:", n_distinct(master_visit$RID, master_visit$VISCODE2), "\n\n")

# Coverage of key summary scores across MRI visits
key_scores <- c("MMSE_MMSCORE", "CDR_CDRSB", "CDR_CDGLOBAL",
                "ADAS_TOTSCORE", "ADAS_TOTAL13", "FAQ_FAQTOTAL",
                "NB_LDELTOTAL", "NB_AVDEL30MIN", "NB_TRAASCOR",
                "NB_TRABSCOR", "NB_CATANIMSC", "NB_DIGITSCOR")

coverage <- tibble(
  variable    = key_scores,
  n_present   = map_int(key_scores, ~ sum(!is.na(master_visit[[.x]]))),
  pct_present = round(map_dbl(key_scores, ~ mean(!is.na(master_visit[[.x]]))) * 100, 1)
)

print(coverage, n = Inf)


#Merge APOE data
apoe <- apoe_raw %>%
  select(RID, GENOTYPE) %>%
  separate(GENOTYPE, into = c("APOE_A1", "APOE_A2"), sep = "/", convert = TRUE) %>%
  mutate(
    APOE4_count   = (APOE_A1 == 4) + (APOE_A2 == 4),
    APOE4_carrier = APOE4_count > 0
  )

cat("APOE: Rows:", nrow(apoe), "| Unique RIDs:", n_distinct(apoe$RID), "\n")
apoe %>% count(APOE4_count) %>% print()

# Merge into master_visit
pre_rows <- nrow(master_visit)
master_visit <- master_visit %>% left_join(apoe, by = "RID")

cat("\nAfter APOE merge:", nrow(master_visit), "rows (expected:", pre_rows, ")\n")
cat("Participants with APOE:", sum(!is.na(master_visit$APOE_A1)), "/", nrow(master_visit), "\n\n")

# UPENNBIOMK — use MEDIAN batch, check VISCODE values
biomarker <- fread(file.path(data_dir, "UPENNBIOMK_MASTER.csv"), showProgress = FALSE) %>%
  as_tibble() %>%
  filter(BATCH == "MEDIAN") %>%
  select(RID, VISCODE, ABETA, TAU, PTAU) %>%
  rename(CSF_ABETA = ABETA, CSF_TAU = TAU, CSF_PTAU = PTAU)

cat("Biomarker (MEDIAN batch):\n")
cat("  Rows:", nrow(biomarker), "| Unique RIDs:", n_distinct(biomarker$RID), "\n")
cat("  VISCODE values:\n")
biomarker %>% count(VISCODE, sort = TRUE) %>% print(n = 20)

# Join biomarker using VISCODE → VISCODE2 mapping
pre_rows <- nrow(master_visit)

master_visit <- master_visit %>%
  left_join(biomarker, by = c("RID", "VISCODE2" = "VISCODE"))

cat("After biomarker merge:", nrow(master_visit), "rows (expected:", pre_rows, ")\n")
cat("MRI visits with CSF ABETA:", sum(!is.na(master_visit$CSF_ABETA)), "/", nrow(master_visit), "\n")
cat("Unique participants with CSF:", n_distinct(master_visit$RID[!is.na(master_visit$CSF_ABETA)]), "/ 602\n\n")

# Final master_visit summary
cat("FINAL master_visit:\n")
cat("  Rows:", nrow(master_visit), "\n")
cat("  Columns:", ncol(master_visit), "\n")
cat("  Unique RIDs:", n_distinct(master_visit$RID), "\n")
cat("  Unique RID+VISCODE2:", n_distinct(master_visit$RID, master_visit$VISCODE2), "\n\n")

# Check no row multiplication occurred at any step
cat("Row integrity check:", nrow(master_visit) == 2954, "\n")


# Check what columns exist in fs_clinical_90day_linked
names(fs_clinical_90day_linked)

# Select all columns, rename MRI_VISCODE2 to match the join key
dx_link <- fs_clinical_90day_linked %>%
  rename(VISCODE2 = MRI_VISCODE2) %>%
  rename_with(~ paste0("DX_", .), -c(RID, VISCODE2)) %>%
  # DX_DIAGNOSIS is redundant naming — simplify
  rename(DIAGNOSIS = DX_DIAGNOSIS)

pre_rows <- nrow(master_visit)

master_visit <- master_visit %>%
  left_join(dx_link, by = c("RID", "VISCODE2"))

cat("After diagnosis merge:", nrow(master_visit), "rows (expected:", pre_rows, ")\n")
cat("Row integrity:", nrow(master_visit) == 2954, "\n")
cat("MRI visits with diagnosis:", sum(!is.na(master_visit$DIAGNOSIS)), "/", nrow(master_visit), "\n")
cat("MRI visits without diagnosis:", sum(is.na(master_visit$DIAGNOSIS)), "\n\n")

# Diagnosis distribution
master_visit %>%
  count(DIAGNOSIS, sort = TRUE) %>%
  mutate(pct = round(n / sum(n) * 100, 1)) %>%
  print() #2954 rows & 520 columns

cat("\nFinal master_visit:", nrow(master_visit), "rows x", ncol(master_visit), "cols\n")
cat("Unique RIDs:", n_distinct(master_visit$RID), "\n")

missing_dx_context <- missing_dx %>%
  left_join(
    master_visit %>%
      filter(!is.na(DIAGNOSIS)) %>%
      select(RID, VISCODE2, EXAMDATE, DIAGNOSIS),
    by = "RID",
    suffix = c("_miss", "_dx"),
    relationship = "many-to-many"
  ) %>%
  mutate(day_gap = as.numeric(EXAMDATE_dx - EXAMDATE_miss)) %>%
  group_by(RID, VISCODE2_miss) %>%
  summarise(
    dx_before = {
      neg <- which(day_gap < 0)
      if (length(neg) > 0) DIAGNOSIS[neg[which.max(day_gap[neg])]] else NA_character_
    },
    dx_after = {
      pos <- which(day_gap > 0)
      if (length(pos) > 0) DIAGNOSIS[pos[which.min(day_gap[pos])]] else NA_character_
    },
    .groups = "drop"
  ) %>%
  mutate(
    flanking_match = case_when(
      is.na(dx_before) | is.na(dx_after) ~ "one_side_missing",
      dx_before == dx_after              ~ "same",
      TRUE                               ~ "different"
    )
  )

missing_dx_context %>% count(flanking_match, sort = TRUE)

# Build the imputation lookup
impute_safe <- missing_dx_context %>%
  mutate(
    DIAGNOSIS_IMPUTED = case_when(
      # Same flanking: safe
      flanking_match == "same" ~ dx_before,
      # No-after but Dementia: safe (doesn't revert)
      flanking_match == "one_side_missing" & !is.na(dx_before) & dx_before == "Dementia" ~ "Dementia",
      # No-before: carry backward from next dx
      flanking_match == "one_side_missing" & is.na(dx_before) & !is.na(dx_after) ~ dx_after,
      # Everything else: leave as NA
      TRUE ~ NA_character_
    ),
    DIAGNOSIS_SOURCE = case_when(
      !is.na(DIAGNOSIS_IMPUTED) & flanking_match == "same" ~ "imputed_flanking",
      !is.na(DIAGNOSIS_IMPUTED) & is.na(dx_after) ~ "imputed_carry_forward",
      !is.na(DIAGNOSIS_IMPUTED) & is.na(dx_before) ~ "imputed_carry_backward",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(DIAGNOSIS_IMPUTED)) %>%
  select(RID, VISCODE2_miss, DIAGNOSIS_IMPUTED, DIAGNOSIS_SOURCE)

cat("Visits imputed:", nrow(impute_safe), "\n")
impute_safe %>% count(DIAGNOSIS_SOURCE, DIAGNOSIS_IMPUTED) %>% print()

# Apply imputation to master_visit while preserving provenance
master_visit <- master_visit %>%
  left_join(impute_safe, by = c("RID", "VISCODE2" = "VISCODE2_miss")) %>%
  mutate(
    DIAGNOSIS_FINAL  = coalesce(DIAGNOSIS, DIAGNOSIS_IMPUTED),
    DIAGNOSIS_SOURCE = case_when(
      !is.na(DIAGNOSIS)         ~ "observed",
      !is.na(DIAGNOSIS_IMPUTED) ~ DIAGNOSIS_SOURCE,
      TRUE                      ~ "missing"
    )
  )

# Validate
cat("\nDiagnosis coverage after imputation:\n")
master_visit %>% count(DIAGNOSIS_SOURCE, sort = TRUE) %>% print()

cat("\nDiagnosis distribution (DIAGNOSIS_FINAL):\n")
master_visit %>%
  count(DIAGNOSIS_FINAL, sort = TRUE) %>%
  mutate(pct = round(n / sum(n) * 100, 1)) %>%
  print()

cat("\nRemaining missing:", sum(is.na(master_visit$DIAGNOSIS_FINAL)), "visits\n")

#############Data visualization plots###############
library(ggplot2)
library(patchwork)

baseline_dx <- master_visit %>%
  filter(!is.na(EXAMDATE)) %>%
  arrange(RID, EXAMDATE) %>%
  distinct(RID, .keep_all = TRUE)

cat("Baseline participants:", n_distinct(baseline_dx$RID), "\n")
baseline_dx %>%
  count(DIAGNOSIS_FINAL, sort = TRUE) %>%
  mutate(pct = round(n / sum(n) * 100, 1)) %>%
  print()

cat("Last known diagnosis:\n")
last_known_dx %>% count(LAST_DX, sort = TRUE) %>% print()

theme_pub <- theme(
  plot.title = element_text(face = "bold", size = 14),
  plot.subtitle = element_text(color = "grey40", size = 11),
  axis.title = element_text(face = "bold", size = 12),
  axis.text = element_text(size = 11, color = "black"),
  axis.line = element_line(color = "black", linewidth = 0.5),
  axis.ticks = element_line(color = "black", linewidth = 0.4),
  panel.background = element_blank(),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  plot.background = element_blank(),
  legend.position = "bottom",
  plot.margin = margin(10, 15, 10, 15)
)

#Baseline diagnosis distribution
p1 <- baseline_dx %>%
  filter(!is.na(DIAGNOSIS_FINAL)) %>%
  mutate(dx = factor(DIAGNOSIS_FINAL, levels = c("CN", "MCI", "Dementia"))) %>%
  count(dx) %>%
  ggplot(aes(x = dx, y = n, fill = dx)) +
  geom_col(width = 0.6, show.legend = FALSE) +
  geom_text(aes(label = n), vjust = -0.5, fontface = "bold", size = 4.5) +
  scale_fill_manual(values = c("CN" = "#4CAF50", "MCI" = "#FF9800", "Dementia" = "#F44336")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(title = "Baseline Diagnosis Distribution",
       subtitle = paste0("N = ", sum(!is.na(baseline_dx$DIAGNOSIS_FINAL)), " participants"),
       x = NULL, y = "Number of Participants") +
  theme_pub

ggsave("baseline_diagnosis.png", p1, width = 6, height = 5, dpi = 300, bg = "white")
p1

# Last visit diagnosis distribution
last_known_dx <- master_visit %>%
  filter(!is.na(DIAGNOSIS_FINAL), !is.na(EXAMDATE)) %>%
  arrange(RID, desc(EXAMDATE)) %>%
  distinct(RID, .keep_all = TRUE) %>%
  select(RID, LAST_DX = DIAGNOSIS_FINAL)

p2 <- last_known_dx %>%
  mutate(dx = factor(LAST_DX, levels = c("CN", "MCI", "Dementia"))) %>%
  count(dx) %>%
  ggplot(aes(x = dx, y = n, fill = dx)) +
  geom_col(width = 0.6, show.legend = FALSE) +
  geom_text(aes(label = n), vjust = -0.5, fontface = "bold", size = 4.5) +
  scale_fill_manual(values = c("CN" = "#4CAF50", "MCI" = "#FF9800", "Dementia" = "#F44336")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(title = "Last Known Diagnosis Distribution",
       subtitle = paste0("N = ", n_distinct(last_known_dx$RID),
                         " participants (last available clinical diagnosis)"),
       x = NULL, y = "Number of Participants") +
  theme_pub

ggsave("last_known_diagnosis.png", p2, width = 6, height = 5, dpi = 300, bg = "white")
p2


#Number of MRI visits per participant
visits_per_pt <- master_visit %>% count(RID, name = "n_visits")
visits_summary <- visits_per_pt %>% count(n_visits, name = "count")

p3 <- visits_summary %>%
  ggplot(aes(x = n_visits, y = count)) +
  geom_col(fill = "#2196F3", width = 0.7) +
  geom_text(aes(label = count), vjust = -0.5, fontface = "bold", size = 3.8) +
  scale_x_continuous(breaks = seq(min(visits_summary$n_visits),
                                  max(visits_summary$n_visits), 1)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(title = "MRI Visits per Participant",
       subtitle = paste0("Median: ", median(visits_per_pt$n_visits),
                         " (Range: ", min(visits_per_pt$n_visits),
                         "–", max(visits_per_pt$n_visits), ")"),
       x = "Number of MRI Visits", y = "Number of Participants") +
  theme_pub

ggsave("mri_visits_per_participant.png", p3, width = 6, height = 5, dpi = 300, bg = "white")
p3

#Follow-up duration
followup <- master_visit %>%
  group_by(RID) %>%
  summarise(fu_months = as.numeric(max(EXAMDATE) - min(EXAMDATE)) / 30.44,
            .groups = "drop")

max_fu <- ceiling(max(followup$fu_months) / 6) * 6

p4 <- followup %>%
  ggplot(aes(x = fu_months)) +
  geom_histogram(binwidth = 6, fill = "#9C27B0", color = "white", alpha = 0.85) +
  stat_bin(binwidth = 6, geom = "text", aes(label = after_stat(count)),
           vjust = -0.5, fontface = "bold", size = 3.8) +
  scale_x_continuous(breaks = seq(0, max_fu, 6)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(title = "Follow-up Duration",
       subtitle = paste0("Median: ", round(median(followup$fu_months), 1),
                         " months (Range: ", round(min(followup$fu_months), 1),
                         "–", round(max(followup$fu_months), 1), ")"),
       x = "Follow-up (Months)", y = "Number of Participants") +
  theme_pub

ggsave("followup_duration.png", p4, width = 7, height = 5, dpi = 300, bg = "white")
p4


###data missingness plots#######
# Select only the essential summary variables
key_vars <- c(
  "DIAGNOSIS_FINAL",
  "MMSE_MMSCORE",
  "CDR_CDRSB", "CDR_CDGLOBAL",
  "ADAS_TOTSCORE", "ADAS_TOTAL13",
  "FAQ_FAQTOTAL",
  "NB_LDELTOTAL",
  "NB_AVDEL30MIN",
  "NB_TRABSCOR",
  "CSF_ABETA", "CSF_TAU",
  "PTGENDER", "PTEDUCAT", "APOE_A1"
)

# Clean readable labels
var_labels <- c(
  "Clinical Diagnosis",
  "MMSE Total Score",
  "CDR Sum of Boxes", "CDR Global Score",
  "ADAS-Cog 11", "ADAS-Cog 13",
  "FAQ Total",
  "Logical Memory (Delayed Recall)",
  "RAVLT (30-min Delayed Recall)",
  "Trail Making Test B",
  "CSF Amyloid-\u03B242", "CSF Total Tau",
  "Sex", "Education (Years)", "APOE Genotype"
)

viscode_order <- c("scmri", "sc", "bl", "m03", "m06", "m12",
                   "m24", "m36", "m48", "m60", "m72", "m84")

# Clean visit labels for x-axis
visit_labels <- c("Screening\nMRI", "Screening", "Baseline", "Month 3",
                  "Month 6", "Month 12", "Month 24", "Month 36",
                  "Month 48", "Month 60", "Month 72", "Month 84")

missingness <- master_visit %>%
  filter(VISCODE2 %in% viscode_order) %>%
  mutate(VISCODE2 = factor(VISCODE2, levels = viscode_order)) %>%
  group_by(VISCODE2) %>%
  summarise(
    n_visits = n(),
    across(all_of(key_vars), ~ round(mean(is.na(.)) * 100, 1)),
    .groups = "drop"
  )

miss_long <- missingness %>%
  select(-n_visits) %>%
  pivot_longer(-VISCODE2, names_to = "variable", values_to = "pct_missing") %>%
  mutate(
    variable = factor(variable, levels = rev(key_vars), labels = rev(var_labels))
  )

p_miss <- miss_long %>%
  ggplot(aes(x = VISCODE2, y = variable, fill = pct_missing)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = round(pct_missing)), size = 3, color = "black") +
  scale_x_discrete(labels = visit_labels) +
  scale_fill_gradient2(low = "#4CAF50", mid = "#FFF9C4", high = "#F44336",
                       midpoint = 50, limits = c(0, 100),
                       name = "% Missing") +
  labs(title = "Missingness Pattern Across Visits and Variables",
       subtitle = paste0("N = ", n_distinct(master_visit$RID),
                         " participants | Values show % missing per visit"),
       x = NULL, y = NULL) +
  theme_pub +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
    axis.text.y = element_text(size = 10),
    legend.position = "right"
  )

ggsave("missingness_heatmap.png", p_miss, width = 12, height = 6, dpi = 300, bg = "white")
p_miss



































































