source(here::here("config.R"))


# ==============================================================
# STEP: Deep Structural Audit — Cognitive & Supporting Datasets
# Purpose: Determine unit of observation for each file
# Risk:   Unknown granularity → unsafe joins → silent row duplication
# ==============================================================

library(data.table)
library(dplyr)
library(purrr)

data_dir <- ADNIMERGE_DIR

# Files to audit (excluding FreeSurfer files we've already examined + excluded)
audit_files <- c(
  # Cognitive
  "MMSE.csv", "CDR.csv", "ADAS.csv", "MOCA.csv",
  "FAQ.csv", "GDSCALE.csv", "NEUROBAT.csv",
  # Supporting
  "PTDEMOG.csv", "APOERES.csv", "UPENNBIOMK_MASTER.csv",
  "ADSL.csv", "REGISTRY.csv", "MEDHIST.csv",
  # MRI metadata / QC
  "MRI3META.csv", "MRIMETA.csv", "MRIQC.csv",
  # Diagnosis (already used, but include for completeness)
  "DXSUM.csv"
)

granularity_audit <- map_dfr(audit_files, function(f) {
  
  x <- fread(file.path(data_dir, f), showProgress = FALSE)
  
  has_rid <- "RID" %in% names(x)
  has_vc  <- "VISCODE" %in% names(x)
  has_vc2 <- "VISCODE2" %in% names(x)
  has_ed  <- "EXAMDATE" %in% names(x)
  
  n_rows <- nrow(x)
  n_rid  <- if (has_rid) n_distinct(x$RID) else NA_integer_
  
  # Determine candidate key and check uniqueness
  if (has_rid & has_vc2) {
    key_type     <- "RID + VISCODE2"
    n_combos     <- n_distinct(x[, .(RID, VISCODE2)])
    dup_tbl      <- x[, .N, by = .(RID, VISCODE2)][N > 1]
    n_dup_combos <- nrow(dup_tbl)
    max_per_key  <- if (n_dup_combos > 0) max(dup_tbl$N) else 1L
  } else if (has_rid & has_vc & !has_vc2) {
    key_type     <- "RID + VISCODE"
    n_combos     <- n_distinct(x[, .(RID, VISCODE)])
    dup_tbl      <- x[, .N, by = .(RID, VISCODE)][N > 1]
    n_dup_combos <- nrow(dup_tbl)
    max_per_key  <- if (n_dup_combos > 0) max(dup_tbl$N) else 1L
  } else if (has_rid & !has_vc & !has_vc2) {
    key_type     <- "RID only"
    n_combos     <- n_rid
    dup_tbl      <- x[, .N, by = .(RID)][N > 1]
    n_dup_combos <- nrow(dup_tbl)
    max_per_key  <- if (n_dup_combos > 0) max(dup_tbl$N) else 1L
  } else {
    key_type     <- "NO RID"
    n_combos     <- NA_integer_
    n_dup_combos <- NA_integer_
    max_per_key  <- NA_integer_
  }
  
  tibble(
    file             = f,
    n_rows           = n_rows,
    n_participants   = n_rid,
    key_tested       = key_type,
    n_unique_keys    = n_combos,
    n_duplicate_keys = n_dup_combos,
    max_rows_per_key = max_per_key,
    key_is_unique    = if (!is.na(n_dup_combos)) n_dup_combos == 0 else NA,
    has_EXAMDATE     = has_ed,
    n_examdate       = if (has_ed) sum(!is.na(x$EXAMDATE)) else NA_integer_
  )
})

# Print full audit table
print(granularity_audit, n = Inf, width = Inf)

library(data.table)

data_dir <- ADNIMERGE_DIR

cog_files <- c("MMSE.csv", "CDR.csv", "ADAS.csv",
               "FAQ.csv", "NEUROBAT.csv")

for (f in cog_files) {
  x <- fread(file.path(data_dir, f), showProgress = FALSE)
  dups <- x[, .N, by = .(RID, VISCODE2)][N > 1]
  
  cat("\n", f, "— Duplicated RID+VISCODE2:", nrow(dups), "\n")
  
  if (nrow(dups) > 0) {
    dup_rows <- x[dups, on = .(RID, VISCODE2)]
    id_cols <- intersect(
      c("RID", "VISCODE", "VISCODE2",
        "ORIGPROT", "COLPROT", "Phase", "update_stamp"),
      names(dup_rows)
    )
    cat("Key columns:\n")
    print(as.data.frame(dup_rows[, ..id_cols]))
  }
}

mmse <- fread(file.path(data_dir, "MMSE.csv"), showProgress = FALSE) %>% as_tibble()

rid830 <- mmse %>% filter(RID == 830, VISCODE2 == "m60")

cat("VISCODE:", rid830$VISCODE, "\n")
cat("COLPROT:", rid830$COLPROT, "\n\n")

meta_cols <- c("Phase", "ID", "SITEID", "ORIGPROT", "COLPROT",
               "VISCODE", "PTID", "update_stamp")

score_cols <- setdiff(names(rid830), meta_cols)

comparison <- tibble(
  column = score_cols,
  row1   = unlist(rid830[1, score_cols]),
  row2   = unlist(rid830[2, score_cols])
) %>%
  mutate(match = row1 == row2 | (is.na(row1) & is.na(row2)))

cat("Columns that DIFFER:\n")
comparison %>% filter(!match) %>% print(n = Inf)

cat("\nTotal differing:", sum(!comparison$match), "of", nrow(comparison), "\n")

####
cog_files <- c("MMSE.csv", "CDR.csv", "ADAS.csv", "FAQ.csv", "NEUROBAT.csv")

meta_cols <- c("Phase", "ID", "SITEID", "ORIGPROT", "COLPROT",
               "VISCODE", "PTID", "update_stamp")

walk(cog_files, function(f) {
  
  x <- fread(file.path(data_dir, f), showProgress = FALSE) %>% as_tibble()
  
  dup <- x %>% filter(RID == 830, VISCODE2 == "m60")
  
  if (nrow(dup) != 2) return(invisible(NULL))
  
  score_cols <- setdiff(names(dup), meta_cols)
  
  diffs <- map_chr(score_cols, function(col) {
    v1 <- as.character(dup[[col]][1])
    v2 <- as.character(dup[[col]][2])
    if (identical(v1, v2)) "SAME" else paste(v1, "vs", v2)
  }) %>%
    set_names(score_cols)
  
  diff_cols <- diffs[diffs != "SAME"]
  
  cat("\n", f, "\n")
  cat("  VISCODE:", dup$VISCODE, " | COLPROT:", dup$COLPROT, "\n")
  cat("  Columns differing:", length(diff_cols), "of", length(score_cols), "\n")
  
  if (length(diff_cols) > 0) {
    walk2(names(diff_cols), diff_cols, ~ cat("   ", .x, ":", .y, "\n"))
  }
})

library(tidyverse)
library(data.table)

data_dir <- ADNIMERGE_DIR

dxsum <- fread(file.path(data_dir, "DXSUM.csv"), showProgress = FALSE) %>% as_tibble()

dups <- dxsum %>%
  count(RID, VISCODE2) %>%
  filter(n > 1)

cat("Duplicate RID+VISCODE2 in DXSUM:", nrow(dups), "\n\n")
print(dups)

# Show the duplicate rows
dxsum %>%
  semi_join(dups, by = c("RID", "VISCODE2")) %>%
  select(RID, VISCODE, VISCODE2, EXAMDATE, ORIGPROT, COLPROT,
         DIAGNOSIS, update_stamp) %>%
  arrange(RID, VISCODE2) %>%
  print(n = Inf)

cat("RID 830 in fs_qc_longitudinal:",
    830 %in% fs_qc_longitudinal$RID, "\n")

cat("RID 6785 in fs_qc_longitudinal:",
    6785 %in% fs_qc_longitudinal$RID, "\n")

ptdemog <- fread(file.path(data_dir, "PTDEMOG.csv"), showProgress = FALSE) %>% as_tibble()

# How many rows per RID+VISCODE2 combination?
ptdemog %>%
  count(RID, VISCODE2) %>%
  count(n, name = "n_combos")

ptdemog_dups <- ptdemog %>%
  count(RID, VISCODE2) %>%
  filter(n > 1)

# Show the COLPROT and VISCODE patterns in the duplicated rows
ptdemog %>%
  semi_join(ptdemog_dups, by = c("RID", "VISCODE2")) %>%
  count(VISCODE, VISCODE2, ORIGPROT, COLPROT, sort = TRUE) %>%
  print(n = 20)

names(ptdemog)

demog <- ptdemog %>%
  # Prefer screening visit, then baseline
  arrange(RID, match(VISCODE2, c("sc", "bl", "m06", "m12"))) %>%
  distinct(RID, .keep_all = TRUE) %>%
  select(RID, PTGENDER, PTDOBYY, PTEDUCAT, PTETHCAT, PTRACCAT,
         PTMARRY, PTHAND)

# Verify: exactly one row per participant
cat("Rows:", nrow(demog), "\n")
cat("Unique RIDs:", n_distinct(demog$RID), "\n\n")

# Quick look at the values
demog %>% count(PTGENDER)
demog %>% count(PTRACCAT, sort = TRUE)
demog %>% count(PTETHCAT, sort = TRUE)
summary(demog$PTEDUCAT)
summary(demog$PTDOBYY)

cat("fs_qc_longitudinal:\n")
cat("  Rows:", nrow(fs_qc_longitudinal), "\n")
cat("  Unique RIDs:", n_distinct(fs_qc_longitudinal$RID), "\n")
cat("  Unique RID+VISCODE2:", n_distinct(fs_qc_longitudinal$RID, fs_qc_longitudinal$VISCODE2), "\n")



