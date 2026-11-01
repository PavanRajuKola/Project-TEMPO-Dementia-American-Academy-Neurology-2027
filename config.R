# ==============================================================
# TEMPO — project configuration
# Edit these paths to point at your local copy of the ADNI data.
# ADNI data are NOT included in this repository (see README).
# ==============================================================

# ADNIMERGE CSV exports (UCSFFSL51.csv, DXSUM.csv, ADAS.csv, FAQ.csv, ...)
ADNIMERGE_DIR  <- Sys.getenv("ADNIMERGE_DIR",  "data/ADNIMERGE_Output")

# ADNIMERGE2 R package data folder (.rda files: UWNPSYCHSUM, FAQ, ...)
ADNIMERGE2_DIR <- Sys.getenv("ADNIMERGE2_DIR", "data/ADNIMERGE2/data")

# Output folders
TABLE_DIR  <- here::here("tables")
FIGURE_DIR <- here::here("figures")
dir.create(TABLE_DIR,  showWarnings = FALSE, recursive = TRUE)
dir.create(FIGURE_DIR, showWarnings = FALSE, recursive = TRUE)
