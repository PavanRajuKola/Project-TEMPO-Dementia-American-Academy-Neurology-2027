# Install the R packages used by the TEMPO analysis scripts
pkgs <- c("tidyverse", "data.table", "dplyr", "purrr", "ggplot2", "patchwork",
          "lme4", "lmerTest", "corrplot", "ggrepel", "ggseg", "here")
install.packages(setdiff(pkgs, rownames(installed.packages())))
