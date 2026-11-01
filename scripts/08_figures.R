source(here::here("config.R"))

#####plots for MRI feature selection
wscore_features_12m <- c("InfLatVentricle", "Hippocampus")
wscore_features_18m <- c("InfLatVentricle", "Hippocampus", "TemporoOccipital", "TotalCortVol")

status_cols <- c("W-score feature" = "#B2182B", "Reliable (prediction set)" = "#2166AC", "Excluded" = "grey60")
dx_cols     <- c(CN = "#4DAF4A", MCI = "#E41A1C")

#reliability at 12 vs 18 months for every feature, with a dashed line at the 0.20 cut-off
p_dumbbell <- feature_table_comp %>%
  mutate(Region = fct_reorder(recode(region, !!!region_labels), reliability_18m)) %>%
  ggplot() +
  geom_segment(aes(x = reliability_12m, xend = reliability_18m, y = Region, yend = Region), colour = "grey60") +
  geom_point(aes(reliability_12m, Region, colour = "12 months"), size = 3) +
  geom_point(aes(reliability_18m, Region, colour = "18 months"), size = 3) +
  geom_vline(xintercept = 0.20, linetype = "dashed") +
  scale_colour_manual(values = c("12 months" = "#92C5DE", "18 months" = "#2166AC"), name = "Window") +
  labs(title = "Slope reliability by follow-up window", x = "Slope reliability", y = NULL) +
  theme_pub
p_dumbbell
ggsave("Fig1_reliability_dumbbell.png", p_dumbbell, width = 7, height = 6.5, dpi = 300)

#annual % change in CN vs MCI for each window's reliable features
plot_apc <- function(tbl, win) {
  tbl %>% filter(.data[[paste0("keep_", win)]]) %>%
    transmute(Region = recode(region, !!!region_labels),
              CN = .data[[paste0("APC_CN_", win)]], MCI = .data[[paste0("APC_MCI_", win)]]) %>%
    pivot_longer(c(CN, MCI), names_to = "Baseline diagnosis", values_to = "APC") %>%
    mutate(Region = fct_reorder(Region, APC)) %>%
    ggplot(aes(APC, Region, fill = `Baseline diagnosis`)) +
    geom_col(position = position_dodge(0.7), width = 0.65) +
    geom_vline(xintercept = 0) +
    scale_fill_manual(values = dx_cols) +
    labs(title = paste0("Annual % change by baseline diagnosis: ", sub("m", "-month", win), " window"),
         x = "Annual percentage change (%/year)", y = NULL) +
    theme_pub
}
p_apc_12m <- plot_apc(feature_table_comp, "12m"); p_apc_18m <- plot_apc(feature_table_comp, "18m")
p_apc_12m; p_apc_18m
ggsave("Fig2a_APC_12m.png", p_apc_12m, width = 7, height = 4,   dpi = 300)
ggsave("Fig2b_APC_18m.png", p_apc_18m, width = 7, height = 7.5, dpi = 300)


##reliability vs rate of change for each window, with the W-score features highlighted
library(ggrepel)   # install.packages("ggrepel") if needed

plot_rel_apc <- function(tbl, win, wfeat) {
  tbl %>%
    mutate(Region = recode(region, !!!region_labels),
           Status = case_when(region %in% wfeat & .data[[paste0("keep_", win)]] ~ "W-score feature",
                              .data[[paste0("keep_", win)]] ~ "Reliable (prediction set)",
                              TRUE ~ "Excluded"),
           apc = abs(.data[[paste0("APC_MCI_", win)]]),
           rel = .data[[paste0("reliability_", win)]]) %>%
    ggplot(aes(apc, rel, colour = Status)) +
    geom_hline(yintercept = 0.20, linetype = "dashed") +
    geom_point(size = 3) +
    geom_text_repel(aes(label = Region), size = 3, show.legend = FALSE, max.overlaps = 30) +
    scale_x_log10() +
    scale_colour_manual(values = status_cols) +
    labs(title = paste0("Reliability vs rate of change: ", sub("m", "-month", win), " window"),
         x = "|Annual % change| in MCI (log scale)", y = "Slope reliability") +
    theme_pub
}
p_sc_12m <- plot_rel_apc(feature_table_comp, "12m", wscore_features_12m)
p_sc_18m <- plot_rel_apc(feature_table_comp, "18m", wscore_features_18m)
p_sc_12m; p_sc_18m
ggsave("Fig3a_reliability_vs_APC_12m.png", p_sc_12m, width = 8, height = 6, dpi = 300)
ggsave("Fig3b_reliability_vs_APC_18m.png", p_sc_18m, width = 8, height = 6, dpi = 300)

##nstall ggseg and list the exact cortical and subcortical region names
install.packages("ggseg")
library(ggseg)
sort(unique(dk$data$region))
sort(unique(aseg$data$region))

dk_atlas   <- dk()
aseg_atlas <- aseg()

# Region names (works whether the atlas stores them in $data or at the top level)
# Search every part of the atlas for a column called "region" or "label"

# Search every part of the atlas for a column called "region" or "label"
find_regions <- function(atlas) {
  parts <- if (is.list(atlas)) atlas else list(atlas)
  hits <- map(parts, function(p) {
    df <- tryCatch(as.data.frame(p), error = function(e) NULL)
    if (is.null(df)) return(NULL)
    col <- intersect(c("region", "label"), names(df))
    if (length(col)) sort(unique(na.omit(df[[col[1]]]))) else NULL
  })
  compact(hits)
}

find_regions(dk_atlas)
find_regions(aseg_atlas)

#Corrected brain maps: reliability maps for each window and the 18-month W-score feature map

cort_map <- c(MiddleTemporal = "middle temporal", Fusiform = "fusiform",
              LateralOccipital = "lateral occipital", Precuneus = "precuneus")
sub_map  <- list(Hippocampus = "Hippocampus", Amygdala = "Amygdala", Caudate = "Caudate",
                 Thalamus = c("Thalamus", "Thalamus Proper"))   # ventricles not in aseg

rel_fill  <- scale_fill_gradient(low = "#FDD49E", high = "#B30000", limits = c(0.2, 0.7),
                                 na.value = "grey90", name = "Reliability")
map_theme <- theme_void() + theme(plot.title = element_text(face = "bold", hjust = 0.5))

# Turn feature names into atlas region names (one feature can map to several rows)
to_atlas <- function(df, map) {
  df %>% filter(region %in% names(map)) %>%
    mutate(region = map[region]) %>% unnest(region)
}

brain_maps <- function(win) {
  kept <- feature_table_comp %>% filter(.data[[paste0("keep_", win)]]) %>%
    select(region, reliability = all_of(paste0("reliability_", win)))
  lab <- sub("m", "-month", win)
  cort_df <- to_atlas(kept, as.list(cort_map))
  sub_df  <- to_atlas(kept, sub_map)
  list(
    cortical = if (nrow(cort_df) > 0)
      ggseg(cort_df, atlas = dk_atlas, mapping = aes(fill = reliability), colour = "grey40") +
      rel_fill + labs(title = paste("Reliable cortical features:", lab, "window")) + map_theme,
    subcortical =
      ggseg(sub_df, atlas = aseg_atlas, mapping = aes(fill = reliability), colour = "grey40") +
      rel_fill + labs(title = paste("Reliable subcortical features:", lab, "window")) + map_theme
  )
}

# One brain map, in ggseg 2 syntax; empty data gives NULL (e.g. no reliable cortex at 12m)

brain_plot <- function(df, atlas, fill_var, fill_scale, title) {
  if (nrow(df) == 0) return(NULL)                       # e.g. no reliable cortex at 12 months
  ggplot() +
    geom_brain(data = df, atlas = atlas,
               mapping = aes(fill = !!sym(fill_var)), colour = "grey40") +
    fill_scale + labs(title = title) + map_theme
}

maps_12m <- brain_maps("12m"); maps_18m <- brain_maps("18m")

# Check the data reaching the plots has the expected columns
print(to_atlas(feature_table_comp %>% filter(keep_18m) %>%
                 select(region, reliability = reliability_18m), sub_map))

maps_12m$subcortical; maps_18m$cortical; maps_18m$subcortical
ggsave("Fig5a_brain_subcortical_12m.png", maps_12m$subcortical, width = 7, height = 4, dpi = 300)
ggsave("Fig5b_brain_cortical_18m.png",    maps_18m$cortical,    width = 9, height = 4, dpi = 300)
ggsave("Fig5c_brain_subcortical_18m.png", maps_18m$subcortical, width = 7, height = 4, dpi = 300)

# One colour per W-score feature; regions not in the data stay grey
w_fill <- scale_fill_manual(values = c("Temporo-occipital composite" = "#2166AC", "Hippocampus" = "#B2182B"),
                            na.value = "grey90", na.translate = FALSE, name = "W-score feature")

p_w_cort <- brain_plot(tibble(region  = c("middle temporal", "fusiform", "lateral occipital"),
                              feature = "Temporo-occipital composite"),
                       dk_atlas, "feature", w_fill, "W-score MRI features (18-month): cortex")
p_w_sub  <- brain_plot(tibble(region = "Hippocampus", feature = "Hippocampus"),
                       aseg_atlas, "feature", w_fill, "W-score MRI features (18-month): subcortex")
p_w_cort; p_w_sub
ggsave("Fig6a_Wscore_features_cortex_18m.png",    p_w_cort, width = 9, height = 4, dpi = 300)
ggsave("Fig6b_Wscore_features_subcortex_18m.png", p_w_sub,  width = 7, height = 4, dpi = 300)
p_w_cort <- brain_plot(tibble(region  = c("middle temporal", "fusiform", "lateral occipital"),
                              feature = "Temporo-occipital composite"),
                       dk_atlas, "feature", w_fill, "W-score MRI features (18-month): cortex")
p_w_sub  <- brain_plot(tibble(region = "Hippocampus", feature = "Hippocampus"),
                       aseg_atlas, "feature", w_fill, "W-score MRI features (18-month): subcortex")
p_w_cort; p_w_sub
ggsave("Fig6a_Wscore_features_cortex_18m.png",    p_w_cort, width = 9, height = 4, dpi = 300)
ggsave("Fig6b_Wscore_features_subcortex_18m.png", p_w_sub,  width = 7, height = 4, dpi = 300)

##########
audit <- tribble(
  ~Step, ~Rule, ~`12-month`, ~`18-month`,
  "1. Candidate regions",
  "Literature-based FreeSurfer volumes (UCSFFSL51), ICV-normalised: 15 bilateral pairs + CorticalGM + TotalCortVol + ThirdVentricle",
  "33", "33",
  "2. Bilateral composites",
  "Left + right summed to reduce measurement noise; temporo-occipital composite (middle temporal + fusiform + lateral occipital) added from slope correlations",
  as.character(nrow(feature_table_comp)), as.character(nrow(feature_table_comp)),
  "3. LME per feature",
  "Volume ~ time × baseline diagnosis + (1 + time | RID), REML; at least 2 visits within the window",
  "N = 460", "N = 459",
  "4. Reliability rule",
  "Reliability ≥ 0.20, non-singular fit, tissue shrinking and ventricles expanding in both CN and MCI",
  as.character(sum(feature_table_comp$keep_12m)), as.character(sum(feature_table_comp$keep_18m)),
  "5. Redundancy check",
  "Feature pairs with |r| ≥ 0.80 flagged",
  paste(nrow(redund_12m), "pair (lateral vs inferior lateral ventricle)"),
  paste(nrow(redund_18m), "pair (lateral vs inferior lateral ventricle)"),
  "6. W-score selection",
  "Following van Loenhoud 2019: one ventricle, hippocampus, regional cortex, global cortex; one regression per feature",
  as.character(ncol(wscore_mri_12m) - 1), as.character(ncol(wscore_mri_18m) - 1)
)
print(audit, width = Inf)
write_csv(audit, "Audit_MRI_feature_selection.csv")

# Excluded features in each window, with the first rule each one failed
excluded <- feature_table_comp %>%
  pivot_longer(-region, names_to = c(".value", "window"), names_pattern = "(.*)_(12m|18m)") %>%
  filter(!keep) %>%
  mutate(Reason = case_when(singular            ~ "Singular fit",
                            reliability < 0.20 ~ "Reliability < 0.20",
                            TRUE               ~ "Implausible direction of change"),
         Region = recode(region, !!!region_labels)) %>%
  select(Window = window, Region, Reliability = reliability, Reason) %>%
  arrange(Window, desc(Reliability))
print(excluded, n = Inf)
write_csv(excluded, "Audit_MRI_excluded_features.csv")

















