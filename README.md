# TEMPO: Temporal Evaluation of Multimodal Progression Outcomes in Mild Cognitive Impairment and Cognitively Normal Individuals

**(AAN 2027 submission)** · Abstract submitted to the American Academy of Neurology (AAN) Annual Meeting 2027

![R](https://img.shields.io/badge/language-R-276DC3) ![Data](https://img.shields.io/badge/data-ADNI-informational) ![License](https://img.shields.io/badge/license-MIT-green)

---

## Overview

TEMPO looks at how brain structure, cognition and everyday functioning **change over time** in people with **mild cognitive impairment (MCI)** compared with **cognitively normal (CN)** individuals. Using longitudinal data from the Alzheimer's Disease Neuroimaging Initiative (ADNI), it estimates rates of change over **12-month** and **18-month** observation windows and tests which features worsen faster in MCI.

## Objective

To identify MRI, cognitive and functional features that show greater longitudinal worsening in individuals with MCI than in CN individuals over 12- and 18-month observation windows.

## Background

Describing the multimodal changes that separate CN and MCI groups may improve understanding of the longitudinal patterns associated with early neurodegenerative change.

## Design / Methods

| Component | Details |
|---|---|
| **Data source** | Alzheimer's Disease Neuroimaging Initiative (ADNI), ADNIMERGE / ADNIMERGE2 |
| **Imaging** | 3-T MRI processed with the FreeSurfer longitudinal pipeline (UCSF FreeSurfer 5.1, `UCSFFSL51`) |
| **Groups** | Baseline CN vs baseline MCI. Excluded: baseline dementia, MCI→CN reverters and converters within the observation window |
| **Windows** | 12 months and 18 months of follow-up (≥ 2 scans per participant) |
| **Statistical model** | Linear mixed-effects models (`lme4` / `lmerTest`): `outcome ~ time_years * BL_DIAGNOSIS + (1 + time_years | RID)`. These handle repeated measures within individuals and unequal follow-up |
| **MRI rate of change** | Annual percentage change (APC) relative to diagnosis-specific baseline means; volumes normalised to intracranial volume (ICV) |
| **Between-group test** | Time × diagnosis interaction, with Benjamini–Hochberg FDR correction within each window |

**18 MRI features** (bilateral regions averaged):

- **Ventricular (3):** lateral ventricles, inferior lateral ventricles, third ventricle
- **Cortical (6):** entorhinal, fusiform, middle temporal, precuneus, insula, lateral occipital
- **Subcortical (7):** hippocampus, amygdala, thalamus, putamen, caudate, nucleus accumbens, pallidum
- **Global (2):** total gray matter volume, total cortical volume (cortical GM + cortical WM)

**Cognitive and functional outcomes:**

- ADAS-Cog 13
- ADNI composite scores: memory (ADNI-MEM), executive function (ADNI-EF), language (ADNI-LAN), visuospatial (ADNI-VS)
- Functional Activities Questionnaire (FAQ) for everyday functioning

## Key results

**12-month window**

- MCI showed greater relative volume loss in the **entorhinal cortex (−2.80 %/yr)** and **insula (−1.20 %/yr)**, plus greater ventricular expansion than CN.
- The **insula** was the only region with a significant between-group APC difference after correction: **−1.03 %/yr (95 % CI −1.67 to −0.39; adjusted p = 0.032)**.
- Cognitive and functional rates of decline were similar between groups after correction.

**18-month window**

| Feature | MCI vs CN difference in APC | Adjusted p |
|---|---|---|
| Lateral ventricles | +1.61 %/yr | ≤ 0.032 |
| Fusiform gyrus | −0.90 %/yr | ≤ 0.032 |
| Entorhinal cortex | −1.17 %/yr | ≤ 0.032 |
| Middle temporal gyrus | −0.73 %/yr | ≤ 0.032 |
| Insula | −0.53 %/yr | ≤ 0.032 |

- MCI showed faster decline in **memory (adjusted p = 0.0061)** and **everyday functioning, FAQ (adjusted p = 0.0062)**.

## Conclusions

MCI was associated with greater relative changes in selected ventricular, temporal and limbic structures, with more regional differences detectable over 18 months than over 12 months. Greater memory decline and worsening everyday functioning were also evident in MCI over 18 months. These findings support longitudinal multimodal assessment to characterise the features that distinguish MCI from cognitively normal status.

---

## Repository structure

```
.
├── config.R                 # Data and output paths (edit this first)
├── install_packages.R       # Installs the required R packages
├── scripts/
│   ├── 01_freesurfer_data_audit_and_merge.R    # FreeSurfer table audit (version, field strength, duplicates,
│   │                                           # longitudinal coverage); merge with diagnosis and cognitive data
│   ├── 02_adni_clinical_data_audit.R           # Unit-of-observation checks for MMSE, ADAS, FAQ, DXSUM, PTDEMOG
│   ├── 03_cohort_definition_and_baseline.R     # Baseline diagnosis, trajectory groups, early converters,
│   │                                           # baseline characteristics
│   ├── 04_mri_feature_selection.R              # 18 bilateral/global features, total cortical volume, cohort filters
│   ├── 05_mri_trajectories_lme_apc.R           # ICV normalisation, 12/18-month cohorts, LME slopes, APC,
│   │                                           # slope reliability, final MRI metrics table
│   ├── 06_cognitive_trajectories_adas13.R      # ADAS-Cog 13 LME timed from first MRI; visit-timing figures;
│   │                                           # sensitivity analyses
│   ├── 07_cognitive_domains_and_function.R     # ADNI-MEM/EF/LAN/VS composites and FAQ; extended-window check
│   └── 08_figures.R                            # APC, reliability and ggseg brain-map figures
├── figures/                 # Generated figures (not tracked)
└── tables/                  # Generated tables (not tracked)
```

## How to run

1. **Get ADNI data access.** Apply at [adni.loni.usc.edu](https://adni.loni.usc.edu/data-samples/access-data/), then download the ADNIMERGE tables (including `UCSFFSL51.csv`, `DXSUM.csv`, `ADAS.csv`, `FAQ.csv`, `PTDEMOG.csv`) and the ADNIMERGE2 data.
2. **Set the paths** in `config.R`, or set the environment variables `ADNIMERGE_DIR` and `ADNIMERGE2_DIR`.
3. **Install the packages:**
   ```r
   source("install_packages.R")
   ```
4. **Run the scripts in order** (`01` → `08`) from the project root, in the same R session. Later scripts reuse objects created by earlier ones (for example `master_visit`, `trajectory_cohort`, `lme_12m`, `lme_18m_data`).

> **Note:** These scripts were written as an interactive, step-by-step analysis workflow, with audit and integrity checks kept inline. Run them one after another in a single session rather than as standalone programs.

**Software:** R (≥ 4.2) with tidyverse, data.table, lme4, lmerTest, corrplot, ggrepel, ggseg, patchwork, here.

## Data availability

The data used here come from the **Alzheimer's Disease Neuroimaging Initiative (ADNI)** database (adni.loni.usc.edu). Under the ADNI Data Use Agreement, **no ADNI data are included in this repository**. `.gitignore` excludes all data files. Researchers can apply for access directly through ADNI.

## Acknowledgement

Data collection and sharing for this project were funded by the Alzheimer's Disease Neuroimaging Initiative (ADNI). As such, the investigators within ADNI contributed to the design and implementation of ADNI and/or provided data but did not participate in the analysis or writing of this work. A complete listing of ADNI investigators can be found on the ADNI website.

## Citation

If you use this code, please cite:

```
Kola PR. TEMPO: Temporal Evaluation of Multimodal Progression Outcomes in Mild Cognitive
Impairment and Cognitively Normal Individuals. Abstract submitted to the American Academy
of Neurology Annual Meeting, 2027.
```

## License

Code released under the [MIT License](LICENSE).
