library(dplyr)
library(ggplot2)
library(lme4)
library(stringr)
library(emmeans)
library(tibble)
library(tidyr)
library(nlme)
library(stringr)

#### =============================== arguments =============================== ####

protein = "VCAM1"
synapse_type = "VGLUT1-PSD95"

#### =============================== get input =============================== ####
if (protein == "VCAM1" && synapse_type == "VGLUT1-PSD95"){
    data <- read.csv("/mnt/image-analysis/synapse-counting/IHC_Exp9_VCAM1_VGLUT1-PSD95_output_data_20251201_194547/metric_results.csv")
} else if (protein == "VCAM1" && synapse_type == "VGAT-GEPH") {
   data <- read.csv("/mnt/image-analysis/synapse-counting/IHC_Exp12_VCAM1_VGAT-GEPH_output_data_20251201_194539/metric_results.csv")
} else if (protein == "GPR37L1" && synapse_type == "VGLUT1-PSD95") {
    data <- read.csv("/mnt/image-analysis/synapse-counting/IHC_Exp10_GPR37L1_VGLUT1-PSD95_output_data_20251201_194537/metric_results.csv")
} else if (protein == "GPR37L1" && synapse_type == "VGAT-GEPH") {
    data <- read.csv("/mnt/image-analysis/synapse-counting/IHC_Exp11_GPR37L1_VGAT-GEPH_output_data_20251201_194548/metric_results.csv")
}

#### =============================== data wrangling =============================== ####

metrics <- c(
  "local_peak_colocalized_spots",
  "overlap_coeff",
  "overlap_um2",
  "pearson_cor",
  "presynapse_image_mfi",
  "postsynapse_image_mfi",
  "pre_puncta_density_per_100_um2",
  "post_puncta_density_per_100_um2",
  "pre_staining_area_um2",
  "post_staining_area_um2",
  "pre_mean_puncta_size_um2",
  "post_mean_puncta_size_um2"
)

data <- data %>%
  mutate(sample_ID = str_c(gRNA, Brain, sep = "_"))

grouped_brain <- data %>%
  group_by(sample_ID, hippocampal_layer) %>%
  summarise(
    across(all_of(metrics), ~ mean(.x, na.rm = TRUE)),
    .groups = "drop"
  )

df_melted <- grouped_brain %>%
  pivot_longer(
    cols = all_of(metrics),
    names_to = "metric",
    values_to = "value"
  )

df_melted <- df_melted %>%
  mutate(
    synapse_hippocampal_layer_metric =
      str_c(synapse_type, hippocampal_layer, metric, sep = "_")
  )

df_pivot <- df_melted %>%
  select(sample_ID, synapse_hippocampal_layer_metric, value) %>%
  pivot_wider(
    names_from = synapse_hippocampal_layer_metric,
    values_from = value
  )

## in the end you need to aggregated the vglut1 and vgat data


process_data_for_metric <- function(data, metric){
    data_metric_assessed <- data %>%
        mutate(
            section = str_extract(img_filename, "section-\\d+")
        ) %>%
        select(
            !!sym(metric),
            gRNA,
            hippocampal_layer,
            section,
            Brain
        )
    data_metric_assessed$gRNA <- factor(data_metric_assessed$gRNA)
    data_metric_assessed$hippocampal_layer <- factor(data_metric_assessed$hippocampal_layer)
    data_metric_assessed$section <- factor(data_metric_assessed$section)
    data_metric_assessed$Brain <- factor(data_metric_assessed$Brain)
    
    return(data_metric_assessed)
}

data_2 <- process_data_for_metric(data, metric = "local_peak_colocalized_spots")

model_for_pca <- lmer(
  formula = as.formula(paste0(metric_assessed, " ~ 1 + (1 | Brain)")),
  data = data_metric,
  control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000))
)