library(dplyr)
library(ggplot2)
library(lme4)
library(stringr)
library(emmeans)
library(tibble)
library(tidyr)
library(nlme)

remotes::install_github("nx10/httpgd")
library(httpgd)

# For VLGUT1-PSD95
metrics_to_assess = c("local_peak_colocalized_spots", # model_1 works, CA3 and DG layers are sign
                      "overlap_coeff", # model_1: works
                      "overlap_um2", # model_1: boundary (singular) fit: see help('isSingular'), model_2 works
                      "pearson_cor", # model_1: works
                      "presynapse_image_mfi", # model_1: works, CA3 and DG layers are sign
                      "postsynapse_image_mfi", # model_1: works, CA3 SL is sign
                      "pre_puncta_density_per_100_um2", # model_1: boundary (singular) fit: see help('isSingular'), model_2 works
                      "post_puncta_density_per_100_um2", # model_1: works
                      "post_puncta_density_per_100_um2", # model_1: works
                      "pre_staining_area_um2", # model_1: works
                      "post_staining_area_um2", # model_1: works
                      "pre_mean_puncta_size_um2", # model_1: works
                      "post_mean_puncta_size_um2") # model_1: works
# For VGAT-GEPH
metrics_to_assess = c("local_peak_colocalized_spots", # model_1: boundary (singular) fit: see help('isSingular'), model_2 works but not quite intuitive
                      "overlap_coeff", # model_1: boundary (singular) fit: see help('isSingular'), model_2 works
                      "overlap_um2", # model_1: boundary (singular) fit: see help('isSingular'), model_2 works
                      "pearson_cor", # model_1: works
                      "presynapse_image_mfi", # model_1: works, CA1 SR, CA3 SL, CA3 SO and CA3 SR are sign
                      "postsynapse_image_mfi", # model_1: boundary (singular) fit: see help('isSingular'), CA3 SR is sign, model_2 works, CA3 SR is sign
                      "pre_puncta_density_per_100_um2", # model_1: works, CA1 SO is sign
                      "post_puncta_density_per_100_um2", # model_1: boundary (singular) fit: see help('isSingular'), model_2 works
                      "pre_staining_area_um2", # model_1: works
                      "post_staining_area_um2", # model_1: boundary (singular) fit: see help('isSingular'), model_2 works
                      "pre_mean_puncta_size_um2", # model_1: works
                      "post_mean_puncta_size_um2") # model_1: works


metric_assessed = "local_peak_colocalized_spots"

# loading the data
vcam1_data <- read.csv("/mnt/image-analysis/synapse-counting/VCAM1/VCAM1-LacZ_VGAT-GEPH_output_data/metric_results.csv")
dim(vcam1_data)

# data preprocessing for only one dependent variable (one synaptic metric, here it is presynaptic MFI)
vcam1_data_metric_assessed <- vcam1_data %>%
    mutate(
    section = str_extract(img_filename, "section-\\d+")
    ) %>%
    select(
    metric_assessed,
    gRNA,
    hippocampal_layer,
    section,
    Brain
    )


# OPTION 1:

# calculating first the log2 ratio on section level & then taking the mean per brain
log2_ratio_section <- vcam1_data_metric_assessed %>%
  pivot_wider(
    names_from = gRNA,
    values_from = local_peak_colocalized_spots
  ) %>% 
  group_by(Brain, section, hippocampal_layer) %>%
  mutate(log2_ratio = log2(`VCAM1-gRNA`/`LacZ-gRNA`)) %>%
  group_by(Brain, hippocampal_layer) %>%
  mutate(!!metric_assessed := mean(log2_ratio)) %>%
  select(hippocampal_layer, Brain, !!metric_assessed) %>%
  distinct()

# dynamically new column name
new_metric_col_name = paste0(metric_assessed, "_log2_ratio")

log2_ratio_section <- log2_ratio_section %>%
  rename(!!new_metric_col_name := metric_assessed)



# OPTION 2:

# calculating the mean per brain
mean_per_brain_layer <- vcam1_data_metric_assessed %>%
  group_by(hippocampal_layer, Brain, gRNA) %>%
  mutate(mean_metric = mean(!!rlang::sym(metric_assessed))) %>%
  select(gRNA, hippocampal_layer, Brain, mean_metric) %>%
  distinct()

# calculating the log2 ratio
log2_ratio <- mean_per_brain_layer %>%
  pivot_wider(
    names_from = gRNA,
    values_from = mean_metric
  ) %>%
  group_by(hippocampal_layer, Brain) %>%
  mutate(log2_ratio = log2(`VCAM1-gRNA`/`LacZ-gRNA`))

# dynamically new column name
new_metric_col_name = paste0(metric_assessed, "_log2_ratio")

log2_ratio <- log2_ratio %>%
  rename(!!new_metric_col_name := log2_ratio)



# # making the model
# model <- lmer(local_peak_colocalized_spots_log2_ratio ~ hippocampal_layer +    # fixed effects with hippocampal layer
#               (1|Brain),                          # between-brain variation because of perfusion/viral injection
#               data = log2_ratio,
#               REML=TRUE)

# for assessing all metrics
model_1 <- lmer(formula = as.formula(paste0(metric_assessed, "_log2_ratio ~ hippocampal_layer + 
                                                            (1|Brain)")),                             # between-brain variation because of perfusion/viral injection
              data = log2_ratio,
              REML=TRUE)

model_2 <- lm(formula = as.formula(paste0(metric_assessed, "_log2_ratio ~ hippocampal_layer")),
               data = log2_ratio)


# first nlme model
model_3 <- lme(
  fixed = local_peak_colocalized_spots_log2_ratio ~ hippocampal_layer,
  data = log2_ratio_section,
  random = ~1 | Brain,
  weights = varIdent(form = ~1 | hippocampal_layer)
)

# checking the model
summary(model_3)
plot(model_3)
qqnorm(resid(model_3))
qqline(resid(model_3))

# extracting the estimated mean per layer (without reference)
em_means <- emmeans(model_3, ~hippocampal_layer)
summary(em_means)

# producing the output table with reference layer
test_output <- test(em_means)
results_df <- as.data.frame(test_output)
results_df$p.value.adj = p.adjust(results_df$p.value, method = "fdr")
