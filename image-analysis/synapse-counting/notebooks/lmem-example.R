# Singular fit warnings:
# For VGLUT1-PSD95 with least singular fit warnings
metrics_to_assess = c("local_peak_colocalized_spots", # ok with model_2, model_3 and improved_model
                      "overlap_coeff", # ok with model_2 and model_3
                      "overlap_um2", # ok with model_2 and model_3
                      "pearson_cor", # ok with model_1, model_2 and improved_model
                      "presynapse_image_mfi", # ok with model_2, model_3 and improved_model
                      "postsynapse_image_mfi", # ok with model_1, model_2 and model_3 and improved_model
                      "pre_puncta_density_per_100_um2", # only model_2 and model_3
                      "post_puncta_density_per_100_um2", # only model_2 and model_1
                      "pre_staining_area_um2", # ok with model_2, model_3 and model_4
                      "post_staining_area_um2", # ok with model_1, model_2 and improved_model
                      "pre_mean_puncta_size_um2", # ok with model_2, model_3 and improved_model
                      "post_mean_puncta_size_um2") # ok with model_2, model_3 and improved_model
# For VGAT-GEPH
metrics_to_assess = c("local_peak_colocalized_spots", # ok with improved_model
                      "overlap_coeff", # no model
                      "overlap_um2", # ok with model_2 and model_3
                      "pearson_cor", # ok with model_2, model_3 and improved_model
                      "presynapse_image_mfi", # ok with model_1 and improved_model
                      "postsynapse_image_mfi", # no model
                      "pre_puncta_density_per_100_um2", # ok with model_2, model_3 and improved_model
                      "post_puncta_density_per_100_um2", # ok with model_2
                      "pre_staining_area_um2", # ok with model_2, model_3 and improved_model
                      "post_staining_area_um2", # ok with model_2
                      "pre_mean_puncta_size_um2", # ok with model_2 and improved_model
                      "post_mean_puncta_size_um2") # ok with model_2, model_3 and improved_model


# getting the libraries
library(dplyr)
library(ggplot2)
library(lme4)
library(stringr)
library(emmeans)

# loading the data
vcam1_data <- read.csv("/Volumes/KINGSTON/data/phd/image-analysis/synapse-counting/VCAM1/VCAM1-LacZ_VGLUT1-PSD95_output_data/metric_results.csv")

# preprocessing the data
vcam1_data_metric_assessed <- vcam1_data %>%
    mutate(
    section = str_extract(img_filename, "section-\\d+")
    ) %>%
    select(
    local_peak_colocalized_spots,
    gRNA,
    hippocampal_layer,
    section,
    Brain
    )

# check the distribution
hist(vcam1_data_metric_assessed$local_peak_colocalized_spots)

# check the distribution
vcam1_data_metric_assessed$local_peak_colocalized_spots_log <- log2(vcam1_data_metric_assessed$local_peak_colocalized_spots)

hist(vcam1_data_metric_assessed$local_peak_colocalized_spots_log)

# scaling the data
vcam1_data_metric_assessed$local_peak_colocalized_spots_log_scaled <- scale(vcam1_data_metric_assessed$local_peak_colocalized_spots_log)
hist(vcam1_data_metric_assessed$local_peak_colocalized_spots_log_scaled)

vcam1_data_metric_assessed$local_peak_colocalized_spots_scaled <- scale(vcam1_data_metric_assessed$local_peak_colocalized_spots)



# modifying the gRNA and hippocampal_layer column from character to factor for one hot encoding
vcam1_data_metric_assessed$gRNA <- as.factor(vcam1_data_metric_assessed$gRNA) 
vcam1_data_metric_assessed$hippocampal_layer <- as.factor(vcam1_data_metric_assessed$hippocampal_layer) 
class(vcam1_data_metric_assessed$gRNA)
class(vcam1_data_metric_assessed$hippocampal_layer)

# model 1
model_1 <- lmer(local_peak_colocalized_spots_log_scaled ~ gRNA * hippocampal_layer + 
                (1 + gRNA | Brain), 
                data = vcam1_data_metric_assessed)

# model 2
model_2 <- lmer(local_peak_colocalized_spots_log_scaled ~ gRNA * hippocampal_layer + 
                     (1|Brain) +                     # between-brain variation because of perfusion/viral injection                      # paired design of hemispheres
                     (1|Brain:hippocampal_layer),    # captures layer specific variation within each brain
                     data = vcam1_data_metric_assessed,
                     REML=TRUE)

# model 3
model_3 <- lmer(local_peak_colocalized_spots_log_scaled ~ gRNA * hippocampal_layer + 
                     (1|Brain) +                                   # between-brain variation because of perfusion/viral injection                      # paired design of hemispheres
                     (1|Brain:hippocampal_layer) +              # captures layer specific variation within each brain
                     (1|Brain:section:hippocampal_layer),                                   # captures layer specific variation within each brain
                     data = vcam1_data_metric_assessed,
                     REML=TRUE)

# improved_model
improved_model <- lmer(local_peak_colocalized_spots_log_scaled ~ gRNA + hippocampal_layer + 
                     (1|Brain) +                             # between-brain variation because of perfusion/viral injection
                     (1|Brain:gRNA) +                        # paired design of hemispheres
                     (1|Brain:section:hippocampal_layer),    # nested measurements
                     data = vcam1_data_metric_assessed,
                     REML=TRUE)

improved_model_wolog <- lmer(local_peak_colocalized_spots_scaled ~ gRNA + hippocampal_layer + 
                     (1|Brain) +                             # between-brain variation because of perfusion/viral injection
                     (1|Brain:gRNA) +                        # paired design of hemispheres
                     (1|Brain:section:hippocampal_layer),    # nested measurements
                     data = vcam1_data_metric_assessed,
                     REML=TRUE)


summary(improved_model)

plot(model_3)
qqnorm(resid(model_3))
qqline(resid(model_3))

# comparing models
anova(model_1, model_2, model_3, improved_model, improved_model_wolog)

# pairwise comparison
emm_results <- emmeans(glm_model, pairwise ~ gRNA | hippocampal_layer)

contrasts_data <- data.frame(emm_results$contrasts)
means_data <- data.frame(emm_results$emmeans)

results_df <- data.frame(
  Layer = contrasts_data$hippocampal_layer,
  LacZ_mean = means_data$emmean[seq(1, nrow(means_data), 2)],  # Get LacZ means
  VCAM1_mean = means_data$emmean[seq(2, nrow(means_data), 2)], # Get VCAM1 means
  Estimate = contrasts_data$estimate,
  SE = contrasts_data$SE,
  t_ratio = contrasts_data$t.ratio,
  p_value = contrasts_data$p.value,
  p_adj_FDR = p.adjust(contrasts_data$p.value, method = "fdr")
)

glm_model <- glmer(local_peak_colocalized_spots ~ gRNA * hippocampal_layer + 
                     (1|Brain) +                             # between-brain variation
                     (1|Brain:gRNA) +                        # paired design of hemispheres
                     (1|Brain:section:hippocampal_layer),    # nested measurements
                     data = vcam1_data_metric_assessed,
                     family = poisson)


