# test script for linear mixed effect models
# now just trying out a lme for a single synaptic metric
# getting the libraries
library(dplyr)
library(ggplot2)
library(lme4)
library(stringr)
library(emmeans)
library(tibble)
library(tidyr)

#remotes::install_github("nx10/httpgd")
library(httpgd)
hgd()

metrics_to_assess = c("local_peak_colocalized_spots",
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
                      "post_mean_puncta_size_um2")

metric_assessed <- "overlap_coeff"

# loading the data
vcam1_data <- read.csv("/mnt/image-analysis/synapse-counting/IHC_Exp9_VCAM1_VGLUT1-PSD95_output_data_20251201_145813/metric_results.csv")
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

# Log-transform the raw data
vcam1_data_metric_assessed[[paste0(metric_assessed, "_log")]] <- log(vcam1_data_metric_assessed[[metric_assessed]] + 1)

# standardizing the values
vcam1_data_metric_assessed[[paste0(metric_assessed, "_scaled")]] <- scale(vcam1_data_metric_assessed[paste0(metric_assessed, "_log")])

# checking distribution
# boxplot(vcam1_data_metric_assessed$presynapse_image_mfi)

# modifying the gRNA and hippocampal_layer column from character to factor for one hot encoding
vcam1_data_metric_assessed$gRNA <- as.factor(vcam1_data_metric_assessed$gRNA) 
vcam1_data_metric_assessed$hippocampal_layer <- as.factor(vcam1_data_metric_assessed$hippocampal_layer) 
class(vcam1_data_metric_assessed$gRNA)
class(vcam1_data_metric_assessed$hippocampal_layer)

# checking all the values in a boxplot
plot(vcam1_data_metric_assessed[[metric_assessed]] ~ vcam1_data_metric_assessed$gRNA, col="skyblue", xlab="gRNA", ylab="variable")

# model
improved_model <- lmer(
    formula = as.formula(paste0(metric_assessed, "_scaled ~ gRNA * hippocampal_layer + 
                                (1|Brain) + (1|Brain:gRNA) + (1|Brain:section:hippocampal_layer)")),
    data = vcam1_data_metric_assessed,
    REML = TRUE)

summary(improved_model)

# The best model to fit my data
improved_model <- lmer(metric_assessed ~ gRNA * hippocampal_layer +   # fixed effects with interaction between gRNA & hippocampal layer
                     (1|Brain) +                             # between-brain variation because of perfusion/viral injection
                     (1|Brain:gRNA) +                        # paired design of hemispheres
                     (1|Brain:section:hippocampal_layer),    # nested measurements
                     data = vcam1_data_metric_assessed,
                     REML=TRUE)
                     
summary(improved_model)

plot(improved_model)
qqnorm(resid(improved_model))
qqline(resid(improved_model))

# checking pairwise comparisons
emm_results <- emmeans(improved_model, pairwise ~ gRNA | hippocampal_layer)

contrasts_data <- data.frame(emm_results$contrasts)
means_data <- data.frame(emm_results$emmeans)

results_df <- data.frame(
  Layer = contrasts_data$hippocampal_layer,
  LacZ_mean = means_data$emmean[seq(1, nrow(means_data), 2)],  # Get LacZ means
  VCAM1_mean = means_data$emmean[seq(2, nrow(means_data), 2)], # Get VCAM1 means
  Estimate = contrasts_data$estimate,
  SE = contrasts_data$SE,
  t_ratio = contrasts_data$z.ratio,
  p_value = contrasts_data$p.value,
  p_adj_FDR = p.adjust(contrasts_data$p.value, method = "fdr")
  )


# Other models tested:

# trying to model with a simple linear regression
basic_lm <- lm(presynapse_image_mfi_scaled ~ gRNA, data = vcam1_data_metric_assessed)
summary(basic_lm)

# nested mixed effects model
mixed_lm <- lmer(presynapse_image_mfi_scaled ~ gRNA + 
                (1|Brain/section/hippocampal_layer), 
                data = vcam1_data_metric_assessed)
summary(mixed_lm)

# model: interaction between gRNA and layer
mixed_lm2 <- lmer(presynapse_image_mfi_scaled ~ gRNA * hippocampal_layer +
                  (1|Brain/section/hippocampal_layer),
                  data = vcam1_data_metric_assessed)
summary(mixed_lm2)

## Generalized models
hist(vcam1_data_metric_assessed$presynapse_image_mfi)

# taking into accound the non-normality of the data WITH ORIGINAL DATA
glm_model <- glmer(presynapse_image_mfi ~ gRNA * hippocampal_layer + 
                     (1|Brain) +                             # between-brain variation
                     (1|Brain:gRNA) +                        # paired design of hemispheres
                     (1|Brain:section:hippocampal_layer),    # nested measurements
                     data = vcam1_data_metric_assessed,
                     family = Gamma(link = "log"),
                     control = glmerControl(optimizer = "bobyqa",
                                            optCtrl = list(maxfun = 200000)))
summary(glm_model)

# QC
plot(glm_model)
qqnorm(resid(glm_model))
qqline(resid(glm_model))



###############
################ Function

library(dplyr)
library(ggplot2)
library(lme4)
library(stringr)
library(emmeans)

synapse_analysis <- function(data_file, metric_to_assess) {
  # Load the data
  data <- read.csv(data_file)
  
  # Preprocess the data
  data <- data %>%
    mutate(
      section = str_extract(img_filename, "section-\\d+")
    ) %>%
    select(
      metric_to_assess,
      gRNA,
      hippocampal_layer,
      section,
      Brain
    )
  
  # Standardize the values
  data[[paste0(metric_to_assess, "_scaled")]] <- scale(data[[metric_to_assess]])
  
  # Convert factors
  data$gRNA <- as.factor(data$gRNA)
  data$hippocampal_layer <- as.factor(data$hippocampal_layer)
  
  # Fit the model
  model <- lmer(
    formula = as.formula(paste0(metric_to_assess, "_scaled ~ gRNA * hippocampal_layer + 
                               (1|Brain) + (1|Brain:gRNA) + (1|Brain:section:hippocampal_layer)")),
    data = data,
    REML = TRUE
  )
  
  # Summarize the model
  summary(model)
  
  # Diagnostic plots
  plot(model)
  qqnorm(resid(model))
  qqline(resid(model))
  
  # Pairwise comparisons
  emm_results <- emmeans(model, pairwise ~ gRNA | hippocampal_layer)
  
  contrasts_data <- data.frame(emm_results$contrasts)
  means_data <- data.frame(emm_results$emmeans)
  
  results_df <- data.frame(
    Layer = contrasts_data$hippocampal_layer,
    LacZ_mean = means_data$emmean[seq(1, nrow(means_data), 2)],
    VCAM1_mean = means_data$emmean[seq(2, nrow(means_data), 2)],
    Estimate = contrasts_data$estimate,
    SE = contrasts_data$SE,
    t_ratio = contrasts_data$t.ratio,
    p_value = contrasts_data$p.value,
    p_adj_FDR = p.adjust(contrasts_data$p.value, method = "fdr")
  )
  
  return(results_df)
}