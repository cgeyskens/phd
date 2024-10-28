# test script for linear mixed effect models
# now just trying out a lme for a single synaptic metric

# loading the data
vcam1_data <- read.csv("/Volumes/KINGSTON/data/phd/image-analysis/synapse-counting/VCAM1/VCAM1-LacZ_VGAT-GEPH_output_data/metric_results.csv")

# getting the libraries
library(dplyr)
library(ggplot2)
library(lme4)

# data preprocessing for only one dependent variable (one synaptic metric, here it is presynaptic MFI)
vcam1_data_vglut_mfi <- vcam1_data %>%
    mutate(
    section = str_extract(img_filename, "section-\\d+")
    ) %>%
    select(
    local_peak_colocalized_spots,
    gRNA,
    hippocampal_layer,
    section,
    Brain
    ) %>%
    filter(hippocampal_layer =="DG ML")

# standardizing the values
vcam1_data_vglut_mfi$local_peak_colocalized_spots_standardized <- standardize(vcam1_data_vglut_mfi$local_peak_colocalized_spots)

# checking distribution
hist(vcam1_data_vglut_mfi$local_peak_colocalized_spots_standardized)

# modifying the gRNA and hippocampal_layer column from character to factor for one hot encoding
vcam1_data_vglut_mfi$gRNA <- as.factor(vcam1_data_vglut_mfi$gRNA) 
vcam1_data_vglut_mfi$hippocampal_layer <- as.factor(vcam1_data_vglut_mfi$hippocampal_layer) 
class(vcam1_data_vglut_mfi$gRNA)
class(vcam1_data_vglut_mfi$hippocampal_layer)

# checking all the values in a boxplot
boxplot(vcam1_data_vglut_mfi$local_peak_colocalized_spots_standardized ~ vcam1_data_vglut_mfi$gRNA, col="skyblue", xlab="gRNA", ylab="variable")

# trying to model with a simple linear regression
basic.lm <- lm(local_peak_colocalized_spots_standardized ~ gRNA, data = vcam1_data_vglut_mfi)
summary(basic.lm)

# nested mixed effects model
mixed.lm <- lmer(local_peak_colocalized_spots_standardized ~ gRNA + 
                (1|Brain/section/hippocampal_layer), 
                data = vcam1_data_vglut_mfi)
summary(mixed.lm)

# model: interaction between gRNA and layer
mixed.lm2 <- lmer(local_peak_colocalized_spots_standardized ~ gRNA * hippocampal_layer +
                  (1|Brain/section/hippocampal_layer),
                  data = vcam1_data_vglut_mfi)
summary(mixed.lm2)

# improved model specification
improved_model <- lmer(local_peak_colocalized_spots_standardized ~ gRNA * hippocampal_layer + 
                     (1|Brain) +                             # between-brain variation
                     (1|Brain:gRNA) +                        # paired design of hemispheres
                     (1|Brain:section:hippocampal_layer),    # nested measurements
                     data = vcam1_data_vglut_mfi,
                     REML=TRUE)
summary(improved_model)

plot(improved_model)
qqnorm(resid(improved_model))
qqline(resid(improved_model))

anova(improved_model, mixed.lm2)

# Extract p-values and coefficient names
coef_table <- summary(improved_model)$coefficients
p_values <- coef_table[,"Pr(>|t|)"]
coef_names <- rownames(coef_table)

# Apply Holm-Bonferroni correction
p_adjusted_holm <- p.adjust(p_values, method = "fdr")

# Create a comparison table
results_table <- data.frame(
  Coefficient = coef_names,
  Original_P = p_values,
  Adjusted_P = p_adjusted_holm,
  Significant = p_adjusted_holm < 0.05
)

# Print results with all p-values
print(data.frame(
  Layer = results_df$layer,
  Estimate = round(results_df$estimate, 3),
  T_value = round(results_df$t_value, 3),
  P_value_uncorrected = round(results_df$p_value, 4),
  P_value_bonferroni = round(results_df$p_value_bonferroni, 4),
  P_value_FDR = round(results_df$p_value_fdr, 4)
))