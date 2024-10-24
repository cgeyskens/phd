# test script for linear mixed effect models
# now just trying out a lme for a single synaptic metric

# loading the data
vcam1_data <- read.csv("/Volumes/KINGSTON/data/phd/image-analysis/synapse-counting/VCAM1/VCAM1-LacZ_VGLUT1-PSD95_output_data/metric_results.csv")

# getting the libraries
library(dplyr)
library(ggplot2)
library(lme4)
library(stringr)
library(stargazer)
library(broom)
library(datawizard)
install.packages("emmeans")
library(emmeans)

# data preprocessing for only one dependent variable (one synaptic metric)
vcam1_data_vglut_mfi <- vcam1_data %>%
    mutate(
    section = str_extract(img_filename, "section-\\d+")
    ) %>%
    select(
    presynapse_image_mfi,
    gRNA,
    hippocampal_layer,
    section,
    Brain
    ) 

# checking the histogram
hist(vcam1_data_vglut_mfi$presynapse_image_mfi)

# standardizing the values
vcam1_data_vglut_mfi$presynapse_mfi_standardized <- standardize(vcam1_data_vglut_mfi$presynapse_image_mfi)
hist(vcam1_data_vglut_mfi$presynapse_mfi_standardized)

# modifying the gRNA and hippocampal_layer column from character to factor for one hot encoding
vcam1_data_vglut_mfi$gRNA <- as.factor(vcam1_data_vglut_mfi$gRNA) 
vcam1_data_vglut_mfi$hippocampal_layer <- as.factor(vcam1_data_vglut_mfi$hippocampal_layer) 
class(vcam1_data_vglut_mfi$gRNA)
class(vcam1_data_vglut_mfi$hippocampal_layer)

# checking all the values
boxplot(vcam1_data_vglut_mfi$presynapse_mfi_standardized ~ vcam1_data_vglut_mfi$gRNA, col="skyblue", xlab="pre_mfo", ylab="gRNA")

# trying to model with a simple linear regression
basic.lm <- lm(presynapse_mfi_standardized ~ gRNA, data = vcam1_data_vglut_mfi)
summary(basic.lm)

# nested histogram
mixed.lm <- lmer(presynapse_mfi_standardized ~ gRNA + (1|Brain/section/hippocampal_layer), data = vcam1_data_vglut_mfi)
summary(mixed.lm)

# improved model:

# Model: interaction between gRNA and layer
mixed.lm2 <- lmer(presynapse_mfi_standardized ~ gRNA * hippocampal_layer +
                  (1|Brain/section),
                  data = vcam1_data_vglut_mfi)
summary(mixed.lm2)
plot(mixed.lm2)

qqnorm(resid(mixed.lm2))
qqline(resid(mixed.lm2))

anova(mixed.lm, mixed.lm2)

library(ggplot2)

ggplot(vcam1_data_vglut_mfi, aes(x = hippocampal_layer, y = presynapse_mfi_standardized, color = gRNA)) +
  geom_point() +
  stat_summary(fun.data = mean_cl_boot, geom = "errorbar", width = 0.2) +
  stat_summary(fun = mean, geom = "line") +
  theme_minimal() +
  labs(title = "Effect of gRNA on Presynaptic Intensity Across Hippocampal Layers",
       x = "Hippocampal Layer", y = "Presynapse MFI (Standardized)")

install.packages("emmeans")
library(emmeans)
R.version

library(ggplot2)
library(dplyr)

# Create data frame with estimates, t-values and calculate p-values
results_df <- data.frame(
    layer = c("CA3 SO", "CA3 SL", "CA3 SR", "DG Hilus", "DG ML", "CA1 SO", "CA1 SR"),
    estimate = c(-1.27561, -1.03171, -0.67940, -0.63545, -0.40098, -0.28282, -0.20889),
    t_value = c(-4.234, -3.425, -2.255, -2.109, -1.331, -0.939, -0.693),
    std_error = c(0.30126, 0.30126, 0.30126, 0.30126, 0.30126, 0.30126, 0.30126)
)

# Calculate p-values (using degrees of freedom from model)
# df = number of observations - number of fixed effects parameters
df <- 176 - 16  # from your model summary
results_df$p_value <- 2 * (1 - pt(abs(results_df$t_value), df))

# Add significance stars
results_df$stars <- ifelse(results_df$p_value < 0.001, "***",
                          ifelse(results_df$p_value < 0.01, "**",
                                 ifelse(results_df$p_value < 0.05, "*", "ns")))

# Print results
print(results_df)

# Plot p-values
ggplot(results_df, aes(x = layer, y = -log10(p_value))) +
    geom_bar(stat = "identity") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
    theme_classic() +
    labs(x = "Hippocampal Layer", 
         y = "-log10(p-value)",
         title = "Statistical Significance by Layer") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

# For effect sizes (estimates) plot
ggplot(results_df, aes(x = layer, y = estimate)) +
    geom_bar(stat = "identity") +
    geom_errorbar(aes(ymin = estimate - std_error, 
                      ymax = estimate + std_error), 
                  width = 0.2) +
    theme_classic() +
    labs(x = "Hippocampal Layer", 
         y = "Effect Size (Estimate)",
         title = "Effect Sizes by Layer") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))


# Calculate adjusted p-values using different methods
results_df$p_value_bonferroni <- p.adjust(results_df$p_value, method = "bonferroni")
results_df$p_value_fdr <- p.adjust(results_df$p_value, method = "fdr")  # Benjamini-Hochberg

# Add adjusted significance stars
results_df$stars_bonferroni <- ifelse(results_df$p_value_bonferroni < 0.001, "***",
                                     ifelse(results_df$p_value_bonferroni < 0.01, "**",
                                           ifelse(results_df$p_value_bonferroni < 0.05, "*", "ns")))

results_df$stars_fdr <- ifelse(results_df$p_value_fdr < 0.001, "***",
                              ifelse(results_df$p_value_fdr < 0.01, "**",
                                    ifelse(results_df$p_value_fdr < 0.05, "*", "ns")))

# Print results with all p-values
print(data.frame(
  Layer = results_df$layer,
  Estimate = round(results_df$estimate, 3),
  T_value = round(results_df$t_value, 3),
  P_value_uncorrected = round(results_df$p_value, 4),
  P_value_bonferroni = round(results_df$p_value_bonferroni, 4),
  P_value_FDR = round(results_df$p_value_fdr, 4)
))