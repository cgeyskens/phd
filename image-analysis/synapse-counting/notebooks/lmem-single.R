# test script for linear mixed effect models
# now just trying out a lme for a single synaptic metric

# loading the data
vcam1_data <- read.csv("/Volumes/KINGSTON/data/phd/image-analysis/synapse-counting/VCAM1/VCAM1-LacZ_VGLUT1-PSD95_output_data/metric_results.csv")

# getting the libraries
library(dplyr)
library(ggplot2)
library(lme4)
library(stringr)

# data preprocessing for only one dependent variable (one synaptic metric, here it is presynaptic MFI)
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

# standardizing the values
vcam1_data_vglut_mfi$presynapse_image_mfi_scaled <- scale(vcam1_data_vglut_mfi$presynapse_image_mfi)

# checking distribution
hist(vcam1_data_vglut_mfi$presynapse_image_mfi_scaled)

# modifying the gRNA and hippocampal_layer column from character to factor for one hot encoding
vcam1_data_vglut_mfi$gRNA <- as.factor(vcam1_data_vglut_mfi$gRNA) 
vcam1_data_vglut_mfi$hippocampal_layer <- as.factor(vcam1_data_vglut_mfi$hippocampal_layer) 
class(vcam1_data_vglut_mfi$gRNA)
class(vcam1_data_vglut_mfi$hippocampal_layer)

# checking all the values in a boxplot
boxplot(vcam1_data_vglut_mfi$presynapse_image_mfi_scaled ~ vcam1_data_vglut_mfi$gRNA, col="skyblue", xlab="gRNA", ylab="variable")

# trying to model with a simple linear regression
basic_lm <- lm(presynapse_image_mfi_scaled ~ gRNA, data = vcam1_data_vglut_mfi)
summary(basic_lm)

# nested mixed effects model
mixed_lm <- lmer(presynapse_image_mfi_scaled ~ gRNA + 
                (1|Brain/section/hippocampal_layer), 
                data = vcam1_data_vglut_mfi)
summary(mixed_lm)

# model: interaction between gRNA and layer
mixed_lm2 <- lmer(presynapse_image_mfi_scaled ~ gRNA * hippocampal_layer +
                  (1|Brain/section/hippocampal_layer),
                  data = vcam1_data_vglut_mfi)
summary(mixed_lm2)

# improved model
improved_model <- lmer(presynapse_image_mfi_scaled ~ gRNA * hippocampal_layer + 
                     (1|Brain) +                             # between-brain variation because of perfusion/viral injection
                     (1|Brain:gRNA) +                        # paired design of hemispheres
                     (1|Brain:section:hippocampal_layer),    # nested measurements
                     data = vcam1_data_vglut_mfi,
                     REML=TRUE)
summary(improved_model)

plot(improved_model)
qqnorm(resid(improved_model))
qqline(resid(improved_model))

anova(improved_model, mixed_lm2)


## Generalized models
hist(vcam1_data_vglut_mfi$presynapse_image_mfi)

# taking into accound the non-normality of the data WITH ORIGINAL DATA
glm_model <- glmer(presynapse_image_mfi ~ gRNA * hippocampal_layer + 
                     (1|Brain) +                             # between-brain variation
                     (1|Brain:gRNA) +                        # paired design of hemispheres
                     (1|Brain:section:hippocampal_layer),    # nested measurements
                     data = vcam1_data_vglut_mfi,
                     family = Gamma(link = "log"),
                     control = glmerControl(optimizer = "bobyqa",
                                            optCtrl = list(maxfun = 200000)))
summary(glm_model)

# QC
plot(glm_model)
qqnorm(resid(glm_model))
qqline(resid(glm_model))

anova(glm_model, mixed_lm2) # not comparable because input data is different

# Extract p-values and coefficient names
coef_table <- summary(improved_model)$coefficients
p_values <- coef_table[,"Pr(>|t|)"]
coef_names <- rownames(coef_table)

# Apply Holm-Bonferroni correction
p_adjusted_fdr <- p.adjust(p_values, method = "fdr")

# Create a comparison table
results_table <- data.frame(
  Coefficient = coef_names,
  Original_P = p_values,
  Adjusted_P = p_adjusted_fdr,
  Significant = p_adjusted_fdr < 0.05
)
