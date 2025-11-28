library(dplyr)
library(ggplot2)
library(lme4)
library(stringr)
library(emmeans)
library(tibble)
library(tidyr)
library(nlme)
library(stringr)

remotes::install_github("nx10/httpgd")
library(httpgd)
hgd() # open the server for plotting

metric_assessed = "VGLUT1.PSD95_CA1.SLM_local_peak_colocalized_spots"  

# load in already log2 transformed mean data from PCA.ipynb
vcam1_data <- read.csv("/mnt/image-analysis/synapse-counting/vcam1_vglut1_vgat_wide_mean_log2_data.csv")
dim(vcam1_data)

colnames(vcam1_data)

# filter on the feature and add brain & gRNA
data_metric <- vcam1_data %>%
    mutate(Brain = word(sample_ID, 2, sep = "_")) %>%
    mutate(gRNA = word(sample_ID, 1, sep = "_")) %>%
    select(Brain, gRNA, metric_assessed)

# ensure the gRNA is a factor
data_metric$gRNA <- factor(data_metric$gRNA, levels = c("LacZ-gRNA", "VCAM1-gRNA"))

# making the model for the PCA without the gRNA fixed effect but taking out the Brain effect
model_for_pca <- lme(
  fixed = as.formula(paste0(metric_assessed, " ~ 1")),
  data = data_metric,
  random = ~1 | Brain,
  control = lmeControl(maxIter = 100, msMaxIter = 100, opt = "optim", singular.ok = TRUE)
)

summary(model_for_pca)
plot(model_for_pca)
qqnorm(resid(model_for_pca))
qqline(resid(model_for_pca))

# get the residuals for each hemisphere
# residual = observed value - (overall mean + brain-specific deviation)
print(residuals(model_for_pca))



############## Function to get all the residuals for each metric

# loading in the data
vcam1_data <- read.csv("/mnt/image-analysis/synapse-counting/vcam1_vglut1_vgat_wide_mean_log2_data.csv")

# get all the metrics for looping
column_names <- as.list(colnames(vcam1_data))
metrics <- column_names[column_names != "sample_ID"]

# add the gRNA and Brain column to the original data
vcam1_data_with_Brain_gRNA <- vcam1_data %>%
    mutate(Brain = word(sample_ID, 2, sep = "_")) %>%
    mutate(gRNA = word(sample_ID, 1, sep = "_"))

# loop through the metrics, subset data, model it and get the residuals
residuals_list <- list()

for (metric_assessed in metrics) {
    # subset data
    data_metric <- vcam1_data_with_Brain_gRNA %>%
        select(Brain, gRNA, metric_assessed)

    # ensure the gRNA is a factor
    data_metric$gRNA <- factor(data_metric$gRNA, levels = c("LacZ-gRNA", "VCAM1-gRNA"))

    # model it
    model_for_pca <- lme(
        fixed = as.formula(paste0(metric_assessed, " ~ 1")),
        data = data_metric,
        random = ~1 | Brain,
        control = lmeControl(maxIter = 100, msMaxIter = 100, opt = "optim", singular.ok = TRUE)
        )
    residuals_list[[metric_assessed]] <- residuals(model_for_pca)
}

residuals_for_pca <- as.data.frame(residuals_list)

prefix <- "residual_"

residuals_for_pca <- as.data.frame(residuals_list) %>%
  rename_with(~ paste0(prefix, .x)) %>%
  bind_cols(vcam1_data_with_Brain_gRNA %>% select(sample_ID)) %>%
  column_to_rownames(., var = "sample_ID")


#### Residuals for PCA
# saving 
output_directory <- "/mnt/image-analysis/synapse-counting/VCAM1"
output_filepath <- file.path(output_directory, "vcam1_lmem_residuals_results.csv")
write.csv(residuals_for_pca, file = output_filepath, row.names = TRUE)


#### Do the PCA
pca_results <- prcomp(residuals_for_pca, scale. = TRUE)
summary(pca_results) 

# Prepare data for ggplot2
plot_data_pca <- as.data.frame(pca_results$x) %>%
  # Add back the original 'Brain' and 'gRNA' columns from the processed data
  # to use for coloring/shaping in the plot
  bind_cols(vcam1_data_with_Brain_gRNA %>% select(Brain, gRNA))

# Plotting the Scree plot
plot(pca_results, type = "l", main = "Scree Plot of LMM Residuals (Brain Effect Removed)")

# Plot the PCA plot
condition_colors <- c("ip" = "#21a0e2", "igg" = "#c6c6c6")

# 6. Create the PCA scatter plot
ggplot(plot_data_pca, aes(x = PC1, y = PC2, color = gRNA, label = Brain)) +
  geom_point(size = 4, alpha = 0.8) +
  #geom_text(hjust = -0.1, vjust = 0.1, size = 5) + # Add sample labels
  # scale_color_manual(values = condition_colors) +
  labs(
    title = "PCA Plot of Samples",
    x = paste0("PC1 (", round(summary(pca_results)$importance[2,1]*100, 2), "%)"),
    y = paste0("PC2 (", round(summary(pca_results)$importance[2,2]*100, 2), "%)")
  ) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    axis.line = element_line(linewidth = 0.75, color = "black"),
    legend.position = "right",
    axis.ticks = element_line(color = "black", linewidth = 0.75, size = 1),
    axis.text.y = element_text(size = 16),
    axis.title.y = element_text(size = 18, family = "Arial"), 
    axis.title.x = element_text(size = 18, family = "Arial"), 
    plot.title = element_text(size = 20, family = "Arial"), 
  ) + 
  xlim(-10, 10) + 
  ylim(-20, 20)
