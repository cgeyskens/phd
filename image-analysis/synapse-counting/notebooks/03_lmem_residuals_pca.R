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

protein = "GPR37L1"

#### =============================== get input =============================== ####
if (protein == "VCAM1"){
    data_vglut1 <- read.csv("/mnt/image-analysis/synapse-counting/IHC_Exp9_VCAM1_VGLUT1-PSD95_output_data_20251201_194547/metric_results.csv")
    data_vgat <- read.csv("/mnt/image-analysis/synapse-counting/IHC_Exp12_VCAM1_VGAT-GEPH_output_data_20251201_194539/metric_results.csv")
} else if (protein == "GPR37L1") {
    data_vglut1 <- read.csv("/mnt/image-analysis/synapse-counting/IHC_Exp10_GPR37L1_VGLUT1-PSD95_output_data_20251201_194537/metric_results.csv")
    data_vgat <- read.csv("/mnt/image-analysis/synapse-counting/IHC_Exp11_GPR37L1_VGAT-GEPH_output_data_20251201_194548/metric_results.csv")
}

#### =============================== data processing =============================== ####

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

process_data <- function (data, synapse_type){
  data_formatted <- data %>%
    mutate( # add column sample_ID
      sample_ID = str_c(gRNA, Brain, sep = "_")
    ) %>% 
    group_by( # first log2 on section level, then mean per brain
      sample_ID, 
      hippocampal_layer
    ) %>% 
    summarise(
      across(all_of(metrics), ~ mean(log2(.x + 1), na.rm = FALSE)),
      .groups = "drop"
    ) %>%
    pivot_longer( # into long format
      cols = all_of(metrics),
      names_to = "metric",
      values_to = "value"
    ) %>%
    mutate( # adding a new column (synapse_type + hippocampal_layer + metric)
      synapse_hippocampal_layer_metric =
        str_c(synapse_type, hippocampal_layer, metric, sep = "_")
    ) %>%
    select( # pivoting back into wide format
      sample_ID, 
      synapse_hippocampal_layer_metric, 
      value) %>% 
    pivot_wider(
      names_from = synapse_hippocampal_layer_metric,
      values_from = value
    ) %>%
    rename_with( # fixing column names (no spaces)
      ~ gsub("[ -]+", "_", .x)
    ) 

    return(data_formatted)
}

# applying the function to VGLUT1 & VGAT data
vglut1_data_format <- process_data(data = data_vglut1, synapse_type = "VGLUT1-PSD95")
vgat_data_format <- process_data(data = data_vgat, synapse_type = "VGAT-GEPH")

# merging the two data frames
df_merged <- merge(vglut1_data_format, vgat_data_format, by = "sample_ID")


#### =============================== Linear Mixed Effects Modelling for PCA =============================== ####

# get all the metrics for looping
column_names <- as.list(colnames(df_merged))
metrics <- column_names[column_names != "sample_ID"]

# add the gRNA and Brain column to the original data
df_merged_with_Brain_gRNA <- df_merged %>%
    mutate(Brain = word(sample_ID, 2, sep = "_")) %>%
    mutate(gRNA = word(sample_ID, 1, sep = "_"))

# loop through the metrics, subset data, model it and get the residuals
residuals_list <- list()

for (metric_assessed in metrics) {
    data_metric <- df_merged_with_Brain_gRNA %>%
        select(Brain, gRNA, metric_assessed)

    # ensure the gRNA is a factor
    data_metric$gRNA <- factor(
      data_metric$gRNA, 
      levels = c("LacZ-gRNA", paste0(protein, "-gRNA"))
    )


    # model it
    model_for_pca <- lmer(
      formula = as.formula(paste0(metric_assessed, " ~ 1 + (1 | Brain)")),
      data = data_metric,
      control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000))
    )
    # get the residuals
    residuals_list[[metric_assessed]] <- residuals(model_for_pca)
}

residuals_for_pca <- as.data.frame(residuals_list) %>%
  tibble::as_tibble() %>%
  rename_with(~ paste0("residual_", .x)) %>%
  bind_cols(df_merged_with_Brain_gRNA %>% 
  select(sample_ID)) %>%
  column_to_rownames(., var = "sample_ID")

#### ============================= saving the residuals =============================== ####

output_directory <- "/mnt/image-analysis/synapse-counting/results_202512"
output_filepath <- file.path(output_directory, paste0(protein, "_lmem_residuals_results.csv"))
write.csv(residuals_for_pca, file = output_filepath, row.names = TRUE)


#### =============================== PCA on residuals =============================== ####

#### perform PCA
pca_results <- prcomp(residuals_for_pca, scale. = TRUE)
summary(pca_results) 

pca_scores <- as.data.frame(pca_results$x)
pca_variance <- pca_results$sdev^2
pca_percentage_variance <- round(100 * pca_variance / sum(pca_variance), 1)

# prepare data for pca plot
plot_data_pca <- as.data.frame(pca_results$x) %>%
  bind_cols(df_merged_with_Brain_gRNA %>% select(Brain, gRNA))

# scree plot
plot(pca_results, type = "l", main = "Scree Plot")

#### plot the PCA plot
# color setting
if (protein == "VCAM1") {
  gRNA_color <- "#21a0e2"   
} else if (protein == "GPR37L1") {
  gRNA_color <- "#e28d21"   
} else {
  warning("You didn’t set arguments correctly")
  gRNA_color <- "#000000"   
}
condition_colors <- c(
  setNames(gRNA_color, paste0(protein, "-gRNA")),
  "LacZ-gRNA" = "#c6c6c6"
)

# create the PCA plot
p <- ggplot(plot_data_pca, aes(x = PC1, y = PC2, color = gRNA, label = Brain)) +
            geom_point(size = 4) +
            #geom_text(hjust = -0.1, vjust = 0.1, size = 5) + # Add sample labels
            scale_color_manual(values = condition_colors) +
            labs(
                title = "PCA Plot of Hemispheres",
                x = paste0("PC1 (", pca_percentage_variance[1], "% variance explained)"),
                y = paste0("PC2 (", pca_percentage_variance[2], "% variance explained)")
            ) +
            theme_minimal() +
            theme(
                panel.background = element_blank(),
                panel.grid = element_blank(),
                axis.line = element_line(size = 1.2, color = "black"),   
                axis.ticks = element_line(size = 1.2, color = "black"), 
                axis.ticks.length = unit(0.2, "cm"),                     
                axis.text.x = element_text(angle = 45, hjust = 1, size = 16), 
                axis.text.y = element_text(size = 16),
                axis.title.y = element_text(size = 18, family = "Arial"),
                axis.title.x = element_text(size = 18, family = "Arial"), 
                plot.title = element_text(size = 20, family = "Arial"),
                #legend.position = "right",
                text = element_text(family = "Arial") 
            ) +
            scale_x_continuous(
                expand = c(0, 0),
                limits = c(-12, 12),
                breaks = seq(-10, 10, by = 5)  
            ) +
            scale_y_continuous(
                expand = c(0, 0),
                limits = c(-12, 12),
                breaks = seq(-10, 10, by = 5) 
            )
p

# saving
ggsave(paste0(protein, "_pca_paper.svg"), 
    plot = p, 
    device = cairo_pdf,
    width = 23, height = 20, units = "cm", dpi=300)
