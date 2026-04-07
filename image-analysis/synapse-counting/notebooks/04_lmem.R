# ============================================================================
# Title: Testing for differences in synaptic metrics between LacZ-gRNA- and VCAM1-gRNA- or GPR37L1-gRNA-injected hemispheres
# Input: Output of the pipeline with the all metrics analyzed
# Description: 
# (1) data processing
# (2) linear mixed effects modelling to extract residuals and account for brain-brain variability
# (3) PCA on residuals to visualize differences between hemispheres based synaptic metrics
# Date: Nov 2025
# ============================================================================


library(dplyr)
library(ggplot2)
library(lme4)
library(stringr)
library(emmeans)
library(tibble)
library(tidyr)

library(httpgd)
hgd()


#### -------------------------------- arguments --------------------------------- #####

# experiment type (crispr or oe)
experiment = "oe"
# protein
protein = "VCAM1"
# synapse type
synapse_type = "VGAT-GEPH"

# loading the data [for each protein (VCAM1 vs GPR37L1) and synapse type (VGLUT1-PSD95 vs VGAT-GEPH)]
data <- read.csv("/mnt/image-analysis/synapse-counting-oe/results_202604/VGAT-GEPH_output_data_20260403_171952/metric_results.csv")
dim(data)

# metrics to assess
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


#### -------------------------------- functions --------------------------------- #####

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

log2_metric <- function(data, metric){
    metric_sym <- rlang::sym(metric)
    log2_col_name <- paste0(metric, "_log2")

    data_log2 <- data %>%
        mutate(
            !!log2_col_name := log2(!!metric_sym + 1),
        ) 
    return (data_log2)
}

model_1 <- function(data, metric){
    model <- lme4::lmer(
      formula = as.formula(paste0(metric, "_log2 ~ gRNA * hippocampal_layer + 
                                  (1|Brain)"
                                  )
                            ),
      data = data,
      REML = TRUE)
}

posthoc_analysis <- function(model){
    emm_results <- emmeans(
        model, pairwise ~ gRNA | hippocampal_layer
    )
    contrasts_data <- as.data.frame(emm_results$contrasts)
    means_data <- as.data.frame(emm_results$emmeans)

    means_wide <- means_data %>%
        select(gRNA, hippocampal_layer, emmean) %>%
        pivot_wider(
            names_from = gRNA, 
            values_from = emmean,
            names_prefix = "Mean_"
        )

    results_df <- merge(
        contrasts_data,
        means_wide,
        by = c("hippocampal_layer")
    )
    results_df$p.value.adj <- p.adjust(results_df$p.value, method = "fdr")

    return(results_df)
}

check_diagnostics <- function(model, metric_name, model_type) {
    
    model_residuals <- residuals(model)
    par(mfrow = c(1, 2))
    
    qqnorm(model_residuals, 
           main = paste("Q-Q Plot for", metric_name, " (", model_type, ")"))
    qqline(model_residuals, col = "red")
    
    plot(fitted(model), model_residuals,
         xlab = "Fitted Values",
         ylab = "Residuals",
         main = paste("Residuals vs. Fitted (", model_type, ")"))
    abline(h = 0, col = "red")
    
    par(mfrow = c(1, 1))
}

run_analysis <- function(data, metrics, protein){
    all_results_list <- list() 
    for (metric in metrics){
        # preprocessing
        preprocessed_metric <- process_data_for_metric(data = data, metric = metric)
        preprocessed_metric_log2 <- log2_metric(data = preprocessed_metric, metric = metric)

        # modelling
        model <- model_1(data = preprocessed_metric_log2, metric = metric)
        
        # diagnostics check
        check_diagnostics(
            model = model, 
            metric_name = metric, 
            model_type = "LMEM using lme4"
        )

        # posthoc analysis
        results <- posthoc_analysis(model=model)

        results$analyzed_metric <- metric
        results$analyzed_protein <- protein
    
        all_results_list[[metric]] <- results
    }

    final_combined_df <- bind_rows(all_results_list)

    return (final_combined_df)
}


#### -------------------------------- actual analysis --------------------------------- ####
data_analyzed <- run_analysis(
    data = data,
    metrics = metrics_to_assess,
    protein = protein
)


#### -------------------------------- saving --------------------------------- ####
if (experiment == "crispr"){
    output_directory <- "/mnt/image-analysis/synapse-counting/results_20251203"
    output_filename <- paste0(protein, "_", synapse_type, "_crispr_lmem_analysis_results.csv")
    output_filepath <- file.path(output_directory, output_filename)
} else if (experiment == "oe") {
    output_directory <- "/mnt/image-analysis/synapse-counting-oe/results_202604"
    output_filename <- paste0(protein, "_", synapse_type, "_oe_lmem_analysis_results.csv")
    output_filepath <- file.path(output_directory, output_filename)
}
write.csv(data_analyzed, file = output_filepath, row.names = FALSE)
