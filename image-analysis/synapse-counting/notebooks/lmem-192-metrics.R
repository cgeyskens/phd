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

# loading in the data
vcam1_data <- read.csv("/mnt/image-analysis/synapse-counting/vcam1_vglut1_vgat_wide_mean_log2_data.csv")

# add the gRNA and Brain column to the original data
vcam1_data_with_Brain_gRNA <- vcam1_data %>%
    mutate(Brain = word(sample_ID, 2, sep = "_")) %>%
    mutate(gRNA = word(sample_ID, 1, sep = "_"))

# ensure the gRNA is a factor
vcam1_data_with_Brain_gRNA$gRNA <- factor(vcam1_data_with_Brain_gRNA$gRNA, levels = c("LacZ-gRNA", "VCAM1-gRNA"))

# making the model for the PCA without the gRNA fixed effect but taking out the Brain effect
model_for_192_metrics <- lme(
  fixed = as.formula(paste0(metric_assessed, " ~ gRNA")),
  data = vcam1_data_with_Brain_gRNA,
  random = ~1 | Brain,
  control = lmeControl(maxIter = 100, msMaxIter = 100, opt = "optim", singular.ok = TRUE)
)

summary(model_for_192_metrics)

plot(model_for_192_metrics)
qqnorm(resid(model_for_192_metrics))
qqline(resid(model_for_192_metrics))

emm_results <- emmeans(model_for_192_metrics, ~ gRNA, contr = "pairwise")

contrasts_data <- data.frame(emm_results$contrasts)
means_data <- data.frame(emm_results$emmeans)

results_df <- data.frame(
  metric_assessed = metric_assessed,
  LacZ_mean = means_data$emmean[seq(0, nrow(means_data), 2)],  # Get LacZ means
  VCAM1_mean = means_data$emmean[seq(1, nrow(means_data), 2)], # Get VCAM1 means
  Estimate = contrasts_data$estimate,
  SE = contrasts_data$SE,
  t_ratio = contrasts_data$t.ratio,
  p_value = contrasts_data$p.value
  )

############################# Functions ###############################

# loading in the data
vcam1_data <- read.csv("/mnt/image-analysis/synapse-counting/vcam1_vglut1_vgat_wide_mean_log2_data.csv")

# subset dataxs
process_data_for_metric <- function(data, metric_assessed){
    data_metric <- data %>%
        select(Brain, gRNA, metric_assessed)
    return(data_metric)
}

# model function
modelling <- function(data, metric_assessed){
    model_for_each_192_metrics <- lme(
        fixed = as.formula(paste0(metric_assessed, " ~ gRNA")),
        data = data,
        random = ~1 | Brain,
        control = lmeControl(maxIter = 100, msMaxIter = 100, opt = "optim", singular.ok = TRUE)
        )
    return(model_for_each_192_metrics)
}

# extract model params
model_extract_params <- function(model_to_extract, metric_assessed){
    
    emm_results <- emmeans(model_to_extract, ~ gRNA, contr = "pairwise")
    
    contrasts_data <- data.frame(emm_results$contrasts)
    means_data <- data.frame(emm_results$emmeans)
    
    results_df <- data.frame(
        metric_assessed = metric_assessed,
        LacZ_mean = means_data$emmean[seq(0, nrow(means_data), 2)],  # Get LacZ means
        VCAM1_mean = means_data$emmean[seq(1, nrow(means_data), 2)], # Get VCAM1 means
        Estimate = contrasts_data$estimate,
        SE = contrasts_data$SE,
        t_ratio = contrasts_data$t.ratio,
        p_value = contrasts_data$p.value
        )
    return(results_df)
}

# main function for full analysis
run_full_analysis <- function(data){

    # get all the metrics for looping
    column_names <- as.list(colnames(data))
    metrics <- column_names[column_names != "sample_ID"]

    # add the gRNA and Brain column to the original data
    vcam1_data_with_Brain_gRNA <- data %>%
        mutate(Brain = word(sample_ID, 2, sep = "_")) %>%
        mutate(gRNA = word(sample_ID, 1, sep = "_"))

    # ensure the gRNA is a factor
    vcam1_data_with_Brain_gRNA$gRNA <- factor(vcam1_data_with_Brain_gRNA$gRNA, levels = c("LacZ-gRNA", "VCAM1-gRNA"))

    
    all_results_list <- list() # empty list
    
    for (metric_assessed in metrics){
        preprocessed_metric <- process_data_for_metric(data = vcam1_data_with_Brain_gRNA, metric_assessed =  metric_assessed)
        model <- modelling(data = preprocessed_metric, metric = metric_assessed)
        results <- model_extract_params(model_to_extract = model, metric_assessed = metric_assessed)
    
        all_results_list[[metric_assessed]] <- results
    }

    final_combined_df <- bind_rows(all_results_list)

    return (final_combined_df)
}

# running the pipeline
data_analyzed <- run_full_analysis(
    data = vcam1_data
)

# adjusting the pvalues
data_analyzed_adj_pval <- data_analyzed %>%
    mutate(adj_pval = p.adjust(p_value, method = "fdr")
    )
