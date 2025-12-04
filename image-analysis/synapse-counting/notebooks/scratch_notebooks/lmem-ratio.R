library(dplyr)
library(ggplot2)
library(lme4)
library(stringr)
library(emmeans)
library(tibble)
library(tidyr)
library(nlme)

# remotes::install_github("nx10/httpgd")
library(httpgd)
hgd()

####################### For VCAM1 ###
# For VLGUT1-PSD95
metrics_to_assess = c("local_peak_colocalized_spots", # model_1 works, CA3 and DG layers are sign. model_3 works, CA3 SO & SR are borderline sign, DG Hilus is sign. model_5 works, DG Hilus sign, CA3-SL CA3-SO CA3 CA3-SR <0.1
                      "overlap_coeff", # model_1: works. model_3 works, model_5 works.
                      "overlap_um2", # model_1: boundary (singular) fit: see help('isSingular'), model_2 works. model_3 works. DG Hilus borderline sign. model_5 works.
                      "pearson_cor", # model_1: works. model_3 works. model_5 works.
                      "presynapse_image_mfi", # model_1: works, CA3 and DG layers are sign. model_3 works, CA3 layers is sign. model_5 works, CA3 layers are sign, DG-Hilus <0.1
                      "postsynapse_image_mfi", # model_1: works, CA3 SL is sign. model_3 works. model_5 works.
                      "pre_puncta_density_per_100_um2", # model_1: boundary (singular) fit: see help('isSingular'), model_2 works. model_3 works. model_5 works.
                      "post_puncta_density_per_100_um2", # model_1: works, model_3 works. model_5 works.
                      "pre_staining_area_um2", # model_1: works, model_3 works. model_5 works
                      "post_staining_area_um2", # model_1: works, model_3 works. model_5 works.
                      "pre_mean_puncta_size_um2", # model_1: works, model_3 works. model_5 works.
                      "post_mean_puncta_size_um2") # model_1: works, model_3 works. model_5 works.
# For VGAT-GEPH
metrics_to_assess = c("local_peak_colocalized_spots", # model_1: boundary (singular) fit: see help('isSingular'), model_2 works but not quite intuitive. model_3 works, CA3 SL & SR sign. model_5 works, CA3-SL CA3-SR sign. CA3 SO <0.1
                      "overlap_coeff", # model_1: boundary (singular) fit: see help('isSingular'), model_2 works. model_3 works. model_5 works.
                      "overlap_um2", # model_1: boundary (singular) fit: see help('isSingular'), model_2 works, model_3 doesnt work. model_4 works. model_5 works.
                      "pearson_cor", # model_1: works. model_3 works. model_5 works.
                      "presynapse_image_mfi", # model_1: works, CA1 SR, CA3 SL, CA3 SO and CA3 SR are sign. model_3 works. model_5 works.
                      "postsynapse_image_mfi", # model_1: boundary (singular) fit: see help('isSingular'), CA3 SR is sign, model_2 works, CA3 SR is sign. model_3 works. model_5 works, CA1-SLM CA1-SR and CA3 SR are <0.1.
                      "pre_puncta_density_per_100_um2", # model_1: works, CA1 SO is sign. model_3 works. model_5 works, CA3-SR <0.1
                      "post_puncta_density_per_100_um2", # model_1: boundary (singular) fit: see help('isSingular'), model_2 works. model_3 does not work. model_5 works.
                      "pre_staining_area_um2", # model_1: works. model_3 works. model_5 works, CA3 SO <0.1
                      "post_staining_area_um2", # model_1: boundary (singular) fit: see help('isSingular'), model_2 works. model_3 does not work. model_5 works.
                      "pre_mean_puncta_size_um2", # model_1: works. model_3 works. model_5 works.
                      "post_mean_puncta_size_um2") # model_1: works. model_3 works. model_5 works.

##################### For GPR37L1
# For VGLUT1-PSD95
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

# FOR VGAT-GEPH
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

protein = "VCAM1"
metric_assessed = "local_peak_colocalized_spots"

# loading the data
vcam1_data <- read.csv("/mnt/image-analysis/synapse-counting/IHC_Exp9_VCAM1_VGLUT1-PSD95_output_data_20251128_135431/metric_results.csv")

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


# OPTION 1: this is the better option

# calculating first the log2 ratio on section level & then taking the mean per brain
log2_ratio_section <- vcam1_data_metric_assessed %>%
  pivot_wider(
    names_from = gRNA,
    values_from = metric_assessed
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
  fixed = as.formula(paste0(metric_assessed, "_log2_ratio ~ hippocampal_layer")),
  data = log2_ratio_section,
  random = ~1 | Brain,
  weights = varIdent(form = ~1 | hippocampal_layer)
)

# more iterations then model_3
model_4 <- lme(
  fixed = as.formula(paste0(metric_assessed, "_log2_ratio ~ hippocampal_layer")),
  data = log2_ratio_section,
  random = ~1 | Brain,
  weights = varIdent(form = ~1 | hippocampal_layer),
  control = lmeControl(maxIter = 100, msMaxIter = 100, opt = "nlminb", singular.ok = TRUE)
)

# more iterations and optim method
model_5 <- lme(
  fixed = as.formula(paste0(metric_assessed, "_log2_ratio ~ hippocampal_layer")),
  data = log2_ratio_section,
  random = ~1 | Brain,
  weights = varIdent(form = ~1 | hippocampal_layer),
  control = lmeControl(maxIter = 100, msMaxIter = 100, opt = "optim", singular.ok = TRUE)
)

# checking the model
summary(model_5)
plot(model_5)
qqnorm(resid(model_5))
qqline(resid(model_5))

# extracting the estimated mean per layer (without reference)
em_means <- emmeans(model_5, ~hippocampal_layer)
summary(em_means)

# producing the output table with reference layer
test_output <- test(em_means)
results_df <- as.data.frame(test_output)
results_df$p.value.adj = p.adjust(results_df$p.value, method = "fdr")









################################## Functions ################################
# loading the data
vcam1_data <- read.csv("/mnt/image-analysis/synapse-counting/IHC_Exp12_VCAM1_VGAT-GEPH_output_data_20251201_194539/metric_results.csv")
dim(vcam1_data)

# protein
protein = "VCAM1"

# synapse type
synapse_type = "VGAT-GEPH"

# metrics
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

# preprocess the data for current metric
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

# take the log2 ratios at the section level
data_log2 <- function(data_to_format, metric, protein){
    protein_gRNA_col_name <- paste0(protein, "-gRNA")
    lacZ_gRNA_col_name <- "LacZ-gRNA"
    
    # add pseudocount in case of zero's
    pseudocount <- 1e-6

    log2_ratio_section <- data_to_format %>%
        mutate(
          !!sym(metric) := !!sym(metric) + pseudocount
        ) %>%
        pivot_wider(
            names_from = gRNA,
            values_from = !!sym(metric)
        ) %>% 
        group_by(Brain, section, hippocampal_layer) %>%
        mutate(log2_ratio = log2(!!sym(protein_gRNA_col_name) / !!sym(lacZ_gRNA_col_name))) %>%
        group_by(Brain, hippocampal_layer) %>%
        mutate(!!sym(metric) := mean(log2_ratio)) %>%
        select(hippocampal_layer, Brain, !!sym(metric)) %>%
        distinct()
    
    new_metric_col_name = paste0(metric, "_log2_ratio")

    log2_ratio <- log2_ratio_section %>%
        rename(!!sym(new_metric_col_name) := !!sym(metric))

    return(log2_ratio)
}

# run the linear mixed effects model - nlme on ratio
modelling_1 <- function(data, metric){
    model <- nlme::lme(
        fixed = as.formula(paste0(metric, "_log2_ratio ~ hippocampal_layer")),
        data = data,
        random = ~1 | Brain,
        weights = varIdent(form = ~1 | hippocampal_layer),
        control = lmeControl(maxIter = 100, msMaxIter = 100, opt = "optim", singular.ok = TRUE)
    )
    return(model)
}

# extract model params
model_extract_params <- function(model_to_extract){
    em_means <- emmeans(model_to_extract, ~hippocampal_layer)
    results_df <- as.data.frame(test(em_means, null = 0, adjust = "fdr"))
    return(results_df)
}

# Function to generate diagnostic plots for a fitted nlme::lme model
check_diagnostics_nlme <- function(model, metric_name) {

    model_residuals <- residuals(model, type = "normalized")
    
    fitted_values <- fitted(model)

    opar <- par(mfrow = c(1, 2)) # Save original parameters

    qqnorm(model_residuals, 
           main = paste("Q-Q Plot:", metric_name))
    qqline(model_residuals, col = "red")
    
    plot(fitted_values, model_residuals,
         xlab = "Fitted Values",
         ylab = "Normalized Residuals",
         main = paste("Residuals vs. Fitted:", metric_name))
    abline(h = 0, col = "red", lty = 2)
    
    par(opar) 
}

# main function for full analysis
run_full_analysis <- function(data, metrics_to_assess, protein){
    all_results_list <- list() # empty list
    for (metric_assessed in metrics_to_assess){
        preprocessed_metric <- process_data_for_metric(data = data, metric = metric_assessed)
        preprocessed_log2 <- data_log2(data_to_format = preprocessed_metric, metric = metric_assessed, protein = protein)
        
        model <- modelling_1(data = preprocessed_log2, metric = metric_assessed)

        check_diagnostics_nlme(model = model, metric_name = metric_assessed)

        results <- model_extract_params(model_to_extract = model)
        results$analyzed_metric <- metric_assessed
        results$analyzed_protein <- protein
    
        all_results_list[[metric_assessed]] <- results
    }

    final_combined_df <- bind_rows(all_results_list)

    return (final_combined_df)
}

# running the pipeline
data_analyzed_nlme_ratio <- run_full_analysis(
    data = vcam1_data,
    metrics_to_assess = metrics_to_assess,
    protein = protein
)





# saving 
output_directory <- "/mnt/image-analysis/synapse-counting/results_202512"
output_filename <- paste0(protein, "_", synapse_type, "_lmem_analysis_results.csv")
output_filepath <- file.path(output_directory, output_filename)
write.csv(data_analyzed, file = output_filepath, row.names = FALSE)
