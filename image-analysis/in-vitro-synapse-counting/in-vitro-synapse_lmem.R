# ============================================================================
# Title: Analysis of in vitro synapse data
# Description: First prep the data, then normalizes & log2 transform, then
# linear mixed model fitting and visualizing for puncta density & size.
# Author: Cydric Geyskens
# Date: Oct 2025
# ============================================================================

#### ====================== loading packages ======================== ####


# (Dev Container from synapse-counting)
library(tidyr)
library(tibble)
library(dplyr)
library(stringr)
library(emmeans)
library(stringr)
library(lme4)
library(RColorBrewer)
library(ggplot2)

# install.packages("readxl")
library(readxl)

# install.packages("ggbeeswarm")
library(ggbeeswarm)

# remotes::install_github("nx10/httpgd")
library(httpgd)
hgd() # open the server for plotting


#### ====================== Data prep & cleaning ======================== ####

# loading in the data
pre_post_in_vitro <- read_excel(
  "/mnt/image-analysis/in-vitro-synapse-counting/all_data/exp17_20_23_24_data.xlsx", "raw_data"
)
dim(pre_post_in_vitro)

# filter out FOVs that were not included in the analysis (had NA values)
data_clean <- pre_post_in_vitro[!is.na(pre_post_in_vitro$raw_synapse_puncta_nr), ]
dim(data_clean)

# ensure that treatment is a factor for the downstream model
data_clean$treatment <- factor(data_clean$treatment , levels = c("Fc", "VCAM1-Fc"))


#### ============== Calculate the synapse density =================== ####

data_clean <- data_clean %>%
  mutate(
    raw_synapse_count_per_100um = (raw_synapse_puncta_nr / total_dendritic_length_um)*100,  
    raw_presynapse_count_per_100um = (raw_presynapse_puncta_nr / total_dendritic_length_um)*100,
    raw_postsynapse_count_per_100um = (raw_postsynapse_puncta_nr / total_dendritic_length_um)*100
  )

#### ================= Normalization & Log2 Transformation =========== ####

# get the Fc control means for each metric across the experiments
fc_means <- data_clean %>%
  filter(treatment == "Fc") %>%
  summarize(
    across(
      c(raw_synapse_count_per_100um,
        raw_presynapse_count_per_100um,
        raw_postsynapse_count_per_100um,
        raw_synapse_puncta_size,
        raw_presynapse_puncta_size,
        raw_postsynapse_puncta_size
      ),
      mean,
      .names = "global_fc_{.col}_mean"
    )
  )

# make new columns with the normalized & log2 transformed metrics
data_clean <- data_clean %>%
  mutate(
    norm_synapse_count_per_100um = raw_synapse_count_per_100um / fc_means$global_fc_raw_synapse_count_per_100um_mean,
    norm_presynapse_count_per_100um = raw_presynapse_count_per_100um / fc_means$global_fc_raw_presynapse_count_per_100um_mean,
    norm_postsynapse_count_per_100um = raw_postsynapse_count_per_100um / fc_means$global_fc_raw_postsynapse_count_per_100um_mean,
    norm_synapse_puncta_size = raw_synapse_puncta_size / fc_means$global_fc_raw_synapse_puncta_size_mean,
    norm_presynapse_puncta_size = raw_presynapse_puncta_size / fc_means$global_fc_raw_presynapse_puncta_size_mean,
    norm_postsynapse_puncta_size = raw_postsynapse_puncta_size / fc_means$global_fc_raw_postsynapse_puncta_size_mean
  ) %>%
  mutate(across(c(norm_synapse_count_per_100um,
                  norm_presynapse_count_per_100um,
                  norm_postsynapse_count_per_100um,
                  norm_synapse_puncta_size,
                  norm_presynapse_puncta_size,
                  norm_postsynapse_puncta_size
                ),
                ~ log2(1 + .), # to handle zero values
                .names = "{.col}_log2"
              )
  )


#### =========== Linear Mixed Model fitting, puncta density ========== #####

synapse_count_model <- lmer(norm_synapse_count_per_100um_log2 ~ treatment + (1 | experiment), data = data_clean)
summary(synapse_count_model)

plot(synapse_count_model)
qqnorm(resid(synapse_count_model))
qqline(resid(synapse_count_model))

synapse_count_emm <- emmeans(synapse_count_model, ~ treatment)
summary(synapse_count_emm)

synapse_count_pairwise_comparisons <- contrast(synapse_count_emm, method = "pairwise")
summary(synapse_count_pairwise_comparisons)

# presynapse count model
presynapse_count_model <- lmer(norm_presynapse_count_per_100um_log2 ~ treatment + (1 | experiment), data = data_clean)
summary(presynapse_count_model)

plot(presynapse_count_model)
qqnorm(resid(presynapse_count_model))
qqline(resid(presynapse_count_model))

presynapse_count_emm <- emmeans(presynapse_count_model, ~ treatment)
summary(presynapse_count_emm)

presynapse_count_pairwise_comparisons <- contrast(presynapse_count_emm, method = "pairwise")
summary(presynapse_count_pairwise_comparisons)

# postsynapse count model
postsynapse_count_model <- lmer(norm_postsynapse_count_per_100um_log2 ~ treatment + (1 | experiment), data = data_clean)
summary(postsynapse_count_model)

plot(postsynapse_count_model)
qqnorm(resid(postsynapse_count_model))
qqline(resid(postsynapse_count_model))

postsynapse_count_emm <- emmeans(postsynapse_count_model, ~ treatment)
summary(postsynapse_count_emm)

postsynapse_count_pairwise_comparisons <- contrast(postsynapse_count_emm, method = "pairwise")
summary(postsynapse_count_pairwise_comparisons)


#### =============== Data Visualization, puncta density =============== ####

y_value_to_visualize = "norm_postsynapse_count_per_100um_log2"

# calculate the mean of each experiment to also plot the means
experiment_means <- data_clean %>%
  group_by(experiment, treatment) %>%
  summarize(mean_norm_synapse_count_per_100um_log2 = mean(norm_synapse_count_per_100um_log2),
            mean_norm_presynapse_count_per_100um_log2 = mean(norm_presynapse_count_per_100um_log2),
            mean_norm_postsynapse_count_per_100um_log2 = mean(norm_postsynapse_count_per_100um_log2)
  )

# create a column combining experiment and treatment for the coloring of the points
data_clean <- data_clean %>%
  mutate(exp_treat = paste0(experiment, "_", treatment))
experiment_means <- experiment_means %>%
  mutate(exp_treat = paste0(experiment, "_", treatment))

# get unique experiments and treatments
experiments <- sort(unique(as.character(data_clean$experiment)))
n_exp <- length(experiments)

# create color mapping: greys for Fc, PuBu for VCAM1-Fc
experiment_colors <- c()

for (i in seq_along(experiments)) {
  exp <- experiments[i]
  # Fc: Greys palette 
  grey_colors <- brewer.pal(9, "Greys")[5:8]
  experiment_colors[paste0(exp, "_Fc")] <- grey_colors[i]
  # VCAM1-Fc: use PuBu palette
  pu_colors <- brewer.pal(9, "PuBu")[5:8]
  experiment_colors[paste0(exp, "_VCAM1-Fc")] <- pu_colors[i]
}

# actual plotting
p <- ggplot(data_clean, aes(x = treatment, 
                            y = !!sym(y_value_to_visualize), 
                            color = exp_treat)) +
    geom_beeswarm(cex = 4, size = 4, alpha = 0.8, priority = "ascending") + 
    geom_point(data = experiment_means, 
              aes(x = as.numeric(as.factor(treatment)) + 0.45, 
              y = !!sym(paste0("mean_", y_value_to_visualize)), 
              color = exp_treat),
              size = 7, shape = 19, alpha = 0.6, show.legend = FALSE) +  
    labs(title = "FOVs by Experiment and Treatment with Mean Synapse Counts",
        y = paste("Mean", gsub("_", " ", y_value_to_visualize)),
        color = "Experiment") +
    scale_color_manual(values = experiment_colors) +
    theme(
      legend.position = "top",  
      panel.background = element_blank(),
      panel.grid = element_blank(), 
      axis.line = element_line(color = "black", size = 1),  
      axis.ticks.x = element_blank(),
      axis.ticks.y = element_line(color = "black", size = 1),
      axis.ticks.length.y = unit(.25, "cm"),
      axis.text = element_text(size = 14), 
      axis.title.x = element_blank(),
      axis.title.y = element_text(size = 16)
    ) + 
    ylim(0, 3) + 
    coord_fixed(ratio = 1.8) 

p



# visualize the data
# p <- ggplot(data_clean, aes(x = treatment, 
#                             y = !!sym(y_value_to_visualize), 
#                             color = as.factor(experiment))) +
#     geom_beeswarm(cex = 4, size = 3, alpha = 0.8, priority = "ascending") + 
#     geom_point(data = experiment_means, 
#               aes(x = as.numeric(treatment) + 0.45 , 
#               y = !!sym(paste0("mean_", y_value_to_visualize)), 
#               color = as.factor(experiment)),
#               size = 5, shape = 19, alpha = 0.6, show.legend = FALSE) +  
#     labs(title = "FOVs by Experiment and Treatment with Mean Synapse Counts",
#         y = paste("Mean", gsub("_", " ", y_value_to_visualize)),
#         color = "Experiment") +
#     scale_color_brewer(palette = "Dark2") +
#     theme(
#       legend.position = "top",  
#       panel.background = element_blank(),
#       panel.grid = element_blank(), 
#       axis.line = element_line(color = "black", size = 1),  
#       axis.ticks.x = element_blank(),
#       axis.ticks.y = element_line(color = "black", size = 1),
#       axis.ticks.length.y = unit(.25, "cm"),
#       axis.text = element_text(size = 14), 
#       axis.title.x = element_blank(),
#       axis.title.y = element_text(size = 16)
#     ) + 
#     ylim(0, 3) + 
#     coord_fixed(ratio = 3) 
# p

# ggsave("presynapse_inh_norm.png", plot = p,
#        width = 1500, height = 1800, units = "px", bg = "white", dpi = 300)



#### ================== Model fitting, puncta size ======================= ####

# synapse size model
synapse_size_model <- lmer(norm_synapse_puncta_size_log2 ~ treatment + (1 | experiment), data = data_clean)
summary(synapse_size_model)

plot(synapse_size_model)
qqnorm(resid(synapse_size_model))
qqline(resid(synapse_size_model))

synapse_size_emm <- emmeans(synapse_size_model, ~ treatment)
summary(synapse_size_emm)

synapse_size_pairwise_comparisons <- contrast(synapse_size_emm, method = "pairwise")
summary(synapse_size_pairwise_comparisons)

# presynapse size model
presynapse_size_model <- lmer(norm_presynapse_puncta_size_log2 ~ treatment + (1 | experiment), data = data_clean)
summary(presynapse_size_model)

plot(presynapse_size_model)
qqnorm(resid(presynapse_size_model))
qqline(resid(presynapse_size_model))

presynapse_size_emm <- emmeans(presynapse_size_model, ~ treatment)
summary(synapse_size_emm)

presynapse_size_pairwise_comparisons <- contrast(presynapse_size_emm, method = "pairwise")
summary(presynapse_size_pairwise_comparisons)

# postsynapse size model
postsynapse_size_model <- lmer(norm_postsynapse_puncta_size_log2 ~ treatment + (1 | experiment), data = data_clean)
summary(postsynapse_size_model)

plot(postsynapse_size_model)
qqnorm(resid(postsynapse_size_model))
qqline(resid(postsynapse_size_model))

postsynapse_size_emm <- emmeans(postsynapse_size_model, ~ treatment)
summary(synapse_size_emm)

postsynapse_size_pairwise_comparisons <- contrast(postsynapse_size_emm, method = "pairwise")
summary(postsynapse_size_pairwise_comparisons)


#### ================== Data Visualization, puncta size ================== ####

y_value_to_visualize = "norm_postsynapse_puncta_size_log2"

# calculate the mean of each experiment to also plot the means
experiment_means <- data_clean %>%
  group_by(experiment, treatment) %>%
  summarize(mean_norm_synapse_puncta_size_log2 = mean(norm_synapse_puncta_size_log2),
            mean_norm_presynapse_puncta_size_log2 = mean(norm_presynapse_puncta_size_log2),
            mean_norm_postsynapse_puncta_size_log2 = mean(norm_postsynapse_puncta_size_log2)
  )

# create a column combining experiment and treatment for the coloring of the points
data_clean <- data_clean %>%
  mutate(exp_treat = paste0(experiment, "_", treatment))
experiment_means <- experiment_means %>%
  mutate(exp_treat = paste0(experiment, "_", treatment))

# get unique experiments and treatments
experiments <- sort(unique(as.character(data_clean$experiment)))
n_exp <- length(experiments)

# create color mapping: greys for Fc, PuBu for VCAM1-Fc
experiment_colors <- c()

for (i in seq_along(experiments)) {
  exp <- experiments[i]
  # Fc: Greys palette 
  grey_colors <- brewer.pal(9, "Greys")[5:8]
  experiment_colors[paste0(exp, "_Fc")] <- grey_colors[i]
  # VCAM1-Fc: use PuBu palette
  pu_colors <- brewer.pal(9, "PuBu")[5:8]
  experiment_colors[paste0(exp, "_VCAM1-Fc")] <- pu_colors[i]
}

# actual plotting
p <- ggplot(data_clean, aes(x = treatment, 
                            y = !!sym(y_value_to_visualize), 
                            color = exp_treat)) +
    geom_beeswarm(cex = 4, size = 4, alpha = 0.8, priority = "ascending") + 
    geom_point(data = experiment_means, 
              aes(x = as.numeric(as.factor(treatment)) + 0.45, 
              y = !!sym(paste0("mean_", y_value_to_visualize)), 
              color = exp_treat),
              size = 7, shape = 19, alpha = 0.6, show.legend = FALSE) +  
    labs(title = "FOVs by Experiment and Treatment with Mean Synapse Counts",
        y = paste("Mean", gsub("_", " ", y_value_to_visualize)),
        color = "Experiment") +
    scale_color_manual(values = experiment_colors) +
    theme(
      legend.position = "top",  
      panel.background = element_blank(),
      panel.grid = element_blank(), 
      axis.line = element_line(color = "black", size = 1),  
      axis.ticks.x = element_blank(),
      axis.ticks.y = element_line(color = "black", size = 1),
      axis.ticks.length.y = unit(.25, "cm"),
      axis.text = element_text(size = 14), 
      axis.title.x = element_blank(),
      axis.title.y = element_text(size = 16)
    ) + 
    ylim(0, 2) + 
    coord_fixed(ratio = 2.5) 

p

# # visualize the data
# p <- ggplot(data_clean, aes(x = treatment, y = !!sym(y_value_to_visualize), color = as.factor(experiment))) +
#     geom_beeswarm(cex = 4, size = 3, alpha = 0.8, priority = "ascending") + 
#     geom_point(data = experiment_means, 
#               aes(x = as.numeric(treatment) + 0.45 , y = !!sym(paste0("mean_", y_value_to_visualize)), color = as.factor(experiment)),
#               size = 5, shape = 19, alpha = 0.6, show.legend = FALSE) +  
#     labs(title = "FOVs by Experiment and Treatment with Mean Synapse Counts",
#         y = paste("Mean", gsub("_", " ", y_value_to_visualize)),
#         color = "Experiment") +
#     scale_color_brewer(palette = "Dark2") +
#     theme(
#       legend.position = "top",  
#       panel.background = element_blank(),
#       panel.grid = element_blank(), 
#       axis.line = element_line(color = "black", size = 1),  
#       axis.ticks.x = element_blank(),
#       axis.ticks.y = element_line(color = "black", size = 1),
#       axis.ticks.length.y = unit(.25, "cm"),
#       axis.text = element_text(size = 14), 
#       axis.title.x = element_blank(),
#       axis.title.y = element_text(size = 16)
#     ) + 
#     ylim(0, 3.8) + 
#     coord_fixed(ratio = 1.2) 
# p

# ggsave("postsynapse_size_inh_norm.png", plot = p, 
#        width = 1500, height = 1800, units = "px", bg = "white", dpi = 300)

