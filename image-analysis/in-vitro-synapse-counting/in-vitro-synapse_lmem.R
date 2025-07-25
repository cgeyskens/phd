#### Trying to fit a linear mixed model
# (Dev Container from synapse-counting)

library(dplyr)
library(ggplot2)
library(lme4)
library(stringr)
library(emmeans)
library(tibble)
library(tidyr)
library(nlme)
library(stringr)

install.packages("ggbeeswarm")
library(ggbeeswarm)

remotes::install_github("nx10/httpgd")
library(httpgd)
hgd() # open the server for plotting


############################### Data prep & cleaning ################################

# loading in the data
pre_post_in_vitro <- read.csv("/mnt/image-analysis/in-vitro-synapse-counting/all_data/exp10_11_12_13_data_cydric.csv")
dim(pre_post_in_vitro)

# filter out FOVs that were not included in the analysis (had NA values)
data_clean <- pre_post_in_vitro[!is.na(pre_post_in_vitro$raw_synapse_count_per_100um), ]
dim(data_clean)

# ensure that treatment is a factor
data_clean$treatment <- factor(data_clean$treatment , levels = c("Fc", "VCAM1-Fc"))

# normalization by global Fc mean over all experiments
fc_global_synapse_mean <- data_clean %>%
  filter(treatment == "Fc") %>% 
  summarize(global_fc_synapse_mean = mean(raw_synapse_count_per_100um),
            global_fc_presynapse_mean = mean(raw_presynapse_count_per_100um),
            global_fc_postsynapse_mean = mean(raw_postsynapse_count_per_100um),
            global_fc_synapse_size_mean = mean(raw_synapse_puncta_size),
            global_fc_presynapse_size_mean = mean(raw_presynapse_puncta_size),
            global_fc_postsynapse_size_mean = mean(raw_postsynapse_puncta_size)
  )

# for synapse count
global_fc_synapse_mean_value <- fc_global_synapse_mean$global_fc_synapse_mean
global_fc_presynapse_mean_value <- fc_global_synapse_mean$global_fc_presynapse_mean
global_fc_postsynapse_mean_value <- fc_global_synapse_mean$global_fc_postsynapse_mean

data_clean$normalized_synapse_count_per_100um <- data_clean$raw_synapse_count_per_100um / global_fc_synapse_mean_value
data_clean$normalized_presynapse_count_per_100um <- data_clean$raw_presynapse_count_per_100um / global_fc_presynapse_mean_value
data_clean$normalized_postsynapse_count_per_100um <- data_clean$raw_postsynapse_count_per_100um / global_fc_postsynapse_mean_value

# for synapse size
global_fc_synapse_size_mean_value <- fc_global_synapse_mean$global_fc_synapse_size_mean
global_fc_presynapse_size_mean_value <- fc_global_synapse_mean$global_fc_presynapse_size_mean
global_fc_postsynapse_size_mean_value <- fc_global_synapse_mean$global_fc_postsynapse_size_mean

data_clean$normalized_synapse_size <- data_clean$raw_synapse_puncta_size / global_fc_synapse_size_mean_value
data_clean$normalized_presynapse_size <- data_clean$raw_presynapse_puncta_size / global_fc_presynapse_size_mean_value
data_clean$normalized_postsynapse_size <- data_clean$raw_postsynapse_puncta_size / global_fc_postsynapse_size_mean_value


############################### Model fitting, puncta density ################################

# synapse count model
synapse_count_model <- lmer(normalized_synapse_count_per_100um ~ treatment + (1 | experiment), data = data_clean)
summary(synapse_count_model)

plot(synapse_count_model)
qqnorm(resid(synapse_count_model))
qqline(resid(synapse_count_model))

synapse_count_emm <- emmeans(synapse_count_model, ~ treatment)
summary(synapse_count_emm)

synapse_count_pairwise_comparisons <- contrast(synapse_count_emm, method = "pairwise")
summary(synapse_count_pairwise_comparisons)

# presynapse count model
presynapse_count_model <- lmer(normalized_presynapse_count_per_100um ~ treatment + (1 | experiment), data = data_clean)
summary(presynapse_count_model)

plot(presynapse_count_model)
qqnorm(resid(presynapse_count_model))
qqline(resid(presynapse_count_model))

presynapse_count_emm <- emmeans(presynapse_count_model, ~ treatment)
summary(presynapse_count_emm)

presynapse_count_pairwise_comparisons <- contrast(presynapse_count_emm, method = "pairwise")
summary(presynapse_count_pairwise_comparisons)

# postsynapse count model
postsynapse_count_model <- lmer(normalized_postsynapse_count_per_100um ~ treatment + (1 | experiment), data = data_clean)
summary(postsynapse_count_model)

plot(postsynapse_count_model)
qqnorm(resid(postsynapse_count_model))
qqline(resid(postsynapse_count_model))

postsynapse_count_emm <- emmeans(postsynapse_count_model, ~ treatment)
summary(postsynapse_count_emm)

postsynapse_count_pairwise_comparisons <- contrast(postsynapse_count_emm, method = "pairwise")
summary(postsynapse_count_pairwise_comparisons)


############################### Data Visualization, puncta density ################################

y_value_to_visualize = "normalized_postsynapse_count_per_100um"

# calculate the mean of each experiment
experiment_means <- data_clean %>%
  group_by(experiment, treatment) %>%
  summarize(mean_normalized_synapse_count_per_100um = mean(normalized_synapse_count_per_100um),
            mean_normalized_presynapse_count_per_100um = mean(normalized_presynapse_count_per_100um),
            mean_normalized_postsynapse_count_per_100um = mean(normalized_postsynapse_count_per_100um)
  )

# visualize the data
p <- ggplot(data_clean, aes(x = treatment, y = !!sym(y_value_to_visualize), color = as.factor(experiment))) +
    geom_beeswarm(cex = 4, size = 3, alpha = 0.8, priority = "ascending") + 
    geom_point(data = experiment_means, 
              aes(x = as.numeric(treatment) + 0.45 , y = !!sym(paste0("mean_", y_value_to_visualize)), color = as.factor(experiment)),
              size = 5, shape = 19, alpha = 0.6, show.legend = FALSE) +  
    labs(title = "FOVs by Experiment and Treatment with Mean Synapse Counts",
        y = paste("Mean", gsub("_", " ", y_value_to_visualize)),
        color = "Experiment") +
    scale_color_brewer(palette = "Dark2") +
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
    ylim(0, 3.8) + 
    coord_fixed(ratio = 1.2) 
p

ggsave("postsynapse_exc_norm.png", plot = p, 
       width = 1500, height = 1800, units = "px", bg = "white", dpi = 300)



############################### Model fitting, puncta size ################################

# synapse size model
synapse_size_model <- lmer(normalized_synapse_size ~ treatment + (1 | experiment), data = data_clean)
summary(synapse_size_model)

plot(synapse_size_model)
qqnorm(resid(synapse_size_model))
qqline(resid(synapse_size_model))

synapse_size_emm <- emmeans(synapse_size_model, ~ treatment)
summary(synapse_size_emm)

synapse_size_pairwise_comparisons <- contrast(synapse_size_emm, method = "pairwise")
summary(synapse_size_pairwise_comparisons)

# presynapse size model
presynapse_size_model <- lmer(normalized_presynapse_size ~ treatment + (1 | experiment), data = data_clean)
summary(presynapse_size_model)

plot(presynapse_size_model)
qqnorm(resid(presynapse_size_model))
qqline(resid(presynapse_size_model))

presynapse_size_emm <- emmeans(presynapse_size_model, ~ treatment)
summary(synapse_size_emm)

presynapse_size_pairwise_comparisons <- contrast(presynapse_size_emm, method = "pairwise")
summary(presynapse_size_pairwise_comparisons)

# postsynapse size model
postsynapse_size_model <- lmer(normalized_postsynapse_size ~ treatment + (1 | experiment), data = data_clean)
summary(postsynapse_size_model)

plot(postsynapse_size_model)
qqnorm(resid(postsynapse_size_model))
qqline(resid(postsynapse_size_model))

postsynapse_size_emm <- emmeans(postsynapse_size_model, ~ treatment)
summary(synapse_size_emm)

postsynapse_size_pairwise_comparisons <- contrast(postsynapse_size_emm, method = "pairwise")
summary(postsynapse_size_pairwise_comparisons)


############################### Data Visualization, puncta size ################################

y_value_to_visualize = "normalized_postsynapse_size"

# calculate the mean of each experiment
experiment_means <- data_clean %>%
  group_by(experiment, treatment) %>%
  summarize(mean_normalized_synapse_size = mean(normalized_synapse_size),
            mean_normalized_presynapse_size = mean(normalized_presynapse_size),
            mean_normalized_postsynapse_size = mean(normalized_postsynapse_size)
  )

# visualize the data
p <- ggplot(data_clean, aes(x = treatment, y = !!sym(y_value_to_visualize), color = as.factor(experiment))) +
    geom_beeswarm(cex = 4, size = 3, alpha = 0.8, priority = "ascending") + 
    geom_point(data = experiment_means, 
              aes(x = as.numeric(treatment) + 0.45 , y = !!sym(paste0("mean_", y_value_to_visualize)), color = as.factor(experiment)),
              size = 5, shape = 19, alpha = 0.6, show.legend = FALSE) +  
    labs(title = "FOVs by Experiment and Treatment with Mean Synapse Counts",
        y = paste("Mean", gsub("_", " ", y_value_to_visualize)),
        color = "Experiment") +
    scale_color_brewer(palette = "Dark2") +
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
    ylim(0, 3.8) + 
    coord_fixed(ratio = 1.2) 
p

ggsave("postsynapse_size_exc_norm.png", plot = p, 
       width = 1500, height = 1800, units = "px", bg = "white", dpi = 300)


