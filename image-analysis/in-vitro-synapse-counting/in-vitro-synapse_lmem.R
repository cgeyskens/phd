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
library(extrafont)

install.packages("ggbeeswarm")
library(ggbeeswarm)

remotes::install_github("nx10/httpgd")
library(httpgd)
hgd() # open the server for plotting


############################### Data prep & cleaning ################################

# loading in the data
pre_post_in_vitro <- read.csv("/mnt/image-analysis/in-vitro-synapse-counting/exp17_20_data.csv")
dim(pre_post_in_vitro)

# filter out FOVs that were not included in the analysis (had NA values)
data_clean <- pre_post_in_vitro[!is.na(pre_post_in_vitro$raw_synapse_count_per_100um), ]
dim(data_clean)

# ensure that treatment is a factor
data_clean$treatment <- factor(data_clean$treatment , levels = c("Fc", "VCAM1-Fc"))

# normalization by Fc mean over all experiments
fc_global_synapse_mean <- data_clean %>%
  filter(treatment == "Fc") %>% 
  summarize(global_fc_synapse_mean = mean(raw_synapse_count_per_100um),
            global_fc_presynapse_mean = mean(raw_presynapse_count_per_100um),
            fc_global_postsynapse_mean = mean(raw_postsynapse_count_per_100um)
  )

global_fc_synapse_mean_value <- fc_global_synapse_mean$global_fc_synapse_mean
global_fc_presynapse_mean_value <- fc_global_synapse_mean$global_fc_presynapse_mean
fc_global_postsynapse_mean_value <- fc_global_synapse_mean$fc_global_postsynapse_mean

data_clean$normalized_synapse_count_per_100um <- data_clean$raw_synapse_count_per_100um / global_fc_synapse_mean_value
data_clean$normalized_presynapse_count_per_100um <- data_clean$raw_presynapse_count_per_100um / global_fc_presynapse_mean_value
data_clean$normalized_postsynapse_count_per_100um <- data_clean$raw_postsynapse_count_per_100um / fc_global_postsynapse_mean_value


############################### Model fitting ################################

# synapse model
synapse_model <- lmer(normalized_synapse_count_per_100um ~ treatment + (1 | experiment), data = data_clean)
summary(synapse_model)

plot(synapse_model)
qqnorm(resid(synapse_model))
qqline(resid(synapse_model))

synapse_emm <- emmeans(synapse_model, ~ treatment)
summary(synapse_emm)

synapse_pairwise_comparisons <- contrast(synapse_emm, method = "pairwise")
summary(synapse_pairwise_comparisons)

# presynapse model
presynapse_model <- lmer(normalized_presynapse_count_per_100um ~ treatment + (1 | experiment), data = data_clean)
summary(presynapse_model)

plot(presynapse_model)
qqnorm(resid(presynapse_model))
qqline(resid(presynapse_model))

presynapse_emm <- emmeans(presynapse_model, ~ treatment)
summary(presynapse_emm)

presynapse_pairwise_comparisons <- contrast(presynapse_emm, method = "pairwise")
summary(presynapse_pairwise_comparisons)

# postsynapse model
postsynapse_model <- lmer(normalized_postsynapse_count_per_100um ~ treatment + (1 | experiment), data = data_clean)
summary(postsynapse_model)

plot(postsynapse_model)
qqnorm(resid(postsynapse_model))
qqline(resid(postsynapse_model))

postsynapse_emm <- emmeans(postsynapse_model, ~ treatment)
summary(postsynapse_emm)

postsynapse_pairwise_comparisons <- contrast(postsynapse_emm, method = "pairwise")
summary(postsynapse_pairwise_comparisons)


############################### Data Visualization ################################

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

ggsave("postsynapse_inh_norm.png", plot = p, 
       width = 1500, height = 1800, units = "px", bg = "white", dpi = 300)
