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
remotes::install_github("nx10/httpgd")
library(httpgd)
hgd() # open the server for plotting


############################### Data prep & cleaning ################################

# loading in the data
vglut1_psd95_in_vitro <- read.csv("/mnt/image-analysis/in-vitro-synapse-counting/exp10_11_12_13_data.csv")
dim(vglut1_psd95_in_vitro)

# filter out FOVs that were not included in the analysis (had NA values)
data_clean <- vglut1_psd95_in_vitro[!is.na(vglut1_psd95_in_vitro$raw_synapse_count_per_100um), ]
dim(data_clean)

# ensure that treatment is a factor
data_clean$treatment <- factor(data_clean$treatment , levels = c("Fc", "VCAM1-Fc"))


############################### Model fitting ################################

# synapse model
synapse_model <- lmer(raw_synapse_count_per_100um ~ treatment + (1 | experiment), data = data_clean)
summary(synapse_model)

plot(synapse_model)
qqnorm(resid(synapse_model))
qqline(resid(synapse_model))

synapse_emm <- emmeans(synapse_model, ~ treatment)
summary(synapse_emm)

synapse_pairwise_comparisons <- contrast(synapse_emm, method = "pairwise")
summary(synapse_pairwise_comparisons)

# presynapse model
presynapse_model <- lmer(raw_presynapse_count_per_100um ~ treatment + (1 | experiment), data = data_clean)
summary(presynapse_model)

plot(presynapse_model)
qqnorm(resid(presynapse_model))
qqline(resid(presynapse_model))

presynapse_emm <- emmeans(presynapse_model, ~ treatment)
summary(presynapse_emm)

presynapse_pairwise_comparisons <- contrast(presynapse_emm, method = "pairwise")
summary(presynapse_pairwise_comparisons)

# postsynapse model
postsynapse_model <- lmer(raw_postsynapse_count_per_100um ~ treatment + (1 | experiment), data = data_clean)
summary(postsynapse_model)

plot(postsynapse_model)
qqnorm(resid(postsynapse_model))
qqline(resid(postsynapse_model))

postsynapse_emm <- emmeans(postsynapse_model, ~ treatment)
summary(postynapse_emm)

postsynapse_pairwise_comparisons <- contrast(postsynapse_emm, method = "pairwise")
summary(postsynapse_pairwise_comparisons)


############################### Data Visualization ################################

y_value_to_visualize = "raw_synapse_count_per_100um"

# calculate the mean of each experiment
experiment_means <- data_clean %>%
  group_by(experiment, treatment) %>%
  summarize(mean_raw_synapse_count_per_100um = mean(raw_synapse_count_per_100um),
            mean_raw_presynapse_count_per_100um = mean(raw_presynapse_count_per_100um),
            mean_raw_postsynapse_count_per_100um = mean(raw_postsynapse_count_per_100um)
  )

# visualize the data
ggplot(data_clean, aes(x = treatment, y = !!sym(y_value_to_visualize), color = as.factor(experiment))) +
  geom_jitter(width = 0.15, size = 4, alpha = 0.6) + 
  geom_point(data = experiment_means, 
             aes(x = as.numeric(treatment) + 0.35 , y = !!sym(paste0("mean_", y_value_to_visualize)), color = as.factor(experiment)),
             size = 7, shape = 19, show.legend = FALSE) +  
  labs(title = "FOVs by Experiment and Treatment with Mean Synapse Counts",
       y = paste("Mean", gsub("_", " ", y_value_to_visualize), "per 100um"),
       color = "Experiment") +
  scale_color_brewer(palette = "Paired") +
  theme_minimal(base_family = "Noto Sans") + 
  theme(
    legend.position = "top",  
    panel.grid = element_blank(),  
    axis.line = element_line(color = "black", size = 1),  
    axis.ticks = element_blank(),
    axis.text = element_text(size = 14), 
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 16)
  ) + 
  coord_fixed(ratio = 0.15)
