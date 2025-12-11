library(dplyr)
library(ggplot2)
library(readxl)
library(emmeans)
library(RColorBrewer)
library(ggbeeswarm)


#### =============================== get input =============================== ####
data <- read_excel(
  "/mnt/image-analysis/synapse-counting/results_202512/VCAM1_gRNA-QC_WB.xlsx"
)


#### ============= normalization to LacZ-gRNA#1 for each exp ================= ####
lacZ_vals <- data %>% 
  filter(gRNA == "LacZ-gRNA1") %>% 
  select(experiment, lacZ_val = vcam1_norm_protein)

data_norm <- data %>% 
  left_join(lacZ_vals, by = "experiment") %>% 
  mutate(vcam1_rel_to_LacZ = vcam1_norm_protein / lacZ_val)


#### ======================== statistics =================================== ####

# testing normality per gRNA group
data_norm %>%
  group_by(gRNA) %>%
  summarise(p_value = shapiro.test(vcam1_norm_protein)$p.value)

# check visually
ggplot(data_norm, aes(sample = vcam1_rel_to_LacZ)) +
  stat_qq() +
  stat_qq_line() +
  facet_wrap(~ gRNA)

# anova
fit <- lm(vcam1_rel_to_LacZ ~ gRNA + experiment, data = data_norm)
anova(fit)

# posthoc with Dunnet's test
emm <- emmeans(fit, "gRNA")
contrast(emm, method = "trt.vs.ctrl", ref = "LacZ-gRNA1", adjust = "BH")


#### ======================== plotting =================================== ####

data_norm <- data_norm %>%
  mutate(exp_gRNA = paste0(experiment, "_", gRNA))

data_norm <- data_norm %>%
  mutate(
    group2 = ifelse(grepl("VCAM1", gRNA),
                    "VCAM1-gRNA",
                    "GFP/LacZ")  
  )

data_means <- data_norm %>%
  group_by(gRNA, group2) %>%
  summarise(
    mean_vcam1_rel_to_LacZ = mean(vcam1_rel_to_LacZ, na.rm = TRUE),
    .groups = "drop"
  )


p <- ggplot() +
  geom_col(
    data = data_means,
    aes(x = gRNA,
        y = mean_vcam1_rel_to_LacZ,
        fill = group2),
    width = 0.6,
    alpha = 0.4,
    color = NA
  ) +
  geom_beeswarm(
    data = data_norm,
    aes(x = gRNA,
        y = vcam1_rel_to_LacZ,
        color = group2),
    cex = 2,
    size = 5,
    alpha = 0.9,
    priority = "ascending",
    stroke = 0,
    corral = "wrap",
    corral.width = 0.6
  ) +
  labs(
    title = "vcam1_rel_to_LacZ across gRNA groups",
    y = "vcam1_rel_to_LacZ",
    color = "Group",
    fill  = "Group"
  ) +
  scale_color_manual(
    values = c("GFP/LacZ" = "grey20",
               "VCAM1-gRNA" = "#21a0e2")
  ) +
  scale_fill_manual(
    values = c("GFP/LacZ" = "#c6c6c6",
               "VCAM1-gRNA" = "#21a0e2")
  ) +
  theme(
    legend.position = "top",
    panel.background = element_blank(),
    panel.grid = element_blank(),
    axis.line = element_line(color = "black", size = 1.2),
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_line(color = "black", size = 1.2),
    axis.ticks.length.y = unit(.3, "cm"),
    axis.text  = element_text(size = 14),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 16)
  ) +
  scale_y_continuous(
    expand = c(0, 0),
    limits = c(0, 2)   # tweak if needed
  )
p

ggsave(paste0("in-vitro_vcam1-gRNA_QC_plot_paper.svg"), 
    plot = p, 
    device = cairo_pdf,
    width = 25, height = 20, units = "cm", dpi=300)
