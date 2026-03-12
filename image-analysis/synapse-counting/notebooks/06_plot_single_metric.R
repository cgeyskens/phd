# ============================================================================
# Title: Visualization of specific metrics 
# Input: Output of the pipeline with the all metrics analyzed
# Description: 
# (1) data processing
# (2) normalization
# (3) plotting
# Date: Nov 2025
# ============================================================================

library(dplyr)
library(ggplot2)

#### =============================== arguments =============================== ####
protein = "GPR37L1"
synapse = "VGAT-GEPH"
layer = "CA3 SL"
metric = "post_puncta_density_per_100_um2"

#### =============================== get input =============================== ####
if (protein == "VCAM1"){
    data_vglut1 <- read.csv("/mnt/image-analysis/synapse-counting/IHC_Exp9_VCAM1_VGLUT1-PSD95_output_data_20251201_194547/metric_results.csv")
    data_vgat <- read.csv("/mnt/image-analysis/synapse-counting/IHC_Exp12_VCAM1_VGAT-GEPH_output_data_20251201_194539/metric_results.csv")
} else if (protein == "GPR37L1") {
    data_vglut1 <- read.csv("/mnt/image-analysis/synapse-counting/IHC_Exp10_GPR37L1_VGLUT1-PSD95_output_data_20251201_194537/metric_results.csv")
    data_vgat <- read.csv("/mnt/image-analysis/synapse-counting/IHC_Exp11_GPR37L1_VGAT-GEPH_output_data_20251201_194548/metric_results.csv")
}

if (synapse == "VGLUT1-PSD95"){
  data = data_vglut1
} else if (synapse == "VGAT-GEPH") {
  data = data_vgat
}

#### =============================== data processing =============================== ####
# subset data
metric_data <- data %>%
  filter(hippocampal_layer == layer) %>%
  select(all_of(metric), gRNA, Brain)

# # get the LacZ global value
# lacz_global <- ca1_sr_data_local_peak_max %>%
#   filter(gRNA == "LacZ-gRNA") %>%
#   summarise(m = mean(local_peak_colocalized_spots, na.rm = TRUE)) %>%
#   pull(m)

# # normalize
# df_norm <- ca1_sr_data_local_peak_max %>%
#   mutate(local_peak_colocalized_spots_norm = local_peak_colocalized_spots / lacz_global)

# get brain-level means
brain_means <- metric_data %>%
  group_by(gRNA, Brain) %>%
  summarise(brain_mean = mean(.data[[metric]], na.rm = TRUE),
            .groups = "drop")

# global LacZ mean
lacz_global <- brain_means %>%
  filter(gRNA == "LacZ-gRNA") %>%
  summarise(m = mean(brain_mean, na.rm = TRUE)) %>%
  pull(m)

# FOV-level normalized + numeric x
fov_df <- metric_data %>%
  mutate(
    y_norm = .data[[metric]] / lacz_global,
    x = ifelse(gRNA == "LacZ-gRNA", 1, 2)
  )

# brain-level normalized + nudged numeric x (inside)
brain_df <- brain_means %>%
  mutate(
    y_norm = brain_mean / lacz_global,
    x = ifelse(gRNA == "LacZ-gRNA", 1, 2),
    xn = ifelse(gRNA == "LacZ-gRNA", 1 + 0.2, 2 - 0.2)
  )


#### =============================== plotting =============================== ####
# specifying the colors
if (protein == "VCAM1"){
  grna_cols <- c(
    "VCAM1-gRNA" = "#21a0e2",
    "LacZ-gRNA" = "#c6c6c6"
  )
} else if (protein == "GPR37L1") {
  grna_cols <- c(
    "GPR37L1-gRNA" = "#e28d21",
    "LacZ-gRNA" = "#c6c6c6"
  )
}

p <- ggplot() +
  # FOV points (outside-ish via jitter)
  geom_point(
    data = fov_df,
    aes(x = x, y = y_norm, color = gRNA),
    # position = position_jitter(width = 0.10, height = 0),
    alpha = 0.6, size = 5
  ) +
  # paired brain means line (inside)
  geom_line(
    data = brain_df,
    aes(x = xn, y = y_norm, group = Brain),
    linewidth = 0.8, color = "black"
  ) +
  # brain mean points (inside, colored, different shape)
  geom_point(
    data = brain_df,
    aes(x = xn, y = y_norm, fill = gRNA),
    shape = 21, size = 8, stroke = 0, alpha = 1
  ) +
  scale_x_continuous(
  ) +
  labs(
    x = NULL,
    y = metric
  ) +
  scale_color_manual(
    values = grna_cols
  ) + 
  scale_fill_manual(
     values = grna_cols
  ) + 
  theme_classic() +
  coord_cartesian(xlim = c(0.7, 2.2)) +
  theme(
    legend.position = "none",
    panel.background = element_blank(),
    panel.grid = element_blank(), 
    axis.line = element_line(color = "black", size = 1.2),  
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_line(color = "black", size = 1.2),
    axis.ticks.length.y = unit(.3, "cm"),
    axis.text = element_text(size = 14), 
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 16)
    )  +
    scale_y_continuous(
        expand = c(0, 0),
        limits = c(0, 2),
        breaks = c(0, 1, 2) 
    ) +  
    scale_x_continuous(
        expand = c(0, 0),
        limits = c(0.8, 2.2),
        breaks = c(0, 1, 2) 
    ) +
    coord_fixed(ratio = 1.8)  # 0r 1.3
p

ggsave(paste0(synapse, "_", layer, "_", metric, "_paper.svg"), 
    plot = p, 
    device = cairo_pdf,
    width = 8, height = 20, units = "cm", dpi=300)
