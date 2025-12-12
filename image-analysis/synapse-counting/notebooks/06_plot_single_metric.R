library(dplyr)
library(ggplot2)





# import data
file_path_vglut1 <- "/mnt/image-analysis/synapse-counting/IHC_Exp9_VCAM1_VGLUT1-PSD95_output_data_20251201_194547/metric_results.csv"
df_vglut1 = read.csv(file_path_vglut1)

# subset data
ca1_sr_data_local_peak_max <- df_vglut1 %>%
  filter(hippocampal_layer == "CA1 SR") %>%
  select(local_peak_colocalized_spots, gRNA, Brain)

# get the LacZ global value
lacz_global <- ca1_sr_data_local_peak_max %>%
  filter(gRNA == "LacZ-gRNA") %>%
  summarise(m = mean(local_peak_colocalized_spots, na.rm = TRUE)) %>%
  pull(m)

# normalize
df_norm <- ca1_sr_data_local_peak_max %>%
  mutate(local_peak_colocalized_spots_norm = local_peak_colocalized_spots / lacz_global)

# get brain mean
df <- ca1_sr_data_local_peak_max %>%
  group_by(gRNA, Brain) %>%
  summarise(brain_mean = mean(local_peak_colocalized_spots, na.rm = TRUE),
            .groups = "drop")
library(dplyr)
library(ggplot2)

# brain-level means
brain_means <- ca1_sr_data_local_peak_max %>%
  group_by(gRNA, Brain) %>%
  summarise(brain_mean = mean(local_peak_colocalized_spots, na.rm = TRUE),
            .groups = "drop")

# global LacZ mean
lacz_global <- brain_means %>%
  filter(gRNA == "LacZ-gRNA") %>%
  summarise(m = mean(brain_mean, na.rm = TRUE)) %>%
  pull(m)

# FOV-level normalized + numeric x
fov_df <- ca1_sr_data_local_peak_max %>%
  mutate(
    y_norm = local_peak_colocalized_spots / lacz_global,
    x = ifelse(gRNA == "LacZ-gRNA", 1, 2)
  )

# brain-level normalized + nudged numeric x (inside)
brain_df <- brain_means %>%
  mutate(
    y_norm = brain_mean / lacz_global,
    x = ifelse(gRNA == "LacZ-gRNA", 1, 2),
    xn = ifelse(gRNA == "LacZ-gRNA", 1 + 0.12, 2 - 0.12)
  )

p <- ggplot() +
  # FOV points (outside-ish via jitter)
  geom_point(
    data = fov_df,
    aes(x = x, y = y_norm, color = gRNA),
    # position = position_jitter(width = 0.10, height = 0),
    alpha = 0.75, size = 2
  ) +
  # paired brain means line (inside)
  geom_line(
    data = brain_df,
    aes(x = xn, y = y_norm, group = Brain),
    linewidth = 0.6, color = "black"
  ) +
  # brain mean points (inside, colored, different shape)
  geom_point(
    data = brain_df,
    aes(x = xn, y = y_norm, fill = gRNA),
    shape = 21, color = "black", size = 3.2, stroke = 1
  ) +
  scale_x_continuous(
  ) +
  labs(
    x = NULL,
    y = "local_peak_colocalized_spots / global LacZ mean"
  ) +
  theme_classic() +
  coord_cartesian(xlim = c(0.7, 3)) +
  theme(
    legend.position = "none",
    axis.ticks.x = element_blank(),
    ) + 
p

ggsave(paste0("VGLUT1-PSD95_CA1_SR_metric_paper.svg"), 
    plot = p, 
    device = cairo_pdf,
    width = 8, height = 20, units = "cm", dpi=300)
