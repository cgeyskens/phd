# test script for multivariate linear mixed-effect model
# need to have a generalized mixed effects model, because the data is not normally distributed
# loading data
vcam1_data <- read.csv("/mnt/image-analysis/synapse-counting/VCAM1/VCAM1-LacZ_VGLUT1-PSD95_output_data/metric_results.csv")
dim(vcam1_data)

# getting the libraries
library(dplyr)
library(ggplot2)
library(tidyr)

# data preprocessing for plotting
long_df <- gather(vcam1_data, 
                    key="variable", 
                    value="value", 
                    local_peak_colocalized_spots, 
                    overlap_coeff,
                    overlap_um2,
                    pearson_cor,
                    presynapse_image_mfi,
                    postsynapse_image_mfi,
                    pre_puncta_density_per_100_um2,
                    post_puncta_density_per_100_um2,
                    pre_staining_area_um2,
                    post_staining_area_um2,
                    pre_mean_puncta_size_um2,
                    post_mean_puncta_size_um2)
dim(long_df)

# plotting the distribution of each metric
ggplot(long_df, aes(x = value)) +
  geom_histogram(bins = 50, color = "white", fill = "blue", alpha = 1) +  # Adjust binwidth as needed
  facet_wrap(~ variable, scales = "free") +  # Create separate panels for each variable
  theme_bw() +
  labs(title = "Histograms", 
       x = "Value", 
       y = "Count") +
  theme(legend.position = "none") 





