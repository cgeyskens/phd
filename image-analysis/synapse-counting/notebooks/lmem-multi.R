# test script for multivariate linear mixed-effect model
# need to have a generalized mixed effects model, because the data is not normally distributed
# loading data
vcam1_data <- read.csv("/Volumes/KINGSTON/data/phd/image-analysis/synapse-counting/VCAM1/VCAM1-LacZ_VGAT-GEPH_output_data/metric_results.csv")

# getting the libraries
library(dplyr)
library(ggplot2)
library(DHARMa)

ggplot(your_data, aes(x = your_variable)) +
  geom_histogram(fill = "lightblue", 
                color = "black", 
                bins = 30) +
  geom_density(aes(y = ..density..), color = "red") + # Add density curve
  theme_minimal() +
  labs(title = "Distribution of Your Variable",
       x = "Your Variable Name",
       y = "Count") +
  theme(plot.title = element_text(size = 14, hjust = 0.5),
        axis.title = element_text(size = 12))

