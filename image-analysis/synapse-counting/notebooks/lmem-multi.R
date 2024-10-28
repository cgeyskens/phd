# test script for multivariate linear mixed-effect model

# loading data
vcam1_data <- read.csv("/Volumes/KINGSTON/data/phd/image-analysis/synapse-counting/VCAM1/VCAM1-LacZ_VGAT-GEPH_output_data/metric_results.csv")

# getting the libraries
library(dplyr)
library(ggplot2)
library(glmmTMB)

install.packages("glmmTMB")
