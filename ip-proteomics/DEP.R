library(DEP)
library(dplyr)
library(SummarizedExperiment)

# loading the data
exp17_data <- read.delim("/mnt/data/ip-proteomics/exp17-my-diann-run/diann_output.gg_matrix.tsv")
View(exp17_data)

# ////
# checking the raw value of certain proteins
vcam1_data <- exp17_data[exp17_data$Genes == 'Csnk2b',]
View(vcam1_data)
# \\\\

# Rename columns
new_names <- c(
    "gpr37l1_ip" = "X.ip.proteomics.ms.convert.output.2ul_CG_1.mzML",
    "gpr37l1_igg" = "X.ip.proteomics.ms.convert.output.2ul_CG_2.mzML",
    "vcam1_ip_r1" = "X.ip.proteomics.ms.convert.output.2ul_CG_3.mzML",
    "vcam1_igg_r1" = "X.ip.proteomics.ms.convert.output.2ul_CG_4.mzML",
    "vcam1_ip_r2" = "X.ip.proteomics.ms.convert.output.2ul_CG_5.mzML",
    "vcam1_igg_r2" = "X.ip.proteomics.ms.convert.output.2ul_CG_6.mzML",
    "vcam1_ip_r3" = "X.ip.proteomics.ms.convert.output.2ul_CG_7.mzML", 
    "vcam1_igg_r3" = "X.ip.proteomics.ms.convert.output.2ul_CG_8.mzML", 
    "vcam1_ip_r4" = "X.ip.proteomics.ms.convert.output.2ul_CG_9.mzML",
    "vcam1_igg_r4" = "X.ip.proteomics.ms.convert.output.2ul_CG_10.mzML"
)
colnames(exp17_data)[match(new_names, colnames(exp17_data))] <- names(new_names)

# remove gpr37l1 data
df <- exp17_data %>%
    select(-c("gpr37l1_ip", "gpr37l1_igg"))

# constructing the summarizedExperiment object
my_data_unique <- make_unique(
  proteins = df, 
  names = "Genes", 
  ids = "Genes",
  delim = ";")

columns = as.integer(2:9) # specifying which columns

exp_design = as.data.frame(
  list(
  label = c("vcam1_ip_r1", "vcam1_igg_r1", "vcam1_ip_r2", "vcam1_igg_r2", "vcam1_ip_r3", "vcam1_igg_r3", "vcam1_ip_r4", "vcam1_igg_r4"),
  condition = c("ip", "igg", "ip", "igg", "ip", "igg", "ip", "igg"),
  replicate = c("1", "1", "2", "2", "3", "3", "4", "4")
  )
)

se <- make_se(my_data_unique, columns, exp_design)
se

# imputate missing values
se_imp <- impute(se, fun = "MinProb")

# normalization
se_norm <- normalize_vsn(se_imp)

# Visualize normalization
plot_normalization(se_norm, se_imp)

assay(se_imp)["Vcam1",]

# 12. Differential expression analysis
data_diff <- test_diff(se_imp, type = "manual", test = c("ip_vs_igg"))

# 13. Add rejections
dep <- add_rejections(data_diff, alpha = 0.05, lfc = log2(1.5))

# Plot PCA and Volcano
plot_pca(dep, x = 1, y = 2, n = 100, point_size = 4)
plot_volcano(dep, contrast = "ip_vs_igg", label_size = 2, add_names = TRUE)


