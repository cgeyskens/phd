library(DEP)
library(dplyr)
library(SummarizedExperiment)
library(ggplot2)
library(readr)

# loading the data
exp17_data <- read_tsv("/mnt/ip-proteomics/exp17-my-diann-run/diann_output.pg_matrix.tsv")
View(exp17_data)

# ////
# checking the raw value of certain proteins
vcam1_data <- exp17_data[exp17_data$Genes == 'Vcam1',]
View(vcam1_data)
# \\\\

# checking colnames
colnames(exp17_data)

# vcam1 synglio column names
vcam1_synglio_names <- c(
    "vcam1_p14_ip_r1" = "D:\\Proteomics2024\\CBD\\JorisdeWitLab\\PCF000184\\CG_1.wiff",
    "vcam1_p14_igg_r1" = "D:\\Proteomics2024\\CBD\\JorisdeWitLab\\PCF000184\\CG_2.wiff",
    "vcam1_p14_ip_r2" = "D:\\Proteomics2024\\CBD\\JorisdeWitLab\\PCF000184\\CG_3.wiff",
    "vcam1_p14_igg_r2" = "D:\\Proteomics2024\\CBD\\JorisdeWitLab\\PCF000184\\CG_4.wiff",
    "vcam1_p14_ip_r3" = "D:\\Proteomics2024\\CBD\\JorisdeWitLab\\PCF000184\\CG_5.wiff",
    "vcam1_p14_igg_r3" = "D:\\Proteomics2024\\CBD\\JorisdeWitLab\\PCF000184\\CG_6.wiff",
    "vcam1_p14_ip_r4" = "D:\\Proteomics2024\\CBD\\JorisdeWitLab\\PCF000184\\CG_7.wiff", 
    "vcam1_p14_igg_r4" = "D:\\Proteomics2024\\CBD\\JorisdeWitLab\\PCF000184\\CG_8.wiff", 
    "vcam1_p28_ip_r1" = "D:\\Proteomics2024\\CBD\\JorisdeWitLab\\PCF000184\\CG_9.wiff",
    "vcam1_p28_igg_r1" = "D:\\Proteomics2024\\CBD\\JorisdeWitLab\\PCF000184\\CG_10.wiff",
    "vcam1_p28_ip_r2" = "D:\\Proteomics2024\\CBD\\JorisdeWitLab\\PCF000184\\CG_11.wiff",
    "vcam1_p28_igg_r2" = "D:\\Proteomics2024\\CBD\\JorisdeWitLab\\PCF000184\\CG_12.wiff",
    "vcam1_p28_ip_r3" = "D:\\Proteomics2024\\CBD\\JorisdeWitLab\\PCF000184\\CG_13.wiff",
    "vcam1_p28_igg_r3" = "D:\\Proteomics2024\\CBD\\JorisdeWitLab\\PCF000184\\CG_14.wiff",
    "vcam1_p28_ip_r4" = "D:\\Proteomics2024\\CBD\\JorisdeWitLab\\PCF000184\\CG_15.wiff",
    "vcam1_p28_igg_r4" = "D:\\Proteomics2024\\CBD\\JorisdeWitLab\\PCF000184\\CG_16.wiff"
)

colnames(exp17_data)[match(vcam1_synglio_names, colnames(exp17_data))] <- names(vcam1_synglio_names)

# gpr37l1 colnames
gpr37l1_names <- c(
    "gpr37l1_ip_r1" = "D:\\Proteomics2025\\VIB\\CBD\\JorisDeWitLab\\PCF000321\\30SPD\\PCF000321_30spd\\CG_01.wiff",
    "gpr37l1_igg_r1" = "D:\\Proteomics2025\\VIB\\CBD\\JorisDeWitLab\\PCF000321\\30SPD\\PCF000321_30spd\\CG_02.wiff",
    "gpr37l1_ip_r2" = "D:\\Proteomics2025\\VIB\\CBD\\JorisDeWitLab\\PCF000321\\30SPD\\PCF000321_30spd\\CG_03.wiff",
    "gpr37l1_igg_r2" = "D:\\Proteomics2025\\VIB\\CBD\\JorisDeWitLab\\PCF000321\\30SPD\\PCF000321_30spd\\CG_04.wiff",
    "gpr37l1_ip_r3" = "D:\\Proteomics2025\\VIB\\CBD\\JorisDeWitLab\\PCF000321\\30SPD\\PCF000321_30spd\\CG_05.wiff",
    "gpr37l1_igg_r3" = "D:\\Proteomics2025\\VIB\\CBD\\JorisDeWitLab\\PCF000321\\30SPD\\PCF000321_30spd\\CG_06.wiff",
    "gpr37l1_ip_r4" = "D:\\Proteomics2025\\VIB\\CBD\\JorisDeWitLab\\PCF000321\\30SPD\\PCF000321_30spd\\CG_07.wiff",
    "gpr37l1_igg_r4" = "D:\\Proteomics2025\\VIB\\CBD\\JorisDeWitLab\\PCF000321\\30SPD\\PCF000321_30spd\\CG_08.wiff"
)

colnames(exp17_data)[match(gpr37l1_names, colnames(exp17_data))] <- names(gpr37l1_names)

# Rename columns
new_names <- c(
    "gpr37l1_ip" = "/ip-proteomics/ms-convert-output/2ul_CG_1.mzML",
    "gpr37l1_igg" = "/ip-proteomics/ms-convert-output/2ul_CG_2.mzML",
    "vcam1_ip_r1" = "/ip-proteomics/ms-convert-output/2ul_CG_3.mzML",
    "vcam1_igg_r1"= "/ip-proteomics/ms-convert-output/2ul_CG_4.mzML",
    "vcam1_ip_r2" = "/ip-proteomics/ms-convert-output/2ul_CG_5.mzML",
    "vcam1_igg_r2"= "/ip-proteomics/ms-convert-output/2ul_CG_6.mzML",
    "vcam1_ip_r3" = "/ip-proteomics/ms-convert-output/2ul_CG_7.mzML",
    "vcam1_igg_r3"= "/ip-proteomics/ms-convert-output/2ul_CG_8.mzML",
    "vcam1_ip_r4" = "/ip-proteomics/ms-convert-output/2ul_CG_9.mzML",
    "vcam1_igg_r4" = "/ip-proteomics/ms-convert-output/2ul_CG_10.mzML"
)
colnames(exp17_data)[match(new_names, colnames(exp17_data))] <- names(new_names)

# Pedro colnames
pedro_names <- c(
    "gpr37l1_ip" = "D..Proteomics2025.VIB.CBD.JorisDeWit.PCF000270.2ul_CG_1.wiff",
    "gpr37l1_igg" = "D..Proteomics2025.VIB.CBD.JorisDeWit.PCF000270.2ul_CG_2.wiff",
    "vcam1_ip_r1" = "D..Proteomics2025.VIB.CBD.JorisDeWit.PCF000270.2ul_CG_3.wiff",
    "vcam1_igg_r1" = "D..Proteomics2025.VIB.CBD.JorisDeWit.PCF000270.2ul_CG_4.wiff",
    "vcam1_ip_r2" = "D..Proteomics2025.VIB.CBD.JorisDeWit.PCF000270.2ul_CG_5.wiff",
    "vcam1_igg_r2" = "D..Proteomics2025.VIB.CBD.JorisDeWit.PCF000270.2ul_CG_6.wiff",
    "vcam1_ip_r3" = "D..Proteomics2025.VIB.CBD.JorisDeWit.PCF000270.2ul_CG_7.wiff", 
    "vcam1_igg_r3" = "D..Proteomics2025.VIB.CBD.JorisDeWit.PCF000270.2ul_CG_8.wiff", 
    "vcam1_ip_r4" = "D..Proteomics2025.VIB.CBD.JorisDeWit.PCF000270.2ul_CG_9.wiff",
    "vcam1_igg_r4" = "D..Proteomics2025.VIB.CBD.JorisDeWit.PCF000270.2ul_CG_10.wiff"
)
colnames(exp17_data)[match(pedro_names, colnames(exp17_data))] <- names(pedro_names)

# remove gpr37l1 data
df <- exp17_data %>%
    select(-c(
        "gpr37l1_ip",
        "gpr37l1_igg"
      )
    )

# constructing the summarizedExperiment object
my_data_unique <- make_unique(
  proteins = df, 
  names = "Genes", 
  ids = "Protein.Group",
  delim = ";")

my_data_unique[my_data_unique$Genes == 'Vcam1',]

columns = as.integer(5:12) # specifying which columns have the wanted intensities

exp_design = as.data.frame(
  list(
  label = c(
    # "vcam1_p14_ip_r1",
    # "vcam1_p14_igg_r1",
    # "vcam1_p14_ip_r2",
    # "vcam1_p14_igg_r2",
    # "vcam1_p14_ip_r3",
    # "vcam1_p14_igg_r3",
    # "vcam1_p14_ip_r4",
    # "vcam1_p14_igg_r4",
    # "vcam1_p28_ip_r1",
    # "vcam1_p28_igg_r1",
    # "vcam1_p28_ip_r2",
    # "vcam1_p28_igg_r2",
    # "vcam1_p28_ip_r3",
    # "vcam1_p28_igg_r3",
    # "vcam1_p28_ip_r4",
    # "vcam1_p28_igg_r4",
    "vcam1_ip_r1", 
    "vcam1_igg_r1",
    "vcam1_ip_r2" ,
    "vcam1_igg_r2",
    "vcam1_ip_r3" ,
    "vcam1_igg_r3",
    "vcam1_ip_r4" ,
    "vcam1_igg_r4"
  ),
  condition = c(
    # "p14_ip", 
    # "p14_igg", 
    # "p14_ip", 
    # "p14_igg", 
    # "p14_ip", 
    # "p14_igg", 
    # "p14_ip", 
    # "p14_igg", 
    # "p28_ip", 
    # "p28_igg", 
    # "p28_ip", 
    # "p28_igg", 
    # "p28_ip", 
    # "p28_igg", 
    # "p28_ip", 
    # "p28_igg",
    "ip", 
    "igg",
    "ip" ,
    "igg",
    "ip" ,
    "igg",
    "ip" ,
    "igg"        
    ),
  replicate = c("1", "1", "2", "2", "3", "3", "4", "4")
  )
)

se <- make_se(my_data_unique, columns, exp_design)
se

# Plot nr of proteins
proteins_per_sample <- plot_numbers(se)
proteins_per_sample <- proteins_per_sample + geom_text(aes(label = sum, y= sum),  
                                     vjust = 2) 
proteins_per_sample

# imputate missing values
se_imp <- impute(se, fun = "MinProb")

# variance stabilization
se_norm <- normalize_vsn(se_imp)

# Visualize normalization
plot_normalization(se_norm, se_imp, se)

assay(se_norm)["Vcam1", ]


# Differential expression analysis
data_diff <- DEP::test_diff(se_norm, type = "manual", test = c("ip_vs_igg"))

# Add rejections
dep <- add_rejections(data_diff, alpha = 0.05, lfc = log2(1))

# Plot PCA and Volcano
plot_pca(dep, x = 1, y = 2, n = 100, point_size = 4)
plot_volcano(dep, contrast = "ip_vs_igg", label_size = 2, add_names = TRUE)

# Plot the pearson's correlation plot
plot_cor(dep, significant = TRUE, lower = 0, upper = 1, pal = "Reds")

# Plot heatmap
plot_heatmap(dep, type = "centered", #kmeans = TRUE, k = 2,
             col_limit = 4, show_row_names = TRUE,
             indicate = c("condition", "replicate"),
             row_font_size = 4)



# line_plot of some log2 raw intensities
extract_protein_data <- function(se, protein_name) {
  df <- data.frame(assay(se)[protein_name, ])
  colnames(df) <- "raw_log2_intensities"
  df$condition <- rownames(df)
  df$protein <- protein_name
  return(df)
}

protein_list <- c(
  "Vcam1", 
  "Prrt1", 
  "Pmch", 
  "Disp2", 
  "Sparc", 
  "Igdcc4", 
  "Pcdhb14", 
  "Chgb"
  )
df_merge <- bind_rows(lapply(protein_list, extract_protein_data, se = se))

ggplot(data=df_merge, aes(x=condition, y=raw_log2_intensities, group=protein)) +
  geom_line(aes(color=protein))+
  geom_point(aes(color=protein)) + 
  scale_x_discrete() +
  scale_y_continuous(expand=c(0, 0), limits=c(0, 22)) + 
  theme_grey(base_size = 22)

dep_results <- get_results(dep)

dep_ip_proteins <- dep_results %>%
                          filter(ip_vs_igg_p.val <=0.05 & ip_vs_igg_ratio > 1)

write.csv(dep_ip_proteins,"vcam1_dep_ip_proteins.csv", row.names = TRUE)
