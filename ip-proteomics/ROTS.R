library(ROTS)
library(readr)
library(dplyr)
library(tibble)
library(imputeLCMD)

# Loading the data
exp17_data <- read_tsv("/mnt/data/ip-proteomics/exp17-my-diann-run/diann_output.pg_matrix.tsv")
View(exp17_data)

# checking colnames
colnames(exp17_data)

# gpr37l1 colnames
exp17_names <- c(
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

colnames(exp17_data)[match(exp17_names, colnames(exp17_data))] <- names(exp17_names)

### Getting the data in the right format

# removing the gpr37l1 columns
df <- exp17_data %>%
    as_tibble() %>%
    select(-c(
      "gpr37l1_ip" ,
      "gpr37l1_igg"
      ) # remove the gpr37l1 columns
    ) %>%
    mutate(
        Genes = if_else(is.na(Genes), Protein.Group, Genes) # If the Genes column has NA values, replace them with values from Group.Protein
    ) %>%
    column_to_rownames(var = "Genes") %>% # The Genes column as rownames
    select(-Protein.Group, -Protein.Names, -First.Protein.Description) %>%  # remove these columns
    filter(rowSums(is.na(.)) < ncol(.)) # removes the rows if they have NAs in all columns

### log2 transformation
df_log2 <- df %>%
  mutate(across(where(is.numeric), log2))

### imputation
data_matrix <- as.matrix(df_log2)
imputed_matrix <- impute.MinProb(dataSet.mvs = data_matrix,
                                      q = 0.01,      
                                      tune.sigma = 1) 
imputed_df <- as.data.frame(imputed_matrix)

# checking the imputed values for a given protein
row_values <- imputed_df["Prrt1", ]
print(row_values)

### Doing ROTS

# Group design
groups = c(0,1,0,1,0,1,0,1)
groups

# Running ROTS
results = ROTS(data = imputed_df, groups = groups, B = 1000, K = 5000, seed = 1234)
summary(results, fdr = 0.05)

# plotting the results
plot(results, fdr = 0.05, type = "volcano")

plot(results, fdr = 0.05, type = "heatmap")
