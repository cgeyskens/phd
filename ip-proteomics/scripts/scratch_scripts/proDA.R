library(DEP)
library(dplyr)
library(SummarizedExperiment)
library(proDA)
library(readr)

# loading the data
exp17_data <- read_tsv("/mnt/ip-proteomics/exp17-my-diann-run/diann_output.pg_matrix.tsv")
View(exp17_data)


# checking the raw value of certain proteins
vcam1_data <- exp17_data[exp17_data$Protein.Names == 'GSLG1_MOUSE',]
View(vcam1_data)

# checking colnames
colnames(exp17_data)

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

# remove gpr37l1 data
df <- exp17_data %>%
    select(-c(
        "gpr37l1_ip",
        "gpr37l1_igg"
      )
    )

intensity_colnames <- grep("vcam1", colnames(df), value=TRUE)
head(intensity_colnames)

###### Create the intensity matrix #####
abundance_matrix <- as.matrix(df[, intensity_colnames])

rownames(abundance_matrix) <- df$Protein.Names


##### remove rows with NA as rowname #####
# Get the row names
row_names <- rownames(abundance_matrix)

# Identify row names that start with "NA"
rows_to_filter <- is.na(row_names)

# Filter out rows where the row name starts with "NA"
filtered_matrix <- abundance_matrix[!rows_to_filter, , drop = FALSE]

##### remove rows where all columns have NA #####
# remove rows where all columns have NA (from the gpr37l1 conditions)
rows_all_na <- apply(filtered_matrix, 1, function(row) all(is.na(row)))
# Filter out rows where all elements are NA
filtered_matrix <- filtered_matrix[!rows_all_na, , drop = FALSE]



###### log2 normalization #####
abundance_matrix <- log2(filtered_matrix)

barplot(colSums(is.na(abundance_matrix)),
        ylab = "# missing values")

boxplot(abundance_matrix,
        ylab = "Intensity Distribution")


##### normalization #####
normalized_abundance_matrix <- median_normalization(abundance_matrix)

##### getting the sample metadata in a dataframes #####
sample_info_df <- data.frame(name = colnames(filtered_matrix),
                             stringsAsFactors = FALSE)
sample_info_df$condition <- substr(sample_info_df$name, 7, nchar(sample_info_df$name)  - 3)
sample_info_df$replicate <- as.numeric(substr(sample_info_df$name, nchar(sample_info_df$name) - 0, nchar(sample_info_df$name)))
sample_info_df


##### Fitting the proDA model
fit <- proDA(
    normalized_abundance_matrix, 
    design = ~ condition + replicate, 
    data_is_log_transformed = TRUE,
    col_data = sample_info_df, 
    reference_level = "igg"
    )
fit

result_names((fit))
test <- proDA::test_diff(
    fit, 
    contrast = conditionip,
    pval_adjust_method = "fdr"
    )

proda_ip_proteins <- subset(test, pval < 0.05 & diff > 1)

write.csv(proda_ip_proteins,"vcam1_proda_ip_proteins.csv", row.names = TRUE)



#### Do the differential testing in another way, using the summarized experiment

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

view(se$assays)

# fitting proDA
fit <- proDA(se, design = ~ condition)
fit
# DEA
result_names((fit))
test <- proDA::test_diff(
    fit, 
    contrast = conditionip,
    #reduced_model = ~1,
    pval_adjust_method = "BH",
    sort_by = "pval"
    )
head(test)

check <- subset(test, pval < 0.05 & diff >= 1)

result_names((fit))
