library(MSstats)
library(ggplot2)
library(MSstatsConvert)
library(data.table)
library(tidyverse)
library(imputeLCMD)

# setting the path to the main DIANN report
input_file_path = "/mnt/data/ip-proteomics/exp17-my-diann-run/diann_output.tsv"

# constructing the annotation file for experimental design
Run <- c(
    "2ul_CG_3",
    "2ul_CG_4",
    "2ul_CG_5",
    "2ul_CG_6",
    "2ul_CG_7",
    "2ul_CG_8",
    "2ul_CG_9",
    "2ul_CG_10"
)
Condition <- c(
    "VCAM1-IP",
    "VCAM1-IgG",
    "VCAM1-IP",
    "VCAM1-IgG",
    "VCAM1-IP",
    "VCAM1-IgG",
    "VCAM1-IP",
    "VCAM1-IgG"
)

BioReplicate <- c(
    "1",
    "1",
    "2",
    "2",
    "3",
    "3",
    "4",
    "4"
)

df <- data.frame(
    Run = Run,
    Condition = Condition,
    BioReplicate = BioReplicate
)

getwd()

current_dir <- getwd()
output_file <- file.path(current_dir, "annotations_diann.csv")
write.csv(df, file = output_file, row.names = FALSE)

annotation_file_path = "/workspaces/phd/ip-proteomics/annotations_diann.csv"

# importing the diann data
input = data.table::fread(input_file_path)
annot = data.table::fread(annotation_file_path)

# need to change column name because of an issue: see https://groups.google.com/g/msstats/c/5EwcYsqTSc0/m/k-F8D-iVAgAJ
setnames(input, old = "Fragment.Quant.Raw", new = "FragmentQuantCorrected")
check <- head(input)

# remove each row with gpr37l1 data
input_adj <- input[Run != "2ul_CG_1" & Run != "2ul_CG_2"]
check <- head(input_adj)
length(unique(input_adj$Protein.Group))

?DIANNtoMSstatsFormat

raw_data = DIANNtoMSstatsFormat(
    input_adj, 
    annotation = annot, 
    global_qvalue_cutoff = 0.01,
    qvalue_cutoff = 0.01,
    pg_qvalue_cutoff = 0.01,
    useUniquePeptide = FALSE,
    removeFewMeasurements = FALSE,
    removeOxidationMpeptides = FALSE,
    removeProtein_with1Feature = FALSE,
    use_log_file = FALSE,
    append = FALSE,
    verbose = TRUE,
    log_file_path = NULL,
    MBR = TRUE
)
check <- head(raw_data)

length(unique(raw_data$ProteinName))

check <- raw_data %>%
  as_tibble() %>%
  head()

?dataProcess

# Processing the data
summarized <- dataProcess(
    raw_data,
    logTrans = 2,
    normalization = FALSE,
    featureSubset = "topN",
    min_feature_count = 1,
    n_top_feature = 200,
    summaryMethod = "TMP",
    equalFeatureVar = TRUE,
    censoredInt = "NA",
    MBimpute = TRUE, # imputation on peptide level
    verbose = TRUE, 
    log_file_path = NULL
)


# QC plots
dataProcessPlots(
    data=summarized, 
    type="ProfilePlot", 
    address = FALSE, 
    which.Protein = "VCAM1_MOUSE"
)

check <- summarized$FeatureLevel[summarized$FeatureLevel$Protein == "VCAM1_MOUSE", ]

check <- summarized$ProteinLevelData

summarized$FeatureLevelData %>%
  filter(!censored) %>%
  ggplot(aes(x=originalRUN, y=newABUNDANCE, fill=GROUP)) +
  geom_boxplot() +
  coord_flip()

summarized$ProteinLevelData %>%
  ggplot(aes(x=originalRUN, y=LogIntensities, fill=GROUP)) +
  geom_boxplot() +
  coord_flip() 

summarized$ProteinLevelData %>%
  ggplot(aes(x=LogIntensities, colour=GROUP, group=originalRUN)) +
  geom_density(linewidth=1) 


any(is.na(summarized$ProteinLevelData$LogIntensities)) 


## PCA clustering of samples
wide_data <- summarized$ProteinLevelData %>%
  select(Protein, originalRUN, LogIntensities) %>%
  group_by(Protein) %>%
  ungroup() %>%
  pivot_wider(
    names_from = originalRUN,
    values_from = LogIntensities
  ) %>%
  column_to_rownames(var = "Protein")

# Check for missing values in the wide data
any(is.na(wide_data))
sum(is.na(wide_data))

# imputation on protein level
data_matrix <- as.matrix(wide_data)
imputed_matrix <- impute.MinProb(dataSet.mvs = data_matrix,
                                      q = 0.01,      
                                      tune.sigma = 1) 
imputed_df <- as.data.frame(imputed_matrix)

imputed_df_vcam1 <- filter(imputed_df, rownames(imputed_df) == "VCAM1_MOUSE") 

?prcomp

pca_results <- imputed_df %>%
  t() %>%
  prcomp(scale=FALSE, center = TRUE)

view(pca_results$x)

pca_results$x %>%
  # as_tibble(rownames="Raw.file") %>%
  # left_join(annot) %>%
  ggplot(aes(x=PC1, y=PC2, fill=Condition)) +
  geom_point(pch=21, size=6) +
  ggtitle("PCA plot of all samples")

