library(limma)
library(readr)
library(dplyr)
library(tibble)
library(imputeLCMD)

# Loading the data
exp17_data <- read_tsv("/mnt/ip-proteomics/exp17-my-diann-run/diann_output.pg_matrix.tsv")
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
                                      q = 0.00001,      
                                      tune.sigma = 1) 
imputed_df <- as.data.frame(imputed_matrix)

# checking the imputed values for a given protein
row_values <- imputed_df["Vcam1", ]
print(row_values)


##### Limma #####
# Create factors for your experimental design
condition <- factor(c(
    "ip", 
    "igg", 
    "ip", 
    "igg",
    "ip", 
    "igg",
    "ip", 
    "igg"
    )
)
replicate <- factor(c(
    "rep_1", 
    "rep_1", 
    "rep_2", 
    "rep_2",
    "rep_3",
    "rep_3", 
    "rep_4", 
    "rep_4"
    )
)
design <- model.matrix(~replicate + condition)
colnames(design)
# Fit the linear model
fit1 <- limma::lmFit(imputed_df, design)
cont <- makeContrasts(conditionip, levels = design)
fit2 <- contrasts.fit(fit1, contrasts = cont)
fit3 <- eBayes(fit2)

results <- topTable(fit3, adjust = "fdr", sort.by = "P", n = 3070)
limma_ip_proteins <- results[results$P.Val < 0.05 & results$logFC > 1, ]

write.csv(dep_ip_proteins,"vcam1_limma_ip_proteins.csv", row.names = TRUE)



### Enhanced volcano plot
library(EnhancedVolcano)


  # create custom key-value pairs for 'high', 'low', 'mid' expression by fold-change
  # this can be achieved with nested ifelse statements
  keyvals <- ifelse(
    results$logFC > 1 & results$P.Val < 0.05, '#21a0e2',
      ifelse(results$logFC < 1, '#c6c6c6',
        '#c6c6c6'))
  keyvals[is.na(keyvals)] <- '#c6c6c6'
  names(keyvals)[keyvals == '#21a0e2'] <- 'ip'
  names(keyvals)[keyvals == '#c6c6c6'] <- 'mid'
  names(keyvals)[keyvals == '#c6c6c6'] <- 'low'

keyvals.shape <- ifelse(
results$logFC > 1 & results$adj.P.Val < 0.05, 19,
      ifelse(results$logFC < 1, 21,
        21))
  keyvals.shape[is.na(keyvals.shape)] <- 21
  names(keyvals.shape)[keyvals.shape == 19] <- 'pp'
  names(keyvals.shape)[keyvals.shape == 21] <- 'mid'
  names(keyvals.shape)[keyvals.shape == 21] <- 'low'

select_proteins <- as.character(rownames(results[results$logFC > 1 & results$adj.P.Val < 0.05, ]))
select_proteins <- c("Vcam1")
EnhancedVolcano(results,
    lab = rownames(results),
    x = 'logFC',
    y = 'P.Value',
    selectLab = select_proteins,
    #selectLab = rownames(results)[which(names(keyvals) %in% c('ip'))],
    colCustom = keyvals,
    shapeCustom = keyvals.shape,
    drawConnectors = TRUE,
    xlim = c(-10, 10),
    ylim = c(-0.5, 6),
    pCutoff = 0.0007,
    cutoffLineType = 'blank',
    FCcutoff = 1,
    colAlpha = 1,
    shape = 21,
    pointSize = 4.0,
    labSize = 4.0,
    gridlines.major = FALSE,
    gridlines.minor = FALSE,
    arrowheads = FALSE,
    parseLabels = TRUE,
    lengthConnectors = unit(0.1, 'npc'),
    ) + theme(text=element_text(size=4,  family="Arial"))
