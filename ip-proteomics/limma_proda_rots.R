library(readr)
library(dplyr)
library(tibble)
library(ggplot2)
library(tidyr)
library(stringr)
library(UpSetR)
library(EnhancedVolcano)

library(limma)
library(imputeLCMD)
library(ROTS)
library(proDA)

##### Loading the data
raw_data <- read_tsv("/mnt/ip-proteomics/exp10_vcam1_synglio_replicates.pg_matrix.tsv")
View(raw_data)

##### Filter out the antibodies fragments: starting with "Ig"
antibodies_condition1 <- startsWith(raw_data$Protein.Group, "A0A")
antibodies_condition2 <- startsWith(raw_data$Genes, "Ig")
antibodies_to_filter_out <- antibodies_condition1 & antibodies_condition2

data_ab_filtered <- raw_data[!antibodies_to_filter_out, ]
data_ab_filtered_out <- raw_data[antibodies_to_filter_out, ] # to check the filtered out proteins

##### Filter out proteins with no genes charactar (eg they are human or otherwise)
gene_na_rows <- is.na(data_ab_filtered$Genes)

data_na_filtered <- data_ab_filtered[!gene_na_rows, ]
data_na_filtered_out <- data_ab_filtered[gene_na_rows, ]

##### Renaming the columns

# checking colnames
colnames(data_na_filtered)
# new colnames
gpr37l1_col_names <- c(
    "gpr37l1_ip_r1" = "/ip-proteomics/exp19-ms-convert-output/CG_01.mzML",
    "gpr37l1_igg_r1"= "/ip-proteomics/exp19-ms-convert-output/CG_02.mzML",
    "gpr37l1_ip_r2" = "/ip-proteomics/exp19-ms-convert-output/CG_03.mzML",
    "gpr37l1_igg_r2"= "/ip-proteomics/exp19-ms-convert-output/CG_04.mzML",
    "gpr37l1_ip_r3" = "/ip-proteomics/exp19-ms-convert-output/CG_05.mzML",
    "gpr37l1_igg_r3"= "/ip-proteomics/exp19-ms-convert-output/CG_06.mzML",
    "gpr37l1_ip_r4" = "/ip-proteomics/exp19-ms-convert-output/CG_07.mzML",
    "gpr37l1_igg_r4"= "/ip-proteomics/exp19-ms-convert-output/CG_08.mzML"
)

vcam1_exp10_synglio_col_names <- c(
  "p14_vcam1_ip_r1" = "D:\\Proteomics2024\\CBD\\JorisdeWitLab\\PCF000184\\CG_1.wiff",
  "p14_vcam1_igg_r1" = "D:\\Proteomics2024\\CBD\\JorisdeWitLab\\PCF000184\\CG_2.wiff",
  "p14_vcam1_ip_r2" = "D:\\Proteomics2024\\CBD\\JorisdeWitLab\\PCF000184\\CG_3.wiff",
  "p14_vcam1_igg_r2" = "D:\\Proteomics2024\\CBD\\JorisdeWitLab\\PCF000184\\CG_4.wiff",
  "p14_vcam1_ip_r3" = "D:\\Proteomics2024\\CBD\\JorisdeWitLab\\PCF000184\\CG_5.wiff",
  "p14_vcam1_igg_r3" = "D:\\Proteomics2024\\CBD\\JorisdeWitLab\\PCF000184\\CG_6.wiff",
  "p14_vcam1_ip_r4" = "D:\\Proteomics2024\\CBD\\JorisdeWitLab\\PCF000184\\CG_7.wiff",
  "p14_vcam1_igg_r4" = "D:\\Proteomics2024\\CBD\\JorisdeWitLab\\PCF000184\\CG_8.wiff",
  "p28_vcam1_ip_r1" = "D:\\Proteomics2024\\CBD\\JorisdeWitLab\\PCF000184\\CG_9.wiff",
  "p28_vcam1_igg_r1" = "D:\\Proteomics2024\\CBD\\JorisdeWitLab\\PCF000184\\CG_10.wiff",
  "p28_vcam1_ip_r2" = "D:\\Proteomics2024\\CBD\\JorisdeWitLab\\PCF000184\\CG_11.wiff",
  "p28_vcam1_igg_r2" = "D:\\Proteomics2024\\CBD\\JorisdeWitLab\\PCF000184\\CG_12.wiff",
  "p28_vcam1_ip_r3" = "D:\\Proteomics2024\\CBD\\JorisdeWitLab\\PCF000184\\CG_13.wiff",
  "p28_vcam1_igg_r3" = "D:\\Proteomics2024\\CBD\\JorisdeWitLab\\PCF000184\\CG_14.wiff",
  "p28_vcam1_ip_r4" = "D:\\Proteomics2024\\CBD\\JorisdeWitLab\\PCF000184\\CG_15.wiff",
  "p28_vcam1_igg_r4" = "D:\\Proteomics2024\\CBD\\JorisdeWitLab\\PCF000184\\CG_16.wiff"
)

colnames(data_na_filtered)[match(vcam1_exp10_synglio_col_names, colnames(data_na_filtered))] <- names(vcam1_exp10_synglio_col_names)

##### further data wrangling
df <- data_na_filtered %>%
        as_tibble() %>%
        column_to_rownames(var = "Genes") %>% # The Genes column as rownames
        select(-Protein.Group, -Protein.Names, -First.Protein.Description)
check <- df["Ntng1", ]

##### log2 transformation
df_log2 <- df %>%
  mutate(across(where(is.numeric), log2))

##### imputation
data_matrix <- as.matrix(df_log2)
imputed_matrix <- impute.MinProb(dataSet.mvs = data_matrix,
                                      q = 0.001,      
                                      tune.sigma = 1) 
imputed_df <- as.data.frame(imputed_matrix)

##### normalization
# df_norm <- as.data.frame(scale(imputed_matrix, center = TRUE, scale = TRUE))

# checking the imputed values for a given protein
row_values <- imputed_df["Vcam1", ]
print(row_values)


##### QC plots: nr of protein identified per sample

plot_df <- df %>%
  # Count the number of non-NA entries (identified proteins) for each sample
  summarise(across(everything(), ~sum(!is.na(.)))) %>%
  # Transpose the data frame to have samples as a column
  pivot_longer(cols = everything(), names_to = "Sample", values_to = "UniqueProteins") %>%
  # Extract the condition (ip or igg) from the sample name
  mutate(Condition = str_extract(Sample, "(ip|igg)"))


condition_colors <- c("ip" = "#21a0e2", "igg" = "#c6c6c6") # You can choose your desired colors

# Create the bar plot using ggplot2
ggplot(plot_df, aes(x = Sample, y = UniqueProteins, fill = Condition)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = UniqueProteins),
            hjust = 0.5, # Center the text horizontally
            vjust = 2, # Center the text vertically
            color = "white",
            size = 5) + # Adjust text size if needed
  scale_fill_manual(values = condition_colors) + # Apply the defined colors
  labs(
    title = "Number of Unique Proteins Identified in Each Sample",
    x = "Sample",
    y = "Number of Unique Proteins"
  ) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(), # Remove all grid lines
    axis.line = element_line(linewidth = 0.75, color = "black"), # Make axis lines bolder and black
    axis.text.x = element_text(angle = 45, hjust = 1, size = 16),# Rotate x-axis labels
    axis.text.y = element_text(size = 16),
    axis.title.y = element_text(size = 18, family = "Arial"), 
    axis.title.x = element_blank(), 
    plot.title = element_text(size = 20, family = "Arial"),
    legend.position = "right",
    text = element_text(family = "Arial")  # Adjust legend position as needed
  )

##### QC plots: log2 intensity plots

plot_df_long <- imputed_df %>%
  rownames_to_column(var = "Protein") %>%
  pivot_longer(cols = -Protein, names_to = "Sample", values_to = "Intensity") %>%
  # Extract the condition (ip or igg) from the sample name
  mutate(Condition = str_extract(Sample, "(ip|igg)"))

# Print the first few rows of the long format data frame
print(head(plot_df_long))

# Set the colors for the conditions
condition_colors <- c("ip" = "#21a0e2", "igg" = "#c6c6c6")# Choose your desired colors

# Create the box plot
ggplot(plot_df_long, aes(x = Sample, y = Intensity, fill = Condition)) +
  geom_boxplot() +
  scale_fill_manual(values = condition_colors) + # Apply the defined colors
  labs(
    title = "Protein Intensity Distribution per Sample",
    x = "Sample",
    y = "Log2 Intensity",
    fill = "Condition" # Set the legend title
  ) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(), # Remove all grid lines
    axis.line = element_line(linewidth = 0.75, color = "black"), # Make axis lines bolder and black
    axis.text.x = element_text(angle = 45, hjust = 1, size = 16),
    axis.text.y = element_text(size = 16),
    axis.title.y = element_text(size = 18, family = "Arial"), 
    axis.title.x = element_blank(), 
    plot.title = element_text(size = 20, family = "Arial"),
    legend.position = "right"
  ) + scale_y_continuous(limits = c(0, NA))

##### QC plots: PCA plot

# 1. Transpose the data frame so that samples are rows and proteins are columns
df_t <- as.data.frame(t(imputed_df))

# 2. Perform Principal Component Analysis (PCA)
pca_result <- prcomp(df_t, scale. = TRUE) # scale. = TRUE for scaling the data

# 3. Extract the principal components and their variances
pca_scores <- as.data.frame(pca_result$x)
pca_variance <- pca_result$sdev^2
pca_percentage_variance <- round(100 * pca_variance / sum(pca_variance), 1)

# 4. Create a data frame for plotting
pca_plot_df <- pca_scores %>%
  # Add the sample names as a column
  mutate(Sample = rownames(df_t)) %>%
  # Extract the condition (ip or igg) from the sample name
  mutate(Condition = stringr::str_extract(Sample, "(ip|igg)"))

# Print the first few rows of the PCA plotting data frame
print(head(pca_plot_df))

# 5. Define colors for the conditions
condition_colors <- c("ip" = "#21a0e2", "igg" = "#c6c6c6")

# 6. Create the PCA scatter plot
ggplot(pca_plot_df, aes(x = PC1, y = PC2, color = Condition, label = Sample)) +
  geom_point(size = 6) +
  #geom_text(hjust = -0.1, vjust = 0.1, size = 5) + # Add sample labels
  scale_color_manual(values = condition_colors) +
  labs(
    title = "PCA Plot of Samples",
    x = paste0("PC1 (", pca_percentage_variance[1], "% variance explained)"),
    y = paste0("PC2 (", pca_percentage_variance[2], "% variance explained)")
  ) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    axis.line = element_line(linewidth = 0.75, color = "black"),
    legend.position = "right",
    axis.ticks = element_line(color = "black", linewidth = 0.75, size = 1),
    axis.text.y = element_text(size = 16),
    axis.title.y = element_text(size = 18, family = "Arial"), 
    axis.title.x = element_text(size = 18, family = "Arial"), 
    plot.title = element_text(size = 20, family = "Arial"), 
  ) + 
  xlim(-30, 30)


##### DEA Limma #####

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

results_limma <- topTable(fit3, adjust = "fdr", sort.by = "P", n = 3070)
limma_ip_proteins <- results_limma[results_limma$adj.P.Val < 0.05 & results_limma$logFC > 1, ]

write.csv(limma_ip_proteins,"gpr37l1_limma_ip_proteins.csv", row.names = TRUE)


#### Limma volcano plot without annotations

results_limma <- tibble::rownames_to_column(results_limma, "Genes")

# making the labels for volcano plot
labels <- results_limma$Genes
labels_upper <- paste0(toupper(labels))

select_gpr37l1 <- c("Gpr37l1")
select_proteins_upper <- paste0(toupper(select_gpr37l1))

# making the statement color the points
keyvals <- ifelse(
    results_limma$logFC > 1 & results_limma$adj.P.Val < 0.05, '#21a0e2', '#c6c6c6')
keyvals[is.na(keyvals)] <- '#c6c6c6'
names(keyvals)[keyvals == '#21a0e2'] <- 'ip'
names(keyvals)[keyvals == '#c6c6c6'] <- 'unspecific'

# making the volcano plot
EnhancedVolcano(results_limma,
    lab = labels_upper,
    x = 'logFC',
    y = 'P.Value',
    selectLab = select_proteins_upper,
    colCustom = keyvals,
    #shapeCustom = keyvals.shape,
    drawConnectors = TRUE,
    xlim = c(-8, 8),
    ylim = c(-0.5, 6),
    cutoffLineType = 'blank',
    colAlpha = 1,
    shape = 21,
    pointSize = 4.0,
    labSize = 6,
    gridlines.major = FALSE,
    gridlines.minor = FALSE,
    arrowheads = FALSE,
    parseLabels = TRUE,
    lengthConnectors = unit(0.1, 'npc'),
    ) + theme(text=element_text(size=4,  family="Arial"))



### !!! Now go to the annotations.R script to get the annotations

##### Limma volcano plot with annotations

# load in annotations
annotations_data <- read_csv("gpr37l1_limma_ip_proteins_annotations.csv")

# merge annotations with the results_limma
# results_limma <- tibble::rownames_to_column(results_limma, "Genes")
results_limma <- tibble::as_tibble(results_limma)
merged_df_for_volcano <- left_join(results_limma, 
                          select(annotations_data, "Genes", 
                          "in_van_oostrum_2023",
                          "in_sorokina_2021",
                          "sorokina_2021_synaptosome",
                          "sorokina_2021_presynaptic",
                          "sorokina_2021_postsynaptic",
                          "in_syngo",
                          "uniprot_membrane_related"))


# making the statement color the points
keyvals <- ifelse(
    merged_df_for_volcano$logFC > 1 & merged_df_for_volcano$adj.P.Val < 0.05, '#21a0e2', '#c6c6c6')
keyvals[is.na(keyvals)] <- '#c6c6c6'
names(keyvals)[keyvals == '#21a0e2'] <- 'ip'
names(keyvals)[keyvals == '#c6c6c6'] <- 'unspecific'

# making the statement to shape the points 
keyvals.shape <- ifelse(
      merged_df_for_volcano$uniprot_membrane_related == TRUE & merged_df_for_volcano$in_syngo ==TRUE, 
      19, 
      21)
keyvals.shape[is.na(keyvals.shape)] <- 21
names(keyvals.shape)[keyvals.shape == 19] <- 'yes'
names(keyvals.shape)[keyvals.shape == 21] <- 'low'

# checking how many labels fit the statements
table(merged_df_for_volcano$uniprot_membrane_related, merged_df_for_volcano$in_syngo)

# making the labels for volcano plot
labels <- merged_df_for_volcano$Genes
labels_upper <- paste0(toupper(labels))

select_proteins_with_na <- as.character(merged_df_for_volcano$Genes[merged_df_for_volcano$uniprot_membrane_related == TRUE & merged_df_for_volcano$in_syngo ==TRUE])
select_proteins <- select_proteins_with_na[!is.na(select_proteins_with_na)]
select_proteins_upper <- paste0(toupper(select_proteins))


# making the volcano plot
EnhancedVolcano(merged_df_for_volcano,
    lab = labels_upper,
    x = 'logFC',
    y = 'P.Value',
    selectLab = select_proteins_upper,
    colCustom = keyvals,
    shapeCustom = keyvals.shape,
    drawConnectors = TRUE,
    xlim = c(-8, 8),
    ylim = c(-0.5, 6),
    cutoffLineType = 'blank',
    colAlpha = 1,
    shape = 21,
    pointSize = 4.0,
    #parseLabels = TRUE,
    labSize = 6,
    gridlines.major = FALSE,
    gridlines.minor = FALSE,
    arrowheads = FALSE,
    parseLabels = TRUE,
    lengthConnectors = unit(0.1, 'npc'),
    ) + theme(text=element_text(size=4,  family="Arial"))


##### DEA ROTS 

# Group design
groups = c(0,1,0,1,0,1,0,1)
groups

# Running ROTS
results = ROTS(data = imputed_df, groups = groups, B = 1000, K = 5000, seed = 1234)
summary(results, fdr = 0.1)

# plotting the results
plot(results, fdr = 0.05, type = "volcano")

# Getting the data in dataframe
logFC <- results$logfc
pvalue <- results$pvalue
adj.pvalue <- results$FDR

rots_results <- as.data.frame(cbind(row.names(imputed_df), logFC, pvalue, adj.pvalue))
rots_results <- rots_results %>%
  mutate(
    logFC = as.numeric(logFC),
    pvalue = as.numeric(pvalue),
    adj.pvalue = as.numeric(adj.pvalue)
    # Add more columns as needed
  )

rots_ip_proteins <- rots_results %>%
                          subset(adj.pvalue <= 0.05 & logFC > 1)
write.csv(rots_ip_proteins,"gpr37l1_rots_ip_proteins.csv", row.names = TRUE)

##### ROTS volcano plot 

rots_results <- tibble::rownames_to_column(rots_results, "Genes")

# making the labels for volcano plot
labels <- rots_results$Genes
labels_upper <- paste0(toupper(labels))

select_gpr37l1 <- c("Gpr37l1")
select_proteins_upper <- paste0(toupper(select_gpr37l1))

# create custom key-value pairs for 'high', 'low', 'mid' expression by fold-change
# this can be achieved with nested ifelse statements
keyvals <- ifelse(
rots_results$logFC > 1 & rots_results$adj.pvalue < 0.05, '#21a0e2',
    ifelse(rots_results$logFC < 1, '#c6c6c6',
    '#c6c6c6'))
keyvals[is.na(keyvals)] <- '#c6c6c6'
names(keyvals)[keyvals == '#21a0e2'] <- 'ip'
names(keyvals)[keyvals == '#c6c6c6'] <- 'mid'
names(keyvals)[keyvals == '#c6c6c6'] <- 'low'


EnhancedVolcano(rots_results,
    lab = labels_upper,
    x = 'logFC',
    y = 'pvalue',
    selectLab = select_proteins_upper,
    colCustom = keyvals,
    #shapeCustom = keyvals.shape,
    drawConnectors = TRUE,
    xlim = c(-7, 7),
    ylim = c(-0.2, 5),
    cutoffLineType = 'blank',
    colAlpha = 1,
    shape = 21,
    pointSize = 4.0,
    labSize = 6,
    gridlines.major = FALSE,
    gridlines.minor = FALSE,
    arrowheads = FALSE,
    parseLabels = TRUE,
    lengthConnectors = unit(0.1, 'npc'),
    ) + theme(text=element_text(size=4,  family="Arial"))

check <- rots_results[rots_results$Genes == "Gpr37l1", ]

##### DEA proDA

# needs to be in matrix format
log2_matrix <- as.matrix(df_log2)

# plots from proDA package
barplot(colSums(is.na(log2_matrix)),
        ylab = "# missing values")

boxplot(log2_matrix,
        ylab = "Intensity Distribution")

#normalization 
normalized_matrix <- scale(log2_matrix)

# getting the sample metadata in a dataframes 
sample_info_df <- data.frame(name = colnames(log2_matrix),
                             stringsAsFactors = FALSE)
sample_info_df$condition <- substr(sample_info_df$name, 9, nchar(sample_info_df$name)  - 3)
sample_info_df$replicate <- as.numeric(substr(sample_info_df$name, nchar(sample_info_df$name) - 0, nchar(sample_info_df$name)))
sample_info_df

##### Fitting the proDA model
fit <- proDA(
    log2_matrix, 
    design = ~ condition + replicate, 
    data_is_log_transformed = TRUE,
    col_data = sample_info_df, 
    reference_level = "igg"
    )
fit

results_proda <- proDA::test_diff(
    fit, 
    contrast = conditionip,
    pval_adjust_method = "fdr"
    )
proda_ip_proteins <- subset(results_proda, adj_pval < 0.05 & diff > 1)



##### proDA volcano plot 

names(results_proda)[names(results_proda) == "name"] <- "Genes"
31431
# making the labels for volcano plot
labels <- results_proda$Genes
labels_upper <- paste0(toupper(labels))

select_gpr37l1 <- c("Gpr37l1")
select_proteins_upper <- paste0(toupper(select_gpr37l1))

# create custom key-value pairs for 'high', 'low', 'mid' expression by fold-change
# this can be achieved with nested ifelse statements
keyvals <- ifelse(
results_proda$diff > 1 & results_proda$adj_pval < 0.05, '#21a0e2',
    ifelse(results_proda$logFC < 1, '#c6c6c6',
    '#c6c6c6'))
keyvals[is.na(keyvals)] <- '#c6c6c6'
names(keyvals)[keyvals == '#21a0e2'] <- 'ip'
names(keyvals)[keyvals == '#c6c6c6'] <- 'mid'
names(keyvals)[keyvals == '#c6c6c6'] <- 'low'


EnhancedVolcano(results_proda,
    lab = labels_upper,
    x = 'diff',
    y = 'pval',
    selectLab = select_proteins_upper,
    colCustom = keyvals,
    drawConnectors = TRUE,
    xlim = c(-7, 7),
    ylim = c(-0.5, 7.5),
    cutoffLineType = 'blank',
    colAlpha = 1,
    shape = 21,
    pointSize = 4.0,
    labSize = 6,
    gridlines.major = FALSE,
    gridlines.minor = FALSE,
    arrowheads = FALSE,
    parseLabels = TRUE,
    lengthConnectors = unit(0.1, 'npc'),
    ) + theme(text=element_text(size=4,  family="Arial"))


##### Upsetplot of Limma, ROTS & proDA
limma_set <- rownames(limma_ip_proteins)
rots_set <- rownames(rots_ip_proteins)
proda_set <- rownames(proda_ip_proteins)

list_of_sets <- list(
  Limma = limma_set,
  ROTS = rots_set,
  proDA = proda_set
)

upset(fromList(list_of_sets),
    text.scale = 2)


##### Line plots of selected proteins

extract_protein_data <- function(df, protein_name) {
  # Check if protein exists in rownames
  if(!(protein_name %in% rownames(df))) {
    return(NULL)
  }
  
  # Extract the row for the specified protein
  protein_row <- df[protein_name, , drop = FALSE]
  
  # Convert to long format
  protein_data <- data.frame(
    condition = colnames(protein_row),
    raw_log2_intensities = as.numeric(as.matrix(protein_row)),
    protein = protein_name
  )
  
  return(protein_data)
}

# Assuming your dataframe is called 'df_proteins'
protein_list <- c(
  "Vcam1"
  # "Tenm3", 
  # "Tenm4", 
  # "Ntng1", 
  # "Vangl2",
  # "Rgs7",
  # "Gabrb1",
  # "Sst",
  # "Nf1",
  # "Trappc4",
  # "Stxbp5",
  # "Cdh9"
)

# Create an empty list to store results
results_list <- list()

# Process each protein
for(protein in protein_list) {
  result <- extract_protein_data(df_log2, protein)
  if(!is.null(result)) {
    results_list[[protein]] <- result
  }
}

# Combine all results
df_merge <- bind_rows(results_list)

# Plot with NA handling
ggplot(data=df_merge, aes(x=condition, y=raw_log2_intensities, group=protein)) +
  geom_line(aes(color=protein), na.rm = TRUE, size = 1) +  # Skip NA values when drawing lines
  geom_point(aes(color=protein), na.rm = TRUE, size = 5) +  # Skip NA values when plotting points
  scale_x_discrete() +
  scale_y_continuous(expand=c(0, 0), limits=c(0, 20)) + 
  theme_minimal(base_size = 22) +
  theme(
    panel.background = element_blank(),
    panel.grid = element_blank(),
    axis.line = element_line(size = 1.2),
    axis.ticks = element_line(size = 1.2),
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  labs(title="Raw Intensity Levels",
       x="Condition",
       y="Log2 Intensity") + theme(text=element_text(size=20,  family="Arial"))
