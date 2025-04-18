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
raw_data <- read_tsv("/mnt/ip-proteomics/exp19-my-diann-run/exp19-diann_output.pg_matrix.tsv")
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
colnames(data_na_filtered)[match(gpr37l1_col_names, colnames(data_na_filtered))] <- names(gpr37l1_col_names)

##### further data wrangling
df <- data_na_filtered %>%
        as_tibble() %>%
        column_to_rownames(var = "Genes") %>% # The Genes column as rownames
        select(-Protein.Group, -Protein.Names, -First.Protein.Description)


##### log2 transformation
df_log2 <- df %>%
  mutate(across(where(is.numeric), log2))

##### imputation
data_matrix <- as.matrix(df_log2)
imputed_matrix <- impute.MinProb(dataSet.mvs = data_matrix,
                                      q = 0.01,      
                                      tune.sigma = 1) 
imputed_df <- as.data.frame(imputed_matrix)

##### normalization
# df_norm <- as.data.frame(scale(imputed_matrix, center = TRUE, scale = TRUE))

# checking the imputed values for a given protein
row_values <- imputed_df["Gpr37l1", ]
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
    axis.text.x = element_text(angle = 45, hjust = 1), # Rotate x-axis labels
    legend.position = "right" # Adjust legend position as needed
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
    axis.text.x = element_text(angle = 45, hjust = 1),
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
  #geom_text(hjust = -0.1, vjust = 0.1, size = 3) + # Add sample labels
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
    axis.ticks = element_line(color = "black", linewidth = 0.75, size = 1) 
  )


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



##### Limma volcano plot 

# # to find the p value cutoff for the volcano plot
# adj.P.Val_cutoff <- 0.05
# sorted_results_limma <- results_limma[order(results_limma$P.Val), ]

# limma_rank_cutoff <- which(sorted_results_limma$adj.P.Val <= adj.P.Val_cutoff)
# limma_dynamic_p_cutoff <- sorted_results_limma$P.Val[max(limma_rank_cutoff)]


# load in annotations
annotations_data <- read_csv("gpr37l1_limma_ip_proteins_annotations.csv")

# create custom key-value pairs for 'high', 'low', 'mid' expression by fold-change
# this can be achieved with nested ifelse statements
keyvals <- ifelse(
results_limma$logFC > 1 & results_limma$adj.P.Val < 0.05, '#21a0e2',
    ifelse(results_limma$logFC < 1, '#c6c6c6',
    '#c6c6c6'))
keyvals[is.na(keyvals)] <- '#c6c6c6'
names(keyvals)[keyvals == '#21a0e2'] <- 'ip'
names(keyvals)[keyvals == '#c6c6c6'] <- 'mid'
names(keyvals)[keyvals == '#c6c6c6'] <- 'low'


# # Create a shape vector based on annotations
# keyvals.shape <- rep(21, nrow(annotations_data)) # Default to open circle
# names(keyvals.shape) <- rownames(annotations_data)


# # Assuming your TRUE/FALSE column in annotations_data is named 'is_target'
# if ("uniprot_membrane_related" %in% colnames(annotations_data)) {
#   keyvals.shape[annotations_data$uniprot_membrane_related == TRUE & !is.na(annotations_data$uniprot_membrane_related)] <- 19 # Solid circle if TRUE
# } else {
#   warning("Column 'is_target' not found in annotations_data. Using default shapes.")
# }





keyvals.shape <- ifelse(
results_limma$logFC > 1 & results_limma$adj.P.Val < 0.05, 19,
      ifelse(results_limma$logFC < 1, 21,
        21))
  keyvals.shape[is.na(keyvals.shape)] <- 21
  names(keyvals.shape)[keyvals.shape == 19] <- 'pp'
  names(keyvals.shape)[keyvals.shape == 21] <- 'mid'
  names(keyvals.shape)[keyvals.shape == 21] <- 'low'

select_proteins <- as.character(rownames(results_limma[results_limma$logFC > 1 & results_limma$adj.P.Val < 0.05, ]))
select_proteins <- c("Gpr37l1")
EnhancedVolcano(results_limma,
    lab = rownames(results_limma),
    x = 'logFC',
    y = 'P.Value',
    selectLab = select_proteins,
    #selectLab = rownames(results_limma)[which(names(keyvals) %in% c('ip'))],
    colCustom = keyvals,
    #shapeCustom = keyvals.shape[rownames(annotations_data)],
    drawConnectors = TRUE,
    xlim = c(-10, 10),
    ylim = c(-0.5, 6),
    pCutoff = 0.0000001,
    cutoffLineType = 'blank',
    FCcutoff = 1,
    colAlpha = 1,
    shape = 21,
    pointSize = 4.0,
    labSize = 6.0,
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
results = ROTS(data = df_norm, groups = groups, B = 1000, K = 5000, seed = 1234)
summary(results, fdr = 0.1)

# plotting the results
plot(results, fdr = 0.05, type = "volcano")

# Getting the data in dataframe
logFC <- results$logfc
pvalue <- results$pvalue
adj.pvalue <- results$FDR

rots_results <- as.data.frame(cbind(row.names(df_norm), logFC, pvalue, adj.pvalue))
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

# to find the p value cutoff for the volcano plot
adj.pvalue_cutoff <- 0.05

sorted_rots_results <- rots_results[order(rots_results$pval), ]
n_sorted_rots <- nrow(sorted_rots_results)

rots_rank_cutoff <- which(sorted_rots_results$adj.pvalue <= adj.pvalue_cutoff)
rots_dynamic_p_cutoff <- sorted_rots_results$pvalue[max(rots_rank_cutoff)]

# create custom key-value pairs for 'high', 'low', 'mid' expression by fold-change
# this can be achieved with nested ifelse statements
keyvals <- ifelse(
rots_results$logFC > 1 & rots_results$pvalue < 0.05, '#21a0e2',
    ifelse(rots_results$logFC < 1, '#c6c6c6',
    '#c6c6c6'))
keyvals[is.na(keyvals)] <- '#c6c6c6'
names(keyvals)[keyvals == '#21a0e2'] <- 'ip'
names(keyvals)[keyvals == '#c6c6c6'] <- 'mid'
names(keyvals)[keyvals == '#c6c6c6'] <- 'low'

keyvals.shape <- ifelse(
rots_results$logFC > 1 & rots_results$adj.pvalue < 0.05, 19,
      ifelse(rots_results$logFC < 1, 21,
        21))
  keyvals.shape[is.na(keyvals.shape)] <- 21
  names(keyvals.shape)[keyvals.shape == 19] <- 'pp'
  names(keyvals.shape)[keyvals.shape == 21] <- 'mid'
  names(keyvals.shape)[keyvals.shape == 21] <- 'low'

select_proteins <- as.character(rownames(rots_results[results$logFC > 1 & rots_results$adj.P.Val < 0.05, ]))
select_proteins <- c("Gpr37l1")
EnhancedVolcano(rots_results,
    lab = rownames(rots_results),
    x = 'logFC',
    y = 'pvalue',
    selectLab = select_proteins,
    #selectLab = rownames(rots_results)[which(names(keyvals) %in% c('ip'))],
    colCustom = keyvals,
    shapeCustom = keyvals.shape,
    drawConnectors = TRUE,
    xlim = c(-5, 5),
    ylim = c(-0.2, 5),
    pCutoff = rots_dynamic_p_cutoff,
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


##### DEA proDA

# needs to be in matrix format
log2_matrix <- as.matrix(df_log2)

# plots from proDA package
barplot(colSums(is.na(log2_matrix)),
        ylab = "# missing values")

boxplot(normalized_matrix,
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
    normalized_matrix, 
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

# to find the p value cutoff for the volcano plot
adj_pval_cutoff <- 0.05

sorted_proda_data <- results_proda[order(results_proda$pval), ]
n_sorted_proda <- nrow(sorted_proda_data)

rank_cutoff <- which(sorted_proda_data$adj_pval <= adj_pval_cutoff)
dynamic_p_cutoff <- sorted_proda_data$pval[max(rank_cutoff)]

# create custom key-value pairs for 'high', 'low', 'mid' expression by fold-change
# this can be achieved with nested ifelse statements
keyvals <- ifelse(
results_proda$diff > 1 & results_proda$pval < 0.05, '#21a0e2',
    ifelse(results_proda$diff < 1, '#c6c6c6',
    '#c6c6c6'))
keyvals[is.na(keyvals)] <- '#c6c6c6'
names(keyvals)[keyvals == '#21a0e2'] <- 'ip'
names(keyvals)[keyvals == '#c6c6c6'] <- 'mid'
names(keyvals)[keyvals == '#c6c6c6'] <- 'low'

keyvals.shape <- ifelse(
results_proda$diff > 1 & results_proda$adj_pval < 0.05, 19,
      ifelse(results_proda$diff < 1, 21,
        21))
  keyvals.shape[is.na(keyvals.shape)] <- 21
  names(keyvals.shape)[keyvals.shape == 19] <- 'pp'
  names(keyvals.shape)[keyvals.shape == 21] <- 'mid'
  names(keyvals.shape)[keyvals.shape == 21] <- 'low'

select_proteins <- as.character(rownames(results_proda[results$diff > 1 & results_proda$adj_pval < 0.05, ]))
select_proteins <- c("Gpr37l1")
EnhancedVolcano(results_proda,
    lab = rownames(results_proda),
    x = 'diff',
    y = 'pval',
    selectLab = select_proteins,
    #selectLab = rownames(results_proda)[which(names(keyvals) %in% c('ip'))],
    colCustom = keyvals,
    shapeCustom = keyvals.shape,
    drawConnectors = TRUE,
    xlim = c(-10, 10),
    ylim = c(-0.5, 10),
    pCutoff = dynamic_p_cutoff,
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
  "Gpr37l1",
  "Pfla4", 
  "Rasgef1a", 
  "Frmpd3", 
  "Brd3os"
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
  "Gpr37l1",
  "Pfla4", 
  "Rasgef1a", 
  "Frmpd3", 
  "Brd3os"
)

# Create an empty list to store results
results_list <- list()

# Process each protein
for(protein in protein_list) {
  result <- extract_protein_data(df_proteins, protein)
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
  scale_y_continuous(expand=c(0, 0), limits=c(0, 25)) + 
  theme_minimal(base_size = 22) +
  theme(
    panel.background = element_blank(),
    panel.grid = element_blank(),
    axis.line = element_line(size = 1.2),
    axis.ticks = element_line(size = 1.2),
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  labs(title="Protein Expression Levels",
       x="Condition",
       y="Log2 Intensity")
