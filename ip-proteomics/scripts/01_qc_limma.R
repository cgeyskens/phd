###############################################################################
# Script: 01_qc_limma.R
# Purpose: formal analysis script to analyze the IP MS-DIA proteomics
# Author: Cydric Geyskens
# Date: 2025-10-24
###############################################################################


#### =============================== setup =============================== ####
library(readr)
library(dplyr)
library(tibble)
library(ggplot2)
library(tidyr)
library(stringr)
library(EnhancedVolcano)
library(svglite)

library(limma)
library(imputeLCMD)

# for data viewing
# remotes::install_github("nx10/httpgd")
library(httpgd)
hgd()



#### =============================== arguments =============================== ####
# only need to change these arguments for full analysis
# vcam1
input_data_filepath <- "/mnt/ip-proteomics/analyses-paper/exp17-zenotof-output-1miscleavage-wo-gpr37l1-20251007/exp17-zenotof-diann-output-1miscleavage-wo-gpr37l1.pg_matrix.tsv"
ip_protein = "Vcam1"
# gpr37l1
input_data_filepath <- "/mnt/ip-proteomics/analyses-paper/exp19-output-1miscleavage-20251008/exp19-diann-output-1miscleavage.pg_matrix.tsv"
ip_protein = "Gpr37l1"


#### =============================== loading the data =============================== ####
raw_data <- read_tsv(input_data_filepath)
View(raw_data)


#### =============================== filtering out contanimants =============================== ####

# filter out the antibodies fragments: starting with "Ig" & "A0A"
antibodies_condition1 <- startsWith(raw_data$Protein.Group, "A0A")
antibodies_condition2 <- startsWith(raw_data$Genes, "Ig")
antibodies_to_filter_out <- antibodies_condition1 & antibodies_condition2

data_ab_filtered <- raw_data[!antibodies_to_filter_out, ]
data_ab_filtered_out <- raw_data[antibodies_to_filter_out, ] # to check the filtered out proteins

# filter out proteins with no genes charactar (eg they are human or otherwise)
gene_na_rows <- is.na(data_ab_filtered$Genes)

data_na_filtered <- data_ab_filtered[!gene_na_rows, ]
data_na_filtered_out <- data_ab_filtered[gene_na_rows, ]


#### =============================== renaming the column names =============================== ####

# new colnames according to the experiments
exp17_vcam1_col_names <- c(
  "vcam1_ip_r1"= "/ip-proteomics/exp17-ms-zenotof-convert-output/2ul_CG_3.mzML",
  "vcam1_igg_r1"= "/ip-proteomics/exp17-ms-zenotof-convert-output/2ul_CG_4.mzML",
  "vcam1_ip_r2"= "/ip-proteomics/exp17-ms-zenotof-convert-output/2ul_CG_5.mzML",
  "vcam1_igg_r2"= "/ip-proteomics/exp17-ms-zenotof-convert-output/2ul_CG_6.mzML",
  "vcam1_ip_r3" = "/ip-proteomics/exp17-ms-zenotof-convert-output/2ul_CG_7.mzML",
  "vcam1_igg_r3"= "/ip-proteomics/exp17-ms-zenotof-convert-output/2ul_CG_8.mzML",
  "vcam1_ip_r4" = "/ip-proteomics/exp17-ms-zenotof-convert-output/2ul_CG_9.mzML",
  "vcam1_igg_r4"= "/ip-proteomics/exp17-ms-zenotof-convert-output/2ul_CG_10.mzML"
)

exp19_gpr37l1_col_names <- c(
    "gpr37l1_ip_r1" = "/ip-proteomics/exp19-ms-zenotof-convert-output/CG_01.mzML",
    "gpr37l1_igg_r1"= "/ip-proteomics/exp19-ms-zenotof-convert-output/CG_02.mzML",
    "gpr37l1_ip_r2" = "/ip-proteomics/exp19-ms-zenotof-convert-output/CG_03.mzML",
    "gpr37l1_igg_r2"= "/ip-proteomics/exp19-ms-zenotof-convert-output/CG_04.mzML",
    "gpr37l1_ip_r3" = "/ip-proteomics/exp19-ms-zenotof-convert-output/CG_05.mzML",
    "gpr37l1_igg_r3"= "/ip-proteomics/exp19-ms-zenotof-convert-output/CG_06.mzML",
    "gpr37l1_ip_r4" = "/ip-proteomics/exp19-ms-zenotof-convert-output/CG_07.mzML",
    "gpr37l1_igg_r4"= "/ip-proteomics/exp19-ms-zenotof-convert-output/CG_08.mzML"
)

if (ip_protein == "Vcam1") {
  cols <- exp17_vcam1_col_names
} else if (ip_protein == "Gpr37l1") {
  cols <- exp19_gpr37l1_col_names
} else {
  warning("You didn't set arguments correctly - ip_protein must be either 'vcam1' or 'gpr37l1' ") 
}

colnames(data_na_filtered)[match(cols, colnames(data_na_filtered))] <- names(cols)
colnames(data_na_filtered)


#### =============================== data wrangling =============================== ####
df <- data_na_filtered %>%
        as_tibble() %>%
        column_to_rownames(var = "Genes") %>% # genes column as rownames
        select(
          -Protein.Group, 
          -Protein.Names, 
          -First.Protein.Description,
          -N.Sequences,
          -N.Proteotypic.Sequences
        )

# sanity check
check <- df[ip_protein, ]
print(check)


#### =================================== at least 4 times in one condition ===================== ####

# sorting column names
df_sorted <- df[, sort(names(df))]

# setting conditioms
igg_cond <- df_sorted[, 1:4] #igg
ip_cond <- df_sorted[, 5:8] #ip

# check if they are 4 times in each condition
igg_cond_valid <- rowSums(!is.na(igg_cond)) >= 4
ip_cond_valid <- rowSums(!is.na(ip_cond)) >= 4

# the genes with 4nonNA
rows_with_4_nonNA <- igg_cond_valid | ip_cond_valid

# how many genes got through the filter
sum(rows_with_4_nonNA)

# sanity check
df_4times <- df_sorted[rows_with_4_nonNA, ]
df_check_NA <- df_sorted[!rows_with_4_nonNA, ]


#### =========================== log2 transformation & imputation =============================== ####

# set seed for reproducibilty
set.seed(123)

# log2 transformation
df_log2 <- df_4times %>%
  mutate(across(where(is.numeric), log2))

# sanity check
check <- df_log2[ip_protein, ]
print(check)

# save log2 raw values
write.csv(df_log2, paste0(ip_protein, "_raw_log2_paper.csv"), row.names = TRUE)

# imputation, will use MinDet for reproducibility (no MinProb)
data_matrix <- as.matrix(df_log2)
imputed_matrix <- impute.MinDet(dataSet.mvs = data_matrix,
                                      q = 0.001) # protein missing = really not detected
imputed_df <- as.data.frame(imputed_matrix)

# sanity check
row_values <- imputed_df[ip_protein, ]
print(row_values)

# save imputed dataframe
write.csv(imputed_df, paste0(ip_protein, "_imputed_log2_paper.csv"), row.names = TRUE)



#### =============================== Plotting =============================== ####

# setting custom colors for VCAM1 & GPR37L1
if (ip_protein == "Vcam1") {
  ip_color <- "#21a0e2"   
} else if (ip_protein == "Gpr37l1") {
  ip_color <- "#e28d21"   
} else {
  warning("You didn’t set arguments correctly — ip_protein must be either 'vcam1' or 'gpr37l1'")
  ip_color <- "#000000"   
}

condition_colors <- c("ip" = ip_color, "igg" = "#c6c6c6") 

# QC plot 1: ip protein line plot. Start point: filtered raw data

extract_protein_data <- function(df, protein_name) {
  # check if protein exists in rownames
  if(!(protein_name %in% rownames(df))) {
    return(NULL)
  }
  
  # extract the row for the specified protein
  protein_row <- df[protein_name, , drop = FALSE]
  
  # convert to long format
  protein_data <- data.frame(
    condition = colnames(protein_row),
    raw_log2_intensities = as.numeric(as.matrix(protein_row)),
    protein = protein_name
  )
  
  return(protein_data)
}

# Assuming your dataframe is called 'df_proteins'
protein_list <- c(
    ip_protein
#   "Vcam1",
#   "Cd2ap",
#   "Prrt1",
#   "Aak1",
#   "Chgb",
#   "Aqp4",
#   "Igdcc4",
#   "Slc39a5",
#   "Pmch"
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

# create an empty list to store results
results_list <- list()

# process each protein
for(protein in protein_list) {
  result <- extract_protein_data(df_log2, protein)
  if(!is.null(result)) {
    results_list[[protein]] <- result
  }
}

# combine all results
df_merge <- bind_rows(results_list)

# plot with NA handling
p1 <- ggplot(data=df_merge, aes(x=condition, y=raw_log2_intensities, group=protein)) +
        geom_line(aes(color=protein), na.rm = TRUE, size = 1) +  
        geom_point(aes(color=protein), na.rm = TRUE, size = 5) +  
        scale_x_discrete() +
        scale_y_continuous(expand=c(0, 0), limits=c(0, 20)) + 
        scale_color_manual(values = setNames(ip_color, ip_protein)) +
        theme_minimal(base_size = 22) +
        theme(
            panel.background = element_blank(),
            panel.grid = element_blank(),
            axis.line = element_line(size = 1.2),
            axis.ticks.length = unit(0.2, "cm"),    
            axis.ticks = element_line(size = 1.2),
            axis.text.x = element_text(angle = 45, hjust = 1)
        ) +
        labs(title="Raw Intensity Levels",
            x="Condition",
            y="Log2 Intensity") + theme(text=element_text(size=20,  family="Arial"))
p1
ggsave(paste0(ip_protein, "_line_plot_paper.svg"), 
    plot = p1, 
    device = cairo_pdf,
    width = 25, height = 20, units = "cm", dpi=300)


# QC plot 2: nr of protein identified per sample. Start point: filtered raw data

# prepare the data for the barplot
plot_df <- df %>%
  # count the number of non-NA entries (identified proteins) for each sample
  summarise(across(everything(), ~sum(!is.na(.)))) %>%
  # transpose the data frame to have samples as a column
  pivot_longer(cols = everything(), names_to = "Sample", values_to = "UniqueProteins") %>%
  # extract the condition (ip or igg) from the sample name
  mutate(Condition = str_extract(Sample, "(ip|igg)"))

# create barplot
p2 <- ggplot(plot_df, aes(x = Sample, y = UniqueProteins, fill = Condition)) +
        geom_bar(stat = "identity") +
        # geom_text(aes(label = UniqueProteins),
        #             hjust = 0.5, 
        #             vjust = 2, 
        #             color = "white",
        #             size = 5) + 
        scale_fill_manual(values = condition_colors) + 
        scale_y_continuous(
            expand = c(0, 0),
            limits = c(0, 2500),   
            breaks = c(0, 500, 1000, 1500, 2000, 2500)  
        ) +
        labs(
            title = "Number of Unique Proteins Identified in Each Sample",
            x = "Sample",
            y = "Number of Unique Proteins"
        ) +
        theme_minimal(base_size = 22) +
        theme(
            panel.background = element_blank(),
            panel.grid = element_blank(),
            axis.line = element_line(size = 1.2, color = "black"),   
            axis.ticks = element_line(size = 1.2, color = "black"), 
            axis.ticks.length = unit(0.2, "cm"),                     
            axis.text.x = element_text(angle = 45, hjust = 1, size = 16), 
            axis.text.y = element_text(size = 16),
            axis.title.y = element_text(size = 18, family = "Arial"),
            axis.title.x = element_text(size = 18, family = "Arial"), 
            plot.title = element_text(size = 20, family = "Arial"),
            legend.position = "right",
            text = element_text(family = "Arial")  
        )
p2

ggsave(paste0(ip_protein, "_bar_plot_ids_paper.svg"), 
    plot = p2, 
    device = cairo_pdf,
    width = 25, height = 20, units = "cm", dpi=300)


# QC plot 3: log2 intensity plots per sample. Start point: imputed data.

plot_df_long <- imputed_df %>%
  rownames_to_column(var = "Protein") %>%
  pivot_longer(cols = -Protein, names_to = "Sample", values_to = "Intensity") %>%
  # Extract the condition (ip or igg) from the sample name
  mutate(Condition = str_extract(Sample, "(ip|igg)"))

# Print the first few rows of the long format data frame
print(head(plot_df_long))

# Create the box plot
p3 <- ggplot(plot_df_long, aes(x = Sample, y = Intensity, fill = Condition)) +
        geom_boxplot(linewidth = 1.2, outlier.size = 1.5) + 
        scale_fill_manual(values = condition_colors) +
        labs(
            title = "Protein Intensity Distribution per Sample",
            x = "Sample",
            y = "Log2 Intensity",
            fill = "Condition"
        ) +
        theme_minimal() +
        theme(
            panel.background = element_blank(),
            panel.grid = element_blank(),
            axis.line = element_line(size = 1.2, color = "black"),   
            axis.ticks = element_line(size = 1.2, color = "black"), 
            axis.ticks.length = unit(0.2, "cm"),                     
            axis.text.x = element_text(angle = 45, hjust = 1, size = 16), 
            axis.text.y = element_text(size = 16),
            axis.title.y = element_text(size = 18, family = "Arial"),
            axis.title.x = element_text(size = 18, family = "Arial"), 
            plot.title = element_text(size = 20, family = "Arial"),
            legend.position = "right",
            text = element_text(family = "Arial") 
            ) + 
            scale_y_continuous(
                expand = c(0, 0),
                limits = c(0, 25),   
                breaks = c(0, 5, 10, 15, 20, 25)  
            ) 
p3

ggsave(paste0(ip_protein, "_box_plot_intensity_paper.svg"), 
    plot = p3, 
    device = cairo_pdf,
    width = 20, height = 20, units = "cm", dpi=300)


# QC plot 4: PCA plot. Start point: imputed data.

# transpose the data frame so that samples are rows and proteins are columns
df_t <- as.data.frame(t(imputed_df))

# perform PCA
pca_result <- prcomp(df_t, scale. = TRUE) 

# extract the components and their variances
pca_scores <- as.data.frame(pca_result$x)
pca_variance <- pca_result$sdev^2
pca_percentage_variance <- round(100 * pca_variance / sum(pca_variance), 1)

# create df or plotting
pca_plot_df <- pca_scores %>%
  # Add the sample names as a column
  mutate(Sample = rownames(df_t)) %>%
  # Extract the condition (ip or igg) from the sample name
  mutate(Condition = stringr::str_extract(Sample, "(ip|igg)"))

# check the data
print(head(pca_plot_df))

# create the PCA plot
p4 <- ggplot(pca_plot_df, aes(x = PC1, y = PC2, color = Condition, label = Sample)) +
            geom_point(size = 4) +
            #geom_text(hjust = -0.1, vjust = 0.1, size = 5) + # Add sample labels
            scale_color_manual(values = condition_colors) +
            labs(
                title = "PCA Plot of Samples",
                x = paste0("PC1 (", pca_percentage_variance[1], "% variance explained)"),
                y = paste0("PC2 (", pca_percentage_variance[2], "% variance explained)")
            ) +
            theme_minimal() +
            theme(
                panel.background = element_blank(),
                panel.grid = element_blank(),
                axis.line = element_line(size = 1.2, color = "black"),   
                axis.ticks = element_line(size = 1.2, color = "black"), 
                axis.ticks.length = unit(0.2, "cm"),                     
                axis.text.x = element_text(angle = 45, hjust = 1, size = 16), 
                axis.text.y = element_text(size = 16),
                axis.title.y = element_text(size = 18, family = "Arial"),
                axis.title.x = element_text(size = 18, family = "Arial"), 
                plot.title = element_text(size = 20, family = "Arial"),
                #legend.position = "right",
                text = element_text(family = "Arial") 
            ) +
            scale_x_continuous(
                expand = c(0, 0),
                limits = c(-45, 40),
                breaks = seq(-40, 40, by = 20)  
            ) +
            scale_y_continuous(
                expand = c(0, 0),
                limits = c(-40, 40),
                breaks = seq(-40, 40, by = 20) 
            )
p4

ggsave(paste0(ip_protein, "_pca_paper.svg"), 
    plot = p4, 
    device = cairo_pdf,
    width = 22, height = 20, units = "cm", dpi=300)



#### =============================== DEP: Limma =============================== #####

# create factors for experimental design
condition <- factor(c(
    "igg",
    "igg",
    "igg",
    "igg",
    "ip",
    "ip",
    "ip",
    "ip" 
    )
)
replicate <- factor(c(
    "rep_1", 
    "rep_2", 
    "rep_3",
    "rep_4",
    "rep_1", 
    "rep_2", 
    "rep_3",
    "rep_4"
    )
)
# paired sample design
design <- model.matrix(~replicate + condition)
colnames(design)
# fit the linear model
fit1 <- limma::lmFit(imputed_df, design)
cont <- makeContrasts(conditionip, levels = design)
fit2 <- contrasts.fit(fit1, contrasts = cont)
fit3 <- eBayes(fit2)

# get the dep results
results_limma <- topTable(fit3, adjust = "fdr", sort.by = "P", n = 7000)
limma_ip_proteins <- results_limma[results_limma$adj.P.Val < 0.05 & results_limma$logFC > 1, ]

# before writing the file, also include the Protein.Names
limma_ip_proteins2 <- limma_ip_proteins %>%
  mutate(Genes = rownames(.)) %>%
  left_join(data_na_filtered[, c("Genes", "Protein.Names", "Protein.Group")],
            by = c("Genes" = "Genes"))

write.csv(limma_ip_proteins2, paste0(ip_protein, "_limma_ip_proteins_paper.csv"), row.names = TRUE)


#### =============================== Limma volcano plot without annotations =============================== ####

results_limma <- tibble::rownames_to_column(results_limma, "Genes")

# making the labels for volcano plot
labels <- results_limma$Genes
labels_upper <- paste0(toupper(labels))

select_proteins <- c(ip_protein, "Cdh9", "Cdh8", "Adgrb1", "Nrn1", "Gabrb3", "Gabrb1")
select_proteins <- c()
select_proteins_upper <- paste0(toupper(select_proteins))

# making the statement color the points
keyvals <- ifelse(
    results_limma$logFC > 1 & results_limma$adj.P.Val < 0.05, ip_color, '#c6c6c6')
keyvals[is.na(keyvals)] <- '#c6c6c6'
names(keyvals)[keyvals == ip_color] <- 'ip'
names(keyvals)[keyvals == '#c6c6c6'] <- 'unspecific'

# making the volcano plot
p5 <- EnhancedVolcano(results_limma,
        lab = labels_upper,
        x = 'logFC',
        y = 'P.Value',
        selectLab = select_proteins_upper,
        colCustom = keyvals,
        #shapeCustom = keyvals.shape,
        drawConnectors = TRUE,
        xlim = c(-11, 11),
        ylim = c(-0, 7),
        cutoffLineType = 'blank',
        colAlpha = 1,
        shape = 19,
        pointSize = 4.0,
        labSize = 6,
        gridlines.major = FALSE,
        gridlines.minor = FALSE,
        arrowheads = FALSE,
        parseLabels = TRUE,
        lengthConnectors = unit(0.1, 'npc'),
        ) + 
        theme_minimal(base_size = 22) +  
        theme(
        panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.line = element_line(size = 1.2, color = "black"),   
        axis.ticks = element_line(size = 1.2, color = "black"),  
        axis.ticks.length = unit(0.3, "cm"),                     
        axis.text.x = element_text(size = 16, family = "Arial"),
        axis.text.y = element_text(size = 16, family = "Arial"),
        axis.title.x = element_text(size = 18, family = "Arial"),
        axis.title.y = element_text(size = 18, family = "Arial"),
        plot.title = element_text(size = 20, family = "Arial"),
        text = element_text(family = "Arial"),
        legend.position = "none"
        ) +
        scale_x_continuous(
            expand = c(0, 0),
            limits = c(-8, 8),
            breaks = seq(-8, 8, by = 2)  
        ) +
        scale_y_continuous(
            expand = c(0, 0),
            limits = c(0, 6.5),
            breaks = seq(0, 7, by = 2) 
        )
p5
ggsave(paste0(ip_protein, "_volcano_paper.svg"), 
    plot = p5, 
    device = cairo_pdf,
    width = 20, height = 20, units = "cm", dpi=300)
