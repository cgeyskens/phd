###############################################################################
# Script: 03_log2_rank.R
# Purpose: log2_rank plot using raw imputed intensities
# Author: Cydric Geyskens
# Date: 2025-10-27
###############################################################################

library(patchwork)
library(ggplot2)
library(readr)
library(dplyr)
library(tibble)
library(ggrepel)

library(httpgd)
hgd()

#### =============================== arguments =============================== ####
# only need to change these arguments for full analysis
# vcam1
imputed_log2_file_path <- "results/Vcam1_imputed_log2_paper.csv"
ip_limma_file_path <- "results/Vcam1_limma_ip_proteins_paper.csv"
ip_protein = "Vcam1"

# gpr37l1
imputed_log2_file_path <- "results/Gpr37l1_imputed_log2_paper.csv"
ip_limma_file_path <- "results/Gpr37l1_all_limma_ip_proteins_annotations_paper.csv"
ip_protein = "Gpr37l1"


#### ================================ load data  ============================ ####

# load raw log2 imputed intensities
imput_log2 <- read_csv(imputed_log2_file_path, col_names = TRUE)
colnames(imput_log2)[1] <- "Genes"

# loading in the IP limma proteins specific to experiment
ip_limma <- read_csv(ip_limma_file_path)


#### ============================ data wrangling  ============================ ####

cols_to_drop_vcam1 = c(
  "...1",
  "logFC",
  "AveExpr",
  "t",
  "Protein.Names",
  "Protein.Group",
  "P.Value",
  "adj.P.Val",
  "B"
)

cols_to_drop_gpr37l1 = c(
  "...1",
  "logFC",
  "AveExpr",
  "t",
  "Protein.Names",
  "Protein.Group",
  "P.Value",
  "adj.P.Val",
  "B",
  "in_van_oostrum_2023",
  "in_sorokina_2021",
  "sorokina_2021_synaptosome",
  "sorokina_2021_presynaptic",
  "sorokina_2021_postsynaptic",
  "in_syngo",
  "uniprot_membrane_related"
)

# data wrangling and merging
if (ip_protein == "Gpr37l1") {
  cols_to_drop <- cols_to_drop_gpr37l1
  merged_df <- ip_limma %>%
    left_join(imput_log2, by = c("Genes" = "Genes")) %>%
    filter(uniprot_membrane_related == TRUE) %>%
    select(-all_of(cols_to_drop))
} else if (ip_protein == "Vcam1") {
  cols_to_drop <- cols_to_drop_vcam1
  merged_df <- ip_limma %>%
    left_join(imput_log2, by = c("Genes" = "Genes")) %>%
    select(-all_of(cols_to_drop))
}
# check
colnames(merged_df)


#### ========================== rank proteins ============================ ####

# Identify IP and IgG columns by name pattern (adjust if your naming differs)
ip_cols  <- grep("_ip_",  names(merged_df), value = TRUE, ignore.case = TRUE)
igg_cols <- grep("_igg_", names(merged_df), value = TRUE, ignore.case = TRUE)

# Compute mean log2 per group, convert to ratio, rank, and sort
rank_tbl <- merged_df %>%
  mutate(
    mean_ip_log2  = rowMeans(across(all_of(ip_cols)),  na.rm = TRUE),
    mean_igg_log2 = rowMeans(across(all_of(igg_cols)), na.rm = TRUE),
    ratio         = mean_ip_log2/mean_igg_log2,  
    rank          = rank(ratio, ties.method = "first")
  ) %>%
  arrange(rank)

# write out the ranked genes
rank_tbl_to_write_out <- rank_tbl %>% arrange(desc(rank))
write.csv(rank_tbl_to_write_out, paste0(ip_protein, "_co-ip_ranked.csv"), row.names = TRUE)


# Identify top enriched proteins (highest ranks)
top_candidates <- rank_tbl %>%
  arrange(desc(rank)) 

#### ================================ plotting ========================== ####

# setting custom colors for VCAM1 & GPR37L1
if (ip_protein == "Vcam1") {
  ip_color <- "#21a0e2"   
} else if (ip_protein == "Gpr37l1") {
  ip_color <- "#e28d21"   
} else {
  warning("You didn’t set arguments correctly — ip_protein must be either 'vcam1' or 'gpr37l1'")
  ip_color <- "#000000"   
}

# select the protein to visualize
selected_proteins <- c("Gpr37l1", "Cdh9", "Cdh8", "Nrn1", "Adgrb1", "Gabrb3", "Gabrb1")
#selected_proteins <- c()

# Add a flag to your dataset
rank_tbl$highlight <- ifelse(rank_tbl$Genes %in% selected_proteins, "highlight", "other")

p <- ggplot(rank_tbl, aes(x = rank, y = ratio)) +
  geom_point(
    aes(color = highlight),
    size = 6,
    shape = 16,
    alpha = 0.9
  ) +
  # Annotate only selected proteins
  # geom_text_repel(
  #   data = subset(rank_tbl, highlight == "highlight"),
  #   aes(label = Genes),
  #   size = 4,
  #   color = "#21a0e2",
  #   box.padding = 0.5,
  #   segment.color = "#21a0e2",
  #   segment.size = 0.3
  # ) +
  scale_color_manual(
    values = c("highlight" = ip_color, "other" = "#c6c6c6"),
    guide = "none"  
  ) +
  labs(
    x = "Protein rank (by IP / IgG ratio)",
    y = "IP mean / IgG mean ratio",
    title = "Rank plot of proteins by IP enrichment over IgG",
    subtitle = "Selected proteins highlighted and annotated"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid = element_blank(),   
    axis.line = element_line(color = "black"),
    plot.title = element_text(face = "bold"),
    plot.subtitle = element_text(size = 10)
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
      limits = c(0, 81),
      breaks = seq(0, 80, by = 20)  
  ) +
  scale_y_continuous(
      expand = c(0, 0),
      limits = c(1, 1.75),
      breaks = seq(0, 2.5, by = 0.25) 
  )
p

ggsave(paste0(ip_protein, "_rank_paper.svg"), 
    plot = p, 
    device = cairo_pdf,
    width = 22, height = 20, units = "cm", dpi=300)






