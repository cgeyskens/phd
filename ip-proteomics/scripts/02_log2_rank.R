###############################################################################
# Script: 02_log2_rank.R
# Purpose: log2_rank plot using raw imputed intensities
# Author: Cydric Geyskens
# Date: 2025-10-27
###############################################################################

library(patchwork)
library(ggplot2)
library(readr)
library(httpgd)
library(dplyr)
library (tibble)
library(ggrepel)

library(httpgd)
hgd()

#### =============================== arguments =============================== ####
# only need to change these arguments for full analysis
# vcam1
imputed_log2_file_path <- "Vcam1_imputed_log2_paper.csv"
ip_limma_file_path <- "Vcam1_limma_ip_proteins_paper.csv"
ip_protein = "Vcam1"

# gpr37l1
imputed_log2_file_path <- "Gpr37l1_imputed_log2_paper.csv"
ip_limma_file_path <- "Gpr37l1_limma_ip_proteins_paper.csv"
ip_protein = "Gpr37l1"


#### ======= rank plot, with imputated values, of Log2 values (mean IP/mean IgG) ==== ####

# load raw log2 imputed intensities
imput_log2 <- read_csv(imputed_log2_file_path, col_names = TRUE)
colnames(imput_log2)[1] <- "Genes"

# loading in the IP limma proteins specific to experiment
ip_limma <- read_csv(ip_limma_file_path)

cols_to_drop = c(
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

# wrangling
merged_df <- ip_limma %>%
  left_join(imput_log2, by = c("Genes" = "Genes")) %>%
  select(-all_of(cols_to_drop))

colnames(merged_df)

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

# Identify top enriched proteins (highest ranks)
top_candidates <- rank_tbl %>%
  arrange(desc(rank)) 

# select the protein to visualize
selected_proteins <- c("Vcam1", "Prrt1", "Aqp4", "Igdcc4", "Slc39a5", "Pmch", "Chgb")
#selected_proteins <- c()

# Add a flag to your dataset
rank_tbl$highlight <- ifelse(rank_tbl$Genes %in% selected_proteins, "highlight", "other")

p <- ggplot(rank_tbl, aes(x = rank, y = ratio)) +
  # Plot all points, coloring based on highlight status
  geom_point(
    aes(color = highlight),
    size = 6,
    shape = 16,
    alpha = 0.9
  ) +
  # Annotate only selected proteins
  geom_text_repel(
    data = subset(rank_tbl, highlight == "highlight"),
    aes(label = Genes),
    size = 4,
    color = "#21a0e2",
    box.padding = 0.5,
    segment.color = "#21a0e2",
    segment.size = 0.3
  ) +
  # Customize colors manually
  scale_color_manual(
    values = c("highlight" = "#21a0e2", "other" = "#c6c6c6"),
    guide = "none"  # hides legend
  ) +
  labs(
    x = "Protein rank (by IP / IgG ratio)",
    y = "IP mean / IgG mean ratio",
    title = "Rank plot of proteins by IP enrichment over IgG",
    subtitle = "Selected proteins highlighted and annotated"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid = element_blank(),   # remove gridlines
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
      limits = c(0, 26),
      breaks = seq(0, 25, by = 5)  
  ) +
  scale_y_continuous(
      expand = c(0, 0),
      limits = c(1, 2.5),
      breaks = seq(0, 2.5, by = 0.5) 
  )
p

ggsave(paste0(ip_protein, "_rank.svg"), 
    plot = p, 
    device = cairo_pdf,
    width = 22, height = 20, units = "cm", dpi=300)






