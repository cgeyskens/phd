library(patchwork)
library(ggplot2)
library(readr)
library(httpgd)
library(dplyr)
library (tibble)
hgd()
library(ComplexHeatmap)

#### ==================== GPR37L1 ================== ####


#### ======= Heatmap of raw Log2 values, no imputation & ranked according to median of IP

# load file with raw log2 values
Gpr37l1_1misc_log2 <- read_csv("Gpr37l1_1miscleav_raw_log2.csv", col_names = TRUE)
colnames(Gpr37l1_1misc_log2)[1] <- "Genes"

# loading in the only the IP proteins
Gpr37l1_all_annotations <- read_csv("Gpr37l1_all_limma_ip_proteins_annotations.csv")

# filter on uniprot membrane
Gpr37l1_all_annotations_filtered <- Gpr37l1_all_annotations %>% 
  filter(uniprot_membrane_related == 1)

cols_to_drop = c(
  "in_van_oostrum_2023", 
  "in_sorokina_2021", 
  "sorokina_2021_synaptosome", 
  "sorokina_2021_presynaptic", 
  "sorokina_2021_postsynaptic", 
  "in_syngo", 
  "uniprot_membrane_related",
  "Protein.Names",
  "Protein.Group",
  "...1"
)
# wrangling
merged_df <- Gpr37l1_all_annotations_filtered %>%
  left_join(Gpr37l1_1misc_log2, by = c("Genes" = "Genes")) %>%
  select(-all_of(cols_to_drop))

col_order <- c("gpr37l1_ip_r1", "gpr37l1_ip_r2", "gpr37l1_ip_r3", "gpr37l1_ip_r4",
              "gpr37l1_igg_r1", "gpr37l1_igg_r2", "gpr37l1_igg_r3", "gpr37l1_igg_r4")

mat <- as.data.frame(merged_df)
rownames(mat) <- mat$Genes
mat <- mat[, col_order]

# ranking
ip_mean <- rowMeans(mat[, 1:4], na.rm = TRUE)
ord <- order(-ip_mean)

# visualization
mat_plot <- as.matrix(mat[ord, , drop = FALSE])

pheatmap(
  mat_plot,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  show_rownames = TRUE,
  show_colnames = TRUE,
  scale = "none",
  color = colorRampPalette(c("white", "navy"))(100),
  na_col = "grey90",
  border_color = NA
)


#### ======= Rank plot, with imputated values, of Log2 values (mean IP/mean IgG) ==== ####

# already merged_df from before
Gpr37l1_1misc_imput_log2 <- read_csv("Gpr37l1_1miscleav_imput_log2.csv", col_names = TRUE)
colnames(Gpr37l1_1misc_imput_log2)[1] <- "Genes"

# loading in  the IP proteins annotations for Uniprot filtering
Gpr37l1_all_annotations <- read_csv("Gpr37l1_all_limma_ip_proteins_annotations.csv")

# loading in the IP limma proteins specific to experiment (1misc)
Gpr37l1_ip_limma <- read_csv("Gpr37l1_limma_ip_proteins_zenotof_1miscleav.csv")

# filter on uniprot membrane
Gpr37l1_all_annotations_filtered <- Gpr37l1_all_annotations %>% 
  filter(uniprot_membrane_related == 1)

cols_to_drop = c(
  "in_van_oostrum_2023", 
  "in_sorokina_2021", 
  "sorokina_2021_synaptosome", 
  "sorokina_2021_presynaptic", 
  "sorokina_2021_postsynaptic", 
  "in_syngo", 
  "uniprot_membrane_related",
  "...1.x",
  "logFC",
  "AveExpr",
  "t",
  "Protein.Names.x",
  "Protein.Group.x",
  "...1.y",
  "Protein.Names.y",
  "Protein.Group.y",
  "P.Value",
  "adj.P.Val",
  "B"
)
# wrangling
merged_df <- Gpr37l1_ip_limma %>%
  inner_join(Gpr37l1_all_annotations_filtered, by = c("Genes" = "Genes")) %>%
  left_join(Gpr37l1_1misc_imput_log2, by = c("Genes" = "Genes")) %>%
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

# Identify top 50 enriched proteins (highest ranks)
top50 <- rank_tbl %>%
  arrange(desc(rank)) %>%
  slice(1:50)

# Basic rank plot with top 50 labeled
library(ggrepel)  # For non-overlapping text labels

p <- ggplot(rank_tbl, aes(x = rank, y = ratio)) +
  geom_point(color = "gray30", size = 1.5) +
  geom_point(data = top50, aes(x = rank, y = ratio), 
             color = "red", size = 2) +
  geom_text_repel(data = top50, 
                  aes(x = rank, y = ratio, label = Genes),  # Change 'Protein' to your protein name column
                  size = 3,
                  max.overlaps = 50,
                  box.padding = 0.5,
                  segment.color = "red",
                  segment.size = 0.3) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  labs(
    x = "Protein rank (by IP / IgG ratio)",
    y = "IP mean / IgG mean ratio",
    title = "Rank plot of proteins by IP enrichment over IgG",
    subtitle = "Top 50 enriched proteins labeled in red"
  ) +
  theme_minimal(base_size = 12)

print(p)



#### =========================== VCAM1 ============================ ####


#### ======= Rank plot, with imputated values, of Log2 values (mean IP/mean IgG) ==== ####

# already merged_df from before
Vcam1_imput_log2 <- read_csv("Vcam1_imputed_log2_paper.csv", col_names = TRUE)
colnames(Vcam1_imput_log2)[1] <- "Genes"

# loading in  the IP proteins annotations for Uniprot filtering
Vcam1_all_annotations <- read_csv("Vcam1_all_limma_ip_proteins_annotations_paper.csv")

# loading in the IP limma proteins specific to experiment
Vcam1_ip_limma <- read_csv("Vcam1_limma_ip_proteins_paper.csv")


cols_to_drop = c(
  "in_van_oostrum_2023", 
  "in_sorokina_2021", 
  "sorokina_2021_synaptosome", 
  "sorokina_2021_presynaptic", 
  "sorokina_2021_postsynaptic", 
  "in_syngo", 
  "uniprot_membrane_related",
  "...1.x",
  "logFC.x",
  "AveExpr.x",
  "t.x",
  "Protein.Names.x",
  "Protein.Group.x",
  "...1.y",
  "Protein.Names.y",
  "Protein.Group.y",
  "P.Value.x",
  "adj.P.Val.x",
  "B.x",
  "logFC.y",
  "AveExpr.y",
  "t.y",
  "Protein.Names.y",
  "Protein.Group.y",
  "...1.y",
  "Protein.Names.y",
  "Protein.Group.y",
  "P.Value.y",
  "adj.P.Val.y",
  "B.y"
)
# wrangling
merged_df <- Vcam1_ip_limma %>%
  inner_join(Vcam1_all_annotations, by = c("Genes" = "Genes")) %>%
  left_join(Vcam1_imput_log2, by = c("Genes" = "Genes")) %>%
  select(-all_of(cols_to_drop))

# merged_df <- Vcam1_ip_limma %>%
#   select(-all_of(cols_to_drop))

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

# Identify top 25 enriched proteins (highest ranks)
top_candidates <- rank_tbl %>%
  arrange(desc(rank)) %>%
  slice(1:25)

# Basic rank plot with top 50 labeled
library(ggrepel)  # For non-overlapping text labels

p <- ggplot(rank_tbl, aes(x = rank, y = ratio)) +
  geom_point(color = "gray30", size = 1.5) +
  geom_point(data = top_candidates, aes(x = rank, y = ratio), 
             color = "red", size = 2) +
  geom_text_repel(data = top_candidates, 
                  aes(x = rank, y = ratio, label = Genes),  # Change 'Protein' to your protein name column
                  size = 3,
                  max.overlaps = 50,
                  box.padding = 0.5,
                  segment.color = "red",
                  segment.size = 0.3) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  labs(
    x = "Protein rank (by IP / IgG ratio)",
    y = "IP mean / IgG mean ratio",
    title = "Rank plot of proteins by IP enrichment over IgG",
    subtitle = "Top 50 enriched proteins labeled in red"
  ) +
  theme_minimal(base_size = 12) +
  theme(
  panel.grid = element_blank(),       # removes all gridlines
  panel.border = element_blank(),
  axis.line = element_line(color = "black"),  # adds clean axis lines
  plot.title = element_text(face = "bold"),
  plot.subtitle = element_text(size = 10)
)

p


selected_proteins <- c("Vcam1", "Prrt1", "Aqp4", "Igdcc4", "Slc39a", "Pmch", "Chgb")
selected_proteins <- c()


# Add a flag to your dataset
rank_tbl$highlight <- ifelse(rank_tbl$Genes %in% selected_proteins, "highlight", "other")

p <- ggplot(rank_tbl, aes(x = rank, y = ratio)) +
  # Plot all points, coloring based on highlight status
  geom_point(
    aes(color = highlight),
    size = 3,
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
  # Customize colors manually
  scale_color_manual(
    values = c("highlight" = "#21a0e2", "other" = "black"),
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






