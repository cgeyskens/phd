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
Vcam1_1misc_imput_log2 <- read_csv("Vcam1_1miscleav_imput_log2.csv", col_names = TRUE)
colnames(Vcam1_1misc_imput_log2)[1] <- "Genes"

# loading in  the IP proteins annotations for Uniprot filtering
Vcam1_all_annotations <- read_csv("Vcam1_all_limma_ip_proteins_annotations.csv")

# loading in the IP limma proteins specific to experiment (1misc)
Vcam1_ip_limma <- read_csv("Vcam1_limma_ip_proteins_zenotof_1miscleav_wo_gpr37l1.csv")


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
merged_df <- Vcam1_ip_limma %>%
  inner_join(Vcam1_all_annotations, by = c("Genes" = "Genes")) %>%
  left_join(Vcam1_1misc_imput_log2, by = c("Genes" = "Genes")) %>%
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
top19 <- rank_tbl %>%
  arrange(desc(rank)) %>%
  slice(1:19)

# Basic rank plot with top 50 labeled
library(ggrepel)  # For non-overlapping text labels

p <- ggplot(rank_tbl, aes(x = rank, y = ratio)) +
  geom_point(color = "gray30", size = 1.5) +
  geom_point(data = top19, aes(x = rank, y = ratio), 
             color = "red", size = 2) +
  geom_text_repel(data = top19, 
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








