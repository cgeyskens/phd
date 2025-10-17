library(ggplot2)
library(reshape2)
library(httpgd)
library(patchwork)
library(readr)

hgd()

if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("ComplexHeatmap")
library(ComplexHeatmap)

#### =========================== VCAM1 data ============================ ####

#### Loading in the files
# define files
files <- c(
  Vcam1_astral_wo_1miscleavage   = "Vcam1_limma_ip_proteins_astral_1miscleav_wo_gpr37l1.csv",
  Vcam1_astral_with_1miscleavage = "Vcam1_limma_ip_proteins_astral_1miscleav.csv",
  Vcam1_astral_wo_2miscleavage   = "Vcam1_limma_ip_proteins_astral_2miscleav_wo_gpr37l1.csv",
  Vcam1_astral_with_2miscleavage = "Vcam1_limma_ip_proteins_astral_2miscleav.csv",
  Vcam1_zenotof_wo_1miscleavage  = "Vcam1_limma_ip_proteins_zenotof_1miscleav_wo_gpr37l1.csv",
  Vcam1_zenotof_with_1miscleavage= "Vcam1_limma_ip_proteins_zenotof_1miscleav.csv",
  Vcam1_zenotof_wo_2miscleavage  = "Vcam1_limma_ip_proteins_zenotof_2miscleav_wo_gpr37l1.csv",
  Vcam1_zenotof_with_2miscleavage= "Vcam1_limma_ip_proteins_zenotof_2miscleav.csv"
)

# load file with annotations for all proteins
Vcam1_all_annotations <- read_csv("vcam1_all_limma_ip_proteins_annotations.csv")
cols_to_convert = c(
  "in_van_oostrum_2023", 
  "in_sorokina_2021", 
  "sorokina_2021_synaptosome", 
  "sorokina_2021_presynaptic",
  "sorokina_2021_postsynaptic",
  "in_syngo",
  "uniprot_membrane_related"
)


# convert from boolean to 1/0
colnames(Vcam1_all_annotations)
Vcam1_all_annotations[cols_to_convert] <- lapply(Vcam1_all_annotations[cols_to_convert], function(x) {
  x <- as.integer(x)        
  x[is.na(x)] <- 0          
  return(x)
})


# Helper to pull proteins from column Genes
clean_prots <- function(df) {
  x <- as.character(df$Genes)
  x <- trimws(x)
  unique(na.omit(x[nzchar(x)]))
}

# Read all, produce a *named list* of protein vectors per condition
proteins_by_condition <- lapply(files, function(f) clean_prots(read.csv(f)))
# names(proteins_by_condition) already taken from `files` names

# Build presence/absence matrix (rows = proteins, cols = conditions)
all_proteins <- sort(unique(unlist(proteins_by_condition)))
presence_mat <- sapply(proteins_by_condition, function(p) as.integer(all_proteins %in% p))
rownames(presence_mat) <- all_proteins
mode(presence_mat) <- "integer"


# Merging the annotations and the presence matrix of run conditions
# convert presence matrix to df
presence_df <- presence_mat |>
  as.data.frame() |>
  rownames_to_column(var = "Genes")

# filter out columns in annotation df
cols_to_drop <- c("Protein.Names", "Protein.Group", "...1")
Vcam1_all_annotations_filtered <- Vcam1_all_annotations %>% 
  select(-all_of(cols_to_drop))

# ensuring character type
presence_df$Genes <- as.character(presence_df$Genes)
Vcam1_all_annotations_filtered$Genes <- as.character(Vcam1_all_annotations_filtered$Genes)

# merge
merged_df <- merge(presence_df, Vcam1_all_annotations_filtered, by = "Genes", all.x = TRUE)

# back to matrix for visualization & Genes as rownames
rownames(merged_df) <- merged_df$Genes
merged_df$Genes <- NULL
merged_mat <- as.matrix(merged_df)


#### ======================= Visualization ==================== ####

db_cols <- c(
  "uniprot_membrane_related",
  "in_syngo",
  "in_van_oostrum_2023",
  "in_sorokina_2021",
  "sorokina_2021_synaptosome",
  "sorokina_2021_presynaptic",
  "sorokina_2021_postsynaptic"
)

db_colors <- c(
  "in_syngo"                   = "#577590",
  "uniprot_membrane_related"   = "#277DA1",
  "in_van_oostrum_2023"        = "#F94144",
  "in_sorokina_2021"           = "#F3722C",
  "sorokina_2021_synaptosome"  = "#F9C74F",
  "sorokina_2021_presynaptic"  = "#90BE6D",
  "sorokina_2021_postsynaptic" = "#43AA8B"
)

# 2) Order proteins by frequency using ONLY run columns (exclude DB cols)
run_cols <- setdiff(colnames(merged_mat), db_cols)
freq <- rowSums(merged_mat[, run_cols, drop = FALSE], na.rm = TRUE)
ord  <- order(-freq, rownames(merged_mat))

# 3) Build plotting matrix with DB columns moved to the end
other_cols <- setdiff(colnames(merged_mat), db_cols)   # preserves existing order of run columns
pmat <- merged_mat[ord, c(other_cols, db_cols), drop = FALSE]

# 4) Annotations
conds <- colnames(pmat)

Platform <- ifelse(grepl("astral", conds, ignore.case = TRUE), "Astral", "Zenotof")
WithWO   <- ifelse(grepl("_wo_",   conds, ignore.case = TRUE), "wo gpr37l1", "with gpr37l1")
Miscl    <- ifelse(grepl("_1miscleavage", conds, ignore.case = TRUE), "1",
             ifelse(grepl("_2miscleavage", conds, ignore.case = TRUE), "2", NA))

# Database annotation: label DB columns by their exact name, others = NA (will render blank)
Database <- ifelse(conds %in% db_cols, conds, NA)

# 5) Colors
bin_colors <- c("0" = "white", "1" = "grey")

ann_colors <- list(
  Platform    = c(Astral = "#21B4B7", Zenotof = "#B681F2"),
  `GPR37L1`   = c(`with gpr37l1` = "#7FB069", `wo gpr37l1` = "#E67E80"),
  Miscleavage = c("1" = "#F3A683", "2" = "#778BEB"),
  Database    = db_colors                # distinct color per database
)

top_ha <- HeatmapAnnotation(
  Platform    = Platform,
  `GPR37L1`   = WithWO,
  Miscleavage = Miscl,
  Database    = Database,                
  col = ann_colors,
  annotation_name_side = "left",
  annotation_name_gp = gpar(fontsize = 9),
  simple_anno_size = unit(3.5, "mm")
)

# Optional factors if you keep column_split (not required; you can omit column_split entirely)
Platform_f <- factor(Platform, levels = c("Astral", "Zenotof"))
GPR_f      <- factor(WithWO,   levels = c("with gpr37l1", "wo gpr37l1"))
Miscl_f    <- factor(Miscl,    levels = c("1", "2"))

# 6) Heatmap (no re-overwrites; pmat is final, DBs are last)
ht <- Heatmap(
  pmat,
  name = "Presence",
  col = bin_colors,
  heatmap_legend_param = list(at = c(0,1), labels = c("Absent","Present"), title = "Protein"),
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_names = TRUE,
  show_column_names = TRUE,
  column_names_rot = 45,
  row_names_side = "left",
  column_title = "Condition",
  row_title = "Protein",
  top_annotation = top_ha,
  column_split = data.frame(
    Platform    = Platform_f,
    GPR37L1     = GPR_f,
    Miscleavage = Miscl_f
  ),
  na_col = "white"
)

draw(ht)








#### =========================== GPR37L1 data ============================ ####

#### ====================== loading the data =================== ####
files <- c(
  Gpr37l1_1miscleavage   = "Gpr37l1_limma_ip_proteins_zenotof_1miscleav.csv",
  Gpr37l1_2miscleavage   = "Gpr37l1_limma_ip_proteins_zenotof_2miscleav.csv"
)

# load file with annotations for all proteins
Gpr37l1_all_annotations <- read_csv("Gpr37l1_all_limma_ip_proteins_annotations.csv")
cols_to_convert = c(
  "in_van_oostrum_2023", 
  "in_sorokina_2021", 
  "sorokina_2021_synaptosome", 
  "sorokina_2021_presynaptic",
  "sorokina_2021_postsynaptic",
  "in_syngo",
  "uniprot_membrane_related"
)

# load file with raw log2 values
Gpr37l1_1misc_log2 <- read_csv("Gpr37l1_1miscleav_raw_log2.csv", col_names = TRUE)
colnames(Gpr37l1_1misc_log2)[1] <- "Genes"


# convert from boolean to 1/0
colnames(Vcam1_all_annotations)
Gpr37l1_all_annotations[cols_to_convert] <- lapply(Gpr37l1_all_annotations[cols_to_convert], function(x) {
  x <- as.integer(x)        
  x[is.na(x)] <- 0          
  return(x)
})


# Helper to pull proteins from column Genes
clean_prots <- function(df) {
  x <- as.character(df$Genes)
  x <- trimws(x)
  unique(na.omit(x[nzchar(x)]))
}

# Read all, produce a *named list* of protein vectors per condition
proteins_by_condition <- lapply(files, function(f) clean_prots(read.csv(f)))
# names(proteins_by_condition) already taken from `files` names

# Build presence/absence matrix (rows = proteins, cols = conditions)
all_proteins <- sort(unique(unlist(proteins_by_condition)))
presence_mat <- sapply(proteins_by_condition, function(p) as.integer(all_proteins %in% p))
rownames(presence_mat) <- all_proteins
mode(presence_mat) <- "integer"




# Merging the annotations and the presence matrix of run conditions
# convert presence matrix to df
presence_df <- presence_mat |>
  as.data.frame() |>
  rownames_to_column(var = "Genes")

# filter out columns in annotation df
cols_to_drop <- c("Protein.Names", "Protein.Group", "...1")
Gpr37l1_all_annotations_filtered <- Gpr37l1_all_annotations %>% 
  select(-all_of(cols_to_drop)) %>%
  filter(uniprot_membrane_related == 1)

# ensuring character type
presence_df$Genes <- as.character(presence_df$Genes)
Gpr37l1_all_annotations_filtered$Genes <- as.character(Gpr37l1_all_annotations_filtered$Genes)

# merge
merged_df <- Gpr37l1_all_annotations_filtered %>%
  left_join(presence_df, by = c("Genes" = "Genes")) %>%
  left_join(Gpr37l1_1misc_log2, by = c("Genes" = "Genes")) %>%
  column_to_rownames(var = "Genes")

# only have the raw log2 columns
cols_to_drop <- c(
  "in_van_oostrum_2023", 
  "in_sorokina_2021", 
  "sorokina_2021_synaptosome", 
  "sorokina_2021_presynaptic", 
  "sorokina_2021_postsynaptic", 
  "in_syngo", 
  "uniprot_membrane_related",
  "Gpr37l1_1miscleavage",
  "Gpr37l1_2miscleavage")
Gpr37l1_1misc_log2_candidates <- merged_df %>% 
  select(-all_of(cols_to_drop)) 
# Gpr37l1_1misc_log2_candidates_mat <- as.matrix(Gpr37l1_1misc_log2_candidates)


# save merged_df
write.csv(merged_df, "Gpr37l1_matrix_for_heatmap.csv", row.names = TRUE)


# back to matrix for visualization & Genes as rownames
merged_mat <- as.matrix(merged_df)



#### ======================= Visualization ==================== ####

db_cols <- c(
  "uniprot_membrane_related",
  "in_syngo",
  "in_van_oostrum_2023",
  "in_sorokina_2021",
  "sorokina_2021_synaptosome",
  "sorokina_2021_presynaptic",
  "sorokina_2021_postsynaptic"
)

db_colors <- c(
  "in_syngo"                   = "#577590",
  "uniprot_membrane_related"   = "#277DA1",
  "in_van_oostrum_2023"        = "#F94144",
  "in_sorokina_2021"           = "#F3722C",
  "sorokina_2021_synaptosome"  = "#F9C74F",
  "sorokina_2021_presynaptic"  = "#90BE6D",
  "sorokina_2021_postsynaptic" = "#43AA8B"
)

# 2) Order proteins by frequency using ONLY run columns (exclude DB cols)
run_cols <- setdiff(colnames(merged_mat), db_cols)
freq <- rowSums(merged_mat[, run_cols, drop = FALSE], na.rm = TRUE)
ord  <- order(-freq, rownames(merged_mat))

# 3) Build plotting matrix with DB columns moved to the end
other_cols <- setdiff(colnames(merged_mat), db_cols)   # preserves existing order of run columns
pmat <- merged_mat[ord, c(other_cols, db_cols), drop = FALSE]


# 4) Annotations
conds <- colnames(pmat)

Miscl    <- ifelse(grepl("_1miscleavage", conds, ignore.case = TRUE), "1",
             ifelse(grepl("_2miscleavage", conds, ignore.case = TRUE), "2", NA))

# Database annotation: label DB columns by their exact name, others = NA (will render blank)
Database <- ifelse(conds %in% db_cols, conds, NA)

# 5) Colors
bin_colors <- c("0" = "white", "1" = "grey")

ann_colors <- list(
  Miscleavage = c("1" = "#F3A683", "2" = "#778BEB"),
  Database    = db_colors                # distinct color per database
)

top_ha <- HeatmapAnnotation(
  Miscleavage = Miscl,
  Database    = Database,                
  col = ann_colors,
  annotation_name_side = "left",
  annotation_name_gp = gpar(fontsize = 9),
  simple_anno_size = unit(3.5, "mm")
)

# Optional factors if you keep column_split (not required; you can omit column_split entirely)
Miscl_f    <- factor(Miscl,    levels = c("1", "2"))

# 6) Heatmap (no re-overwrites; pmat is final, DBs are last)
ht <- Heatmap(
  pmat,
  name = "Presence",
  col = bin_colors,
  heatmap_legend_param = list(at = c(0,1), labels = c("Absent","Present"), title = "Protein"),
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_names = TRUE,
  show_column_names = TRUE,
  column_names_rot = 45,
  row_names_side = "left",
  column_title = "Condition",
  row_title = "Protein",
  top_annotation = top_ha,
  na_col = "white",
  row_names_gp = gpar(fontsize = 8),
  column_names_gp = gpar(fontsize = 8)
)

draw(ht)


# heatmap of raw log2 values
col_order <- c("gpr37l1_ip_r1", "gpr37l1_ip_r2", "gpr37l1_ip_r3", "gpr37l1_ip_r4",
              "gpr37l1_igg_r1", "gpr37l1_igg_r2", "gpr37l1_igg_r3", "gpr37l1_igg_r4")



mat_plot <- as.matrix(Gpr37l1_1misc_log2_candidates_mat[ord, col_order, drop = FALSE])

pheatmap(
  mat_plot,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  show_rownames = TRUE,
  show_colnames = TRUE,
  scale = "none",         
  color = colorRampPalette(c("white","navy"))(100),
  na_col = "grey90",
  border_color = NA
)

