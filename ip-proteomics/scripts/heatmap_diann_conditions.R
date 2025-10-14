library(ggplot2)
library(reshape2)
library(httpgd)
library(patchwork)
hgd()

if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("ComplexHeatmap")
library(ComplexHeatmap)

#### =========================== VCAM1 data ============================ ####
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

# Helper to pull proteins from column X
clean_prots <- function(df) {
  x <- as.character(df$X)
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


freq <- rowSums(presence_mat)
ord  <- order(-freq, rownames(presence_mat))
pmat <- presence_mat[ord, , drop = FALSE]


# annotations
conds <- colnames(pmat)

Platform <- ifelse(grepl("astral",   conds, ignore.case = TRUE), "Astral", "Zenotof")
WithWO   <- ifelse(grepl("_wo_",     conds, ignore.case = TRUE), "wo gpr37l1", "with gpr37l1")
Miscl    <- ifelse(grepl("_1miscleavage", conds, ignore.case = TRUE), "1",
             ifelse(grepl("_2miscleavage", conds, ignore.case = TRUE), "2", NA))

# --- Simple binary heatmap ---
bin_colors <- c("0" = "white", "1" = "grey")

ann_colors <- list(
  Platform   = c(Astral = "#21B4B7", Zenotof = "#B681F2"),
  `GPR37L1`  = c(`with gpr37l1` = "#7FB069", `wo gpr37l1` = "#E67E80"),
  Miscleavage = c("1" = "#F3A683", "2" = "#778BEB")
)

top_ha <- HeatmapAnnotation(
  Platform   = Platform,
  `GPR37L1`  = WithWO,
  Miscleavage = Miscl,
  col = ann_colors,
  annotation_name_side = "left",
  annotation_name_gp = gpar(fontsize = 9),
  simple_anno_size = unit(3.5, "mm")  # strip height
)

Platform_f <- factor(Platform, levels = c("Astral", "Zenotof"))
GPR_f      <- factor(WithWO,   levels = c("with gpr37l1", "wo gpr37l1"))
Miscl_f    <- factor(Miscl,    levels = c("1", "2"))


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
    Platform   = Platform_f,
    GPR37L1    = GPR_f,
    Miscleavage= Miscl_f
  ),
  na_col = "white"
)

draw(ht)





#### =========================== GPR37L1 data ============================ ####

#### ===== loading the data ===== ####

gpr37l1_1miscleavage <- read.csv("Gpr37l1_limma_ip_proteins_zenotof_1miscleav.csv")
gpr37l1_2miscleavage <- read.csv("Gpr37l1_limma_ip_proteins_zenotof_2miscleav.csv")


gpr37l1_1miscleavage$X
gpr37l1_2miscleavage$X

