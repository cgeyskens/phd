# Upload your count matrix in the exported.csv file #import dataset -> from text(readr)
# Upload your location file #import dataset -> from text(readr)

# 1. Convert the count matrix from a .csv file to sparse .mtx file and filter out genes ------------------------------------------------

# Load libraries
library(Matrix)
library(dplyr)
library(ggplot2)

# Convert into dataframe
df_A1 <- as.data.frame(A1_DAPI_tiff_measurements)
df_A2 <- as.data.frame(A2_DAPI_tiff_measurements)
df_B2 <- as.data.frame(B2_DAPI_tiff_measurements)
df_C2 <- as.data.frame(C2_DAPI_tiff_measurements)
df_D1 <- as.data.frame(D1_DAPI_tiff_measurements)

# Convert first column as rownames
df_A1_with_row_names <- data.frame(df_A1[,-1], row.names=df_A1[,1])
df_A2_with_row_names <- data.frame(df_A2[,-1], row.names=df_A2[,1])
df_B2_with_row_names <- data.frame(df_B2[,-1], row.names=df_B2[,1])
df_C2_with_row_names <- data.frame(df_C2[,-1], row.names=df_C2[,1])
df_D1_with_row_names <- data.frame(df_D1[,-1], row.names=df_D1[,1])

# Merge all the files together
install.packages("DescTools")
library(DescTools)
df_all_with_row_names <- MultMerge(df_A1_with_row_names, df_A2_with_row_names,
                                   df_B2_with_row_names, df_C2_with_row_names,
                                   df_D1_with_row_names)


# Filter out genes that are not necessary for this analysis
# the genes that need to be filtered out
# the genes of the 2020 Nat Com Holt paper: Agt, Gfap, Unc13, Frzb, Fam107a, Ascl1
# cell type markers that created a lot of noise: Snap25, Aldoc, Slc1a2, Plp1
# genes with low nr of counts: Elfn2, Hapln1, Plnxa3

rowSums(df_all_with_row_names)


genes_to_filter_out_all <- c("Agt", "Gfap", "Unc13c", "Frzb","Fam107a", "Ascl1",
                             "Hapln1", "Elfn2", "Plxna3")

#the operation to filter out the genes
df_all_with_row_names_filter <- df_all_with_row_names[!(row.names(df_all_with_row_names) %in% genes_to_filter_out_all),]


# Convert into sparse matrix format (the file format needed for Seurat)
matrix_all <- data.matrix(df_all_with_row_names_filter)
sparse_matrix_all <- Matrix(matrix_all, sparse=TRUE)

head(sparse_matrix_all)


# 2. Creating Seurat object using F Pestana function -----------------------------------
library(Seurat)
library(clusteringR)
seurat_mtx_all_Pestana <- clusteringSeurat(sparse_matrix_all, "all", metadataAvailable = FALSE,
                                           normalizationMethod = "SCT", mapTypeValue = "umap",
                                           resolutionValue = 0.4)

# Optimization with dotplot
# Best seems to be SCT and resolution of 0.12 for all samples


# 3. Filter the regular Seurat object and spatial seurat object---------------------------
# delete cells that have less then 2 counts and 2 unique genes
seurat_mtx_all_Pestana_filter <- subset(seurat_mtx_all_Pestana, 
                                        subset = nCount_RNA > 2 & nFeature_RNA > 2)


# 4. Draw plots ---------------------------------------------------------------------------

str(seurat_mtx_all_Pestana_filter)

install.packages("visualisR")
library(visualisR)
install.packages("viridisLite")
library(viridis)
library(Seurat)
library(ggpubr)


# Cell type markers list
features_c <- c("Aldh1l1","Aqp4", "Slc1a3", "Slc1a2", "Aldoc","Pecam1", "Flt1","Mog", "Plp1","Csf1r", "C1qa","Gad1", "Gad2", "Slc17a7", "Slc17a6", "Rbfox3","Syt1", "Snap25")


# Candidates list
candidates_c <- c("Gpr37l1", "Hepacam","Lsamp", "Plxnb1", "Vcam1","Alcam","Cadm4", "Ntm", "Neo1", 
                  "Lrp1", "Slc3a2", "Ptprf", "Ncam1","Nlgn3", "Negr1")



# Draw umap plot -----------------------------------------------------------------
dim_plot_Pestana <- drawDimPlot(seurat_mtx_all_Pestana_filter_cell_names, "all", showLegend = TRUE,
                                heightValue = 15, widthValue = 25)


# Rename clusters to cell type names
seurat_mtx_all_Pestana_filter_cell_names <- RenameIdents(object=seurat_mtx_all_Pestana_filter, 
                                                         "0"="0. Dentate Gyrus, Granule Cell Layer ",
                                                         "1"="1. Protoplasmic Astrocytes",
                                                         "2"="2. Oligodendrocytes",
                                                         "3"="3. Cortical Neurons",
                                                         "4"="4. Endothelial Cells",
                                                         "5"="5. Medial Habenula Neurons",
                                                         "6"="6. CA1 Neurons",
                                                         "7"="7. Oligodendrocyte Precursor Cells",
                                                         "8"="8. CA2 & CA3 Neurons",
                                                         "9"="9. Thalamic Neurons",
                                                         "10"="10. Inhibitory Neurons",
                                                         "11"="11. Microglia",
                                                         "12"="12. Dentate Gyrus, Subgranular Cell Layer",
                                                         "13"="13. Choroid Plexus Cells",
                                                         "14"="14. Fibrous Astrocytes",
                                                         "15"="15. CA1, CA2 & CA3 Neurons")

                                                             
## Custom Dot Plot of candidates
dotPlotClusters <- DotPlot(seurat_mtx_all_Pestana_filter, 
                           features = candidates_c,
                           col.min = -1,
                           col.max = 1,
                           split.by = NULL,
                           idents = c(1,0,8,15,10,11,7),
                           dot.min = 0,
                           dot.scale = 15,
                           scale = T,group.by = NULL) + 
  theme(axis.text.x = element_text(angle = 90, size = 8),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid = element_blank(),
        #axis.ticks = element_blank(),
        #axis.line = element_blank(),
        axis.text.x.top = element_text(face = "italic"),
        axis.title.x.bottom = element_text(vjust = 10),
        axis.text.y = element_text(size = 8,face = "italic"),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8),
        legend.position = "right",
        legend.box = "horizontal",
        legend.key.height = unit(0.3,"cm"))+
  labs(x="",y="") +
  scale_x_discrete(position = "top") +
  scale_y_discrete(position = "left") +
  scale_colour_gradient2(low = "#266bae", mid= "grey", high = "#bc373f")

dotPlotClusters

ggsave(filename = paste0("../Results/DotPlots/", format(Sys.time(), "%Y%m%d_%H%M%S"), "DotPlot.png"), 
       plot = dotPlotClusters, width = 40, height = 15,
       path = "../Results/DotPlots/",
       units = "cm",limitsize = FALSE)



## Custom Dot Plot of cell type markers and DEGs

markers_DEGs_c <- c("Rbfox3","Syt1", "Snap25", "Slc17a7", "Slc17a6", "Gad1", "Gad2", "Aldh1l1","Aqp4", "Slc1a3", "Slc1a2", "Aldoc","Pecam1", "Flt1","Mog", "Plp1","Csf1r", "C1qa", "Tnr", "Vcan")


dotPlotClusters_markers <- DotPlot(seurat_mtx_all_Pestana_filter, 
                                   features = markers_DEGs_c,
                                   col.min = -1,
                                   col.max = 1,
                                   split.by = NULL,
                                   dot.min = 0,
                                   dot.scale = 12,
                                   scale = T,group.by = "seurat_clusters",
                                   cluster.idents = FALSE) + 
  theme(axis.text.x = element_text(angle = 90, size = 8),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid = element_blank(),
        #axis.ticks = element_blank(),
        #axis.line = element_blank(),
        axis.text.x.top = element_text(face = "italic"),
        axis.title.x.bottom = element_text(vjust = 10),
        axis.text.y = element_text(size = 8,face = "italic"),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8),
        legend.position = "right",
        legend.box = "horizontal",
        legend.key.height = unit(0.3,"cm"))+
  labs(x="",y="") +
  scale_x_discrete(position = "top") +
  scale_y_discrete(position = "left") +
  scale_colour_gradient2(low = "#266bae", mid= "grey", high = "#bc373f")

dotPlotClusters_markers

dotPlotClusters_final <- dotPlotClusters_markers + scale_y_discrete(limits = rev)

ggsave(filename = paste0("../Results/DotPlots/", format(Sys.time(), "%Y%m%d_%H%M%S"), "DotPlot.png"), 
       plot = dotPlotClusters_final, width = 40, height = 20, dpi = 700,
       path = "../Results/DotPlots/",
       units = "cm",limitsize = FALSE)



# 5. DEGs ------------------------------------------------------------------------------

# identify DEGs between clusters
DEGs_all <- FindAllMarkers(seurat_mtx_all_Pestana_filter, 
                           logfc.threshold = 0.05,
                           min.pct = 0.01,
                           only.pos = TRUE)

# Writing the DEGs file to an excel format
write.csv(DEGs_all, "C:/Users/cydri/Documents/KU Leuven - PhD/Resolve spatial transcriptomics/Analysis - cell segmentation/R file directory/DEGs_all.csv")


# identify DEGs from Protoplasmic Astrocytes (Cluster 1) 
# vs DG Granule cell layer (Cluster 0), CA2 & CA3 Neurons (Cluster 8), 
DEGsAsNe <- FindMarkers(seurat_mtx_all_Pestana_filter, ident.1 = c(0,8), ident.2 = 1,
                          logfc.threshold = 0.001, min.pct = 0.001, only.pos = FALSE)


# identify DEGs from Protoplasmic Astrocytes (Cluster 1) 
# vs DG Granule cell layer (Cluster 0), OPCs (Cluster 7), CA2 & CA3 Neurons (Cluster 8), 
# Inhibitory Neurons (Cluster 10), Microglia (Cluster 11), CA1, CA2, CA3 Neurons (Cluster 15)
DEGsAsNeMg <- FindMarkers(seurat_mtx_all_Pestana_filter, ident.1 = c(0,7,8,10,11,15), ident.2 = 1,
                          logfc.threshold = 0.001, min.pct = 0.001, only.pos = FALSE)




# Volcano plot of DEGSs between certain clusters -------------------------------------------------------------------------
if (!requireNamespace('BiocManager', quietly = TRUE))
  install.packages('BiocManager')

BiocManager::install('EnhancedVolcano')

library(EnhancedVolcano)
library(ggrepel)

# Make volcano of astrocytes vs Dentate gyrus granule cell layer and CA2/CA3 neurons
# first include the gene names as a column named "gene"
DEGsAsNe_gene <- tibble::rownames_to_column(DEGsAsNe, "gene")


# select the gene you want to display onto the volcano plot
genes_display_vol1 <- c("Aldoc", "Aldh1l1", "Aqp4", "Slc1a2", "Slc1a3", 
                        "Gpr37l1", "Hepacam", "Lsamp", "Alcam", "Vcam1", "Plxnb1",
                        "Gad1", "Gad2", "Csf1r", "C1qa", "Slc17a7", "Slc17a6", "Rbfox3","Syt1", "Snap25", "Vcan", "Tnr",
                        "Plp1", "Flt1", "Ntm", "Sirpa", "Mog", "Lrp1", "Cadm4")

keyvals <- ifelse(
  DEGsAsNe_gene$avg_log2FC < -0.15, "#e68613",
  ifelse(DEGsAsNe_gene$avg_log2FC > 0.15, "royalblue",
         "grey"))
keyvals[is.na(keyvals)] <- "black"
names(keyvals)[keyvals=="royalblue"] <- "high"
names(keyvals)[keyvals=="black"] <- "mid"
names(keyvals)[keyvals=="#e68613"] <- "low"


vol1 <- EnhancedVolcano(DEGsAsNe_gene, x="avg_log2FC", y="p_val_adj", 
                        lab = DEGsAsNe_gene$gene,
                        title = "Astrocytes (Cluster 0) vs other (Cluster 1, 2, 6, 7, 8, 9, 10, 11)",
                        cutoffLineType = "longdash",
                        cutoffLineWidth = 1,
                        labSize = 5,
                        colCustom = keyvals,
                        colAlpha = 0.9,
                        pointSize = 5,
                        max.overlaps = 30,
                        pCutoff = 0.05,
                        FCcutoff = 0.15, 
                        selectLab = genes_display_vol1,  
                        gridlines.major = FALSE, 
                        gridlines.minor = FALSE,
                        xlim = c(-4, 4),
                        drawConnectors = TRUE,
                        widthConnectors = 1,
                        typeConnectors = "closed",
                        endsConnectors = "last",
                        arrowheads = FALSE,
                        boxedLabels = TRUE,
                        borderWidth = 1.2,
                        border = "full",
                        borderColour = "black",
                        legendPosition = "none")
vol1 

ggsave(filename = paste0("../Results/VolcanPlots/", format(Sys.time(), "%Y%m%d_%H%M%S"), "VolcanoPlot.png"), 
       plot = vol1, width = 25, height = 20,
       path = "../Results/VolcanoPlots/",
       units = "cm",limitsize = FALSE)



# Make volcano of Astrocytes vs Dentate gyrus granule cell layer, CA2/CA3 neurons, 
# CA1, CA2, CA3 neurons, Inhibitory neurons, Microglia and OPCs

# first include the gene names as a column named "gene"
DEGsAsNeMg_gene <- tibble::rownames_to_column(DEGsAsNeMg, "gene")


# select the gene you want to display onto the volcano plot
genes_display_vol3 <- c("Aldoc", "Aldh1l1", "Aqp4", "Slc1a2", "Slc1a3", 
                        "Gpr37l1", "Hepacam", "Lsamp", "Alcam", "Vcam1", "Plxnb1",
                        "Gad1", "Gad2", "Csf1r", "C1qa", "Slc17a7", "Slc17a6", "Rbfox3","Syt1", "Snap25", "Vcan", "Tnr")

keyvals <- ifelse(
  DEGsAsNeMg_gene$avg_log2FC < -0.15, "#e68613",
  ifelse(DEGsAsNeMg_gene$avg_log2FC > 0.15, "royalblue",
         "grey"))
keyvals[is.na(keyvals)] <- "black"
names(keyvals)[keyvals=="royalblue"] <- "high"
names(keyvals)[keyvals=="black"] <- "mid"
names(keyvals)[keyvals=="#e68613"] <- "low"


vol3 <- EnhancedVolcano(DEGsAsNeMg_gene, x="avg_log2FC", y="p_val_adj", 
                        lab = DEGsAsNeMg_gene$gene,
                        title = "Astrocytes (Cluster 0) vs other (Cluster 1, 2, 6, 7, 8, 9, 10, 11)",
                        cutoffLineType = "longdash",
                        cutoffLineWidth = 1,
                        labSize = 5,
                        colCustom = keyvals,
                        colAlpha = 0.9,
                        pointSize = 5,
                        max.overlaps = 30,
                        pCutoff = 0.05,
                        FCcutoff = 0.15, 
                        selectLab = candidates_c,  
                        gridlines.major = FALSE, 
                        gridlines.minor = FALSE,
                        xlim = c(-3, 3),
                        drawConnectors = TRUE,
                        widthConnectors = 1,
                        typeConnectors = "closed",
                        endsConnectors = "last",
                        arrowheads = FALSE,
                        boxedLabels = TRUE,
                        borderWidth = 1.2,
                        border = "full",
                        borderColour = "black",
                        legendPosition = "none")
vol3 

ggsave(filename = paste0("../Results/VolcanPlots/", format(Sys.time(), "%Y%m%d_%H%M%S"), "VolcanoPlot.png"), 
       plot = vol3, width = 25, height = 20,
       path = "../Results/VolcanoPlots/",
       units = "cm",limitsize = FALSE)



# QC: nr of transcripts/cell over each sample, unique genes/cell over each sample
# in violin plot-------------------------------------------------------------------------------------

#Violin plot of nr transcripts/cell over each sample

VlnPlot1 <- VlnPlot(object = seurat_mtx_all_Pestana_filter,
                    features = "nCount_RNA",
                    group.by = "orig.ident",
                    pt.size = 0,
                    split.by = NULL,
                    cols = NULL,
                    y.max = 2500) + 
  scale_fill_manual(values=c("#818284", "#818284","#818284","#818284","#818284")) +
  theme(axis.line = element_line(colour = 'black', size = 1.5),
        axis.ticks = element_line(colour = "black", size = 1.5),
        axis.ticks.length = unit(0.25, "cm"),
        axis.text.y = element_text(angle = 0,size = 10,family = "Arial")) +
  stat_summary(fun.data = mean_sdl, geom = "point", size = 15, color = "black", shape = 95)

VlnPlot1

ggsave(filename = paste0("../Results/VlnPlots/", format(Sys.time(), "%Y%m%d_%H%M%S"), "VlnPlot.png"), 
       plot = VlnPlot1, width = 15, height = 10,
       path = "../Results/BarPlots/",
       units = "cm",limitsize = FALSE)


# Violin Plot of Unique genes/cell across samples

VlnPlot2 <- VlnPlot(object = seurat_mtx_all_Pestana_filter,
                    features = "nFeature_RNA",
                    group.by = "orig.ident",
                    pt.size = 0,
                    split.by = NULL,
                    cols = NULL,
                    y.max = 100) + 
  scale_fill_manual(values=c("#818284", "#818284","#818284","#818284","#818284")) +
  theme(axis.line = element_line(colour = 'black', size = 1.5),
        axis.ticks = element_line(colour = "black", size = 1.5),
        axis.ticks.length = unit(0.25, "cm"),
        axis.text.y = element_text(angle = 0,size = 10,family = "Arial")) +
  stat_summary(fun.data = mean_sdl, geom = "point", size = 15, color = "black", shape = 95)

VlnPlot2

ggsave(filename = paste0("../Results/VlnPlots/", format(Sys.time(), "%Y%m%d_%H%M%S"), "VlnPlot.png"), 
       plot = VlnPlot2, width = 15, height = 10,
       path = "../Results/BarPlots/",
       units = "cm",limitsize = FALSE)


# Violin plot of Nr of transcripts/cell over each cluster


VlnPlot3 <- VlnPlot(object = seurat_mtx_all_Pestana_filter,
                    features = "nCount_RNA",
                    idents = c(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15),
                    pt.size = 0,
                    split.by = NULL,
                    cols = NULL,
                    y.max = 2500) + 
  scale_fill_manual(values=c("#f8766d", "#e68613","#cd9600","#aba300","#7cae00", "#0cb702", "#00be67", "#00c19a", "#00bfc4", "#00b8e7", "#00a9ff", "#8494ff", "#c77cff", "#ed68ed", "#ff61cc", "#ff68a1")) +
  theme(axis.line = element_line(colour = 'black', size = 1.5),
        axis.ticks = element_line(colour = "black", size = 1.5),
        axis.ticks.length = unit(0.25, "cm"),
        axis.text.y = element_text(angle = 0,size = 10,family = "Arial")) +
  stat_summary(fun.data = mean_sdl, geom = "point", size = 7, color = "black", shape = 95)

VlnPlot3

ggsave(filename = paste0("../Results/VlnPlots/", format(Sys.time(), "%Y%m%d_%H%M%S"), "VlnPlot.png"), 
       plot = VlnPlot3, width = 20, height = 10,
       path = "../Results/BarPlots/",
       units = "cm",limitsize = FALSE)

# Violin plot of unique genes/cell over each cluster


VlnPlot4 <- VlnPlot(object = seurat_mtx_all_Pestana_filter,
                    features = "nFeature_RNA",
                    idents = c(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15),
                    pt.size = 0,
                    split.by = NULL,
                    cols = NULL,
                    y.max = 100) + 
  scale_fill_manual(values=c("#f8766d", "#e68613","#cd9600","#aba300","#7cae00", "#0cb702", "#00be67", "#00c19a", "#00bfc4", "#00b8e7", "#00a9ff", "#8494ff", "#c77cff", "#ed68ed", "#ff61cc", "#ff68a1")) +
  theme(axis.line = element_line(colour = 'black', size = 1.5),
        axis.ticks = element_line(colour = "black", size = 1.5),
        axis.ticks.length = unit(0.25, "cm"),
        axis.text.y = element_text(angle = 0,size = 10,family = "Arial")) +
  stat_summary(fun.data = mean_sdl, geom = "point", size = 7, color = "black", shape = 95)

VlnPlot4

ggsave(filename = paste0("../Results/VlnPlots/", format(Sys.time(), "%Y%m%d_%H%M%S"), "VlnPlot.png"), 
       plot = VlnPlot4, width = 20, height = 10,
       path = "../Results/BarPlots/",
       units = "cm",limitsize = FALSE)


# Barplot for nr of cells across samples

bar_plot_data <- seurat_mtx_all_Pestana_filter@meta.data

library(dplyr)

sum_bar_plot_data_sample<- bar_plot_data %>%
  group_by(orig.ident) %>%
  summarise(n=n())

library(ggplot2)


bar1 <- ggplot(sum_bar_plot_data_sample, aes(x=orig.ident, y=n, fill=orig.ident)) +
  geom_bar(stat = "identity", color="black", size=1) + NoLegend() +
  scale_fill_manual(values=c("#818284", "#818284","#818284","#818284","#818284")) +
  ylim(0,10000) +
  theme(panel.background = element_blank(),
        axis.line = element_line(colour = 'black', size = 2),
        axis.ticks = element_line(colour = "black", size = 2),
        axis.ticks.length = unit(0.25, "cm"),
        axis.text.y = element_text(angle = 0,size = 10,family = "Arial")) +
  labs(y = "Number of cells")

bar1

ggsave(filename = paste0("../Results/BarPlots/", format(Sys.time(), "%Y%m%d_%H%M%S"), "BarPlot.png"), 
       plot = bar1, width = 15, height = 10,
       path = "../Results/BarPlots/",
       units = "cm",limitsize = FALSE)


# Barplot for nr of cells across clusters
library(dplyr)

sum_bar_plot_data_cluster<- bar_plot_data %>%
  group_by(seurat_clusters) %>%
  summarise(n=n())

library(ggplot2)

bar2 <- ggplot(sum_bar_plot_data_cluster, aes(x=seurat_clusters, y=n, fill=seurat_clusters)) +
  geom_bar(stat = "identity", color="black", size=1) + NoLegend() +
  scale_fill_manual(values=c("#f8766d", "#e68613","#cd9600","#aba300","#7cae00", "#0cb702", "#00be67", "#00c19a", "#00bfc4", "#00b8e7", "#00a9ff", "#8494ff", "#c77cff", "#ed68ed", "#ff61cc", "#ff68a1")) +
  ylim(0,10000) +
  theme(panel.background = element_blank(),
        axis.line = element_line(colour = 'black', size = 2),
        axis.ticks = element_line(colour = "black", size = 2),
        axis.ticks.length = unit(0.25, "cm"),
        axis.text.y = element_text(angle = 0,size = 10,family = "Arial")) +
  labs(y = "Number of cells")

bar2

ggsave(filename = paste0("../Results/BarPlots/", format(Sys.time(), "%Y%m%d_%H%M%S"), "BarPlot.png"), 
       plot = bar2, width = 20, height = 9,
       path = "../Results/BarPlots/",
       units = "cm",limitsize = FALSE)





# 10. Get Data for spatial mapping into one sample ----------------------------------

# 10.1 Get the data cluster from one sample ---------------------------------------

# Sample A1
seurat_mtx_all_Pestana_filter_subset_sample_A1 <- subset(seurat_mtx_all_Pestana_filter, orig.ident == "A1")

#Sample A2
seurat_mtx_all_Pestana_filter_subset_sample_A2 <- subset(seurat_mtx_all_Pestana_filter, orig.ident == "A2")

#Sample B2
seurat_mtx_all_Pestana_filter_subset_sample_B2 <- subset(seurat_mtx_all_Pestana_filter, orig.ident == "B2")

# Sample C2
seurat_mtx_all_Pestana_filter_subset_sample_C2 <- subset(seurat_mtx_all_Pestana_filter, orig.ident == "C2")

# Sample C2
seurat_mtx_all_Pestana_filter_subset_sample_D1 <- subset(seurat_mtx_all_Pestana_filter, orig.ident == "D1")


# 10.2. Produce .TXT file for spatial mapping into one sample ----------------------

#Sample A1
Cluster_ID_sample_A1 <- seurat_mtx_all_Pestana_filter_subset_sample_A1@active.ident

Cell_Cluster_ID_sample_A1_df <- as.data.frame(Cluster_ID_sample_A1)

Cell_Cluster_ID_sample_A1_df <- tibble::rownames_to_column(Cell_Cluster_ID_sample_A1_df, "Cell_ID")

write.table(Cell_Cluster_ID_sample_A1_df, file="cell_cluster_sample_A1.txt", sep = "\t",
            row.names = FALSE, col.names = FALSE, quote = FALSE)

#Sample A2
Cluster_ID_sample_A2 <- seurat_mtx_all_Pestana_filter_subset_sample_A2@active.ident

Cell_Cluster_ID_sample_A2_df <- as.data.frame(Cluster_ID_sample_A2)

Cell_Cluster_ID_sample_A2_df <- tibble::rownames_to_column(Cell_Cluster_ID_sample_A2_df, "Cell_ID")

write.table(Cell_Cluster_ID_sample_A2_df, file="cell_cluster_sample_A2.txt", sep = "\t",
            row.names = FALSE, col.names = FALSE, quote = FALSE)

#Sample B2
Cluster_ID_sample_B2 <- seurat_mtx_all_Pestana_filter_subset_sample_B2@active.ident

Cell_Cluster_ID_sample_B2_df <- as.data.frame(Cluster_ID_sample_B2)

Cell_Cluster_ID_sample_B2_df <- tibble::rownames_to_column(Cell_Cluster_ID_sample_B2_df, "Cell_ID")

write.table(Cell_Cluster_ID_sample_B2_df, file="cell_cluster_sample_B2.txt", sep = "\t",
            row.names = FALSE, col.names = FALSE, quote = FALSE)

# Sample C2
Cluster_ID_sample_C2 <- seurat_mtx_all_Pestana_filter_subset_sample_C2@active.ident

Cell_Cluster_ID_sample_C2_df <- as.data.frame(Cluster_ID_sample_C2)

Cell_Cluster_ID_sample_C2_df <- tibble::rownames_to_column(Cell_Cluster_ID_sample_C2_df, "Cell_ID")

write.table(Cell_Cluster_ID_sample_C2_df, file="cell_cluster_sample_C2.txt", sep = "\t",
            row.names = FALSE, col.names = FALSE, quote = FALSE)


# Sample D1
Cluster_ID_sample_D1 <- seurat_mtx_all_Pestana_filter_subset_sample_D1@active.ident

Cell_Cluster_ID_sample_D1_df <- as.data.frame(Cluster_ID_sample_D1)

Cell_Cluster_ID_sample_D1_df <- tibble::rownames_to_column(Cell_Cluster_ID_sample_D1_df, "Cell_ID")

write.table(Cell_Cluster_ID_sample_D1_df, file="cell_cluster_sample_D1.txt", sep = "\t",
            row.names = FALSE, col.names = FALSE, quote = FALSE)



