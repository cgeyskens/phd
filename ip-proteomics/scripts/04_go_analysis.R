###############################################################################
# Script: 04_go_analysis.R
# Purpose: GO analysis for GPR37L1 co-IP hits
# Author: Cydric Geyskens
# Date: 2025-10-27
###############################################################################

library(readr)
library(clusterProfiler)
library(dplyr)
library(ggplot2)
library(org.Mm.eg.db)

library(httpgd)
hgd()


#### ======================= loading the GPR37L1 data ======================= ####
ip_data <- read.csv("results/Gpr37l1_limma_ip_proteins_paper.csv", row.names = 1)
View(ip_data)


#### ============================ GO analysis ============================== ####
values_uniprot <- ip_data[["Protein.Group"]]
length(values_uniprot)

go_bp <- enrichGO(gene = values_uniprot,
                OrgDb = org.Mm.eg.db,
                keyType = "UNIPROT",
                ont = "BP", # or CC or MF
                pvalueCutoff = 0.05,
                pAdjustMethod = "BH",
                readable = TRUE)
head(go_bp, 5)


#### ============================== plotting ============================== ####
p <- dotplot(go_bp, showCategory = 5) +
        scale_size(range = c(2, 15),
                  breaks = c(5, 8, 11, 15)) +
        scale_colour_gradient(low = "#ffb454", high = "#734710") +
        scale_fill_gradient(low = "#ffb454", high = "#734710") +
        theme(
          panel.border = element_rect(colour = "black", fill = NA, size = 1.5),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.ticks = element_line(size = 1.2, color = "black"),
          axis.ticks.length = unit(0.2, "cm")
        )
p

ggsave("gpr37l1_go_analysis.svg", 
    plot = p, 
    device = cairo_pdf,
    width = 28, height = 20, units = "cm", dpi=300)


