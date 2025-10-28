library(readr)
library(clusterProfiler)
library(dplyr)
library(ggplot2)
library(grid) 

### ============================ GO analysis of GPR37L1 data ============================= ###

#BiocManager::install("org.Mm.eg.db")
library(org.Mm.eg.db)

# loading in the data
ip_data <- read.csv("results/Gpr37l1_limma_ip_proteins_paper.csv", row.names = 1)
View(ip_data)

# Clusterprofiler analysis 
values_uniprot <- ip_data[["Protein.Group"]]
length(values_uniprot)

go_bp <- enrichGO(gene = values_uniprot,
                OrgDb = org.Mm.eg.db,
                keyType = "UNIPROT",
                ont = "BP", # or CC or MF
                pvalueCutoff = 0.05,
                pAdjustMethod = "BH",
                #universe = values_background_uniprot, # specifiying background list, if not just the entire genome
                readable = TRUE)
head(go_bp, 8)
as.data.frame(go_bp)

go_cc <- enrichGO(gene = values_uniprot,
                OrgDb = org.Mm.eg.db,
                keyType = "UNIPROT",
                ont = "CC", # or CC or MF
                pvalueCutoff = 0.05,
                pAdjustMethod = "BH",
                #universe = values_background_uniprot, # specifiying background list, if not just the entire genome
                readable = TRUE)
head(summary(go_cc))

go_mf <- enrichGO(gene = values_uniprot,
                OrgDb = org.Mm.eg.db,
                keyType = "UNIPROT",
                ont = "MF", # or CC or MF
                pvalueCutoff = 0.05,
                pAdjustMethod = "BH",
                #universe = values_background_uniprot, # specifiying background list, if not just the entire genome
                readable = TRUE)
head(summary(go_mf))

p <- dotplot(go_bp, showCategory=8)
p
p <- barplot(go_bp, 
        x = "Count",
        drop = TRUE, 
        showCategory = 10, 
        title = "GO Biological Pathways",
        font.size = 8)
p

p <- dotplot(go_bp, showCategory = 5) +
        scale_size(range = c(2, 15),
                  breaks = c(5, 8, 11, 15)) +
        scale_colour_gradient(low = "#ffb454", high = "#734710") +
        scale_fill_gradient(low = "#ffb454", high = "#734710") +
        guides(
          size = guide_legend(
            override.aes = list(
              colour = "black",  # legend dots in black
              fill = "black",
              stroke = 0,        # remove border around legend dots
              shape = 16,
              alpha = 1
            )
          ),
          colour = guide_colorbar()
        ) +
        theme(
          panel.border = element_rect(colour = "black", fill = NA, size = 1.5),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.ticks = element_line(size = 1.2, color = "black"),
          axis.ticks.length = unit(0.2, "cm")
        )

p
