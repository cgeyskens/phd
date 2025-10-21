library(readr)
library(clusterProfiler)
library(dplyr)


BiocManager::install("org.Mm.eg.db")
library(org.Mm.eg.db)

# loading in the data
ip_data <- read.csv("Gpr37l1_limma_ip_proteins_zenotof_1miscleav.csv", row.names = 1)
View(ip_data)

##### Clusterprofiler analysis 
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
head(summary(go_bp))

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

dotplot(go_mf, showCategory=10)


genes <- bitr(
  values_uniprot,
  fromType = "UNIPROT",
  toType = c("ENTREZID", "SYMBOL", "ENSEMBL"),
  OrgDb = org.Mm.eg.db
)


ggo <- groupGO(gene     = values_uniprot,
               OrgDb    = org.Mm.eg.db,
               ont      = "CC",
               level    = 3,
               readable = TRUE)

head(ggo)
