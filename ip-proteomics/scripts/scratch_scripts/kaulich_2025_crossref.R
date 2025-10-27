library(readxl)
library(dplyr)
library(UpSetR)
install.packages("VennDiagram")
library(VennDiagram)

# load excel files from Kaulich et al. 2025
ca1_proteins <- read_excel("/mnt/ip-proteomics/for_annotations/kaulich_2025_fig1g_CA1_protein.xlsx")
ca3_proteins <- read_excel("/mnt/ip-proteomics/for_annotations/kaulich_2025_fig1g_CA3_protein.xlsx")
dg_proteins <- read_excel("/mnt/ip-proteomics/for_annotations/kaulich_2025_fig1g_DG_protein.xlsx")

# loading ip_data
ip_data <- read.csv("Vcam1_limma_ip_proteins_zenotof_1miscleav_wo_gpr37l1.csv", row.names = 1)
View(ip_data)

# filter proteins
ca1_proteins_filtered <- ca1_proteins %>%
  filter(sig_CA1vCA3 == TRUE & sig_CA1vDG == TRUE)

ca3_proteins_filtered <- ca3_proteins %>%
  filter(sig_CA1vCA3 == TRUE & sig_CA3vDG == TRUE)

dg_proteins_filtered <- dg_proteins %>%
  filter(sig_CA3vDG == TRUE & sig_CA1vDG == TRUE)


# plotting
ca1_list <- pull(ca1_proteins_filtered, Genes)
ca3_list <- pull(ca3_proteins_filtered, Genes)
dg_list <- pull(dg_proteins_filtered, Genes)

ip_list <- pull(ip_data, Genes)
ip_list <- toupper(ip_list)

list_of_sets <- list(
  CA1 = ca1_list,
  CA3 = ca3_list,
  DG = dg_list,
  IP = ip_list
)

upset(fromList(list_of_sets),
    text.scale = 2,
    nsets = 4, 
    nintersects = NA, 
    order.by = "freq")

venn.plot <- venn.diagram(
  x = list_of_sets,
  filename = NULL,
  fill = c("skyblue", "lightgreen", "salmon", "gold"),
  alpha = 0.5,
  cex = 1.5,
  cat.cex = 1.5
)
grid::grid.draw(venn.plot)

intersect(ip_list, ca1_list) 
intersect(ip_list, ca3_list)
intersect(ip_list, dg_list)
