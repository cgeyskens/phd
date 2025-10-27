install.packages("VennDiagram")
library(VennDiagram)

# load the ip data
ip_data_vcam1 <- read.csv("Vcam1_limma_ip_proteins_zenotof_1miscleav_wo_gpr37l1.csv", row.names = 1)
ip_data_gpr37l1 <- read.csv("Gpr37l1_limma_ip_proteins_zenotof_1miscleav.csv", row.names = 1)

ip_list_vcam1 <- pull(ip_data_vcam1, Genes)
ip_list_gpr37l1 <- pull(ip_data_gpr37l1, Genes)

list_of_sets <- list(
  VCAM1 = ip_list_vcam1,
  GPR37L1 = ip_list_gpr37l1
)

venn <- venn.diagram(
  x = list_of_sets,
  filename = "vcam1_gpr37l1_venn_plot.png",
  fill = c("skyblue", "lightgreen"),
  alpha = 0.5,
  cex = 1.5,
  cat.cex = 1.5
)
intersect(ip_list_vcam1, ip_list_gpr37l1) 
