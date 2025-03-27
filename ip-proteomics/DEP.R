library(DEP)
library(dplyr)
library(SummarizedExperiment)

exp17_data <- read.delim("/mnt/data/ip-proteomics/exp17-my-diann-run/diann_output.gg_matrix.tsv")
View(exp17_data)

vcam1_data <- exp17_data[exp17_data$Genes == 'Vcam1',]
View(vcam1_data)

data <- as_tibble(exp17_data)

df <- data %>%
    rename(
        gpr37l1_ip = X.ip.proteomics.ms.convert.output.2ul_CG_1.mzML,
        gpr37l1_igg = X.ip.proteomics.ms.convert.output.2ul_CG_2.mzML,
        vcam1_ip_r1 = X.ip.proteomics.ms.convert.output.2ul_CG_3.mzML,
        vcam1_igg_r1 = X.ip.proteomics.ms.convert.output.2ul_CG_4.mzML,
        vcam1_ip_r2 = X.ip.proteomics.ms.convert.output.2ul_CG_5.mzML,
        vcam1_igg_2 = X.ip.proteomics.ms.convert.output.2ul_CG_6.mzML,
        vcam1_ip_r3 = X.ip.proteomics.ms.convert.output.2ul_CG_7.mzML, 
        vcam1_igg_r3 = X.ip.proteomics.ms.convert.output.2ul_CG_8.mzML, 
        vcam1_ip_r4 = X.ip.proteomics.ms.convert.output.2ul_CG_9.mzML,
        vcam1_igg_r4 = X.ip.proteomics.ms.convert.output.2ul_CG_10.mzML
        ) %>%
        select(-c(gpr37l1_ip, gpr37l1_igg)) %>% # remove the gpr37l1 columns 
        filter(!if_all(-Genes, is.na)) 

dim(df)
View(df)

colnames(df)
rownames(df)
data$Genes %>% duplicated() %>% any()


se <- SummarizedExperiment(assays = list(counts = df[, -1]))
rowData(se) <- DataFrame(Gene = df$Genes)
se

plot_numbers(se)
colData(se)$ID <- colnames(se)
colData(se)


# setting the colname
colData(se)$condition <- c("ip", "igg", "ip", "igg","ip", "igg", "ip", "igg")
colData(se)

# plotting
plot_numbers(se)
plot_coverage(se)

# data formatting for normalization, the values are non-numeric
counts_data <- assay(se, "counts")
View(counts_data)

counts_matrix <- as.matrix(counts_data)
rownames(counts_matrix) <- rownames(counts_data)
colnames(counts_matrix) <- colnames(counts_data)

assay(se, "counts", withDimnames = FALSE) <- counts_matrix
class(assay(se, "counts"))

# normalization
data_norm <- normalize_vsn(se)


non_finite_values <- is.finite(assay(se, "counts"))
non_finite_rows <- which(!rowSums(non_finite_values))
non_finite_cols <- which(!colSums(non_finite_values))

# Print the number of non-finite values
cat("Number of non-finite rows:", length(non_finite_rows), "\n")
cat("Number of non-finite columns:", length(non_finite_cols), "\n")

plot_missval(se)
plot_detect(se)
