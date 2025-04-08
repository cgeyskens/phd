library(samr)
library(dplyr)
library(ggplot2)


# Loading the data
exp17_data <- read_tsv("/mnt/data/ip-proteomics/exp17-my-diann-run/diann_output.pg_matrix.tsv")
View(exp17_data)

# checking colnames
colnames(exp17_data)

# gpr37l1 colnames
exp17_names <- c(
    "gpr37l1_ip" = "/ip-proteomics/ms-convert-output/2ul_CG_1.mzML",
    "gpr37l1_igg" = "/ip-proteomics/ms-convert-output/2ul_CG_2.mzML",
    "vcam1_ip_r1" = "/ip-proteomics/ms-convert-output/2ul_CG_3.mzML",
    "vcam1_igg_r1"= "/ip-proteomics/ms-convert-output/2ul_CG_4.mzML",
    "vcam1_ip_r2" = "/ip-proteomics/ms-convert-output/2ul_CG_5.mzML",
    "vcam1_igg_r2"= "/ip-proteomics/ms-convert-output/2ul_CG_6.mzML",
    "vcam1_ip_r3" = "/ip-proteomics/ms-convert-output/2ul_CG_7.mzML",
    "vcam1_igg_r3"= "/ip-proteomics/ms-convert-output/2ul_CG_8.mzML",
    "vcam1_ip_r4" = "/ip-proteomics/ms-convert-output/2ul_CG_9.mzML",
    "vcam1_igg_r4" = "/ip-proteomics/ms-convert-output/2ul_CG_10.mzML"
)

colnames(exp17_data)[match(exp17_names, colnames(exp17_data))] <- names(exp17_names)


# removing the gpr37l1 columns
df <- exp17_data %>%
    as_tibble() %>%
    select(-c(
      "gpr37l1_ip" ,
      "gpr37l1_igg"
      ) # remove the gpr37l1 columns
    ) %>%
    mutate(
        Genes = if_else(is.na(Genes), Protein.Group, Genes) # If the Genes column has NA values, replace them with values from Group.Protein
    ) %>%
    column_to_rownames(var = "Genes") %>% # The Genes column as rownames
    select(-Protein.Group, -Protein.Names, -First.Protein.Description) %>%  # remove these columns
    filter(rowSums(is.na(.)) < ncol(.)) # removes the rows if they have NAs in all columns

### log2 transformation
df_log2 <- df %>%
  mutate(across(where(is.numeric), log2))

### imputation
data_matrix <- as.matrix(df_log2)
imputed_matrix <- impute.MinProb(dataSet.mvs = data_matrix,
                                      q = 0.01,      
                                      tune.sigma = 1) 

### SAMR

# Group design
groups = c("ip","igg","ip","igg","ip","igg","ip","igg")
groups = c(1, 2, 1, 2, 1, 2, 1, 2)
groups

# Running samr
sam_result <- SAM(x = imputed_matrix,
                y = groups,
                resp.type = "Two class unpaired",
                geneid = rownames(imputed_matrix),
                genenames = rownames(imputed_matrix),
                #s0 = 0.1,
                s0.perc = TRUE,
                nperms = 1000,
                center.arrays = FALSE,
                testStatistic = "standard",
                return.x = TRUE,
                random.seed = 1234,
                logged2 = TRUE,
                fdr.output = 0.1,
                eigengene.number = 1
                )

names(sam_result$samr.obj)


## qq plot of samr
samr::samr.plot(sam_result$samr.obj,
            del = 0.45, 
             min.foldchange = 1.5)


## according to script of Pang et al. 
logFC<-log2(sam_result$samr.obj$foldchange)
pvalue=data.frame(samr.pvalues.from.perms(sam_result$samr.obj$tt, sam_result$samr.obj$ttstar))

pvalue["Vcam1", ]
adj.pvalue=p.adjust(pvalue$samr.pvalues.from.perms.sam_result.samr.obj.tt..sam_result.samr.obj.ttstar., method = "BH")

test_results <- cbind(row.names(imputed_matrix), LogFC, pvalue, adj.pvalue)

logFC <- log2(sam_result$samr.obj$foldchange)
pvalue <- data.frame(p.value = samr.pvalues.from.perms(sam_result$samr.obj$tt, sam_result$samr.obj$ttstar))
pvalue["Vcam1", ]
adj.pvalue <- data.frame(adj.p.value = p.adjust(pvalue$p.value, method = "BH"))
test_results <- cbind(row.names(imputed_matrix), logFC, pvalue, adj.pvalue)
colnames(test_results)[1] <- "Gene"

test_results["Chgb", ]

test_results <- test_results %>%
  mutate(neg.log10.p.value = -log10(p.value))


# make the volcano plot

genes_to_highlight <- c("Vcam1", "Prrt1", "Chgb") # Replace with your gene names
highlight_colors <- c("blue", "green", "purple") # Corresponding colors

volcano_plot <- ggplot(test_results, aes(x = logFC, y = neg.log10.p.value, label = Gene)) +
  geom_point(alpha = 0.6) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") + # 
  labs(
    title = "Volcano Plot from SAM Analysis",
    x = "Log2 Fold Change",
    y = "-log10(p-value)"
  ) +
  theme_minimal()

for (i in seq_along(genes_to_highlight)) {
  gene <- genes_to_highlight[i]
  color <- highlight_colors[i]

  volcano_plot <- volcano_plot +
    geom_point(data = subset(test_results, Gene == gene),
               color = color,
               size = 3) +
    geom_text(data = subset(test_results, Gene == gene),
              nudge_y = 0.015 * i, # Adjust nudge for each label
              nudge_x = 0.015 * i,
              size = 3,
              color = color)
}

print(volcano_plot)

