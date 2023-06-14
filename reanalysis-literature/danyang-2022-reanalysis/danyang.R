install.packages("tximport")
install.packages("readr")
library(tximport)
library(readr)

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("edgeR")
library(edgeR)

# Need to use EdgeR instead DESeq2, look at the tutorial from VIB bulk RNAseq online. 
# The reads I got from Danyang are probably normalized by readlength, so the analysis 
# don't need any readlength normalization.
# EdgeR will also get more differentially expressed genes compared to DESeq2
# the publication: https://pubmed.ncbi.nlm.nih.gov/34982959/

# setting wd, loading the data and filtering the data ####

# set wd
setwd("D:/code/reanalysis-literature/Danyang-He-2022")

# loading the counts
counts <- read_delim("GSE127449_rsem.genes.counts.matrix", 
                      delim = "\t", escape_double = FALSE, 
                      trim_ws = TRUE)

# loading the metadata of the columns
metadata_tau_gfp <- read_table("metadata_tau-gfp.txt")

# converting to counts table to a df 
counts_df <- as.data.frame(counts)

# select only the relevant samples (column names)
df <- counts_df[,c("...1", "E1_S1","E2_S2", "E3_S3", "E4_S4", "E5_S5", "E6_S6", "E7_S7",
                   "E8_S8", "E9_S9", "E10_S10", "E11_S11", "E12_S12", "F1_S13",
                   "F2_S14", "F3_S15", "F4_S16")]

#set gene names as rownames
#rownames contains duplicates (because of excel Mar-1 & Mar-2)
length(counts_df[,1])
length(unique(counts_df[,1]))
df<-df[!duplicated(counts_df[,1]),]
rownames(df) <- df[,1] 
head(df)
df<-df[,-1]

head(df)
class(df)
class(metadata_tau_gfp)



# EdgeR ####

# loading data into edgeR object
dge <- DGEList(counts=df, samples = metadata_tau_gfp, group = metadata_tau_gfp$microglia)
dim(dge)
head(dge$counts)
dge$samples

#filter out genes with less then 10 counts in total
dge <- dge[(rowSums(dge$counts)) >= 10,]
dim(dge)

# library size normalization
dge$samples$lib.size
dge <- calcNormFactors(dge)

## DGE analysis
design <- model.matrix(~mice + microglia, dge$samples) 
design

# Estimate the dispersions
dge <- estimateDisp(dge,design=design)

# fit the generalized linear model and perform the DGE test
fit <- glmFit(dge,design)

design

# retrieve results from dge test
lrt <- glmLRT(fit,coef=9)

# retrieving the DGE genes
tt <- topTags(lrt,n=nrow(dge),p.value=0.05)

# retrieving all genes into a table
tt.all <- topTags(lrt,n=nrow(dge))

# checking the genes shown in the manuscripts at Figure S2A
tt$table["Igf1",]
tt$table["Slc2a5",]
tt$table["Cx3cr1",]
tt$table["Spp1",]
tt$table["Anxa5",]

# setting the log fold change threshold
tt.String <- tt$table[abs(tt$table$logFC) >= 1,]




# plotting ####

# create two boolean vectors with all the results for the plotting
up <- (tt.all$table$logFC > 1 & -log10(tt.all$table$FDR) > 2)
down <- (tt.all$table$logFC < -1 & -log10(tt.all$table$FDR) > 2)

# create a new column called DE in the tt.all table that says whether the gene up or down regulated
tt.all$table$DE <- ifelse(up,'up',ifelse(down,'down','not DE'))

# create the volcano plot
library(ggplot2)
ggplot(tt.all$table,aes(logFC,-log10(FDR))) + 
  geom_point(shape=1,aes(color=DE)) +
  scale_color_manual(name="differential expression",
                     values=c("green","black","red"),
                     labels=c(">2 fold DOWN and p<0.01",
                              "|lfc|<1 or p>0.01",
                              ">2 fold UP and p<0.01")) +
  ggtitle("p-value versus fold change") + 
  geom_hline(yintercept=2,color="red",linetype=2) +
  geom_vline(xintercept=-1,color="red",linetype=2) +
  geom_vline(xintercept=1,color="red",linetype=2)












