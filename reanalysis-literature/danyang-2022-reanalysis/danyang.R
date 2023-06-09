install.packages("tximport")
install.packages("readr")
library(tximport)
library(readr)
library(DESeq2)
library(tximeta)

# Need to use EdgeR instead DESeq2, look at the vignette online

# set wd
setwd("D:/code/reanalysis-literature/Danyang-He-2022")

# loading the counts
counts <- read_delim("GSE127449_rsem.genes.counts.matrix", 
                      delim = "\t", escape_double = FALSE, 
                      trim_ws = TRUE)

# loading the metadata of the columns
metadata_tau_gfp <- read_table("metadata_tau-gfp.txt")
View(metadata_tau_gfp)

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

# loading into DeSeq object
dd4 <- DESeqDataSetFromMatrix(countData = round(df),
                              colData = metadata_tau_gfp,
                              design = ~mice + microglia) # error: because of genenames in the first column

# filter genes out with less the 10 counts
dds <- dd4[rowSums(counts(dd4)) >= 10, ]
nrow(counts(dds))

# Load into the DESeq object
# we need the not the wald test but the LRT test because the microglia comes from the brains
dds1 <- DESeq(dds, test="LRT", full = ~mice+microglia, reduced = ~mice) 

# check the dispersion plot
plotDispEsts(dds1)

# Checking the comparison
resultsNames(dds1)

res <- results(dds1,name="microglia_Tau.GFP.Plus_vs_Tau.GFP.Neg")
head(res,n=10)

head(res[order(res$padj),])
head(res[order(res$log2FoldChange),])


#shrink the logfoldchange
res.shr <- lfcShrink(dds1,coef=9)
dim(res.shr)

#removing genes with missing padj values, there were no missing padj values
resFix <- res.shr[!is.na(res.shr$padj),]
dim(resFix)

# shrinking the log2foldchange
resFix[order(resFix$log2FoldChange),]
res["Itga9",]
resFix["Itga9",]
head(res)
head(resFix)
summary(res)
summary(resFix)

# checking the DEG genes
up <- (resFix$log2FoldChange > 0.5 & resFix$padj < 0.05)
down <- (resFix$log2FoldChange < -0.5 & resFix$padj < 0.05)

sum(down)
sum(up)

# change to df for volcano plot
res.shr.df <- as.data.frame(resFix)

# volcano plot
ggplot(res.shr.df,aes(log2FoldChange,-log10(padj))) + 
  geom_point(shape=1,aes(color=DE)) +
  scale_color_manual(name="differential expression",
                     values=c("black","red","green"),
                     labels=c("|lfc|<1 or p>0.01",
                              ">2 fold UP and p<0.01",
                              ">2 fold DOWN and p<0.01")) +
  ggtitle("p-value versus fold change") + 
  geom_hline(yintercept=2,color="red",linetype=2) +
  geom_vline(xintercept=-1,color="red",linetype=2) +
  geom_vline(xintercept=1,color="red",linetype=2)

