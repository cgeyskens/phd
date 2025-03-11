#Join data and calculate average value for heatmap------------------------------------------
# inner_join all tables
head(A1_data4)
dim(A1_data4)
A1_data4

head(A2_data4)
dim(A2_data4)

head(B2_data4)
dim(B2_data4)

head(C2_data4)
dim(C2_data4)

head(D1_data4)
dim(D1_data4)

A1_A2_data4 <- inner_join(A1_data4, A2_data4, by=c("rowname", "colname"))
A1_A2_data4
dim(A1_A2_data4)

A1_A2_B2_data4 <- inner_join(A1_A2_data4, B2_data4, by=c("rowname", "colname"))
dim(A1_A2_B2_data4)
head(A1_A2_B2_data4)

A1_A2_B2_C2_data4 <- inner_join(A1_A2_B2_data4, C2_data4, by=c("rowname", "colname"))
dim(A1_A2_B2_C2_data4)

all_data4 <- inner_join(A1_A2_B2_C2_data4, D1_data4, by=c("rowname", "colname"))
head(all_data4)
dim(all_data4)

# Performing the pearsons correlation test between each samples (10 times) ------------------
cor.test(all_data4$A1_value, all_data4$A2_value)

cor.test(all_data4$A1_value, all_data4$B2_value)

cor.test(all_data4$A1_value, all_data4$C2_value)

cor.test(all_data4$A1_value, all_data4$D1_value)

cor.test(all_data4$A2_value, all_data4$B2_value)

cor.test(all_data4$A2_value, all_data4$C2_value)

cor.test(all_data4$A2_value, all_data4$D1_value)

cor.test(all_data4$B2_value, all_data4$C2_value)

cor.test(all_data4$B2_value, all_data4$D1_value)

cor.test(all_data4$C2_value, all_data4$D1_value)

# Make the heatmap plot across all samples -------------------------------------------

install.packages("ggcorrplot")
library(ggcorrplot)

head(all_data4)
dim(all_data4)

#combine and convert the first two columns of all_data4 as rownames
rownames(all_data4) <- do.call(paste,c(all_data4[c("rowname","colname")], sep="_"))
all_data4_col <- all_data4[,!names(all_data4) %in% c("rowname", "colname")]
head(all_data4_col)
dim(all_data4_col)

class(all_data4_col)

# Compute correlation matrix

corr <- cor(all_data4_col)
corr

pval <- cor_pmat(all_data4_col)
head(pval)

corrplot <- ggcorrplot(corr, hc.order=FALSE, type="lower", 
            colors = c("#266bae","grey","#bc373f"), lab=TRUE)

corrplot

ggsave(filename = paste0("../Results/CorrPlots/", format(Sys.time(), "%Y%m%d_%H%M%S"), "CorrPlot.png"), 
       plot = corrplot, width = 40, height = 40, dpi = 700,
       path = "../Results/CorrPlots/",
       units = "cm",limitsize = FALSE)


# Make the scatter plots with regression line ---------------------------------------

library("ggpubr")

# Sample A1 - A2
scatter_plot_A1_A2 <- ggscatter(all_data4, x="A2_value", y="A1_value",
                                add="reg.line", conf.int=FALSE,
                                cor.coef=FALSE, cor.method="pearson",
                                xlab="Sample A2", ylab="Sample A1",
                                color = "#929292", shape=1, size=1.5,
                                add.params = list(color="black")) +
                              theme(axis.line = element_line(size=0),
                                    panel.border = element_rect(colour = "black", fill = NA, size=2),
                                    axis.ticks = element_line(colour = "black", size = 1),
                                    axis.ticks.length = unit(0.25, "cm"))
scatter_plot_A1_A2

ggsave(filename = paste0("../Results/Scatterplots/", format(Sys.time(), "%Y%m%d_%H%M%S"), "ScatterPlot_A1_A2.png"), 
       plot = scatter_plot_A1_A2, width = 10, height = 10,
       path = "../Results/Scatterplots/",
       units = "cm",limitsize = FALSE)


# Sample A1 - B2
scatter_plot_A1_B2 <- ggscatter(all_data4, x="B2_value", y="A1_value",
                                add="reg.line", conf.int=FALSE,
                                cor.coef=FALSE, cor.method="pearson",
                                xlab="Sample B2", ylab="Sample A1",
                                color = "#929292", shape=1, size=1.5,
                                add.params = list(color="black")) +
  theme(axis.line = element_line(size=0),
        panel.border = element_rect(colour = "black", fill = NA, size=2),
        axis.ticks = element_line(colour = "black", size = 1),
        axis.ticks.length = unit(0.25, "cm"))
scatter_plot_A1_B2

ggsave(filename = paste0("../Results/Scatterplots/", format(Sys.time(), "%Y%m%d_%H%M%S"), "ScatterPlot_A1_B2.png"), 
       plot = scatter_plot_A1_B2, width = 10, height = 10,
       path = "../Results/Scatterplots/",
       units = "cm",limitsize = FALSE)


# Sample A1 - C2
scatter_plot_A1_C2 <- ggscatter(all_data4, x="C2_value", y="A1_value",
                                add="reg.line", conf.int=FALSE,
                                cor.coef=FALSE, cor.method="pearson",
                                xlab="Sample C2", ylab="Sample A1",
                                color = "#929292", shape=1, size=1.5,
                                add.params = list(color="black")) +
  theme(axis.line = element_line(size=0),
        panel.border = element_rect(colour = "black", fill = NA, size=2),
        axis.ticks = element_line(colour = "black", size = 1),
        axis.ticks.length = unit(0.25, "cm"))
scatter_plot_A1_C2

ggsave(filename = paste0("../Results/Scatterplots/", format(Sys.time(), "%Y%m%d_%H%M%S"), "ScatterPlot_A1_C2.png"), 
       plot = scatter_plot_A1_C2, width = 10, height = 10,
       path = "../Results/Scatterplots/",
       units = "cm",limitsize = FALSE)

# Sample A1 - D1
scatter_plot_A1_D1 <- ggscatter(all_data4, x="D1_value", y="A1_value",
                                add="reg.line", conf.int=FALSE,
                                cor.coef=FALSE, cor.method="pearson",
                                xlab="Sample D1", ylab="Sample A1",
                                color = "#929292", shape=1, size=1.5,
                                add.params = list(color="black")) +
  theme(axis.line = element_line(size=0),
        panel.border = element_rect(colour = "black", fill = NA, size=2),
        axis.ticks = element_line(colour = "black", size = 1),
        axis.ticks.length = unit(0.25, "cm"))
scatter_plot_A1_D1

ggsave(filename = paste0("../Results/Scatterplots/", format(Sys.time(), "%Y%m%d_%H%M%S"), "ScatterPlot_A1_D1.png"), 
       plot = scatter_plot_A1_D1, width = 10, height = 10,
       path = "../Results/Scatterplots/",
       units = "cm",limitsize = FALSE)


# Sample A2 - B2
scatter_plot_A2_B2 <- ggscatter(all_data4, x="B2_value", y="A2_value",
                                add="reg.line", conf.int=FALSE,
                                cor.coef=FALSE, cor.method="pearson",
                                xlab="Sample B2", ylab="Sample A2",
                                color = "#929292", shape=1, size=1.5,
                                add.params = list(color="black")) +
  theme(axis.line = element_line(size=0),
        panel.border = element_rect(colour = "black", fill = NA, size=2),
        axis.ticks = element_line(colour = "black", size = 1),
        axis.ticks.length = unit(0.25, "cm"))
scatter_plot_A2_B2

ggsave(filename = paste0("../Results/Scatterplots/", format(Sys.time(), "%Y%m%d_%H%M%S"), "ScatterPlot_A2_B2.png"), 
       plot = scatter_plot_A2_B2, width = 10, height = 10,
       path = "../Results/Scatterplots/",
       units = "cm",limitsize = FALSE)


# Sample A2 - C2
scatter_plot_A2_C2 <- ggscatter(all_data4, x="C2_value", y="A2_value",
                                add="reg.line", conf.int=FALSE,
                                cor.coef=FALSE, cor.method="pearson",
                                xlab="Sample C2", ylab="Sample A2",
                                color = "#929292", shape=1, size=1.5,
                                add.params = list(color="black")) +
  theme(axis.line = element_line(size=0),
        panel.border = element_rect(colour = "black", fill = NA, size=2),
        axis.ticks = element_line(colour = "black", size = 1),
        axis.ticks.length = unit(0.25, "cm"))
scatter_plot_A2_C2

ggsave(filename = paste0("../Results/Scatterplots/", format(Sys.time(), "%Y%m%d_%H%M%S"), "ScatterPlot_A2_C2.png"), 
       plot = scatter_plot_A2_C2, width = 10, height = 10,
       path = "../Results/Scatterplots/",
       units = "cm",limitsize = FALSE)


# Sample A2 - D1
scatter_plot_A2_D1 <- ggscatter(all_data4, x="D1_value", y="A2_value",
                                add="reg.line", conf.int=FALSE,
                                cor.coef=FALSE, cor.method="pearson",
                                xlab="Sample D1", ylab="Sample A2",
                                color = "#929292", shape=1, size=1.5,
                                add.params = list(color="black")) +
  theme(axis.line = element_line(size=0),
        panel.border = element_rect(colour = "black", fill = NA, size=2),
        axis.ticks = element_line(colour = "black", size = 1),
        axis.ticks.length = unit(0.25, "cm"))
scatter_plot_A2_D1

ggsave(filename = paste0("../Results/Scatterplots/", format(Sys.time(), "%Y%m%d_%H%M%S"), "ScatterPlot_A2_D1.png"), 
       plot = scatter_plot_A2_D1, width = 10, height = 10,
       path = "../Results/Scatterplots/",
       units = "cm",limitsize = FALSE)

# Sample B2 - C2
scatter_plot_B2_C2 <- ggscatter(all_data4, x="C2_value", y="B2_value",
                                add="reg.line", conf.int=FALSE,
                                cor.coef=FALSE, cor.method="pearson",
                                xlab="Sample C2", ylab="Sample B2",
                                color = "#929292", shape=1, size=1.5,
                                add.params = list(color="black")) +
  theme(axis.line = element_line(size=0),
        panel.border = element_rect(colour = "black", fill = NA, size=2),
        axis.ticks = element_line(colour = "black", size = 1),
        axis.ticks.length = unit(0.25, "cm"))
scatter_plot_B2_C2

ggsave(filename = paste0("../Results/Scatterplots/", format(Sys.time(), "%Y%m%d_%H%M%S"), "ScatterPlot_B2_C2.png"), 
       plot = scatter_plot_B2_C2, width = 10, height = 10,
       path = "../Results/Scatterplots/",
       units = "cm",limitsize = FALSE)

# Sample B2 - D1
scatter_plot_B2_D1 <- ggscatter(all_data4, x="D1_value", y="B2_value",
                                add="reg.line", conf.int=FALSE,
                                cor.coef=FALSE, cor.method="pearson",
                                xlab="Sample D1", ylab="Sample B2",
                                color = "#929292", shape=1, size=1.5,
                                add.params = list(color="black")) +
  theme(axis.line = element_line(size=0),
        panel.border = element_rect(colour = "black", fill = NA, size=2),
        axis.ticks = element_line(colour = "black", size = 1),
        axis.ticks.length = unit(0.25, "cm"))
scatter_plot_B2_D1

ggsave(filename = paste0("../Results/Scatterplots/", format(Sys.time(), "%Y%m%d_%H%M%S"), "ScatterPlot_B2_D1.png"), 
       plot = scatter_plot_B2_D1, width = 10, height = 10,
       path = "../Results/Scatterplots/",
       units = "cm",limitsize = FALSE)

# Sample C2 - D1
scatter_plot_C2_D1 <- ggscatter(all_data4, x="D1_value", y="C2_value",
                                add="reg.line", conf.int=FALSE,
                                cor.coef=FALSE, cor.method="pearson",
                                xlab="Sample D1", ylab="Sample C2",
                                color = "#929292", shape=1, size=1.5,
                                add.params = list(color="black")) +
  theme(axis.line = element_line(size=0),
        panel.border = element_rect(colour = "black", fill = NA, size=2),
        axis.ticks = element_line(colour = "black", size = 1),
        axis.ticks.length = unit(0.25, "cm"))
scatter_plot_C2_D1

ggsave(filename = paste0("../Results/Scatterplots/", format(Sys.time(), "%Y%m%d_%H%M%S"), "ScatterPlot_C2_D1.png"), 
       plot = scatter_plot_C2_D1, width = 10, height = 10,
       path = "../Results/Scatterplots/",
       units = "cm",limitsize = FALSE)
