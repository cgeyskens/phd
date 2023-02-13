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

#calculate mean of values into new "value" column
all_data5 <- mutate(all_data4, sum_value=A1_value+A2_value+B2_value+C2_value+D1_value)
all_data6 <-mutate(all_data5, value=sum_value/5)
head(all_data6)

#only include rowname and colname
all_data7 <- select(all_data6, rowname, colname, value)
dim(all_data7)
head(all_data7)

#pivor wider format
all_data8 <- pivot_wider(all_data7, id_cols=1, names_from=colname, values_from=value)
all_data8
all_data9 <- column_to_rownames(all_data8, "rowname")
head(all_data9)



#make heatmap--------------------------------------------------------------------------
#The order of you transcripts
#x axis --> first astrocyte markers, then neuron markers, ...
level_order_x <- rev(celltype_markers)

#y-axis --> from high to low astrocyte coloc scores
sort_y <- t(all_data9)
head(sort_y)
sort_y_2 <- as_tibble(sort_y, rownames="rowname")
new <- mutate(sort_y_2, astrocyte_expression=Aldh1l1+Aldoc+Aqp4+Slc1a3 + Slc1a2)
new
sort_y_3 <- arrange(new, desc(astrocyte_expression))

level_order_y <- sort_y_3$rowname
level_order_y


#Specifying Fonts
windowsFonts(A=windowsFont("Arial"))

#make the actual heatmap
g1 <- ggplot(all_data7, aes(x=factor(colname, level=level_order_y), 
                        y=factor(rowname, level=level_order_x), 
                        fill=value)) + geom_tile() +
  coord_fixed() +
  scale_fill_gradientn(name="Colocalization score",
                       limits=c(-1,1), 
                       na.value="yellow",
                       colors=c("purple", "black", "yellow"),
                       breaks=c(-1, 0, 1),
                       labels=c("low", 0, "high")) +
  xlab("") + ylab("") +
  scale_x_discrete(position="top", expand=c(0,0)) +
  scale_y_discrete(expand=c(0,0)) +
  theme(axis.ticks = element_blank(),
        axis.text.x=element_text(family="A", face="italic", color="black", size="2", angle=45, hjust=0, vjust=0),
        axis.text.y=element_text(family="A", face="italic", color="black", size="2.5"),
        axis.title = element_text(vjust=-0.5),
        legend.key.height = unit(5,"mm"),
        legend.key.width=unit(2, "mm"),
        legend.text= element_text(family="A", size="6"),
        legend.title = element_text(family="A", size="6", angle=90),
        legend.box.just="center",
        panel.background = element_blank(),
        panel.border=element_blank())

g1

ggsave("SpatialTranscriptomics_Heatmap_summary_CA3.tiff", width=10, height=10, units = c("cm"), dpi=600)







