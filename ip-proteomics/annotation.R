library(gprofiler2)
library(clusterProfiler)
library(readr)
library(RSQLite)
library(DBI)

# loading in the raw data
raw_data <- read_tsv("/mnt/ip-proteomics/exp19-my-diann-run/exp19-diann_output.pg_matrix.tsv")
View(raw_data)

# Loading in the results
ip_data <- read.csv("gpr37l1_limma_ip_proteins.csv", row.names = 1)
ip_data$Genes <- rownames(ip_data)
view(ip_data)

# merging to get protein IDs
merged_df <- merge(ip_data,
                    raw_data[, c("Genes", "Protein.Group")],
                    by.x = "Genes",
                    by.y = "Genes",
                    all.x = TRUE)





##### Crossreference with van Oostrum et al. 2023 dataset
van_oostrum_2023 <- read_tsv("/mnt/ip-proteomics/for_annotations/van_oostrum_2023.tsv")

# Get the values from the specified column in dataset A
values_a <- merged_df[["Protein.Group"]]

# Get the unique values from the specified column in dataset B
unique_values_b <- unique(van_oostrum_2023[["protein"]])

# Find the matches
matches <- values_a %in% unique_values_b

# Count the number of matches
num_matches <- sum(matches)

# Calculate the total number of rows in dataset A
total_rows_a <- nrow(merged_df)

# Calculate the percentage of matching rows
oostrum_percentage_matches <- (num_matches / total_rows_a) * 100




##### Crossreference with Sorokina et al. 2021

# connect to sql data
connection <- dbConnect(
    drv = SQLite(),
    dbname = "/mnt/ip-proteomics/for_annotations/synaptic.proteome_SR_20210408.db.sqlite"
)
tbl(connection, "FullGenefullPaper") %>% collect() -> sorokina_data
head(sorokina_data)
dbListTables(connection)

# Get the values from the specified column in dataset A
values_a <- merged_df[["Genes"]]

# Get the unique values from the specified column in dataset B
unique_values_b <- unique(sorokina_data[["MouseName"]])

# Find the matches
matches <- values_a %in% unique_values_b

# Count the number of matches
num_matches <- sum(matches)

# Calculate the total number of rows in dataset A
total_rows_a <- nrow(merged_df)

# Calculate the percentage of matching rows
sorokina_percentage_matches <- (num_matches / total_rows_a) * 100

# now calculate how many pre
pre_synaptic_sorokina <- sorokina_data %>% filter(Localisation == "Presynaptic")
unique_values_b_pre <- unique(pre_synaptic_sorokina[["MouseName"]])
matches <- values_a %in% unique_values_b_pre
num_matches <- sum(matches)
total_rows_a <- nrow(merged_df)
sorokina_percentage_presynaptic_matches <- (num_matches / total_rows_a) * 100

# now calculate how many post
post_synaptic_sorokina <- sorokina_data %>% filter(Localisation == "Postsynaptic")
unique_values_b_post <- unique(post_synaptic_sorokina[["MouseName"]])
matches <- values_a %in% unique_values_b_post
num_matches <- sum(matches)
total_rows_a <- nrow(merged_df)
sorokina_percentage_postsynaptic_matches <- (num_matches / total_rows_a) * 100

# now calculate how many synaptosome
synaptosome_synaptic_sorokina <- sorokina_data %>% filter(Localisation == "Synaptosome")
unique_values_b_synaptosome <- unique(synaptosome_synaptic_sorokina[["MouseName"]])
matches <- values_a %in% unique_values_b_synaptosome
num_matches <- sum(matches)
total_rows_a <- nrow(merged_df)
sorokina_percentage_synaptosome_matches <- (num_matches / total_rows_a) * 100

dbDisconnect(connection)



###### SynGO analysis, manually in site
values_a <- merged_df[["Genes"]]
output_string <- paste(values_a, collapse = ", ")

SynGO_percentage_matches = (42 / total_rows_a) * 100


##### make the graph
categories <- c(
    "Van Oostrum 2023",
    "SynGO",
    "Sorokina 2021, postsynaptic",
    "Sorokina 2021, presynaptic", 
    "Sorokina 2021, synaptosome",
    "Sorokina 2021, in database"
    )

percentages <- c(
    oostrum_percentage_matches,
    SynGO_percentage_matches,
    sorokina_percentage_postsynaptic_matches,
    sorokina_percentage_presynaptic_matches,
    sorokina_percentage_synaptosome_matches, 
    sorokina_percentage_matches
    )

# bar colors
bar_colors <- ifelse(categories == "SynGO", "#ffba19",
                   ifelse(categories == "Van Oostrum 2023", "#c2b51c",
                          "#21a0e2"))
barplot(percentages,
        names.arg = NA,
        horiz = TRUE,
        xlim = c(0, 100),
        xlab = "% of co-immunoprecipitated proteins",
        col = bar_colors,  # You can choose a different color
        border = "black",
        las = 1, # Make axis labels horizontal
        axes = FALSE
)
axis(side = 1,
     lwd = 3
)
y_label_y_positions <- c(0.5, 1.7, 2.9, 4.1, 5.3, 6.5) 
text(x = 7, # Adjust the x-coordinate to position labels inside
     y = y_label_y_positions,
     labels = categories,
     cex = 2.5,      # Make the text bigger (adjust as needed)
     col = "black",
     adj = c(0, 0.1),        # Left-align the text
     family = "Arial" # Set the font family to Arial
)

