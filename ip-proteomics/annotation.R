library(gprofiler2)
library(clusterProfiler)
library(readr)
library(RSQLite)
library(DBI)
library(stringr)
library(dplyr)
library(UniprotR)
library(readxl)

# loading in the raw data
raw_data <- read_tsv("/mnt/ip-proteomics/exp19-my-diann-run/exp19-diann_output.pg_matrix.tsv")
View(raw_data)

# Loading in the results
ip_data <- read.csv("gpr37l1_limma_ip_proteins.csv", row.names = 1)
ip_data$Genes <- rownames(ip_data)
View(ip_data)

# merging to get protein IDs
merged_df <- merge(ip_data,
                    raw_data[, c("Genes", "Protein.Group")],
                    by.x = "Genes",
                    by.y = "Genes",
                    all.x = TRUE)


###### function to check presence in datasets
check_presence_in_reference <- function(df_to_modify, column_to_search, reference_values, new_column_name) {
  # Ensure reference_values is a character vector
  reference_values <- as.character(reference_values)

  # Apply the string detection and create the new column
  df_to_modify %>%
    mutate(!!new_column_name := sapply(.data[[column_to_search]], function(text_to_check) {
      any(str_detect(text_to_check, reference_values))
    }))
}

###### function to check sorokina presynaptic/postsynaptic/synaptosome
is_localised <- function(mouse_name, sorokina_df, localisation_type) {
  localised_proteins <- sorokina_df %>%
    filter(Localisation == localisation_type) %>%
    pull(MouseName)
  mouse_name %in% localised_proteins
}


##### Crossreference with van Oostrum et al. 2023 dataset
# loading in the data
van_oostrum_2023 <- read_tsv("/mnt/ip-proteomics/for_annotations/van_oostrum_2023.tsv")

# crossrefence
merged_df_updated_oostrum <- merged_df %>%
  check_presence_in_reference(
    column_to_search = "Protein.Group",
    reference_values = van_oostrum_2023$protein,
    new_column_name = "in_van_oostrum_2023"
  )


##### Crossreference with Sorokina et al. 2021
# connect to sql data
connection <- dbConnect(
    drv = SQLite(),
    dbname = "/mnt/ip-proteomics/for_annotations/synaptic.proteome_SR_20210408.db.sqlite"
)
tbl(connection, "FullGenefullPaper") %>% collect() -> sorokina_data
head(sorokina_data)
dbListTables(connection)
dbDisconnect(connection)

# crossrefence for in.sorokina.database
merged_df_updated_sorokina <- merged_df_updated_oostrum %>%
  check_presence_in_reference(
    column_to_search = "Genes",
    reference_values = sorokina_data$MouseName,
    new_column_name = "in_sorokina_2021"
  )

# crossreference for sorokina_synaptosome
merged_df_updated_sorokina_1 <- merged_df_updated_sorokina %>%
  mutate(sorokina_2021_synaptosome = sapply(Genes, is_localised, sorokina_data, "Synaptosome"))

# crossreference for sorokina_presynaptic
merged_df_updated_sorokina_2 <- merged_df_updated_sorokina_1 %>%
  mutate(sorokina_2021_presynaptic = sapply(Genes, is_localised, sorokina_data, "Presynaptic"))

# crossreference for sorokina_postsynaptic
merged_df_updated_sorokina_3 <- merged_df_updated_sorokina_2 %>%
  mutate(sorokina_2021_postsynaptic = sapply(Genes, is_localised, sorokina_data, "Postsynaptic"))


###### SynGO crossreference
# loading the data
syngo_data <- read_excel("/mnt/ip-proteomics/for_annotations/SynGO_bulk_download_release_20231201/syngo_annotations.xlsx")

# crossreference
merged_df_updated_syngo <- merged_df_updated_sorokina_3 %>%
  check_presence_in_reference(
    column_to_search = "Protein.Group",
    reference_values = syngo_data$uniprot_id,
    new_column_name = "in_syngo"
  )


###### Crossreference with Uniprot subcellular location
values_uniprot <- merged_df[["Protein.Group"]]

## Seperate multiple entries
multiple_accessions <- grep(";", values_uniprot, value = TRUE)  # Identify entries with multiple accessions
separated_accessions <- unlist(strsplit(multiple_accessions, ";")) # Separate the multiple accessions into individual entries
protein_list_cleaned <- values_uniprot[!grepl(";", values_uniprot)] # Create a cleaned protein list by excluding the original combined entries
final_protein_list <- c(protein_list_cleaned, separated_accessions) # Combine the cleaned list with the separated accessions

# Do the Uniprot request request in batches (sometime the connection gets lost)
batch_size <- 50
num_proteins <- length(final_protein_list)
results_list <- list()

# Folder to save intermediate batch results
dir.create("batch_results", showWarnings = FALSE)

# Process the protein list in batches
for (i in 1:ceiling(num_proteins / batch_size)) {
  start_index <- (i - 1) * batch_size + 1
  end_index <- min(i * batch_size, num_proteins)
  current_batch <- final_protein_list[start_index:end_index]

  batch_file <- paste0("batch_results/batch_", i, ".rds")

  # Skip if this batch has already been processed
  if (file.exists(batch_file)) {
    cat(paste("Skipping batch", i, "- already processed.\n"))
    results_list[[i]] <- readRDS(batch_file)
    next
  }

  cat(paste("Processing batch", i, "of", ceiling(num_proteins / batch_size), "...\n"))

  # Try to get subcellular locations for the current batch
  batch_results <- try(GetSubcellular_location(current_batch), silent = TRUE)

  if (inherits(batch_results, "try-error")) {
    cat(paste("Error processing batch", i, ":", batch_results, "\n"))
    # You can save the error details if needed
  } else {
    results_list[[i]] <- batch_results
    saveRDS(batch_results, batch_file)  # Save each batch result to disk
  }

  Sys.sleep(2) # Wait to avoid rate limiting
}

# Optionally combine all results into a single object after processing
subcellular_locations <- do.call(rbind, results_list)

# Define the terms you are looking for (lowercase for case-insensitive matching)
membrane_terms <- tolower(c(
  "cell membrane",
  "cell junction",
  "plasma membrane",
  "secreted",
  "extracellular space",
  "extracellular matrix",
  "cell surface"
))

# Create the new column 'is_membrane_related' initialized to FALSE
subcellular_locations$uniprot_membrane_related <- FALSE

# Get the correct column name
location_column_name <- "Subcellular.location..CC."

# Iterate through each row of the data frame
for (i in 1:nrow(subcellular_locations)) {
  location_text <- tolower(subcellular_locations[[location_column_name]][i])

  # Check if location_text is not NA and has a length greater than 0
  if (!is.na(location_text) && length(location_text) > 0) {
    # Iterate through the membrane-related terms
    for (term in membrane_terms) {
      if (grepl(term, location_text, fixed = TRUE)) {
        subcellular_locations$uniprot_membrane_related[i] <- TRUE
        break # If one term is found, no need to check others for this row
      }
    }
  }
}

subcellular_locations$Protein.Group <- rownames(subcellular_locations)

# merging the dataframe based on the Protein.Group
merged_df_updated <- merged_df_updated_syngo %>%
  left_join(subcellular_locations %>% select(Protein.Group, uniprot_membrane_related), by = "Protein.Group")


write.csv(merged_df_updated,"gpr37l1_limma_ip_proteins_annotations.csv", row.names = TRUE)



##### Clusterprofiler analysis 
BiocManager::install("org.Mm.eg.db")
library(org.Mm.eg.db)

values_uniprot <- merged_df_updated[["Protein.Group"]]
length(values_uniprot)
values_background_uniprot <- raw_data[["Protein.Group"]]
length(values_background_uniprot)

go <- enrichGO(gene = values_uniprot,
                OrgDb = org.Mm.eg.db,
                keyType = "UNIPROT",
                ont = "BP", # or CC or MF
                pvalueCutoff = 0.05,
                pAdjustMethod = "BH",
                #universe = values_background_uniprot, # specifiying background list, if not just the entire genome
                readable = TRUE)
head(summary(go))

dotplot(go, showCategory=30)


##### for gProfiler check through webportal
check <- annotations_data$Genes
print(check)
for (item in check) {
  cat(item, "\n")
}

##### make the percentage match graph

# calculate the percentages for each of the annotations
total_rows <- nrow(merged_df_updated)

# sorokina in db 
num_true_sorokina <- sum(merged_df_updated$in_sorokina_2021, na.rm = TRUE)
sorokina_percentage_matches <- (num_true_sorokina / total_rows) * 100

# sorokina synaptosome
num_true_sorokina_syn <- sum(merged_df_updated$sorokina_2021_synaptosome, na.rm = TRUE)
sorokina_percentage_synaptosome_matches <- (num_true_sorokina_syn / total_rows) * 100

# sorokina presynaptic
num_true_sorokina_pre <- sum(merged_df_updated$sorokina_2021_presynaptic, na.rm = TRUE)
sorokina_percentage_presynaptic_matches <- (num_true_sorokina_pre / total_rows) * 100

# sorokina postsynaptic
num_true_sorokina_post <- sum(merged_df_updated$sorokina_2021_postsynaptic, na.rm = TRUE)
sorokina_percentage_postsynaptic_matches <- (num_true_sorokina_post / total_rows) * 100

# Van Oostrum 2023
num_true_van_oostrunm <- sum(merged_df_updated$in_van_oostrum_2023, na.rm = TRUE)
oostrum_percentage_matches <- (num_true_van_oostrunm/total_rows) * 100

# SynGo
num_true_syngo <- sum(merged_df_updated$in_syngo, na.rm = TRUE)
SynGO_percentage_matches <- (num_true_syngo/total_rows) * 100


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
#ffce82
# bar colors
bar_colors <- ifelse(categories == "SynGO", "#ffadbc",
                   ifelse(categories == "Van Oostrum 2023", "#ffce82",
                          "#e5db99"))

ggplot(data.frame(percentages = percentages, categories = categories),
       aes(x = percentages, y = categories)) +
  geom_bar(stat = "identity", fill = bar_colors, color = "black") +
  scale_x_continuous(limits = c(0, 100),
                    breaks = c(0,20,40,60,80,100)) +
  xlab("% of co-immunoprecipitated proteins") +
  ylab(NULL) + # Remove default y-axis label
  theme_minimal() +
  theme(
    axis.line.x = element_line(linewidth = 0.75, color = "black"),
    axis.ticks.x = element_line(color = "black", linewidth = 0.75, size = 1),
    axis.text.x = element_text(size = 16),
    axis.title.x = element_text(size = 25, family = "Arial"),
    axis.text.y = element_blank(),
    panel.grid.major.y = element_blank(), # Remove horizontal grid lines
    panel.grid.minor = element_blank()    # Remove minor grid lines
  ) +
  geom_text(aes(label = categories),
            x = 7, # Adjust for label positioning
            hjust = 0,
            size = 10,
            color = "black",
            family = "Arial")
