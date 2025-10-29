###############################################################################
# Script: 02_annotations.R
# Purpose: annotate co-IPed hits with relevant databases and datasets: UniProt, 
# Van Oostrum et al. 2023, Sorokina et al. 2021, SynGO
# Author: Cydric Geyskens
# Date: 2025-10-27
###############################################################################

library(readr)
library(RSQLite)
library(stringr)
library(dplyr)
library(UniprotR)
library(readxl)
library(ggplot2)

library(httpgd)
hgd()


#### ============================== arguments =============================== ####
# only need to change these arguments for full analysis

# vcam1
ip_limma_file_path <- "results/Vcam1_limma_ip_proteins_paper.csv"
ip_protein = "Vcam1"

# gpr37l1
ip_limma_file_path <- "results/Gpr37l1_limma_ip_proteins_paper.csv"
ip_protein = "Gpr37l1" 

# loading the ip results
ip_data <- read.csv(ip_limma_file_path, row.names = 1)
View(ip_data)

# loading in the data or setting the path the sql data
van_oostrum_2023 <- read_tsv("/mnt/ip-proteomics/for_annotations/van_oostrum_2023.tsv")
sorokina_db_path <- "/mnt/ip-proteomics/for_annotations/synaptic.proteome_SR_20210408.db.sqlite"
syngo_data <- read_excel("/mnt/ip-proteomics/for_annotations/SynGO_bulk_download_release_20231201/syngo_annotations.xlsx")
uniprot_data_folder <- paste0("uniprot_results_", ip_protein, "_paper_20251027")


#### ============================== functions =============================== ####
# function to check presence in datasets
check_presence_in_reference <- function(df_to_modify, column_to_search, reference_values, new_column_name) {
  # Ensure reference_values is a character vector
  reference_values <- as.character(reference_values)

  # Apply the string detection and create the new column
  df_to_modify %>%
    mutate(!!new_column_name := sapply(.data[[column_to_search]], function(text_to_check) {
      any(str_detect(text_to_check, reference_values))
    }))
}

# function to check sorokina presynaptic/postsynaptic/synaptosome
is_localised <- function(mouse_name, sorokina_df, localisation_type) {
  localised_proteins <- sorokina_df %>%
    filter(Localisation == localisation_type) %>%
    pull(MouseName)
  mouse_name %in% localised_proteins
}


#### ========================== van Oostrum et al. 2023 ====================== ####

merged_df_updated_oostrum <- ip_data %>%
  check_presence_in_reference(
    column_to_search = "Protein.Group",
    reference_values = van_oostrum_2023$protein,
    new_column_name = "in_van_oostrum_2023"
  )


#### ========================== Sorokina et al. 2021 ====================== ####
# connect to sql data
connection <- dbConnect(
    drv = SQLite(),
    dbname = sorokina_db_path
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


#### ========================== SynGO ====================== ####

merged_df_updated_syngo <- merged_df_updated_sorokina_3 %>%
  check_presence_in_reference(
    column_to_search = "Protein.Group",
    reference_values = syngo_data$uniprot_id,
    new_column_name = "in_syngo"
  )

#### ==================== UniProt: cell surface proteins =================== ####

###### crossreference with Uniprot subcellular location
values_uniprot <- ip_data[["Protein.Group"]]

## seperate multiple entries
multiple_accessions <- grep(";", values_uniprot, value = TRUE)  
separated_accessions <- unlist(strsplit(multiple_accessions, ";")) 
protein_list_cleaned <- values_uniprot[!grepl(";", values_uniprot)] 
final_protein_list <- c(protein_list_cleaned, separated_accessions) 

# do the Uniprot request request in batches
batch_size <- 50
num_proteins <- length(final_protein_list)
results_list <- list()

# folder to save intermediate batch results
dir.create(uniprot_data_folder, showWarnings = FALSE)

# process the protein list in batches
for (i in 1:ceiling(num_proteins / batch_size)) {
  start_index <- (i - 1) * batch_size + 1
  end_index <- min(i * batch_size, num_proteins)
  current_batch <- final_protein_list[start_index:end_index]

  batch_file <- paste0(paste0(uniprot_data_folder, "/batch_"), i, ".rds")

  # skip if this batch has already been processed
  if (file.exists(batch_file)) {
    cat(paste("Skipping batch", i, "- already processed.\n"))
    results_list[[i]] <- readRDS(batch_file)
    next
  }

  cat(paste("Processing batch", i, "of", ceiling(num_proteins / batch_size), "...\n"))

  # try to get subcellular locations for the current batch
  batch_results <- try(GetSubcellular_location(current_batch), silent = TRUE)

  if (inherits(batch_results, "try-error")) {
    cat(paste("Error processing batch", i, ":", batch_results, "\n"))
    # You can save the error details if needed
  } else {
    results_list[[i]] <- batch_results
    saveRDS(batch_results, batch_file)  # Save each batch result to disk
  }

  Sys.sleep(2) 
}

# combine all results into a single object after processing
subcellular_locations <- do.call(rbind, results_list)

# membrane related terms
membrane_terms <- tolower(c(
  "cell membrane",
  "secreted",
  "extracellular space",
  "extracellular matrix",
  "cell surface"
))

exclude_terms <- tolower(c(
  "peripheral membrane protein",
  "mitochondrion membrane"
))

# create the new column 'is_membrane_related' initialized to FALSE
subcellular_locations$uniprot_membrane_related <- FALSE

# get the correct column name
location_column_name <- "Subcellular.location..CC."

# iterate through each row of the data frame
for (i in 1:nrow(subcellular_locations)) {
  location_text <- tolower(subcellular_locations[[location_column_name]][i])

  # check if location_text is not NA and has a length greater than 0
  if (!is.na(location_text) && length(location_text) > 0) {

    # first, check if any of the exclusion terms are present
    exclude_match <- any(sapply(exclude_terms, grepl, location_text, fixed = TRUE))

    if (!exclude_match) {
      # only look for membrane-related terms if no exclusion match
      for (term in membrane_terms) {
        if (grepl(term, location_text, fixed = TRUE)) {
          subcellular_locations$uniprot_membrane_related[i] <- TRUE
          break
        }
      }
    } else {
      # explicitly set to FALSE if exclusion term is found
      subcellular_locations$uniprot_membrane_related[i] <- FALSE
    }
  }
}

subcellular_locations$Protein.Group <- rownames(subcellular_locations)

# merging the dataframe based on the Protein.Group
merged_df_updated <- merged_df_updated_syngo %>%
  left_join(subcellular_locations %>% select(Protein.Group, uniprot_membrane_related), by = "Protein.Group")


#### ===================== saving the annotation data ================== ####
write.csv(merged_df_updated, paste0(ip_protein, "_all_limma_ip_proteins_annotations_paper.csv"), row.names = TRUE)


#### ========================== plotting ====================== ####

# setting custom colors for VCAM1 & GPR37L1
if (ip_protein == "Vcam1") {
  ip_color <- "#21a0e2"   
} else if (ip_protein == "Gpr37l1") {
  ip_color <- "#e28d21"   
} else {
  warning("You didn’t set arguments correctly — ip_protein must be either 'vcam1' or 'gpr37l1'")
  ip_color <- "#000000"   
}

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

p <- ggplot(
        data.frame(percentages = percentages, categories = categories),
        aes(x = percentages, y = categories)
        ) +
        geom_bar(stat = "identity", fill = ip_color, color = NA) +
        xlab("% of co-immunoprecipitated proteins") +
        ylab(NULL) +  # No y-axis title
        theme_minimal() +
        theme(
            axis.line = element_line(size = 1.2, color = "black"),
            axis.ticks.x = element_line(size = 1.2, color = "black"),  
            axis.ticks.length = unit(0.2, "cm"),       
            axis.ticks.y = element_blank(),    
            axis.text.y  = element_blank(),    
            axis.text.x  = element_text(size = 16, color = "black"),
            axis.title.x = element_text(size = 25, family = "Arial"),
            panel.grid.major.y = element_blank(),  
            panel.grid.minor = element_blank()     
        ) +
        # geom_text(
        #     aes(label = categories),
        #     x = 7,          
        #     hjust = 0,
        #     size = 10,
        #     color = "black",
        #     family = "Arial"
        # ) +
        scale_x_continuous(
            expand = c(0.03, 0),
            limits = c(0, 100),
            breaks = c(0, 20, 40, 60, 80, 100)
        ) 
p

ggsave(paste0(ip_protein, "_bar_plot_annotations_paper.svg"), 
    plot = p, 
    device = cairo_pdf,
    width = 23, height = 20, units = "cm", dpi=300)






