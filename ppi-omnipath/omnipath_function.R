### With this script I would like to get the interactors for my candidates based on
### the OMNIPATH protein-protein interaction collection of many different databases. 

## Layout
## 1. Installing packages and loading libraries
## 2. Loading OMNIPATH interactors of candidate
## 3. Loading OMNIPATH subcellular annotations of candidate's interactors
## 4. Merging data from step 2 and step 3
## 5. Preprocess data for input into network visualization library (igraph)
## 6. Create the igrap object
## 7. Network visualization with ggraph & saving the output


# 1. Installing required packages and loading the libraries ----
R.version

remove.packages('cli')
install.packages('cli')

remove.packages('rlang')
install.packages("rlang")
install.packages("usethis")
install.packages("devtools")
install.packages("githubinstall")
install.packages("dnet")
install.packages('vctrs')
install.packages("gprofiler2")

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("supraHex")

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("Rgraphviz")

install.packages('igraph')

install_github('saezlab/OmnipathR')

library(usethis)
library(devtools)
require(devtools)
library(tidyr)
library(supraHex)
library(dnet)
library(vctrs)
library(gprofiler2)
library(igraph)
library(ggraph)
library(tidyverse)
library(stringr)

library(OmnipathR)


# Make sure you load the libraries
# For the function, insert the candidate in the format as omnipath_ppi("Vcam1")

# FUNCTION
omnipath_ppi <- function(candidate){


# 2. Getting all interactors of the candidate in one dataframe with annotations---- 

# 2.1 Get the interactors from the 'all interactors' database (also without references).

# 2.1.1. First we get the mouse interactors of the candidate

# Importing all the mouse interactions from the OMNIPATH database
all_interactions_mouse <- import_all_interactions(organism = 10090) 

# Now we filter on the candidate; trim the unnecessary columns; add an extra column of interactors of the candidate;
# filter columns based on the interaction column to exclude duplicate interactions and do this by excluding 
# the duplicate rows with the lowest nr of reference; change the interaction column values to uppercase and add a 
# column 'Species'-mouse; change the mouse_interaction column values to upper case
mouse_candidate_interactors <- all_interactions_mouse %>%
  filter(source_genesymbol == candidate | 
           target_genesymbol == candidate) %>%
  select(-is_directed, 
         -is_stimulation,
         -is_inhibition,
         -consensus_direction,
         -consensus_stimulation,
         -consensus_inhibition,
         -dorothea_level) %>% 
  mutate(interaction=case_when(
    source_genesymbol == candidate ~ target_genesymbol,
    source_genesymbol != candidate ~ source_genesymbol)) %>% 
  arrange(desc(n_references)) %>%
  distinct(interaction, .keep_all = T) %>% 
  setNames(paste0('mouse_', names(.))) %>%
  mutate(mouse_interaction = toupper(mouse_interaction))
                                            

# 2.1.2. Now we get the human interactors of the 'all interactors' database

# Importing all the human interactions from the OMNIPATH database
all_interactions_human <- import_all_interactions(organism = 9606)

# Same steps as for the mouse interactors in 2.1.1. but we include one more step: 
# because of the complexes (eg. ITGA_ITGD), we separate these with "separate_rows" function
human_candidate_interactors <- all_interactions_human %>%
  filter(source_genesymbol == toupper(candidate) | 
           target_genesymbol == toupper(candidate)) %>%
  select(-is_directed, 
         -is_stimulation,
         -is_inhibition,
         -consensus_direction,
         -consensus_stimulation,
         -consensus_inhibition,
         -dorothea_level) %>%
  mutate(interaction=case_when(
    source_genesymbol == toupper(candidate) ~ target_genesymbol,
    source_genesymbol != toupper(candidate) ~ source_genesymbol)) %>%
  separate_rows(11, sep="_") %>%
  arrange(desc(n_references)) %>%
  distinct(interaction, .keep_all = T) %>%
  setNames(paste0('human_', names(.)))
  
  
# 2.2. Merging the dataframes from the mouse and human candidate interactors from step 2.1.

# Merge dataframes
candidate_interactors <- full_join(
  human_candidate_interactors,
  mouse_candidate_interactors,
  by=c("human_interaction"="mouse_interaction"))

# Change the "human_interaction" column to "interactor" & select only necessary columns
candidate_interactors_trim <- candidate_interactors %>%
  mutate(interaction=human_interaction) %>%
  select(interaction,
         human_type, 
         mouse_type, 
         human_n_references,
         mouse_n_references)


# 3. Create dataframe of subcellular location annotations of the candidate interactors ----

# Import subcellular location annotations from the ComPPI database because
# this database collected a localization score for each protein
annotations_ComPPI <- import_omnipath_annotations(resources = "ComPPI")

# Filter the annotations_ComPPI dataset on the interactors of the candidate; 
# pivot the dataframe; get the subcellular location annotation with the highest score
annotations <- annotations_ComPPI %>%
  filter(genesymbol %in% candidate_interactors_trim$interaction) %>%
  pivot_wider(names_from = "label",
              values_from = "value") %>%
  arrange(desc(score)) %>%
  distinct(genesymbol, .keep_all = T)


# 4. Merging the data frames from step 2 & step 3 ----

# Joining thhe "candidate_interactors_trim" with the annotations data frame 
df_merged <- full_join(candidate_interactors_trim, annotations, 
                     by = c("interaction"="genesymbol"))

# Add another column to the dataset called "species"
df_final <- df_merged %>%
  mutate_if(is.numeric, as.character) %>%
  replace(is.na(df_merged),"/") %>%
  mutate(species=case_when(
    human_type != "/" & mouse_type != "/" ~ "human & mouse",
    human_type == "/" & mouse_type != "/" ~ "mouse",
    human_type != "/" & mouse_type == "/" ~ "human" ))

# REMARK: the subcellular location annotations are based on human orthologues,
# because there were no mouse protein annotations in the omnipath database


# 5. Prepare dataframes for input into iGraph object for network analysis----

# For downstream analysis in iGraph, we need two lists: 
# The node list and the edge list
# Attributes of node list = protein names, subcellular locations and species
# Attributes of edge list = protein-protein, nr of references_human,
# nr of references_mouse and type of interaction

# 5.1 Creating the node list

# Filtering out the interactions that don't have a subcellular annotation (eg. tf-miRNAs, miRNA-target, ...);
# Selecting the necessary columns; renaming those columns and annotate the candidate as both human and mouse;
# adding the candidate if the candidate is not there
node_list_intermediate <- df_final %>%
  filter(df_final$location != "/") %>%
  select(interaction,
         location,
         species) %>%
  rename("Protein" = "interaction",
         "Subcellular location" = "location",
         "Species" = "species") 

if(toupper(candidate) %in% node_list_intermediate$Protein) {
  
  node_list <- node_list_intermediate %>%
  mutate(Species = case_when(
    Protein == toupper(candidate) ~ "human & mouse", 
    TRUE ~ Species)) %>%
    mutate(`Subcellular location` = case_when(
      Protein == toupper(candidate) ~ "membrane",
      TRUE ~ `Subcellular location`))

  } else {
  
  node_list <- node_list_intermediate %>%
  add_row("Protein" = case_when(
    !((toupper(candidate)) %in% node_list_intermediate$Protein) ~ toupper(candidate))) %>%
    mutate(Species = case_when(
      Protein == toupper(candidate) ~ "human & mouse", 
      TRUE ~ Species)) %>%
    mutate(`Subcellular location` = case_when(
      Protein == toupper(candidate) ~ "membrane",
      TRUE ~ `Subcellular location`))
}


# 5.2. Creating the edge list ----

# 5.2.1. Join all the interactions from the human and mouse

# We rename the n_references column to n_references_species & add a species column
all_interactions_mouse_n <- all_interactions_mouse  %>%
  rename("n_references_mouse" = "n_references") %>%
  cbind(species = "mouse")
  
all_interactions_human_n <- all_interactions_human %>%
  rename("n_references_human" = "n_references") %>%
  cbind(species = "human")

# Merge these data frames
all_interactions <- bind_rows(all_interactions_human_n, all_interactions_mouse_n)



# 5.2.2. Get all the complexes (only human) that bind to the candidate because these are not included in all interactions
all_complexes_candidate <- all_interactions_human %>%
  filter(source_genesymbol == toupper(candidate) | 
           target_genesymbol == toupper(candidate))  %>%
  select(-is_directed, 
         -is_stimulation,
         -is_inhibition,
         -consensus_direction,
         -consensus_stimulation,
         -consensus_inhibition,
         -dorothea_level) %>%
  mutate(complex = case_when(
    source_genesymbol == toupper(candidate) ~ target_genesymbol,
    source_genesymbol != toupper(candidate) ~ source_genesymbol)) 

if(any(str_detect(all_complexes_candidate$complex, "_"))) {
  
  all_complexes_candidate <- all_complexes_candidate %>%
  filter(str_detect(complex, "_")) %>% # the error comes from here -> it filters based on the "_", when there is none, then the following code would cannot proceed
  mutate(complex_a = ifelse(str_detect(complex, '_'),
                            str_match(complex, '(.*)_')[,2],
                            complex),
         complex_b = str_match(complex, '_(.*)')[,2]) %>%
  arrange(desc(n_references)) %>%
  distinct(complex, .keep_all = T) %>%
  select(source, 
         target, 
         type, 
         sources, 
         references, 
         curation_effort,
         n_references, 
         n_resources, 
         complex_a, 
         complex_b) %>%
  rename("source_genesymbol" = "complex_a",
         "target_genesymbol" = "complex_b",
         "n_references_human" = "n_references")  %>%
  cbind(species = "human")
  
} else {
    all_complexes_candidate <- data.frame(source=character(),
                                          target=character(),
                                          type=character(),
                                          sources=character(),
                                          references=character(),
                                          curation_effort=double(),
                                          n_references_human=double(),
                                          n_resources=integer(),
                                          source_genesymbol=character(),
                                          target_genesymbol=character())
      
  }


# 5.2.3. Merging the complexes dataframe with the all interactors data frame ----

# Create a string with all the interactions of the candidate that we use will to filter
x <- str_to_title(as.vector(node_list$Protein))
y <- toupper(x)
z <- c(x,y)

# Filter all_interaction data frame with the string z
all_interactions_interactors <- filter(all_interactions,
                                        source_genesymbol %in% z,
                                        target_genesymbol %in% z)

# Merge the complex data frame with the filtered all_interactions data frame
interactions_complexes_joined <- full_join(all_interactions_interactors, 
                                            all_complexes_candidate)

# Changing the protein ids for mouse to upper case; make a new column "n_references" 
# which are the total nr of references; select necessary columns; order, filter columns 
# and pivot wider; rename columns
interactions_complexes_intermediate <- interactions_complexes_joined %>%
  mutate(source_genesymbol = toupper(source_genesymbol),
        target_genesymbol = toupper(target_genesymbol)) %>%
  mutate(n_references=coalesce(
    n_references_human, n_references_mouse)) %>%
  select(source_genesymbol, 
         target_genesymbol, 
         type, 
         n_references, 
         species) %>%
  arrange(source_genesymbol, target_genesymbol) %>%
  distinct(across(-n_references), .keep_all = T) %>%
  pivot_wider(names_from = species, values_from = n_references) %>%
  rename("Protein1" = "source_genesymbol",
         "Protein2" = "target_genesymbol",
         "Type of interaction" = "type",
         "Number of references for human" = "human",
         "Number of references for mouse" = "mouse")
  
interactions_complexes_final <- interactions_complexes_intermediate %>%
  replace(is.na(.), 0) %>%
  mutate(`Total nr of references` = rowSums(.[4:5]))
  
edge_list <- interactions_complexes_final



# 6. Create igraph object ----
g <- graph_from_data_frame(d=edge_list, vertices=node_list, directed = FALSE)
degree(g)
plot(g)



# 7. Make a nice network plot ----
set.seed(300) # to get the same output graph every time
g1 <- g %>%
  ggraph(layout = "kk") +
  geom_edge_link(alpha = .5, 
                 aes(width = `Total nr of references`, 
                     color=`Type of interaction`)) +
  geom_node_point(aes(color=`Subcellular location`, shape=Species), 
                  size=10) +
  geom_node_text(aes(label = name),  repel = FALSE, colour="black", size= 3)+
  scale_edge_width(range = c(0, 2),breaks = c(0, 10, 20), limits = c(0,20), labels=c("0", "10", "+20")) +
  scale_shape_manual(values=c(19,18)) +
  scale_color_brewer(palette = "Set2") +
  scale_edge_color_manual(values=c(post_translational = "#1d9cd8", transcriptional="#eb662a"))+
  theme_graph() +
  guides(edge_color = guide_legend(override.aes = list(edge_width =3)),
         edge_width = guide_legend(override.aes = list(size = 3)))
g1



# Saving the plot ----
ggsave(filename = paste0("../Results/Network/", format(Sys.time(), "%Y%m%d_%H%M%S"), "Network.png"), 
       plot = g1, width = 25, height = 20,
       path = "../Results/Network/",
       units = "cm",limitsize = FALSE)

return(g1)

}

# Gene ontology search of all interactors ----

# Exporting proteins to right format
write_excel_csv(node_list, "C:/Users/cydri/Documents/KUL-VIB  - PhD/In silico/Omnipath interactors/Interactors_VCAM1.csv")









