# installation of MSDAP, took a while
install.packages(c("devtools", "tidyverse", "tinytex", "BiocManager"))
tinytex::install_tinytex()

BiocManager::install(c('ProtGenerics', 'MSnbase', 'limma'), update=T, ask=F)

Sys.setenv(R_REMOTES_NO_ERRORS_FROM_WARNINGS="true")

library(ProtGenerics)
library(limma)
packageVersion("limma")
library(MSnbase)
packageVersion("MSnbase")
library(devtools)
install.packages("pdftools")
library(pdftools)
install.packages("archive")
library(archive)
install.packages("arrow", INSTALL_opts = '--no-lock')
library(arrow)
devtools::install_github("ftwkoopmans/msdap", upgrade = "never")
library(msdap)

pandoc_version()
install.packages("pandoc", version = "1.12.3")
library(pandoc)
pandoc_install("3.1.6")
pandoc_version()
sessionInfo()

# Loading the datasets
dataset = import_dataset_diann(filename = "/mnt/ip-proteomics/exp17-my-diann-run/diann_output.tsv")
dataset = import_fasta(dataset, files = c("/mnt/ip-proteomics/UP000000589_10090_15042025.fasta", "/mnt/ip-proteomics/crap_fasta_04032025.fasta"))  

# sample metadata, adjust manually the metadata file
write_template_for_sample_metadata(dataset, "sample_metadata.xlsx")

# load again the metadata
dataset = import_sample_metadata(dataset, filename = "sample_metadata.xlsx")

# set contrast
dataset = setup_contrasts(dataset, contrast_list = list(  c("IP","IgG")  ) )

# quick analysis
dataset = analysis_quickstart(
  dataset,
  #filter_min_detect = 3,            # each peptide must have a good confidence score in at least N samples per group
  #filter_min_quant = 3,             # similarly, the number of reps where the peptide must have a quantitative value
  #filter_fraction_detect = 0.75,    # each peptide must have a good confidence score in at least 75% of samples per group
  filter_fraction_quant = 0.1,     # analogous for quantitative values
  filter_min_peptide_per_prot = 1,  # minimum number of peptides per protein (after applying above filters) required for DEA. Set this to 2 to increase robustness, but note that'll discard approximately 25% of proteins in typical datasets (i.e. that many proteins are only quantified by 1 peptide)
  filter_by_contrast = TRUE,        # only relevant if dataset has 3+ groups. For DEA at each contrast, filters and normalization are applied on the subset of relevant samples within the contrast for efficiency, see further MS-DAP manuscript. Set to FALSE to disable and use traditional "global filtering" (filters are applied to all sample groups, same data table used in all statistics)
  norm_algorithm = c("vsn", "modebetween_protein"), # normalization; first vsn, then modebetween on protein-level (applied sequentially so the MS-DAP modebetween algorithm corrects scaling/balance between-sample-groups). "mwmb" is a good alternative for "vsn"
  dea_algorithm = c("deqms", "msempire", "msqrob"), # statistics; apply multiple methods in parallel/independently
  dea_qvalue_threshold = 0.01,                      # threshold for significance of adjusted p-values in figures and output tables
  dea_log2foldchange_threshold = NA,                # threshold for significance of log2 foldchanges. 0 = disable, NA = automatically infer through bootstrapping
  #diffdetect_min_peptides_observed = 2,             # 'differential detection' only for proteins with at least 2 peptides. The differential detection metric is a niche usecase and mostly serves to identify proteins identified near-exclusively in 1 sample group and not the other
  #diffdetect_min_samples_observed = 3,              # 'differential detection' only for proteins observed in at least 3 samples per group
  #diffdetect_min_fraction_observed = 0.5,           # 'differential detection' only for proteins observed in 50% of samples per group
  output_qc_report = TRUE,                          # optionally, set to FALSE to skip the QC report (not recommended for first-time use)
  output_abundance_tables = TRUE,                   # optionally, set to FALSE to skip the peptide- and protein-abundance table output files
  output_dir = "msdap_results",                     # output directory, here set to "msdap_results" within your working directory. Alternatively provide a full path, eg; output_dir="C:/path/to/myproject",
  output_within_timestamped_subdirectory = TRUE
)
