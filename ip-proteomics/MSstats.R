library(MSstats)
library(ggplot2)
library(MSstatsConvert)
library(data.table)

# setting the path to the main DIANN report
input_file_path = "/mnt/data/ip-proteomics/exp17-my-diann-run/diann_output.tsv"

# constructing the annotation file
Run <- c(
    "/ip-proteomics/ms-convert-output/2ul_CG_1.mzML",
    "/ip-proteomics/ms-convert-output/2ul_CG_2.mzML",
    "/ip-proteomics/ms-convert-output/2ul_CG_3.mzML",
    "/ip-proteomics/ms-convert-output/2ul_CG_4.mzML",
    "/ip-proteomics/ms-convert-output/2ul_CG_5.mzML",
    "/ip-proteomics/ms-convert-output/2ul_CG_6.mzML",
    "/ip-proteomics/ms-convert-output/2ul_CG_7.mzML",
    "/ip-proteomics/ms-convert-output/2ul_CG_8.mzML",
    "/ip-proteomics/ms-convert-output/2ul_CG_9.mzML",
    "/ip-proteomics/ms-convert-output/2ul_CG_10.mzML"
)
Condition <- c(
    "GPR37L1-IP",
    "GPR37L1-IgG",
    "VCAM1-IP",
    "VCAM1-IP",
    "VCAM1-IP",
    "VCAM1-IP",
    "VCAM1-IgG",
    "VCAM1-IgG",
    "VCAM1-IgG",
    "VCAM1-IgG"
)

BioReplicate <- c(
    "1",
    "2",
    "3",
    "4",
    "5",
    "6",
    "7",
    "8",
    "9",
    "10"
)

df <- data.frame(
    Run = Run,
    Condition = Condition,
    BioReplicate = BioReplicate
    )
getwd()

current_dir <- getwd()
output_file <- file.path(current_dir, "annotations_diann.csv")
write.csv(df, file = output_file, row.names = FALSE)

annotation_file_path = "/workspaces/phd/ip-proteomics/annotations_diann.csv"

# importing the diann data
input = data.table::fread(input_file_path)
annot = data.table::fread(annotation_file_path)

# need to change column name because of an issue: see https://groups.google.com/g/msstats/c/5EwcYsqTSc0/m/k-F8D-iVAgAJ
setnames(input, old = "Fragment.Quant.Raw", new = "FragmentQuantCorrected")
head(input)

output = DIANNtoMSstatsFormat(input, annotation = annot, MBR = TRUE,
                                use_log_file = FALSE)

# Processing the data
summarized <- dataProcess(
    output,
    logTrans = 2,
    normalization = FALSE,
    featureSubset = "topN",
    n_top_feature = 3,
    summaryMethod = "TMP",
    equalFeatureVar = TRUE,
    censoredInt = "NA",
    MBimpute = FALSE
)

# QC plots
dataProcessPlots(data=summarized, type="ProfilePlot", 
                 address = FALSE, which.Protein = "Vcam1")
