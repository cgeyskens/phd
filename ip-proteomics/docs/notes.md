## Notes for Co-Ip MS-based proteomics

##### Dev docker container template fro VSCode from:
https://github.com/RamiKrispin/vscode-r

### Input from Pedro

##### Instrument specific parameters:

ZenoTOF range: 350 to 1250 m/z range

##### DIA-NN:

Pedro used version 19.

Mass accuracy: 20 and MS1 accuracy: 12

Always use latest version

MBR: Match-between-runs. If a peptide is abundant in at least one sample, then it will be extracted with high sensitivity across the entire experiment. If you don't set this, it will take the first sample as "bleuprint". The more variable your samples, the more advantagous this MBR is. 

1 missed cleavage, try 2, this will be complexer spectra but could potentially increase the amount of protein identified.

Always change the output folder because it will overwrite it! 

Between windowns-based and linux-based DIA-NN runs, the Aerts lab saw difference in intensities

##### Downstream analysis packages:

Peptide-level analysis: MS Stats, MS-DAP

Protein level analysis: Alphapepstats: doesnt use limma but t-tests & doesnt adjust the p values. But you can make it work with additional code.

Papers: Peng et al. 2024 Nature comms

### Good YT videos:

https://www.youtube.com/watch?v=cGKJzx9IIi4

### Background notes for downstream analysis:

 Peng et al. 2024 Nature Comms.
 - Expression matrix type: 
    - DIA-NN provides directLFQ intensities, you can extract MaxLFQ, top0 and top3 with the 'iq' package from the parquet file.
 - Normalization: 
    - None: as DIA-NN does some normalization
    - Quantile: forces all samples to have the same distribution. Ranks the proteins and replaces them with the average rank value. Not used for DEs studies.
    - Center.mean: substracts the mean intensity from each proteins intensity. This is preferable for IP vs IgG samples.
    - Max: scales all values between 0 and 1. Sensitive to outliers ofcourse.
 - Imputation algo's:
    - Impseq (low mean poisson sequence imputation): assumes missingness is due to low abundance.
    - mindet (mininmal detectable value imputation): replaces values with the lowest detectable intensity in the dataset
    - min (fixed minimal value imputation): replaces values with a fixed small number, often a fraction.
    - nbavg (neighbour-based average imputation): estimates average intensities of similar proteins
 - DEA methods:
    - ROTS (reproducibility-optimized test statistic): permutation-based statistical test. It ranks features and generates a null distribution by random shuffling to estimate significance
    - Limma (linear models for micro-arrays): sensitive to missingness. Applies Bayes shrinkage to improve variance estimation.
    - DEP (Differential enrichment analysis of proteomics data): basically limma with imputation. 
    - SAM (significance analysis of microarrays): a permutation-based statistical test that calculates a modified t-statistics. Uses FDR to control for multiple testing. Introduces a delta parameter to determine significance thresholds based on dataset variability.
    - proDA (probabilistic Drop analysis): works well with missing values, no need for seperate imputation
    - ttest: too simplistic

 From DIA-NN issues:
 https://github.com/vdemichev/DiaNN/issues/1136
 - For IP; it is recommended to not use imputation for missing values (some are not detected in the IgG control condition)
 - For IP; dissabeling normalization & use MBR & no imputation. But if imputation necessary in downstream tools use minimal-value imputation on the protein level. 

From the previous sections and reading:
  - Expression matrix type: directLFQ
  - Normalization: none or else center.mean
  - imputation: if proteins are partially missing in IgG -> use mindet. If proteins are completely absent -> do not impute.
  - DEA: try out ROTS, SAM and ProDA, these will handle missing values better


Nice python-based packages: Gopher
Other gene ontology analysis: clusterProfiler, gProfiler2


### Development notes

V - copy raw data to HPC
V - make apptainer container for DIA-NN (see notes in README: https://github.com/vdemichev/DiaNN/issues 1202#issuecomment-2511182874)

V- try to run DIA-NN with same settings as Pedro on the cluster
  - issue DIA-NN doesnt support .wiff files on linux version. Need to convert from .wiff to .mzML with most popular tool: msconvert.
        V--> solution_1: make apptainer container for msconvert from official dockerhub container (https://hub.docker.com/r/proteowizard/pwiz-skyline-i-agree-to-the-vendor-licenses)
        --> solution_2: make apptainer container from custom dockerhub container (https://github.com/jspaezp/elfragmentador-data#setting-up-msconvert-on-singularity-)
        --> solution_3: perform the file conversion locally with a docker container
        --> solution_4: run DIA-NN with Wine on hpc.
###### 28.02.2025. 
Could run DIA-NN on hpc, one difference with pedro's run: I don't have a contaminants file. Ask from pedro. Also include antibodies in contaminants?
###### 04.03.2025
Could run DIA-NN on hpc, now with the contanimants file from the CRAP website. Check the intensities compared to pedro.
###### 06.03.2025
Check other packages, like MS-DAP and MSstats. https://github.com/wfondrie/msstats-demo/blob/main/msstats-demo.ipynb
###### 20.03.2025
Compared raw and log2 intensities from my DIANN run on cluster vs Pedro's run on windows and found that they were indeed different. Also the DEGs with simple ttest.
###### 01.04.2025
Setup docker environment for R locally with devcontainers.
#### TODO:
- Setup docker environment for R locally with devcontainers.
   - cgeyskens/ip-proteomics:v1 created but could not get the right extensions in the container.
- Try the analysis with DEP (uses limma), SAM (samr), ROTS (rots) for DEA.


- Try out MSstats/MS-DAP in jupyter notebooks. With MSstats there is more support: https://github.com/Vitek-Lab/MSstats/issues/34


### Downstream analyses
#### Protein level:
DEP: analysis done
ROTS: analysis done
SAMR: issue with analysis, package outdated
Limma: analysis done
proDA: analysis done
#### Peptide level:
MSstats: Get stuck when processing data (featureSubset = "all"). If with top_N_features = 300 then I get weird PCA plot, not like the protein-level analysis.
MS-DAP: nice output of pdf report but just deletes rows with too many NAs


V - Also analyze previous synaptogliosomes (Exp10) experiment! Very different log2 intensity values
- Make a script that output the QC from each and test different parameters, to increase the amount of DEPs and have more clustered samples in the PCA.
- make the .config file with the parameters
- Finetune the running parameters (try different parameter)
- Use recent version of DIA-NN


QC plots:
 V- Nr of proteins per sample (barplot)
 - Sample Coefficient of Variation plots
 V- PCA plot of samples
 (- Abundance rank of proteins IP vs IgG)
 V- Intensity plots of raw log2 data per sample (boxplot)
 V- Intensity plots of normalized data (boxplot)
 V- Line plots of interesting candidates or heatmap

 V- Differential expression analysis, Volcano plot, of VCAM1 and membrane proteins (uniprot)
 (- Differential expression analysis, Volcano plot, of SYNGO mentioned proteins)
 V- Differential expression analysis, bar plot, of crossref with synaptic proteome databases, (of Sorokina, Mclean, Croning et al. 2021 synapse proteome database ) 
 V- Differential expression analysis, Volcano plot, of crossref with van Oostrum et al. 2023
 V- Intensity line plots of raw data per sample of selected candidates
 - (SYNGO cellular compartment enrichment analysis plots of synaptic proteins in VCAM1 IP condition)
 - Plot of retention time histogram

 - Use Kaulich et al. 2025 for annotation of hipppocampal synaptic laminae (full tissue) or regions (synaptosomes)
 - Overlap plot of VCAM1 & GPR37L1 IP conditions that are membrane proteins (also include non-mitocochondrial CSPs, see Chan et al. 2025)
 - Gene ontology analysis plots of GPR37L1 IP proteins, use as background the whole mouse proteome
 - Heatmap of GPR37L1 IP membrane proteins that are interesting (extension of the Line plots)
 - Crossreference VCAM1 and GPR37L1 IP datasets with the in silico PPI datasets.