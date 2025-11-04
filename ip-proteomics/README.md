# Synaptosome co-immunoprecipitation coupled with mass spectrometry (MS)

This repository documents the complete workflow for **co-immunoprecipitation (co-IP) coupled with mass spectrometry (MS)** analysis of synaptosomal proteins.  

The workflow covers both **raw MS data processing on the Flemish HPC (Vlaams Supercomputer Centrum, VSC)** and **downstream R-based data analysis** within a reproducible containerized environment.

## MS raw data processing

All raw data were processed on the **Vlaams Supercomputer Centrum (VSC)** using Apptainer containers to ensure reproducibility.

### 1. Raw file conversion

**Input**: `*.mzML` files from the **ZenoTOF 7600**

**Output**: `*.wiff` files

Steps:

1. Pull the ProteoWizard Apptainer container from DockerHub  
   ``` bash 
   hpc/container/pull-proteowizard-container.sh
   ``` 


2. Covert raw files: 

    ```bash
    hpc/job-scripts/ms-convert-job.sh
    ```

### 2. DIA-NN spectral library generation.

**Input**: converted `*.wiff` files

**Output**: DIA-NN spectral library (`.speclib`)

Steps:

1. Build the DIA-NN Apptainer container:
    
    ```bash
    hpc/containers/build-apptainer-diann-2.2.0.sh
    ```
    from `hpc/diann-2.2.0.def`

2. Generate the spectral library:

    ```bash
    hpc/job-script/diann-library-job.sh
    ```


### 3. DIA-NN run of samples against the spectral library

Steps:

1. For GPR37L1: 

```bash
hpc/job-scripts/diann-gpr37l1-job.sh
``` 
with log file: 

```bash
hpc/job-scripts/diann-gpr37l1-job.out
``` 

2. For VCAM1: 

```bash
hpc/job-scripts/diann-vcam1-job.sh
``` 
with log file: 

```bash
hpc/job-scripts/diann-vcam1-job.out
``` 


## MS downstream data analysis

All R-based data analysis was performed inside a **reproducible Docker development container** (.devcontainer) in Visual Studio Code. Please see the [documentation](https://github.com/RamiKrispin/vscode-r) for more info.

Analysis steps:

`scripts/01_qc_limma.R`
    
1. Remove contaminants, i.e. antibodies fragments & human proteins.

2. Retain only proteins that were detected at least 4 times in one condition

3. Data was log2 transformed and missing values imputed using MinDet

4. PCA analysis of samples

5. Differential protein enrichment using Limma

`scripts/02_annotation.R`

1. Annotations of co-IPed proteins with Van Oostrum et al. 2023, Sorokina et al. 2021, SynGO and UniProt (cell surface proteins)

2. Generata the annotation plot

`scripts/03_log2_rank.R`

1. Ranking of the co-IPed proteins according to the ratio = IP (mean) / IgG (mean)

2. In case of GPR37L1, co-IPed are filtered based on UniProt Annotations

3. Produce ranked plots.

`scripts/04_go_analysis.R`

1. Only for GPR37L1 co-IPed proteins; Biological Process annotation

2. Generate GO terms visualization

`scripts/05_interactome_overlap.R`

1. Compute overlap between VCAM1 and GPR37L1 interactomes

2. Plotting of Venn Diagram

