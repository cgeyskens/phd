# Proteome-informed targeted image-based single-cell spatial transcriptomics (scST)

This repository documents the complete workflow for **target image-based single-cell spatial transcriptomics (scST)** analysis of hippocampal WT P28 sections.

All raw data processing and downstream analyses were performed on the **Flemish High-Performance Computing infrastructure**
(Vlaams Supercomputer Centrum, VSC) using Conda-managed Python environments.

## Dependencies

All processing and analyses were performed within a OS-independent Conda environment, defined in: 

- `st-env-cpu.yml` (CPU-only)
- `st-env-gpu.yml` (GPU-enabled; required for Cellpose-SAM & scVI)

These environment files ensure full reproducibility across systems.

## scST raw data processing

Notebook: 

`hpc-scripts/01-raw-data-processing.ipynb`

Processing steps:

1. Transcript coordinates tables and raw DAPI images were imported into a SpatialData Object

2. Image processing: min-max normalization & CLAHE

3. Cell segmentation using Cellpose-SAM

4. Assignment of transcript to segmented cells

## Downstream data analysis

Notebook:

`hpc-scripts/02-downstream-analysis.ipynb`

Analysis steps:

1. Cell & gene-level quality filtering

2. Batch correction using scVI

3. Nearest neighbor graph construction, UMAP dim reduction & Clustering

4. Manual cell type annotation based on spatial localization patterns & cluster-specific marker genes

## Quality control:

Notenook:

`hpc-scripts/03-qc-plots.ipynb`

- Generation of all QC visualizations

## Differential expression analysis

Notebook:

`hpc-scripts/04-diff-expr-analysis.ipynb`

Analysis details:

1. Differential expression analysis using scVI

2. Comparisons performed: 
        - Astro vs CA1
        - Astro vs CA2/CA3
        - Astro vs DG Granule Neurons

3. Visualization of candidate's LogFC, dotplots & spatial expression maps 



