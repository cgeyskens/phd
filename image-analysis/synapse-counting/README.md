# Multi-metric high throughput synapse analysis

This repository documents the complete workflow for **Hippocampal lamina-resolved and multi-metric synapse analysis of CRISPR/Cas9-mediated KO of VCAM1 or GPR37L1**.  

The repo covers both **Custom multi-metric high throughput image processing pipeline** and **downstream data analyses** within a reproducible environment.


### Repository Structure
```
synapse-counting/
├─- docs/
├── hpc/ # all hpc related scripts and outputs
├── notebooks/ # all scripts for setting parameters before pipeline execution & downstream analyses
| ├── 01_in-vitro_VCAM1-gRNA_QC.R  # Analysis of the VCAM1-gRNA QC in astrocyte culture
| ├── 02_lmem_residuals_pca.R  # Residuals of synaptic metrics to account for brain-to-brain variability and PCA
| ├── 03_lmem_residuals_heatmap.ipynb  # Hierachical clustering of residuals and heatmap plotting
| ├── 04_lmem.R  # to formally test for differences in synaptic metrics between gRNA conditions
| ├── 05_heatmap_lmem.ipynb  # plotting the Log2 FC in a heatmap across synaptic metrics x hippocampal laminae
| ├── 06_plot_single_metric.R  # plotting of single synaptic metric
| ├── handcrafted_parameters.ipynb # empirically setting the preprocessing parameters per synapse type and hippocampal lamina
| ├── local-peaks-epoch.ipynb # to check the local peak maxima parameter optimization across N runs
│ └── scratch_noteboooks/ # development scripts, not used in paper
├── pipeline/
| ├── bin/ # scripts used in the pipeline
| ├── nextflow.config # config file for nextflow pipeline
| ├── optimization_param_ranges_local_peaks_202511.json # used as input for the pipeline and empirically setted in `handcrafted_parameters.ipynb`
| ├── preprocess_params_202511.json # used as input for the pipeline and and empirically setted in `handcrafted_parameters.ipynb`
| ├── run-pipeline.nf # main pipeline script
├── synapse_counting/ # Python package used in the pipeline
├── Dockerfile # used for the pipeline
└── synapse-counting-env.yml # used for the pipeline and downstream analysis
```

## Custom multi-metric high throughput image processing pipeline

You can find this in `pipeline/` with detailed documentation.

## Downstream data analysis

