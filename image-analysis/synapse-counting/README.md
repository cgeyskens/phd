# Automated synapse image processing and analysis pipeline

[![PyPI](https://img.shields.io/pypi/v/synapse-counting.svg?color=green)](https://pypi.org/project/synapse-counting/)
![Development Status](https://img.shields.io/badge/development%20status-alpha-red)
[![License](https://img.shields.io/pypi/l/synapse-counting.svg?color=green)](https://github.com/cgeyskens/phd/image-analysis/synapse-counting/synapse_counnting/blob/main/LICENSE)

**synapse-counting** is a Nextflow pipeline for processing and analyzing immunostained synapses from high-zoom images. It takes as input ZEISS Airyscan images (.czi format), preprocessing parameters, parameter ranges for local peak detection, and other parameters necessary for analyzing synapses. When these parameters are set the pipeline will process the imaging data. A custom python package was also developed based on existing scikit-image functions but tailored towards synapse analysis.

<p align="center">
    <img src="image-readme.png">
</p>  

## Pipeline summary
The pipeline will analyze 8 synapse metrics with 5 processes, 4 metrics assess colocalization and 4 metrics assess the pre- and postsynaptic puncta, these metrics include the following:
1. Regular Mander's Coefficient colocalization
2. Pearsons correlation Coefficient colocalization
3. Regular pixel overlap after binarization of the images
4. Local peak maxima colocalization: 
    * For each batch of presynapse and postsynapse images per hippocampal layer, it will try to optimize the  parameters such that the scaled difference between actual and rotated presynapse/postsynapses images is maximized. 
    * It uses Tree-structured Parzen Estimator (TPE) from the optuna library to optimize the parameters.
5. Puncta analysis: pre- and postsynapse mean fluorescence intensity (MFI)
6. Puncta analysis: pre- and postsynaptic puncta density per 100 um2
7. Puncta analysis: pre- and postsynaptic staining area
8. Puncta analysis: pre- and postsynaptic puncta size

Python scripts were parallelized with Nextflow, in each script the image data was parallelized using Dask.

## Usage
In our case the image was gathered as below: 

<p align="center">
    <img src="image-readme-2.png">
</p>  

Running the pipeline example: 
```bash
nextflow run run-pipeline.nf -profile docker
```

## Pipeline input
### A. Image data naming convention
The image files were named as follows: \
`CRISPR-Exp4_IHC-Exp2_Brain-4_section-1_488-VGLUT1_647-PSD95_LacZ-gRNA_63X&3XzoomAiryscan_CA1_SLM.czi` \
for the control hemisphere 

`CRISPR-Exp4_IHC-Exp2_Brain-4_section-1_488-VGLUT1_647-PSD95_GPR37L1-gRNA_63X&3XzoomAiryscan_CA1_SLM.czi` \
for the gRNA hemisphere 

Important segments of the filename, that are needed for the analysis are separated by "_", \
we need from this \
[0] = CRISPR-Exp4, \
[1] = IHC-Exp2, \
[2] = Brain-4, \
[3] = section-1, \
[6] = LacZ-gRNA, \
[8] = CA1, \
[9] = SLM\
you get: CRISPR-Exp4_IHC-Exp2_Brain-4_section-1_CA1_SLM \
in the pipeline we will use certain segments of this filename name, that are hard coded which are the following: \
gRNA column: [6] \
hippocampal layer column: [8] and [9] (last two) \
Brain: [0]

### B. Preprocessing parameters
The processing parameters include: 
1. background substraction with rolling ball
2. gaussian blur
3. CLAHE
4. tophat
5. watershed for binary images
6. minimum puncta size threshold
7. threshold algorithm for binary images

You can set preprocessing parameters inside a JSON file that the pipeline uses. You can customize these parameters before you run the pipeline with a notebook: `notebooks/handcrafted_parameters.ipynb` \
Then you can copy these parameters in: `scripts/preprocess_params.json` 

### C. Parameter ranges for local peak maxima colocalization assessment
For the local peak maxima colocalization, we need to set some ranges for certain parameters to detect local peak maxima's.
1. pre_distance and post_distance: distance between individual (pre- and post) synaptic local peaks
2. pre_threshold & post_threshold: threshold at which a local peak is detected
3. max_distance_um: maximum distance at which a pre and postsynaptic local peak is considered a synapse.

You can set parameter ranges inside a JSON file that the pipeline uses. You can customize these parameters before you run the pipeline with a notebook: `notebooks/local_peaks_maxima.ipynb` \
Then you can copy these parameters in: `scripts/optimization_param_ranges_local_peaks.json` 

### D. Other input parameters
Inside the `nextflow.config` file you can the following pipeline parameters:
1. Whether you would like to run the pipeline inside a conda environment, a docker container or an apptainer container.
2. your path to the input directory
3. your path to the output directory
4. your path to the preprocessing parameters json file
5. your path to the ocal peak optimization parameters json file
6. which name segments to use from the image file name
7. your pre- and postsynapse image channel
8. number of trials for finding the optimal parameters for local peak maxima detection with 

## Pipeline output
The pipeline will produce following files:
1. `all_trial_data.csv` contains exhaustive trial data from local peak maxima optimization
2. `per_trial_data.csv` contains limited data per trial
3. `optimized_param_ranges.csv` contains the optimized parameters for detecting local peaks and assessing synaptic colocalization
4. `metric_results.csv` contains exhaustive results from all the metrics
5. `plotted_data.png` an image that lets you inspect the plotted results
6. `all_statistics.csv` containes the statistics
7. `internal_coloc_controls.png` an image that plots the difference between actual and rotated presynapse/postsynapses colocalization images

## Dependencies
The pipeline can be run inside a Conda environment or a Docker image. Apptainer can also pull the docker image from DockerHub to run it in a HPC system.

## Limitations
The pipeline can only be used for ZEISS image format (.czi).
