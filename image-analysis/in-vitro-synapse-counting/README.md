# Semi-automated image processing

## Dependencies
(1) [Fiji](https://fiji.sc/). Analysis performed with v2.16.0/1.54p; Java 21.0.7 [64-bit]. Downloaded for MacOS Apple Silicon.

(2) [Fiji hiPNAT Plugin](https://github.com/tferr/hIPNAT/tree/1.0.2). Download the .jar file under releases from v1.0.2, drag & drop it under the plugins folder in the Fiji files. Now you should see the `Summarize Skeleton` function under `Analyze/Skeleton`.

## Semi-automatic ImageJ Macro Protocol

Scripts can be found in `/image-processing-scripts/`

### 1.	Calculate mean threshold over all images and channels per experiment

a. Open an image manually.

b. Run `01_macro_z-stack_split_grey.ijm`; it will split the channels, max project, set LUT to grey and delete the original image.

c. Manually find the best threshold (`Image/Adjust/Threshold) for each channel & note it in an excel sheet. Close this image.

d. Do this for every image in the experiment & calculate the mean threshold in the excel sheet. If you have the average mean threshold for all images in an experiment (Fc & Fc-VCAM1), you can now proceed to step 2.


### 2. Clear soma ROIs for unspecific binding

a. Open an image manually

b. Run `01_macro_z-stack_split_grey.ijm`; it will split the channels, max project, set LUT to grey, and delete the original image.

c. Open the ROI Manager (`Analyze/Tools/ROI Manager`), manually draw an ROI around each soma in the image with the polygon tool, add these ROIs to the ROI Manager and save these ROIs in soma_roi folder inside each condition (expX/Fc/).

d. Run `02_macro_clear_rois.ijm` to clear all the information inside the ROIs from all channels, as the antibodies detected unspecific staining in the soma's.

### 3.	Set the mean threshold & create mask

a. In `03_macro_set_threshold.ijm`, copy paste the mean threshold from the excel sheet (from step 1) inside the macro.

b. Run `03_macro_set_threshold.ijm` to set the mean threshold in all channels and create a mask of all channels.

### 4.	Count the synapses

a. Run `04_macro_count_synapses.ijm` to count the pre & post synapses with analyze particles function. Be aware of the channels corresponding to the pre, post or MAP2 marker.

b. Copy paste the count numbers and the average particle size in the excel sheet.

### 5.	Skeletonize to calculate total branch length 

a. Run `05_macro_skeletonize.ijm` to get the total branch length, this will also close all the channels but not the produced tables with numbers.

b. Copy paste the numbers in th excel sheets. Calculate the now the presynapse, postsynapse and synapse density in the excel sheet.

### 6. Close all the images and tables

a. Run `06_macro_closing.ijm` to close all image, tables, and clear our the ROI manager.

b. Restart with Step 2 for the next image.


# Downstream data analysis

## Data

Processed data can be found inside `/data/` folder. Data:

VGLUT1-PSD95 dataset: `exp10_11_12_13_data.xlsx`

VGAT-GEPH dataset: `exp17_20_23_24_data.xlsx`

## Dependencies

A single R script for data analysis was used inside a reproducible Docker development container .devcontainer in Visual Studio Code. Please see the [documentation](https://github.com/RamiKrispin/vscode-r) for more info.

Same Docker Development environment as in `image-analysis/synapse-counting/notebooks/`

## Analysis steps:

`in-vitro-synapse-lmem.R`

1. Filtering out FOVs that were not included in the image processing due to poor image quality (2 for VGLUT1-PSD95 and 2 for VGAT-GEPH)
2. Calculate synapse density = synapse count / dendritic length
3. Normalize values using control (Fc) mean across experiments
4. Log2 transformation for better fit in linear mixed effects model
5. Fitting of a linear mixed effects model to account for inter-experimental variability
6. Post-hoc pairwise comparisons using Tukey's HSD test
7. Do this for VGLUT1-PSD95 synapse puncta density and synapse puncta size
8. Do this for VGAT-GEPH synapse puncta density and synapse puncta size

## Credits
Semi-automatic image processing and Fiji Macro's were based on Efstathia Kotoula's macro's and adjusted to better flow.