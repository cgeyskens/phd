# Semi-automatic in-vitro synapse counting using Fiji Macro's

## Dependencies
(1) [Fiji](https://fiji.sc/). Analysis performed with v2.16.0/1.54p; Java 21.0.7 [64-bit]. Downloaded for MacOS Apple Silicon.

(2) [Fiji hiPNAT Plugin](https://github.com/tferr/hIPNAT/tree/1.0.2). Download the .jar file under releases from v1.0.2, drag & drop it under the plugins folder in the Fiji files. Now you should see the `Summarize Skeleton` function under `Analyze/Skeleton`.

## Analysis protocol

### 1.	Calculate mean threshold over all images and channels per experiment

a. Open an image manually.

b. Run `macro_z-stack_split_grey.ijm`; it will split the channels, max project, set LUT to grey and delete the original image.

c. Manually find the best threshold (`Image/Adjust/Threshold) for each channel & note it in an excel sheet. Close this image.

d. Do this for every image in the experiment & calculate the mean threshold in the excel sheet. If you have the average mean threshold for all images in an experiment (Fc & Fc-VCAM1), you can now proceed to step 2.


### 2. Clear soma ROIs for unspecific binding

a. Open an image manually

b. Run `macro_z-stack_split_grey.ijm`; it will split the channels, max project, set LUT to grey, and delete the original image.

c. Open the ROI Manager (`Analyze/Tools/ROI Manager`), manually draw an ROI around each soma in the image with the polygon tool, add these ROIs to the ROI Manager and save these ROIs in soma_roi folder inside each condition (expX/Fc/).

d. Run `macro_clear_rois.ijm` to clear all the information inside the ROIs from all channels, as the antibodies detected unspecific staining in the soma's.

### 3.	Set the mean threshold & create mask

a. In `macro_set_threshold.ijm`, copy paste the mean threshold from the excel sheet (from step 1) inside the macro.

b. Run `macro_set_threshold.ijm` to set the mean threshold in all channels and create a mask of all channels.

### 4.	Count the synapses

a. Run `macro_count_synapses.ijm` to count the pre & post synapses with analyze particles function. Be aware of the channels corresponding to the pre, post or MAP2 marker.

b. Copy paste the count numbers and the average particle size in the excel sheet.

### 5.	Skeletonize to calculate total branch length 

a. Run `macro_skeletonize.ijm` to get the total branch length, this will also close all the channels but not the produced tables with numbers.

b. Copy paste the numbers in th excel sheets. Calculate the now the presynapse, postsynapse and synapse density in the excel sheet.

### 6. Close all the images and tables

a. Run `macro_closing.ijm` to close all image, tables, and clear our the ROI manager.

b. Restart with Step 2 for the next image.


## Credits
Workflow and Fiji Macro's were based on Efstathia Kotoula's macro's and adjusted to better flow.