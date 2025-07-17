### Semi-automatic in-vitro synapse counting using Fiji Macro's

#### Analysis protocol

1.	Calculate mean threshold over all images and channels per experiment

        a.	Open image manually
        b.	Run “macro_z-stack_split_grey.ijm” for max projection, split the channels and set LUT to grey
        c.	Manually find the best threshold for each channel & write down in excel sheet
        d.	Do this for every image in the experiment & calculate the mean threshold
        e.	Close all the images

2.	Clear soma ROIs for unspecific binding

        a.	Open image manually
        b.	Run “macro_z-stack_split_grey.ijm” to max project, split the channels, set LUT to grey and delete the original image
        c.	Manually draw an ROI around each some in the image with the polygon tool, save these ROIs.
        d.	Run “macro_clear_rois.ijm” to clear all the information inside the ROIs from all channels

3.	Set the mean threshold & create mask

        a.	In “macro_set_threshold”, copy paste the mean threshold from the excel sheet (from step 1) into the macro
        b.	Run “macro_set_threshold” to set the mean threshold in all images and create a mask of all images

4.	Count the synapses

        a.	Run “macro_count_synapses” to count the pre & post synapses with analyze particles function 
        b.	Copy paste the count numbers and the average particle size in the excel sheet

5.	Skeletonize to calculate branch length 

        a.	Run “macro_skeletonize” to get the total branch length, this will also close all the images but not the produced tables with numbers
        b.	Copy paste the numbers in th excel sheets
        c. Close all the tables and start again at Step 2.


##### Credits
Fiji Macro's were based on Efstathia Kotoula's macro's and adjusted.