### If your having troubles connecting the to external hard drive: do the following commands
# first unmount & then mount it again (only in WSL on windows)
sudo umount /mnt/d
sudo mount -t drvfs D: /mnt/d

### You can run the nextflow pipeline with a conda env even if the conda env is not activated in your current directory. Specificy in the nextflow.config file the conda under profiles and then you execute with the following command:
nextflow run test.nf -profile conda # I now specified conda in the config file

### This command prevents the creation of seperate .nextflow.log.x files for each run pipeline
nextflow -log /var/log/nextflow.log run test.nf -profile conda

### After running pipeline: removing the work directory & the nextflow.log.x files
rm -r work/
rm nextflow.log.*

### this command will let you also produce a flowchart
nextflow run scripts/run-pipeline.nf -with-dag nextflow_flowchart.png

### presistent issues
- For some reason nextflow can't find the scripts into the bin directory, temp fix with ${baseDir}/bin/

### resolved
- passing arguments from nextflow to python (commit on 09/08/23) with the argparse library
- passing the outputfolder from nextflow to python (16/08/23) with the argparse libary (multiple parsers)
- specifying the folder in nextflow that needs to be processed in python (16/08/23)
- develop a surogate process in the test pipeline: like threshold, this way you can go along make an entire test-pipeline with each step (22/08/2023)
- make the surrogate function into a batch process that is callable from nextflow (23/08/2023)
- make the 2nd script for the test-pipeline that calculates the mean thresholds (23/08/2023)
- in nextflow, channel the csv to the 2nd script that calculates the mean thresholds (23/08/2023)
- make the 3rd script that actually does the operation: taking the mean thresholds and applying it to the input images (23/08/2023)
- in nextflow, channel the mean thresholds to the 3rd script to actually calculate the overlap in um2 (23/08/2023)
- make a fourth script that makes also nice graphs of: for each hippocampal layer a comparison between actual and rotated overlaps in um2 with the mean thresholds of each hippocampal layer (24/08/2023)
- attach this last script to the pipeline (24/08/2023)

last part: small adjustments & make the code clean: (29/08/2023)
- use "" instead of '', be consistent
- naming of python and nextflow scripts: change the names of python to more intuitive names
- delete other unnessary development scripts (and delete the output of the ipynb files)
- Take the content of this file and make a new .md file that is named "development.md".
- update the ppt file for meeting with Joris
- update the yml file

(03/05/2024)
- include and test suggestions of Benjamin into the original ipynb and py script (CLAHE, tophat filter)
- try to run the different coloc analysis workflows parralell
- implement the pipeline with GPUs on HPC or with Dask on CPUs as scikitimage cannot run on GPUs
- modularize your code: make a /src file for modules and code and then a seperate folder for the scripts
    module for: metadata, preprocessing, calc_synaptic_coloc, calc_synaptic_metrics, visualization
- make a full basic pipeline without parameter optimization so you can get the first results
(16/05/2024)
Adjust existing pipeline
- create seperate directories when running the pipelin    
- Instead of using the names vglut1 and psd95, use presynaptic and postsynaptic (in case of inhibitory markers)
- Specified from the nextflow script: the working folders, can you create the folders from nextflow, this way you only need to specify the the folder with the input images 
- Customize the names from the file input that should be stored for analysis
- adjust/customize existing standard pipeline:
    - plot lines between data points that originates form the same brain.
    - plot also statistics above plots.
    - dynamic with "LacZ-gRNA" vs "candidate-gRNA", sometimes have VCAM1-gRNA, be carefull (line 136 in results.py). 
    - dynamic with hippocampal layers (in case of VGAT-Gephyrin and VGLUT2-PSD95 combinations)
    - add extra parameter in nextflow scripts, if whether in which channel is pre or post synapse. This is coded in python, but needs arguments
    - more defined parameters in nextflow script
    - look into the fluorescence intensity descent when you do binarization of the images: found the error: it was in the argument of the label method
    - include watershed for binarized images
    - adjust preprocessing parameters for every hippocampal layer in every synaptic combination:
        - VGLUT1-PSD95
        - VGLUT2-PSD95 (PSD95 only optimized for Cortex L4)
    - creating some sort of "preprocess profiles" whether VGLUT1-PSD95, VGLUT2-PSD95 and VGAT-Gephyrin
    - filtering out puncta smaller then... because this is background

### Next-up
- adjust preprocessing parameters for every hippocampal layer in every synaptic combination:
    V- VGAT-Gephyrin
- local maxima detection:
    V- try out optuna per image
    V- arrange the output of the optimization
    V- optuna per hippocampal layer; set the ranges for each hippocampal layer
    V- try different optimizers; other ones where having issues
    V- dask-optuna?
    V- test the output of the ranges of the hippocampal layers
    V- recheck your code and make .py files
    V- include it into the nextflow pipeline and check for the optimization efficiency(graph)
- last adjustments pipeline: 
    V- try the nextflow pipeline locally
    V- simplify the pipeline with the directories structure (no intermediate structure)
    V- make the pipeline a bit nicer, with logfile when you run it, leave out the create directories
    V- include two JSON files as input for the pipeline: one for handcrafted parameters ranges for local peaks optimization and one for the handcrafted preprocessing parameters. PREPROCESS PARAMS is all the same except for CLAHE in local peaks! 
    V- go over the scripts and clean it up more and make helper functions, check out the nextflow 
        make helper functions for:
            V- the hippocampal layer extraction from the filename
            V- extract the synaptic marker from the filename
    V- reassess the handcrafted preprocessing parameters
    V- reassess the handcrafted parameter ranges for optimization

    V- ERROR when running with docker. from conda to docker, publish in dockerhub, add profiles to config file 
        V- update nextflow to latest stable version
        V- make a python package in conda
        V- try to put the scripts into a bin
        V- Kris meeting: the problem was with the docker entrypoint inside the dockerfile. For future reference, use just a python environment inside the docker image.
    V- try the pipeline on the HPC
    - was having this issue: https://github.com/conda/conda/issues/7980
        V- just make a regular docker image with pip on local
        V- pull the image from DockHub with Apptainer on HPC  --> PROBLEM: architecture AMD vs ARM, solved; I can shell into the container
        V- copy all the necessary files to the HPC and try out the pipeline
    V- sbatch the pipeline on HPC with a 100 trials for optimization, with apptainer pulling from dockerhub
    - QC: make histogram of you p-values! To see if the distribution is alright.
    - p-value adjustment by Hochberg or Bonferroni; 
        - which p values to group for this? Grouped on hippocampal layer or metric assessed?
    V- How to show effect size? Cohen's d
    - change the heatmap of the prelim results to a dotplot (like Dan suggested)
    V- make the pipeline graph
    V- make documentation
    V- meet with Pavie, to look over the code
        - CLAHE not working for local peaks?
    V- Watchout with CLAHE and local peaks
    - Benjamin Pavie will look at the code
    V- Setup a meeting with the statistics consultancy of KU Leuven
        V- look into multivariate models (mixed models)
        V- Cohen's D is alright for effect size
        V- paired t-test is alright, no p adjustment needs to be made
        V- for the direction analysis: can take the difference between each brains mean
    V- Do Cohen's D for effect size
Downstream analyses to be performed:
    - check for outliers (box plots, Z-scores, robust PCA)
    V- differential or paired PCA: across hemisphere: VCAM1-gRNA, LacZ-gRNA (not on individual images, but on brain means)
    V- Cluster analysis: hierarchical clustering or k-means clustering/gaussian mixture models. First PCA/TSNE then k means.
    - Correlation analysis: between parameters and inter-region correlation
    - Network analysis: to assess the how features co-vary
    - Mixed effect models; can account for multiple regions or replicates
    - Random Forest for feature importance
Adjustments & problems to the pipeline:
    - New features: 
            - for individual images: 
            - for individual puncta: puncta intensity, puncta shape, 
            - colocalization: synaptic distance
    - Problem: preprocessing-induced bias in your dim reduction analysis
            - Adaptive thresholding techniques: Otsu or CLAHE
            - Before applying preprocessing: intensity normalization
            - parameter harmonization, eg for VGLUT1-PSD95 use CA1 SR as reference

## Summary of meeting with statistician Geert Verbeke
- you need to look at the distribution of the metrics in every group (of fixed effects: gRNA and hippocampal layer).
- with linear models you can differ a lot of the normal distribution, they are very robust to this.
- Don't make unnecessary hard with multivariate models (they are a separate field with very advanced statistics).

Conclusion:
- make a univariate linear mixed effects model of every response variable.
- because this is an exploratory experiment, he even suggested to NOT do correction of p values for multiple testing. Mention this in the methods of the paper that you thought about this! But I would like to do Hochberg to get more confidence in the data.
- try the generalized linear mixed effects model, but keep in mind this is not so robust as the linear mixed effects model. This could 
also be good to mention in the methods that I tried another model -> sensitivity analysis with another model.
- if the Q-Q plot looks good, its fine. If its the Q-Q plot is not good, make a transformation with logarithmic scale of the data.
- must do the analysis with the lme4 package because statsmodel won't converge. Make it a function so that it outputs the values for further use in a dataframe.
- add seed for so that it is reproducible
 
## further research into the question (Ask on training multilevels with FLAMES, Feb 21)
- Try to do SVD or LDA instead of PCA
- Really need to apply a multivariate linear mixed effects model, then do a PCA on the residuals ("batch corrected") so that you capture the variance between the samples and not the batch effects
- If not: do an LDA
- It is essentialy a batch effect between samples


For next, after all prelimary data
- do the actual analysis
- make heatmap
- check pylint for right styling of code
- as a control in supplemental figure; can do some clustering with features from the region_props_table

### working test-pipeline
Working nextflow script:
test.nf
Working python scripts:
test-1-threshold.py (threshold-&-different-analysis.ipynb)
test-2-mean-thresholds.py (test-2-mean-thresholds.ipynb)
test-3-overlap.py (test-3-overlap.ipynb)
test-4-figure.py (test-4-figure.ipynb)

### Benchmarking coloc analyses
1. Original: local maxima detection + searching best parameters space with actual vs rotated + distance based colocalization (based on pixel intensity) - need high computing power
2. Pixel classifcation + distance based colocalization (based on ML) - need high computing power

3. Pearsons correlation coef (based on pixel intensities)
4. Manders overlap coef (based on binary images, so thresholded)
5. Overlap area (based on binary images, so thresholded, convert pixels to um2) - also seperate pixel (um2) of vglut1 and psd95 channels

### For synaptic markers seperately (check out Sudhof lab papers)
1. mean fluorescence intensity (MFI) or staining intensity
2. puncta nr or puncta density (puncta nr per 100 um2)
3. staining area
4. punta size (Bosworth et al. 2024)


### When finished:
- Write a good README.md file for future use of the pipeline (from where you should run the pipeline, where the data is, ...). 
- make a class or module out of the analysis? so that the scripts can the class in the module.
- parrallel workflow via nextflow?


### Statiscian meeting with Steffen Fieuws 28.05.2025
- take the mean per brain, do not work on individual sections
- take the log and then the ratio of LacZ-gRNA vs candidate-gRNA & simplify the model 
    (take out the nested measurements (1|Brain:section:hippocampal_layer))
- PCA plot on hemispheres (gRNA) on residuals to correct for the brain batch effect


### Statistician meeting with Steffen Fieuws 09.06.2025
- take log ratio op sectie niveau: 1. ratio log op sectie niveau. 2. mean per brein.
- in lme model, residuele variantie kan verschillen tussen locaties !!!! random intercept.
- SAS code: 
proc mixed data=horizon;
class brain   location;
model logratio= location/solution;
lsmeans location/diff;
repeated/group=location type=simple;
random intercept/subject=brain;
run;
- Code in nlme package:
model_3 <- lme(
  fixed = local_peak_colocalized_spots_log2_ratio ~ hippocampal_layer,
  data = log2_ratio_section,
  random = ~1 | Brain,
  weights = varIdent(form = ~1 | hippocampal_layer)
)


### Figures eventually
Main:
- Dotplot of layers (x-axis) synaptic metrics (y-axis), showing effect size, directions & p values (0.1 FDR)
- PCA plot of samples (hemisphere level on residuals). Look into PLS-DA which is supervised.
Supplemental:
- Heatmap of features over samples
- Cluster analysis on highly variable features on sample level
- Epoch data of local peaks optimization
- Correlation analysis across metrics & between metrics and PCs


### Notes, implementing a linear mixed model on the log2 ratio values 
- most metrics can work with the model. 
    I do see some times conflicting results between the paired t-test and the simplified linear model.
        example: local_peaks_colocalized_spots for VCAM1 & VGAT_GEPH doesnt show any significance, 
                but it is very clear from the data points that there IS a diference.
- The nice thing about this setup is that its a paired design and now we just take the log2 ratio and don't account for the piared design in the model


### Next
Try the models with nlme with log_ratio (VCAM1 vs LacZ) and without. In the last model, put also the treamtment as interaction term like this: 
model_3 <- lme(
  fixed = local_peak_colocalized_spots_log2 ~ gRNA * hippocampal_layer,
  data = log2,
  random = ~1 | Brain,
  weights = varIdent(form = ~1 | hippocampal_layer)
)
- Make the preprocessing parameters for GPR37L1 better.