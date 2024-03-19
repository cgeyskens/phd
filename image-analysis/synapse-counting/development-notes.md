### If your having troubles connecting the to external hard drive: do the following commands
# first unmount & then mount it again
sudo umount /mnt/d
sudo mount -t drvfs D: /mnt/d

### You can run the nextflow pipeline with a conda env even if the conda env is not activated in your current directory. Specificy in the nextflow.config file the conda under profiles and then you execute with the following command:
nextflow run test.nf -profile conda # I now specified conda in the config file

### This command prevents the creation of seperate .nextflow.log.x files for each run pipeline
nextflow -log /var/log/nextflow.log run test.nf -profile conda

### After running pipeline: removing the work directory & the nextflow.log.x files
rm -r work/
rm nextflow.log.*

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

### Next-up
1. Make pipeline with basic coloc workflow in nextflow
last part: small adjustments & make the code clean:
- for the plotting: try to add the mean to the datapoints with R (for this first install R in the conda env)
- customizability: go over the pipeline and scripts try to customize it. Some things that should be specified from the nextflow script (eg the working folders, can you create the folders from nextflow, this way you only need to specify the the folder with the input images). Also customize the nextflow inputs with the params. variable (also go over the nextflow carpentries again that could be useful for this pipeline). Also instead of using the names vglut1 and psd95, use presynaptic and postsynaptic (in case of inhibitory markers).
- Write a good README.md file for future use of the pipeline (from where you should run the pipeline, where the data is, ...). 

### Next-up for after HPC training:
2. include suggestions of Benjamin into the original ipynb and py script (CLAHE, tophat filter)
3. Implement the pipeline onto the HPC using apptainer & conda
4. Implement the pipeline with GPUs on HPC
5. benchmark the different coloc analysis workflows
5. do the actual analysis


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
3. staining area or puncta size
4. punta size (Bosworth et al. 2024)