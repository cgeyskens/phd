#### If your having troubles connecting the to external hard drive: do the following commands
sudo umount /mnt/d
sudo mount -t drvfs D: /mnt/d

#### You can run the nextflow pipeline with a conda env even if the conda env is not activated in your current directory. Specificy in the nextflow.config file the conda under profiles and then you execute with the following command:
nextflow run test.nf -profile conda # I now specified conda in the config file

#### This command prevents the creation of seperate .nextflow.log.x files for each run pipeline
nextflow -log /var/log/nextflow.log run test.nf -profile conda

#### After running pipeline: removing the work directory & the nextflow.log.x files
rm -r work/
rm nextflow.log.*

##### Issues
- For some reason nextflow can't find the scripts into the bin directory, temp fix with ${baseDir}/bin/

#### resolved
- passing arguments from nextflow to python (commit on 09/08/23) with the argparse library
- passing the outputfolder from nextflow to python (16/08/23) with the argparse libary (multiple parsers)
- specifying the folder in nextflow that needs to be processed in python (16/08/23)

##### Next-up
1. with basic coloc workflow
- develop a surogate process in the test pipeline: like threshold, this way you can go along make an entire test-pipeline with each step
- make this into a batch process
- channel the csv to another script to calculate the mean, and channel the mean to another script to actually calculate the coloc

2. include suggestions of Benjamin into the ipynb and py script

3. Try out small pipeline with the surogate coloc analysis

4. benchmark the different coloc analysis workflows


#### test-for-pipeline
Working scripts: test.nf wit test-file-1.py





### Benchmarking coloc analyses
1. Original: local maxima detection + distance based colocalization (based on pixel intensity) - need high computing power
2. Pixel classifcation + distance based colocalization (based on ML) - need high computing power

3. Pearsons correlation coef (based on pixel intensities)
4. Manders overlap coef (based on binary images, so thresholded)
5. Overlap area (based on binary images, so thresholded, convert pixels to um2) - also seperate pixel (um2) of vglut1 and psd95 channels

### For synaptic markers seperately (check out Sudhof lab papers)
1. mean fluorescence intensity (MFI) or staining intensity
2. puncta nr or puncta density (puncta nr per 100 um2)
3. staining area or puncta size