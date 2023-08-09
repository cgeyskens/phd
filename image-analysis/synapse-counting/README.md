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

##### Next-up
1. with basic image parameters (image size)
- specify the input directory with nextflow & batch process all the files inside the specified dir with python & output a single csv file of the input file parameters

- channel the csv to another python script and calculate the mean

2. include suggestions of Benjamin into the ipynb and py script

3. Try out small pipeline with the actual parameters
