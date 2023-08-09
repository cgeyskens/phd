#### If your having troubles connecting the to external hard drive: do the following commands
sudo umount /mnt/d
sudo mount -t drvfs D: /mnt/d

#### You can run the nextflow pipeline with a conda env even if the conda env is not activated in your current directory. Specificy in the nextflow.config file the conda under profiles and then you execute with the following command:
nextflow run test.nf -profile conda

#### This command prevents the creation of seperate .nextflow.log.x files for each run pipeline
nextflow -log /var/log/nextflow.log run test.nf -profile conda

#### After running pipeline: removing the work directory & the nextflow.log.x files
rm -r work/
rm nextflow.log.*

##### Still needs debugging
- For some reason nextflow can't find the scripts into the bin directory, temp fix with ${baseDir}/bin/

##### Next-up

1. Parsing arguments from nextflow to python
- Looking into parsing arguments from nextflow script to python script (argparse library)
- Check also batch processing in python (so that you don't need to import for each image, all the libraries again)
- solution? Defining the directory in nextflow with argparser

2. Parsing arguments from python to nextflow
