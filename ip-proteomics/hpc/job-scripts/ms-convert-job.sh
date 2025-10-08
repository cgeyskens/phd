#!/bin/bash

#SBATCH --cluster=wice
#SBATCH --job-name="ms-convert"
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=10
#SBATCH --time=1:00:00
#SBATCH --account=lp_big_wice_cpu
#SBATCH --partition=dedicated_big_bigmem
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=cydric.geyskens@kuleuven.be

echo JOB STARTED

# set source and destination directories
SOURCE="/scratch/leuven/357/vsc35768/ip-proteomics/exp17-zenotof-raw-data"
DEST="/scratch/leuven/357/vsc35768/ip-proteomics/exp17-ms-zenotof-convert-output"
THREADS=10
CONTAINER="/scratch/leuven/357/vsc35768/ip-proteomics/ms-convert.sif"

# make output dir exists
mkdir -p "$DEST"

# parralell excecution of the file conversion (10 files so 10 threads)
apptainer exec --bind "$SOURCE:/raw-data" --bind "$DEST:/output" "$CONTAINER" bash -c '
    find /raw-data -name "*.wiff" | xargs -I {} -P '"$THREADS"' bash -c "
        # Set output filename and log file
        input_file={}
        output_name=\$(basename {} .wiff).mzML
        log_file=\"/output/\$(basename {} .wiff).log\"

        # logging
        echo \"Processing: \$input_file\"
        echo \"Output: \$output_name\"
        echo \"Log: \$log_file\"

        # run the actual conversion
        wine msconvert --64 --zlib --outdir /output \"\$input_file\" > \"\$log_file\" 2>&1
    "
'
echo "JOB ENDED"