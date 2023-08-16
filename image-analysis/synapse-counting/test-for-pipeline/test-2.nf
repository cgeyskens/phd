#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

// specifying the input and output directories
input_dir = "/mnt/d/code/phd/image-analysis/synapse-counting/test-images-VLGUT1-PSD95/"
output_dir = "/mnt/d/code/phd/image-analysis/synapse-counting/output_data/"

process processFile {
    input:
    path x

    script:
    """
    python ${baseDir}/bin/test-file-2.py ${input_dir} ${output_dir}
    """
}

workflow {
    processFile(
        input_ch = input_dir
    )
}



