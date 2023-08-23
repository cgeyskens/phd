#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

// specifying the input and output directories
input_dir = "/mnt/d/code/phd/image-analysis/synapse-counting/test-images-VLGUT1-PSD95/"
intermediate_dir = "/mnt/d/code/phd/image-analysis/synapse-counting/intermediate_data/"
output_dir = "/mnt/d/code/phd/image-analysis/synapse-counting/output_data/"

process calculateThresholds {
    input:
    path x

    script:
    """
    python ${baseDir}/bin/test-file-2.py ${input_dir} ${intermediate_dir}
    """

    output:
    file x
}
 
process calculateMeanThresholds {
    input:
    path z

    script:
    """
    python ${baseDir}/bin/test-file-3.py ${intermediate_dir} ${output_dir}
    """
}

workflow {
    input_ch = input_dir
    threshold_ch = calculateThresholds(input_ch)
    meanThreshold_ch = calculateMeanThresholds(threshold_ch)
}



