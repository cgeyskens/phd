#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

// specifying the input and output directories
input_dir = "/Volumes/KINGSTON/code/phd/image-analysis/synapse-counting/test-images-VLGUT1-PSD95/"
intermediate_dir = "/Volumes/KINGSTON/code/phd/image-analysis/synapse-counting/intermediate_data/"
output_dir = "/Volumes/KINGSTON/code/phd/image-analysis/synapse-counting/output_data/"

process calculateThresholds {
    input:
    path x

    script:
    """
    python ${baseDir}/bin/test-1-thresholds.py ${input_dir} ${intermediate_dir}
    """

    output:
    file x
}
 
process calculateMeanThresholds {
    input:
    path z

    script:
    """
    python ${baseDir}/bin/test-2-mean-threshold.py ${intermediate_dir} ${intermediate_dir}
    """
    
    output:
    file z
}

process calculateOverlap {
    input:
    path q

    script:
    """
    python ${baseDir}/bin/test-3-overlap.py ${input_dir} ${intermediate_dir} ${output_dir}
    """

    output:
    file q
}

process makeFigure {
    input:
    path w

    script:
    """
    python ${baseDir}/bin/test-4-figure.py ${output_dir} ${output_dir}
    """
}


workflow {
    input_ch = input_dir
    threshold_ch = calculateThresholds(input_ch)
    meanThreshold_ch = calculateMeanThresholds(threshold_ch)
    overlap_ch = calculateOverlap(meanThreshold_ch)
    fifure_ch = makeFigure(overlap_ch)
}



