#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

// specifying the input and output directories
input_dir = "/Volumes/KINGSTON/code/phd/image-analysis/synapse-counting/test-images-VLGUT1-PSD95-A"
output_dir = "/Volumes/KINGSTON/code/phd/image-analysis/synapse-counting/output_data"
name_segments = "0 1 2 3"

process RunPearsonsAnalysis {
    input:
    path x

    script:
    """
    python ${baseDir}/pearsons.py ${input_dir} ${output_dir} ${name_segments}
    """
    output:
    file y
}

process RunMandersAnalysis {
    input:
    path x

    script:
    """
    python ${baseDir}/manders.py ${input_dir} ${output_dir} ${name_segments}
    """
    
    output:
    file z
}

process GatherResults {



}

workflow {
    RunPearsonsAnalysis(
        input_ch = input_dir
        )
    RunMandersAnalysis(
        input_ch = input_dir
    )
    GatherRes(


    )
}