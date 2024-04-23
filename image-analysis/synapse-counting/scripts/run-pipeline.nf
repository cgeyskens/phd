#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

// specifying the input and output directories
input_dir = "/Volumes/KINGSTON/code/phd/image-analysis/synapse-counting/test-images-VLGUT1-PSD95-A"
output_dir = "/Volumes/KINGSTON/code/phd/image-analysis/synapse-counting/output_data"
name_segments = "0 1 2 3 10 11"

process RunPearsonsAnalysis {
    input:
    path x

    script:
    """
    python ${baseDir}/pearsons.py ${input_dir} ${output_dir} ${name_segments}
    """ 
}

process RunMandersAnalysis {
    input:
    path x

    script:
    """
    python ${baseDir}/manders.py ${input_dir} ${output_dir} ${name_segments}
    """ 
}

process OverlapAnalysis {
    input:
    path x

    script:
    """
    python ${baseDir}/overlap.py ${input_dir} ${output_dir} ${name_segments}
    """ 
}

process PunctaAnalysis {
    input:
    path x

    script:
    """
    python ${baseDir}/puncta.py ${input_dir} ${output_dir} ${name_segments}
    """
}

workflow {
    RunPearsonsAnalysis(
        input_ch = input_dir
        )
    RunMandersAnalysis(
        input_ch = input_dir
        )
    OverlapAnalysis(
        input_ch = input_dir
    )
    PunctaAnalysis(
        input_ch = input_dir
    )

}