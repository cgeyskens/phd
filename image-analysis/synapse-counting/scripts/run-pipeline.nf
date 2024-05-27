#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

// the pipeline
process CreateDirectories {
    script:
    """
    mkdir -p ${params.intermediate_dir}
    mkdir -p ${params.output_dir}
    """ 

    output:
    stdout
}

process RunPearsonsAnalysis {
    input:
    path input

    script:
    """
    python ${baseDir}/pearsons.py ${params.input_dir} ${params.intermediate_dir} ${params.name_segments}
    """ 

    output:
    stdout
}

process RunMandersAnalysis {
    input:
    path input

    script:
    """
    python ${baseDir}/manders.py ${params.input_dir} ${params.intermediate_dir} ${params.name_segments}
    """ 

    output:
    stdout
}

process OverlapAnalysis {
    input:
    path input

    script:
    """
    python ${baseDir}/overlap.py ${params.input_dir} ${params.intermediate_dir} ${params.name_segments}
    """ 

    output:
    stdout
}

process PunctaAnalysis {
    input:
    path input

    script:
    """
    python ${baseDir}/puncta.py ${params.input_dir} ${params.intermediate_dir} ${params.name_segments}
    """

    output:
    stdout
}


process ResultsPlots {
    input:
    file from_pearson
    file from_manders
    file from_overlap
    file from_puncta

    script:
    """
    python ${baseDir}/results.py ${params.intermediate_dir} ${params.output_dir}
    """

    output:
    stdout
    }


workflow {
    input_ch = Channel.fromPath(params.input_dir)
    CreateDirectories()
    pearson_ch = RunPearsonsAnalysis(input_ch)
    manders_ch = RunMandersAnalysis(input_ch)
    overlap_ch = OverlapAnalysis(input_ch)
    puncta_ch = PunctaAnalysis(input_ch)
    ResultsPlots(pearson_ch, manders_ch, overlap_ch, puncta_ch)
}
