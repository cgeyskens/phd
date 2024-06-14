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
    python ${baseDir}/pearsons.py \
    --input_dir ${params.input_dir} \
    --output_dir ${params.intermediate_dir} \
    --name_segments ${params.name_segments} \
    --presynapse_channel ${params.presynapse_channel} \
    --postsynapse_channel ${params.postsynapse_channel}
    """ 

    output:
    stdout
}

process RunMandersAnalysis {
    input:
    path input

    script:
    """
    python ${baseDir}/manders.py \
    --input_dir ${params.input_dir} \
    --output_dir ${params.intermediate_dir} \
    --name_segments ${params.name_segments} \
    --presynapse_channel ${params.presynapse_channel} \
    --postsynapse_channel ${params.postsynapse_channel}
    """ 

    output:
    stdout
}

process OverlapAnalysis {
    input:
    path input

    script:
    """
    python ${baseDir}/overlap.py \
    --input_dir ${params.input_dir} \
    --output_dir ${params.intermediate_dir} \
    --name_segments ${params.name_segments} \
    --presynapse_channel ${params.presynapse_channel} \
    --postsynapse_channel ${params.postsynapse_channel}
    """ 

    output:
    stdout
}

process PunctaAnalysis {
    input:
    path input

    script:
    """
    python ${baseDir}/puncta.py \
    --input_dir ${params.input_dir} \
    --output_dir ${params.intermediate_dir} \
    --name_segments ${params.name_segments} \
    --presynapse_channel ${params.presynapse_channel} \
    --postsynapse_channel ${params.postsynapse_channel} \
    --protein_and_synaptic_marker ${params.protein_and_synaptic_marker}
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
    python ${baseDir}/results.py \
    --input_dir ${params.intermediate_dir} \
    --output_dir ${params.output_dir} \
    --protein_and_synaptic_marker ${params.protein_and_synaptic_marker}
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
