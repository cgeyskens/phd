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

process RunLocalPeaksOptimization {
    input:
    path input_images

    script:
    """
    python ${baseDir}/local_peaks_optimize.py \
    --input_dir ${params.input_dir} \
    --output_dir ${params.intermediate_dir} \
    --presynapse_channel ${params.presynapse_channel} \
    --postsynapse_channel ${params.postsynapse_channel} \
    --protein_and_synaptic_marker ${params.protein_and_synaptic_marker} \
    --nr_of_optimization_trials ${params.nr_of_optimization_trials}
    """ 

    output:
    path "local_peak_best_optimization_params.csv", emit: local_peaks_optimized
    stdout
}

process RunLocalPeaksColocalization {
    input:
    path input_images
    path local_peaks_csv

    script:
    """
    python ${baseDir}/local_peaks_coloc.py \
    --input_dir ${params.input_dir} \
    --output_dir ${params.intermediate_dir} \
    --name_segments ${params.name_segments} \
    --presynapse_channel ${params.presynapse_channel} \
    --postsynapse_channel ${params.postsynapse_channel} \
    --protein_and_synaptic_marker ${params.protein_and_synaptic_marker} \
    --optimized_params_df $local_peaks_csv
    """ 
    
    output:
    stdout
}

process RunPearsonsAnalysis {
    input:
    path input_images

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
    path input_images

    script:
    """
    python ${baseDir}/manders.py \
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

process RunOverlapAnalysis {
    input:
    path input_images

    script:
    """
    python ${baseDir}/overlap.py \
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

process RunPunctaAnalysis {
    input:
    path input_images

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
    file from_local_peaks_calc

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
    RunLocalPeaksOptimization(input_ch)
    local_peaks_ch = RunLocalPeaksColocalization(input_ch, RunLocalPeaksOptimization.out.local_peaks_optimized)
    pearson_ch = RunPearsonsAnalysis(input_ch)
    manders_ch = RunMandersAnalysis(input_ch)
    overlap_ch = RunOverlapAnalysis(input_ch)
    puncta_ch = RunPunctaAnalysis(input_ch)
    ResultsPlots(pearson_ch, manders_ch, overlap_ch, puncta_ch, local_peaks_ch)
}
