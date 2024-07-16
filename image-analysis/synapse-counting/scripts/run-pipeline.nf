#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

// helpers for visualization
def input = params.protein_and_synaptic_marker
def parts = input.split('_')
def protein = parts[0] + "_" + parts[1]
def synapse_marker = parts[2] + "_" + parts[3]

// info about run
log.info """\
         =====================================================================
                             SYNAPSE COUNTING PIPELINE
         =====================================================================
         ----------------------------
         Input from         : ${params.input_dir}
         Output to          : ${params.output_dir}
         ----------------------------
         Synaptic markers           : $synapse_marker
         Protein & Control          : $protein
         Nr of optimization trials  : ${params.nr_of_optimization_trials}
         ----------------------------
         Preprocessing params JSON file       : ${params.preprocessing_params}
         Optimization param ranges JSON file  : ${params.optimization_params}
         ----------------------------
         """
         .stripIndent()


// Channels
// protein and synaptic marker
protein_and_synaptic_marker = Channel.from(params.protein_and_synaptic_marker)

// the input and output folder 
input_images = Channel.fromPath(params.input_dir)
output_folder = Channel.fromPath(params.output_dir)

// the preprocessing & optimization ranges parameters
preprocess_params = Channel.fromPath(params.preprocessing_params)
optimization_params = Channel.fromPath(params.optimization_params)

// important name segments
name_segments = Channel.from(params.name_segments)

// pre and postsynapse channels
pre_synapse = Channel.from(params.presynapse_channel)
post_synapse = Channel.from(params.postsynapse_channel)

// number of optimization trials
nr_of_optimization_trials = Channel.from(params.nr_of_optimization_trials)


// the pipeline
process RunLocalPeaksOptimization {
    publishDir path: "${params.output_dir}", mode: 'copy'
    
    input:
    path input_images
    val pre_synapse
    val post_synapse
    val protein_and_synaptic_marker
    val nr_of_optimization_trials
    path preprocess_params
    val optimization_params

    output:
    path "local_peak_best_optimization_params.csv", emit: local_peaks_optimized
    path "local_peak_optimization_trial_data.csv"
    path "local_peak_optimization_all_data.csv"

    script:
    """
    local_peaks_optimize.py \
    --input_dir $input_images \
    --presynapse_channel $pre_synapse \
    --postsynapse_channel $post_synapse \
    --protein_and_synaptic_marker $protein_and_synaptic_marker \
    --nr_of_optimization_trials $nr_of_optimization_trials \
    --preprocessing_params $preprocess_params \
    --opt_param_ranges $optimization_params 
    """ 
}

process RunLocalPeaksColocalization {
    input:
    path input_images
    val name_segments
    val pre_synapse
    val post_synapse
    val protein_and_synaptic_marker
    path preprocess_params
    path local_peaks_opt_csv

    output:
    path "local_peak_coloc.csv", emit: local_peaks_ch

    script:
    """
    local_peaks_coloc.py \
    --input_dir $input_images \
    --name_segments $name_segments \
    --presynapse_channel $pre_synapse \
    --postsynapse_channel $post_synapse \
    --protein_and_synaptic_marker $protein_and_synaptic_marker \
    --preprocessing_params $preprocess_params \
    --optimized_params_df $local_peaks_opt_csv 
    """ 
}

process RunPearsonsAnalysis {
    input:
    path input_images
    val name_segments
    val pre_synapse
    val post_synapse
    val protein_and_synaptic_marker
    path preprocess_params
    
    output:
    path "pearson_results.csv" , emit: pearson_ch

    script:
    """
    pearsons.py \
    --input_dir $input_images \
    --name_segments $name_segments \
    --presynapse_channel $pre_synapse \
    --postsynapse_channel $post_synapse \
    --protein_and_synaptic_marker $protein_and_synaptic_marker  \
    --preprocessing_params $preprocess_params
    """ 
}

process RunMandersAnalysis {
    input:
    path input_images
    val name_segments
    val pre_synapse
    val post_synapse
    val protein_and_synaptic_marker
    path preprocess_params
    
    output:
    path "manders_results.csv", emit: manders_ch

    script:
    """
    manders.py \
    --input_dir $input_images \
    --name_segments $name_segments \
    --presynapse_channel $pre_synapse \
    --postsynapse_channel $post_synapse \
    --protein_and_synaptic_marker $protein_and_synaptic_marker \
    --preprocessing_params $preprocess_params
    """ 

}

process RunOverlapAnalysis {
    input:
    path input_images
    val name_segments
    val pre_synapse
    val post_synapse
    val protein_and_synaptic_marker
    path preprocess_params

    output:
    path "overlap_results.csv" , emit: overlap_ch

    script:
    """
    overlap.py \
    --input_dir $input_images \
    --name_segments $name_segments \
    --presynapse_channel $pre_synapse \
    --postsynapse_channel $post_synapse \
    --protein_and_synaptic_marker $protein_and_synaptic_marker \
    --preprocessing_params $preprocess_params
    """ 
}

process RunPunctaAnalysis {
    input:
    path input_images
    val name_segments
    val pre_synapse
    val post_synapse
    val protein_and_synaptic_marker
    path preprocess_params

    output:
    path "puncta_results.csv" , emit: puncta_ch

    script:
    """
    puncta.py \
    --input_dir $input_images \
    --name_segments $name_segments \
    --presynapse_channel $pre_synapse \
    --postsynapse_channel $post_synapse \
    --protein_and_synaptic_marker $protein_and_synaptic_marker \
    --preprocessing_params $preprocess_params
    """
}

process ResultsPlots {
    publishDir path: "${params.output_dir}", mode: 'copy'
    
    input:
    path from_local_peaks_calc
    path from_manders
    path from_overlap
    path from_pearson
    path from_puncta
    val protein_and_synaptic_marker

    script:
    """
    results.py \
    --local_peaks_df $from_local_peaks_calc \
    --manders_df $from_manders \
    --overlap_df $from_overlap \
    --pearson_df $from_pearson \
    --puncta_df $from_puncta \
    --protein_and_synaptic_marker $protein_and_synaptic_marker
    """

    output:
    path "metric_results.csv"
    path "coloc_metric_internal_controls.png"
    path "all_statistics.csv"
    path "significant_statistics.csv"
    path "combined_plots.png"
}


workflow {
    RunLocalPeaksOptimization(
        input_images,
        pre_synapse,
        post_synapse,
        protein_and_synaptic_marker,
        nr_of_optimization_trials,
        preprocess_params,
        optimization_params
        )
    
    local_peaks_ch = RunLocalPeaksColocalization(
                    input_images,
                    name_segments, 
                    pre_synapse,
                    post_synapse,
                    protein_and_synaptic_marker,
                    preprocess_params,
                    RunLocalPeaksOptimization.out.local_peaks_optimized
                    )
    
    pearson_ch = RunPearsonsAnalysis(
                input_images, 
                name_segments,
                pre_synapse,
                post_synapse,
                protein_and_synaptic_marker,
                preprocess_params
                )
    
    manders_ch = RunMandersAnalysis(
                input_images, 
                name_segments,
                pre_synapse,
                post_synapse,
                protein_and_synaptic_marker,
                preprocess_params
                )
    
    overlap_ch = RunOverlapAnalysis(
                input_images, 
                name_segments,
                pre_synapse,
                post_synapse,
                protein_and_synaptic_marker,
                preprocess_params
                )
    
    puncta_ch = RunPunctaAnalysis(
                input_images, 
                name_segments,
                pre_synapse,
                post_synapse,
                protein_and_synaptic_marker,
                preprocess_params
                )

    ResultsPlots(
        local_peaks_ch, 
        manders_ch, 
        overlap_ch, 
        pearson_ch, 
        puncta_ch,
        protein_and_synaptic_marker
        )
}
