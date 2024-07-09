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


// the pipeline
process RunLocalPeaksOptimization {
    publishDir "${params.output_dir}", mode: 'copy'
    
    input:
    path input_images
    path preprocess_params
    path optimization_params

    script:
    """
    python ${baseDir}/local_peaks_optimize.py \
    --input_dir ${params.input_dir} \
    --presynapse_channel ${params.presynapse_channel} \
    --postsynapse_channel ${params.postsynapse_channel} \
    --protein_and_synaptic_marker ${params.protein_and_synaptic_marker} \
    --nr_of_optimization_trials ${params.nr_of_optimization_trials} \
    --preprocessing_params $preprocess_params \
    --opt_param_ranges $optimization_params 
    """ 

    output:
    path "local_peak_best_optimization_params.csv", emit: local_peaks_optimized
    path "local_peak_optimization_trial_data.csv"
    path "local_peak_optimization_all_data.csv"
}

process RunLocalPeaksColocalization {
    input:
    path input_images
    path local_peaks_opt_csv
    path preprocess_params

    script:
    """
    python ${baseDir}/local_peaks_coloc.py \
    --input_dir ${params.input_dir} \
    --name_segments ${params.name_segments} \
    --presynapse_channel ${params.presynapse_channel} \
    --postsynapse_channel ${params.postsynapse_channel} \
    --protein_and_synaptic_marker ${params.protein_and_synaptic_marker} \
    --preprocessing_params $preprocess_params \
    --optimized_params_df $local_peaks_opt_csv 
    """ 

    output:
    path "local_peak_coloc.csv", emit: local_peaks_ch
}

process RunPearsonsAnalysis {
    input:
    path input_images
    path preprocess_params

    script:
    """
    python ${baseDir}/pearsons.py \
    --input_dir ${params.input_dir} \
    --name_segments ${params.name_segments} \
    --presynapse_channel ${params.presynapse_channel} \
    --postsynapse_channel ${params.postsynapse_channel} \
    --protein_and_synaptic_marker ${params.protein_and_synaptic_marker} \
    --preprocessing_params $preprocess_params
    """ 

    output:
    path "pearson_results.csv", emit: pearson_ch
}

process RunMandersAnalysis {
    input:
    path input_images
    path preprocess_params

    script:
    """
    python ${baseDir}/manders.py \
    --input_dir ${params.input_dir} \
    --name_segments ${params.name_segments} \
    --presynapse_channel ${params.presynapse_channel} \
    --postsynapse_channel ${params.postsynapse_channel} \
    --protein_and_synaptic_marker ${params.protein_and_synaptic_marker} \
    --preprocessing_params $preprocess_params
    """ 

    output:
    path "manders_results.csv", emit: manders_ch
}

process RunOverlapAnalysis {
    input:
    path input_images
    path preprocess_params

    script:
    """
    python ${baseDir}/overlap.py \
    --input_dir ${params.input_dir} \
    --name_segments ${params.name_segments} \
    --presynapse_channel ${params.presynapse_channel} \
    --postsynapse_channel ${params.postsynapse_channel} \
    --protein_and_synaptic_marker ${params.protein_and_synaptic_marker} \
    --preprocessing_params $preprocess_params
    """ 

    output:
    path "overlap_results.csv", emit: overlap_ch
}

process RunPunctaAnalysis {
    input:
    path input_images
    path preprocess_params

    script:
    """
    python ${baseDir}/puncta.py \
    --input_dir ${params.input_dir} \
    --name_segments ${params.name_segments} \
    --presynapse_channel ${params.presynapse_channel} \
    --postsynapse_channel ${params.postsynapse_channel} \
    --protein_and_synaptic_marker ${params.protein_and_synaptic_marker} \
    --preprocessing_params $preprocess_params
    """

    output:
    path "puncta_results.csv", emit: puncta_ch
}

process ResultsPlots {
    publishDir "${params.output_dir}", mode: 'copy'
    
    input:
    path from_local_peaks_calc
    path from_manders
    path from_overlap
    path from_pearson
    path from_puncta

    script:
    """
    python ${baseDir}/results.py \
    --local_peaks_df $from_local_peaks_calc \
    --manders_df $from_manders \
    --overlap_df $from_overlap \
    --pearson_df $from_pearson \
    --puncta_df $from_puncta \
    --protein_and_synaptic_marker ${params.protein_and_synaptic_marker}
    """

    output:
    path "metric_results.csv"
    path "coloc_metric_internal_controls.png"
    path "all_statistics.csv"
    path "significant_statistics.csv"
    path "combined_plots.png"
}


workflow {
    input_images = Channel.fromPath(params.input_dir)
    preprocess_params = Channel.fromPath(params.preprocessing_params)
    optimization_params = Channel.fromPath(params.optimization_params)

    RunLocalPeaksOptimization(input_images, preprocess_params, optimization_params)
    local_peaks_ch = RunLocalPeaksColocalization(
                    input_images, 
                    RunLocalPeaksOptimization.out.local_peaks_optimized, 
                    preprocess_params
                    )

    manders_ch = RunMandersAnalysis(input_images, preprocess_params)
    overlap_ch = RunOverlapAnalysis(input_images, preprocess_params)
    pearson_ch = RunPearsonsAnalysis(input_images, preprocess_params)
    puncta_ch = RunPunctaAnalysis(input_images, preprocess_params)

    ResultsPlots(local_peaks_ch, manders_ch, overlap_ch, pearson_ch, puncta_ch)
}
