#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

// info about run
log.info """\
         =====================================================================
                             SYNAPSE COUNTING PIPELINE
         =====================================================================
         ----------------------------
         Input from         : ${params.input_dir}
         Output to          : ${params.output_dir}
         ----------------------------
         Synaptic marker combo      : ${params.protein_and_synaptic_marker}
         Nr of optimization trials  : ${params.nr_of_optimization_trials}
         ----------------------------
         """
         .stripIndent()


// the pipeline
process RunLocalPeaksOptimization {
    publishDir "${params.output_dir}", mode: 'copy'
    
    input:
    path input_images

    script:
    """
    python ${baseDir}/local_peaks_optimize.py \
    --input_dir ${params.input_dir} \
    --presynapse_channel ${params.presynapse_channel} \
    --postsynapse_channel ${params.postsynapse_channel} \
    --protein_and_synaptic_marker ${params.protein_and_synaptic_marker} \
    --nr_of_optimization_trials ${params.nr_of_optimization_trials}
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

    script:
    """
    python ${baseDir}/local_peaks_coloc.py \
    --input_dir ${params.input_dir} \
    --name_segments ${params.name_segments} \
    --presynapse_channel ${params.presynapse_channel} \
    --postsynapse_channel ${params.postsynapse_channel} \
    --protein_and_synaptic_marker ${params.protein_and_synaptic_marker} \
    --optimized_params_df $local_peaks_opt_csv
    """ 
    
    output:
    path "local_peak_coloc.csv"
}

process RunPearsonsAnalysis {
    input:
    path input_images

    script:
    """
    python ${baseDir}/pearsons.py \
    --input_dir ${params.input_dir} \
    --name_segments ${params.name_segments} \
    --presynapse_channel ${params.presynapse_channel} \
    --postsynapse_channel ${params.postsynapse_channel}
    """ 

    output:
    path "pearson_results.csv"
}

process RunMandersAnalysis {
    input:
    path input_images

    script:
    """
    python ${baseDir}/manders.py \
    --input_dir ${params.input_dir} \
    --name_segments ${params.name_segments} \
    --presynapse_channel ${params.presynapse_channel} \
    --postsynapse_channel ${params.postsynapse_channel} \
    --protein_and_synaptic_marker ${params.protein_and_synaptic_marker}
    """ 

    output:
    path "manders_results.csv"
}

process RunOverlapAnalysis {
    input:
    path input_images

    script:
    """
    python ${baseDir}/overlap.py \
    --input_dir ${params.input_dir} \
    --name_segments ${params.name_segments} \
    --presynapse_channel ${params.presynapse_channel} \
    --postsynapse_channel ${params.postsynapse_channel} \
    --protein_and_synaptic_marker ${params.protein_and_synaptic_marker}
    """ 

    output:
    path "overlap_results.csv"
}

process RunPunctaAnalysis {
    input:
    path input_images

    script:
    """
    python ${baseDir}/puncta.py \
    --input_dir ${params.input_dir} \
    --name_segments ${params.name_segments} \
    --presynapse_channel ${params.presynapse_channel} \
    --postsynapse_channel ${params.postsynapse_channel} \
    --protein_and_synaptic_marker ${params.protein_and_synaptic_marker}
    """

    output:
    path "puncta_results.csv"
}

process ResultsPlots {
    publishDir "${params.output_dir}", mode: 'move'
    
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
}


workflow {
    input_images = Channel.fromPath(params.input_dir)

    RunLocalPeaksOptimization(input_images)
    local_peaks_ch = RunLocalPeaksColocalization(input_images, RunLocalPeaksOptimization.out.local_peaks_optimized)
    manders_ch = RunMandersAnalysis(input_images)
    overlap_ch = RunOverlapAnalysis(input_images)
    pearson_ch = RunPearsonsAnalysis(input_images)
    puncta_ch = RunPunctaAnalysis(input_images)

    ResultsPlots(local_peaks_ch, manders_ch, overlap_ch, pearson_ch, puncta_ch)
}
