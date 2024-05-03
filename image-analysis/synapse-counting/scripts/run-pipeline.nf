#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

// specifying the input and output directories
input_dir = "/Volumes/KINGSTON/code/phd/image-analysis/synapse-counting/images-VCAM1-LacZ"
intermediate_dir = "/Volumes/KINGSTON/code/phd/image-analysis/synapse-counting/intermediate_data"
output_dir = "/Volumes/KINGSTON/code/phd/image-analysis/synapse-counting/output_data"
name_segments = "0 1 2 3 6 8 9"

process RunPearsonsAnalysis {
    input:
    path x

    script:
    """
    python ${baseDir}/pearsons.py ${input_dir} ${intermediate_dir} ${name_segments}
    """ 

    output:
    stdout
}

process RunMandersAnalysis {
    input:
    path x

    script:
    """
    python ${baseDir}/manders.py ${input_dir} ${intermediate_dir} ${name_segments}
    """ 

    output:
    stdout
}

process OverlapAnalysis {
    input:
    path x

    script:
    """
    python ${baseDir}/overlap.py ${input_dir} ${intermediate_dir} ${name_segments}
    """ 

    output:
    stdout
}

process PunctaAnalysis {
    input:
    path x

    script:
    """
    python ${baseDir}/puncta.py ${input_dir} ${intermediate_dir} ${name_segments}
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
    python ${baseDir}/results.py ${intermediate_dir} ${output_dir}
    """

    output:
    stdout
    }


workflow {
    input_ch = input_dir
    pearson_ch = RunPearsonsAnalysis(input_ch)
    manders_ch = RunMandersAnalysis(input_ch)
    overlap_ch = OverlapAnalysis(input_ch)
    puncta_ch = PunctaAnalysis(input_ch)
    ResultsPlots(pearson_ch, manders_ch, overlap_ch, puncta_ch)
}




    // RunPearsonsAnalysis(
    //     input_ch = input_dir
    //     )
    // RunMandersAnalysis(
    //     input_ch = input_dir
    //     )
    // OverlapAnalysis(
    //     input_ch = input_dir
    // )
    // PunctaAnalysis(
    //     input_ch = input_dir
    // )
    // ResultsPlots(
    //     input_ch = intermediate_dir
    // )