#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

// specifying the input and output directories
input_dir = "/mnt/d/code/phd/image-analysis/synapse-counting/test-images-2/OE_Exp1_IHC_Exp1_HA-GPR37L1_555-VGLUT1_647-PSD95_63X_airyscan_1.8zoom_CA1_SO.czi"
output_dir = "/mnt/d/code/phd/image-analysis/synapse-counting/output_data/"

process runPythonScript {
    input:
    path x

    script:
    """
    python ${baseDir}/bin/test-file-1.py ${input_dir} ${output_dir}
    """
}

workflow {
    runPythonScript(
        input_ch = input_dir
        )
}

