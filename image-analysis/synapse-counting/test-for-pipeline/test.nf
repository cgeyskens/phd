#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

params.input_file = "/mnt/d/code/phd/image-analysis/synapse-counting/test-images-VLGUT1-PSD95/OE_Exp1_IHC_Exp1_HA-GPR37L1_555-VGLUT1_647-PSD95_63X_airyscan_1.8zoom_CA1_SO.czi"

process runPythonScript {
    input:
    file input_file_path from params.input_file

    script:
    """
    python ${baseDir}/bin/test-file-1.py input_file_path ${input_file_path}
    """
}

workflow {
    runPythonScript()
}




