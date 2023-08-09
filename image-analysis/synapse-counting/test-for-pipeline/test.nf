#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

input_ch = Channel.fromPath("/mnt/d/code/phd/image-analysis/synapse-counting/test-images-2/*.czi")

process runPythonScript {
    input:
    path input_dir

    script:
    """
    python ${baseDir}/bin/test-file-1.py ${input_dir}
    """
}

workflow {
    runPythonScript(input_ch)
}




