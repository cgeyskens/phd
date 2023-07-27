#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

process runPythonScript {

    script:
    """
    python ${baseDir}/bin/test-file-1.py
    """
}

workflow {
    runPythonScript()
}




