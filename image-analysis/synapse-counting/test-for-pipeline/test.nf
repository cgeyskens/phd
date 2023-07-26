#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

// Define the process that runs the Python script
process runPythonScript {
    
    conda "/home/cgeyskens/anaconda3/envs/image-analysis"
    
    // Input files (you can use glob patterns to match multiple files)
    input:
    file path

    // Define the script to be executed
    script:
    """
    test-file-1.py "$path"
    """

    // Output directory for the generated CSV files
    output:
    file "output_data/test_{path.baseName}.csv"
}


// Create a channel that matches the input files using the provided pattern
input_files = file("/mnt/d/code/phd/image-analysis/synapse-counting/test-images-VLGUT1-PSD95/OE_Exp1_IHC_Exp1_HA-GPR37L1_555-VGLUT1_647-PSD95_63X_airyscan_1.8zoom_CA1_SO.czi")

// Run the process for each input file in parallel
workflow {
    runPythonScript(input_files)
}






