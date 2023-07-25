#!/usr/bin/env nextflow

params.input_dir = "/mnt/d/code/phd/image-analysis/synapse-counting/test-images-VLGUT1-PSD95/"
params.output_dir = "/mnt/d/code/phd/image-analysis/synapse-counting/output_data/" 

// Function to process an image
def processImage(imagePath) {
    script:
    """
    echo "Processing image: $imagePath"
    python test-file-2.py "$imagePath" > output.txt
    """
}

// Define the main workflow
workflow {
    // Get a list of all image files in the input directory
    images = Channel.fromPath(params.input_dir)

    // Process each image in parallel using the 'processImage' function
    processedData = images.map { imagePath ->
        processImage(imagePath)
    }

    // Save the DataFrame to CSV for each image
    process saveToCSV {
        input:
        file imageFile from processedData

        output:
        file "${params.output_dir}/test_${imageFile.baseName}.csv" into outputFiles

        script:
        """
        # Load the DataFrame from 'output.txt' and write to CSV using pandas
        python -c 'import pandas as pd; df = pd.read_csv("output.txt", header=None); df.to_csv("${outputFiles}", index=False)'
        """
    }
}




