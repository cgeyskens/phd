# test-file-2.py

from pathlib import Path
import czifile
import pandas as pd
import os

def process_czi_file(path):
    
    # getting the shape of the image
    image = czifile.imread(path)
    data = image.shape

    # get the filename
    filename = os.path.splitext(os.path.basename(path))[0]
    split_filename = filename.split("_")

    # get only the experimental parameters from the filename
    index_nums = [0, 1, 2, 3, 10, 11] # the indexes of the elements that I would like to extract
    desired_parts = [split_filename[val] for val in index_nums]
    desired_filename = "_".join(desired_parts)

    # converting tuple to df
    df = pd.DataFrame(data, columns = [desired_filename] )

    # pivoting the dataframe
    dic = df.to_dict(orient='dict')
    df_final = pd.DataFrame.from_dict(dic, orient='index')
    
    return df_final

def process_folder(folder_path, output_folder):
    # Process all CZI files in the given folder
    file_paths = Path(folder_path).rglob('*.czi')
    
    for path in file_paths:
        df = process_czi_file(str(path))
        
        # Writing the DataFrame to a CSV file
        output_filename = "test_" + ".csv"
        output_path = os.path.join(output_folder, output_filename)
        df.to_csv(output_path)

if __name__ == "__main__":
    input_folder = "/mnt/d/code/phd/image-analysis/synapse-counting/test-images-2/"
    output_folder = "/mnt/d/code/phd/image-analysis/synapse-counting/output_data/"
    
    process_folder(input_folder, output_folder)

