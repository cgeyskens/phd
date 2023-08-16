# test-file-2.ipynb

from pathlib import Path
import czifile
import pandas as pd
import os
from microfilm.microplot import microshow
import numpy as np
import argparse

# defining the input_folder and output_folder
# input_folder = "/mnt/d/code/phd/image-analysis/synapse-counting/test-images-VLGUT1-PSD95/"
# output_folder = "/mnt/d/code/phd/image-analysis/synapse-counting/output_data/"

# adding parser arguments
parser = argparse.ArgumentParser(description='process input files')
parser.add_argument('input_dir', type=str, help='directory to input files')
parser.add_argument('output_dir', type=str, help='directory to output folder')

args = parser.parse_args()

# assigning the parser arguments
input_folder = args.input_dir
output_folder = args.output_dir

# get a list of files in that input_folder
file_list = os.listdir(input_folder)

# defing the image mean intensity
def process(filename):
    
    # getting the mean intensity of the image
    image = czifile.imread(filename)
    intensity = np.mean(image)
    
    return intensity

# getting the right filenames
def image_filename(filename):
    
    # get the filename
    name_of_file = os.path.splitext(os.path.basename(filename))[0]
    split_filename = name_of_file.split("_")

    # get only the experimental parameters from the filename
    index_nums = [0, 1, 2, 3, 10, 11] # the indexes of the elements that I would like to extract fro; the filename
    desired_parts = [split_filename[val] for val in index_nums]
    desired_filename = "_".join(desired_parts)
    
    return desired_filename

# getting the values
intensities = [process(input_folder + file) for file in file_list]
filenames = [image_filename(input_folder + file) for file in file_list]

# merging into df
d = {"image_names": filenames, "intensities": intensities}
df = pd.DataFrame(d).set_index("image_names")

# writing the df out into a csv
output_filename = "test1" + ".csv"
output_path = output_folder + output_filename
df.to_csv(output_path)