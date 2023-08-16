# test-file-1.py

from pathlib import Path
import czifile
import pandas as pd
import os
import argparse

# input folder file path: "/mnt/d/code/phd/image-analysis/synapse-counting/test-images-2/OE_Exp1_IHC_Exp1_HA-GPR37L1_555-VGLUT1_647-PSD95_63X_airyscan_1.8zoom_CA1_SO.czi"
# output folder file path: "/mnt/d/code/phd/image-analysis/synapse-counting/output_data/"


# adding parser argument for nextflow script
parser = argparse.ArgumentParser(description='process input files')
parser.add_argument('input_dir', type=str, help='directory to input files')
parser.add_argument('output_dir', type=str, help='directory to output folder')

args = parser.parse_args()

# assigning the parser arguments
path = args.input_dir
output_folder = args.output_dir

# define here a bath process function
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

# Writing the dataframe in a csv format into the a specific folder
output_filename = "test_" + desired_filename + ".csv"
output_path = output_folder + output_filename
df_final.to_csv(output_path)
