# test-file-1.py

from pathlib import Path
import czifile
import pandas as pd
import os
import argparse

parser = argparse.ArgumentParser(description='process input file')
parser.add_argument('input_file_path', type=str, help='path to input file')
args = parser.parse_args()

path = args.input_file_path

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
output_path = "/mnt/d/code/phd/image-analysis/synapse-counting/output_data/" + output_filename
df_final.to_csv(output_path)
