# test-2-mean-threshold.py
import csv
import pandas as pd
import argparse


# adding parser arguments
parser = argparse.ArgumentParser(description="calculate mean thresholds")
parser.add_argument("input_dir", type=str, help="directory to input files")
parser.add_argument("output_dir", type=str, help="directory to output folder")

args = parser.parse_args()

# assigning the parser arguments
input_folder = args.input_dir
output_folder = args.output_dir

# input_folder = "/mnt/d/code/phd/image-analysis/synapse-counting/output_data/"
# output_folder = "/mnt/d/code/phd/image-analysis/synapse-counting/output_data/"

# specifying the path to the csv file that needs to be imported
input_file_path = input_folder + "thresholds.csv"

df = pd.read_csv(input_file_path)
        
# calculating the average of the thresholds for vglut1 and psd95
df_threshold_means = df.groupby("image file name").agg({"vglut1_threshold": "mean", "psd95_threshold": "mean"})

# writing to a csv file
output_filename = "mean_threshold_values.csv"
output_path = output_folder + output_filename
df_threshold_means.to_csv(output_path)        

