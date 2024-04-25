import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import glob
import os
import argparse
import sys
import warnings

# with this peice of code, it will recognize the custom modules
project_root = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
sys.path.append(project_root)

# custom modules
from synapse_counting import results_plotting

# adding parser arguments
parser = argparse.ArgumentParser(description="process input files")
parser.add_argument("input_dir", type=str, help="directory to input files")
parser.add_argument("output_dir", type=str, help="directory to output folder")

args = parser.parse_args()

# suppress future warnings
warnings.filterwarnings("ignore", category=FutureWarning)

# assigning the parser arguments
input_folder = args.input_dir
output_folder = args.output_dir

# getting all files in the folder
all_files = glob.glob(f"{input_folder}/*.csv")

# read in the first .csv from the folder
base_df = pd.read_csv(all_files[0])

# setting the index for common column such that we can merge the following dataframes based on the index
base_df.set_index("img_filename", inplace = True)

# loop through the csv and merge with the base_df based on the index
for filename in all_files[1:]:
  df = pd.read_csv(filename)
  print(df.columns)
  df.set_index("img_filename", inplace = True)
  base_df = pd.concat([base_df, df], axis=1)
merged_df = base_df.drop(["Unnamed: 0"], axis = 1)
merged_df = merged_df.reset_index() # this is the final df with all the data

# formats the data in the right format
df_melted = results_plotting.data_formatting(df=merged_df, id_vars=["img_filename"], 
                                            value_vars=["overlap_um", "overlap_um_rot"],
                                            value_name="overlap in um")
# plots the results
p = results_plotting.plot_data(df=df_melted, x = "hippocampal layer", y = "overlap in um", extra_y = 100, title = 'Overlap by Hippocampal Layer and Condition')

# outputs the plot
output_file = os.path.join(output_folder, "synapse_measures.png")
plot = p.get_figure()
plt.savefig(output_file, dpi = 300)
