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

# getting the merged dataframe with all the data
all_files = glob.glob(f"{input_folder}/*.csv") # getting all files in the folder
base_df = pd.read_csv(all_files[0]) # read in the first .csv from the folder
base_df.set_index("img_filename", inplace = True) # setting the index for common column such that we can merge the following dataframes based on the index
for filename in all_files[1:]: # loop through the csv and merge with the base_df based on the index
  df = pd.read_csv(filename)
  print(df.columns)
  df.set_index("img_filename", inplace = True)
  base_df = pd.concat([base_df, df], axis=1)
merged_df = base_df.drop(["Unnamed: 0"], axis = 1)
merged_df = merged_df.reset_index() # this is the final df with all the data

# formats the data in the right format for the plots
df_melted_overlap_um2 = results_plotting.data_formatting(df = merged_df, 
                                                         id_vars = ["img_filename"], 
                                                         value_vars = ["overlap_um2", "overlap_um2_rot"],
                                                         value_name = "overlap in um2")

df_melted_pearson = results_plotting.data_formatting(df = merged_df, 
                                                     id_vars = ["img_filename"], 
                                                     value_vars = ["pearson_cor", "pearson_cor_rot"],
                                                     value_name = "pearsons correlation")

df_melted_manders = results_plotting.data_formatting(df = merged_df, 
                                                     id_vars = ["img_filename"],
                                                     value_vars = ["overlap_coeff", "overlap_coeff_rot"],
                                                     value_name = "Manders overlap coefficient")

# plotting the results
fig, axes = plt.subplots(nrows=3, ncols = 1, figsize=(10, 15))
plot_overlap = results_plotting.plot_data(df = df_melted_overlap_um2, 
                                          x = "hippocampal layer", 
                                          y = "overlap in um2", 
                                          extra_y_upper = 100, 
                                          title = 'Overlap by hippocampal Layer and condition', 
                                          ax=axes[0])
plot_pearson = results_plotting.plot_data(df = df_melted_pearson, 
                                          x = "hippocampal layer", 
                                          y = "pearsons correlation", 
                                          extra_y_upper = 0.1, 
                                          title = 'Pearsons correlation by hippocampal Layer and condition', 
                                          extra_y_lower = -0.05, 
                                          ax=axes[1])
plot_manders = results_plotting.plot_data(df = df_melted_manders, 
                                          x = "hippocampal layer", 
                                          y = "Manders overlap coefficient", 
                                          extra_y_upper = 0.1, 
                                          title = 'Manders overlap coefficient by hippocampal Layer and condition', 
                                          ax=axes[2])

# saving the combined plots into one file
combined_plots_figure = plt.gcf()  
plt.tight_layout()  
output_plot_path = os.path.join(output_folder, "combined_plots.png")
combined_plots_figure.savefig(output_plot_path, dpi = 300)  

