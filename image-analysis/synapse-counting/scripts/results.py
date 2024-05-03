import pandas as pd
import numpy as np
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

### ---------------------------- Plotting the internal rotated controls vs actual --------------------------------- ###
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
                                          hue = "condition",
                                          ax=axes[0])
plot_pearson = results_plotting.plot_data(df = df_melted_pearson, 
                                          x = "hippocampal layer", 
                                          y = "pearsons correlation", 
                                          extra_y_upper = 0.1, 
                                          title = 'Pearsons correlation by hippocampal Layer and condition', 
                                          hue = "condition",
                                          extra_y_lower = -0.05, 
                                          ax=axes[1])
plot_manders = results_plotting.plot_data(df = df_melted_manders, 
                                          x = "hippocampal layer", 
                                          y = "Manders overlap coefficient", 
                                          extra_y_upper = 0.1, 
                                          title = 'Manders overlap coefficient by hippocampal Layer and condition',
                                          hue = "condition", 
                                          ax=axes[2])

# saving the combined plots into one file
internal_control = plt.gcf()  
plt.tight_layout()  
output_internal_control_plot_path = os.path.join(output_folder, "internal_control.png")
internal_control.savefig(output_internal_control_plot_path, dpi = 300)  

### ------------------------------ Plotting the LacZ-gRNA vs the candidate-gRNA --------------------------------- ###
# formatting the data
gRNA_data = merged_df
gRNA_data['gRNA'] = gRNA_data['img_filename'].apply(lambda x: x.split('_')[-3:-2]).apply(lambda x: '_'.join(x))
gRNA_data['hippocampal_layer'] = gRNA_data['img_filename'].apply(lambda x: x.split('_')[-2:]).apply(lambda x: ' '.join(x))

# plotting the results
# Create the figure object
fig, axes = plt.subplots(nrows=6, ncols=2, figsize=(15, 30))

# colocalization metrics
plot_overlap = results_plotting.plot_data(df=gRNA_data, x="hippocampal_layer", y="overlap_um2", extra_y_upper=50, title='Overlap by hippocampal Layer and condition', ax=axes[0,0])
plot_pearson = results_plotting.plot_data(df=gRNA_data, x="hippocampal_layer", y="pearson_cor", extra_y_upper=0.1, title='Pearsons correlation by hippocampal Layer and condition', extra_y_lower=-0.05, ax=axes[0,1])
plot_manders = results_plotting.plot_data(df=gRNA_data, x="hippocampal_layer", y="overlap_coeff", extra_y_upper=0.1, title='Manders overlap coefficient by hippocampal Layer and condition', extra_y_lower=-0.05, ax=axes[1,0])

# single synapse marker metrics
plot_pre_mfi = results_plotting.plot_data(df=gRNA_data, x="hippocampal_layer", y="presynapse_image_mfi", extra_y_upper=200, title='Presynaptic MFI by hippocampal Layer and condition', extra_y_lower=-0.05, ax=axes[2,0])
plot_pre_mfi = results_plotting.plot_data(df=gRNA_data, x="hippocampal_layer", y="postsynapse_image_mfi", extra_y_upper=200, title='Postsynaptic MFI by hippocampal Layer and condition', extra_y_lower=-0.05, ax=axes[2,1])
plot_pre_density = results_plotting.plot_data(df=gRNA_data, x="hippocampal_layer", y="pre_puncta_density_per_100_um2", extra_y_upper=100, title='Presynaptic puncta density by hippocampal Layer and condition', extra_y_lower=-0.05, ax=axes[3,0])
plot_post_density = results_plotting.plot_data(df=gRNA_data, x="hippocampal_layer", y="post_puncta_density_per_100_um2", extra_y_upper=100, title='Postsynaptic puncta density by hippocampal Layer and condition', extra_y_lower=-0.05, ax=axes[3,1])
plot_pre_area = results_plotting.plot_data(df=gRNA_data, x="hippocampal_layer", y="pre_staining_area_um2", extra_y_upper=50, title='Presynaptic staining area by hippocampal Layer and condition', extra_y_lower=-0.05, ax=axes[4,0])
plot_post_area = results_plotting.plot_data(df=gRNA_data, x="hippocampal_layer", y="post_staining_area_um2", extra_y_upper=50, title='Postsynaptic staining area by hippocampal Layer and condition', extra_y_lower=-0.05, ax=axes[4,1])
plot_pre_puncta_size = results_plotting.plot_data(df=gRNA_data, x="hippocampal_layer", y="pre_mean_puncta_size_um2", extra_y_upper=0.05, title='Presynaptic puncta size (um2) by hippocampal Layer and condition', extra_y_lower=-0.05, ax=axes[5,0])
plot_post_puncta_size = results_plotting.plot_data(df=gRNA_data, x="hippocampal_layer", y="post_mean_puncta_size_um2", extra_y_upper=0.05, title='Postsynaptic puncta size (um2) by hippocampal Layer and condition', extra_y_lower=-0.05, ax=axes[5,1])

# Save the figure explicitly
combined_plots_figure = plt.gcf()  # Get the current figure object
plt.tight_layout()  # Adjust layout to prevent overlap
output_combined_plot_path = os.path.join(output_folder, "combined_plots.png")
combined_plots_figure.savefig(output_combined_plot_path)  # Save the figure object

# save the merged_df 
output_csv_path = os.path.join(output_folder, "results.csv")
merged_df.to_csv(output_csv_path)

### -------------------------------------------- calculate statistics ----------------------------------------------- ###
merged_df['Brain'] = merged_df['img_filename'].str.split('_').str[2] # create new column with brain replicates

hippocampal_layers = ["CA1 SLM", "CA1 SR", "CA1 SO", "CA3 SO", "CA3 SL", "CA3 SR", "DG ML", "DG Hilus"]
metrics = ["overlap_um2", "pearson_cor", "overlap_coeff", "presynapse_image_mfi", "postsynapse_image_mfi", "pre_puncta_density_per_100_um2", "post_puncta_density_per_100_um2", "pre_staining_area_um2", "post_staining_area_um2", "pre_mean_puncta_size_um2", "post_mean_puncta_size_um2"]

# calculate the statistics with a t-test
statistics_results = []
for hippocampal_layer in hippocampal_layers:
    for metric in metrics:
        _ , table, _ , p_value = results_plotting.check_statistics(merged_df, hippocampal_layer = hippocampal_layer, metric = metric, candidate_gRNA = "VCAM1-gRNA")
        statistics_results.append({"hippocampal_layer": hippocampal_layer, "metric": metric, "p_value": p_value})

# add another column if the results are significant or not
df_statistics_results = pd.DataFrame(statistics_results) 
df_statistics_results['statistical significant?'] = np.where(df_statistics_results["p_value"] <= 0.05, 'Yes', 'No')

# create a new dataframe with only the significant results
p_value_threshold = 0.05  # significant results
df_significant_statistics_results = df_statistics_results[df_statistics_results["p_value"] <= p_value_threshold]

# save the statistics 
output_all_statistics_path = os.path.join(output_folder, "all_statistics.csv")
df_statistics_results.to_csv(output_all_statistics_path)

output_significant_statistics_path = os.path.join(output_folder, "significant_statistics.csv")
df_significant_statistics_results.to_csv(output_significant_statistics_path)

