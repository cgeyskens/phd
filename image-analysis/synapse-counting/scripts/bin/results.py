#!/usr/bin/env python

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os
import argparse
import sys
import warnings


### --------------------------------- parser arguments and loading custom library --------------------------------- ###

# # Get the parent directory of scripts (parent_directory)
# parent_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))
# # Get the synapse_counting directory
# synapse_counting_dir = os.path.join(parent_dir, 'synapse_counting')
# # Append the synapse_counting directory to sys.path
# sys.path.append(synapse_counting_dir)

from synapse_counting import results_plotting

# adding parser arguments
parser = argparse.ArgumentParser(description="process input files")
parser.add_argument("--local_peaks_df",type=str, help="dataframe with local peaks results for each hippcampal layer")
parser.add_argument("--manders_df",type=str, help="dataframe with manders results for each hippcampal layer")
parser.add_argument("--overlap_df",type=str, help="dataframe with overlap results for each hippcampal layer")
parser.add_argument("--pearson_df",type=str, help="dataframe with pearson results for each hippcampal layer")
parser.add_argument("--puncta_df",type=str, help="dataframe with puncta results for each hippcampal layer")
parser.add_argument("--protein_and_synaptic_marker", type=str, help="protein that you are interested in")
args = parser.parse_args()

# suppress future warnings
warnings.filterwarnings("ignore", category=FutureWarning)

# assigning the parser arguments
protein_and_synaptic_marker = args.protein_and_synaptic_marker

# reading in the dataframes
local_peaks_df = pd.read_csv(args.local_peaks_df)
manders_df = pd.read_csv(args.manders_df)
overlap_df = pd.read_csv(args.overlap_df)
pearson_df = pd.read_csv(args.pearson_df)
puncta_df = pd.read_csv(args.puncta_df)

# same index
local_peaks_df.set_index("img_filename", inplace=True)
manders_df.set_index("img_filename", inplace=True)
overlap_df.set_index("img_filename", inplace=True)
pearson_df.set_index("img_filename", inplace=True)
puncta_df.set_index("img_filename", inplace=True)


### ---------------------------------------------- data wrangling ------------------------------------------------- ###

# getting the merged dataframe with all the data
merged_df = pd.concat([local_peaks_df, manders_df, overlap_df, pearson_df, puncta_df], axis=1)
merged_df = merged_df.drop(["Unnamed: 0"], axis = 1)
merged_df = merged_df.reset_index() # this is the final df with all the data

# formatting the data, to include a column of gRNA and hippocampal layer, important for calculating statistics
merged_df['gRNA'] = merged_df['img_filename'].apply(lambda x: x.split('_')[-3:-2]).apply(lambda x: '_'.join(x))
merged_df['hippocampal_layer'] = merged_df['img_filename'].apply(lambda x: x.split('_')[-2:]).apply(lambda x: ' '.join(x))
merged_df['Brain'] = merged_df['img_filename'].str.split('_').str[2]

# save the merged_df 
merged_df.to_csv("metric_results.csv")


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

df_melted_local_peaks = results_plotting.data_formatting(df = merged_df, 
                                                     id_vars = ["img_filename"],
                                                     value_vars = ["local_peak_colocalized_spots", "local_peak_colocalized_spots_rot"],
                                                     value_name = "Local peaks colozalized spots")

# plotting the results
fig, axes = plt.subplots(nrows=4, ncols = 1, figsize=(10, 25))
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
plot_local_peaks = results_plotting.plot_data(df = df_melted_local_peaks, 
                                          x = "hippocampal layer", 
                                          y = "Local peaks colozalized spots", 
                                          extra_y_upper = 100, 
                                          title = 'Local peaks colozalized spots by hippocampal Layer and condition',
                                          hue = "condition", 
                                          ax=axes[3])

# saving the combined plots into one file
internal_control = plt.gcf()  
plt.tight_layout()
internal_control.savefig("coloc_metric_internal_controls.png", dpi = 300)  

### ------------------------------ Plotting the LacZ-gRNA vs the candidate-gRNA (swarm plot) --------------------------------- ###

# # plotting the results
# # Create the figure object
# fig, axes = plt.subplots(nrows=6, ncols=2, figsize=(15, 30))

# # colocalization metrics
# plot_overlap = results_plotting.plot_data(df=gRNA_data, x="hippocampal_layer", y="overlap_um2", extra_y_upper=50, title='Overlap by hippocampal Layer and condition', ax=axes[0,0])
# plot_pearson = results_plotting.plot_data(df=gRNA_data, x="hippocampal_layer", y="pearson_cor", extra_y_upper=0.1, title='Pearsons correlation by hippocampal Layer and condition', extra_y_lower=-0.05, ax=axes[0,1])
# plot_manders = results_plotting.plot_data(df=gRNA_data, x="hippocampal_layer", y="overlap_coeff", extra_y_upper=0.1, title='Manders overlap coefficient by hippocampal Layer and condition', extra_y_lower=-0.05, ax=axes[1,0])

# # single synapse marker metrics
# plot_pre_mfi = results_plotting.plot_data(df=gRNA_data, x="hippocampal_layer", y="presynapse_image_mfi", extra_y_upper=200, title='Presynaptic MFI by hippocampal Layer and condition', extra_y_lower=-0.05, ax=axes[2,0])
# plot_pre_mfi = results_plotting.plot_data(df=gRNA_data, x="hippocampal_layer", y="postsynapse_image_mfi", extra_y_upper=200, title='Postsynaptic MFI by hippocampal Layer and condition', extra_y_lower=-0.05, ax=axes[2,1])
# plot_pre_density = results_plotting.plot_data(df=gRNA_data, x="hippocampal_layer", y="pre_puncta_density_per_100_um2", extra_y_upper=100, title='Presynaptic puncta density by hippocampal Layer and condition', extra_y_lower=-0.05, ax=axes[3,0])
# plot_post_density = results_plotting.plot_data(df=gRNA_data, x="hippocampal_layer", y="post_puncta_density_per_100_um2", extra_y_upper=100, title='Postsynaptic puncta density by hippocampal Layer and condition', extra_y_lower=-0.05, ax=axes[3,1])
# plot_pre_area = results_plotting.plot_data(df=gRNA_data, x="hippocampal_layer", y="pre_staining_area_um2", extra_y_upper=50, title='Presynaptic staining area by hippocampal Layer and condition', extra_y_lower=-0.05, ax=axes[4,0])
# plot_post_area = results_plotting.plot_data(df=gRNA_data, x="hippocampal_layer", y="post_staining_area_um2", extra_y_upper=50, title='Postsynaptic staining area by hippocampal Layer and condition', extra_y_lower=-0.05, ax=axes[4,1])
# plot_pre_puncta_size = results_plotting.plot_data(df=gRNA_data, x="hippocampal_layer", y="pre_mean_puncta_size_um2", extra_y_upper=0.05, title='Presynaptic puncta size (um2) by hippocampal Layer and condition', extra_y_lower=-0.05, ax=axes[5,0])
# plot_post_puncta_size = results_plotting.plot_data(df=gRNA_data, x="hippocampal_layer", y="post_mean_puncta_size_um2", extra_y_upper=0.05, title='Postsynaptic puncta size (um2) by hippocampal Layer and condition', extra_y_lower=-0.05, ax=axes[5,1])

# # Save the figure explicitly
# combined_plots_figure = plt.gcf()  # Get the current figure object
# plt.tight_layout()  # Adjust layout to prevent overlap
# output_combined_plot_path = os.path.join("combined_plots.png")
# combined_plots_figure.savefig(output_combined_plot_path)  # Save the figure object


### ----------------------- dynamically getting the hippocampal layers and setting the metrics -------------------------------------- ###

# dynamically getting the hippocampal layers into a list from the merged_df & sorting them according to a predefined list
unique_hippocampal_layers = merged_df["hippocampal_layer"].unique().tolist()
desired_order = ["CA1 SO", "CA1 SP", "CA1 SR", "CA1 SLM", 
                 "CA2 SO", "CA2 SP", "CA2 SR", 
                 "CA3 SO", "CA3 SP", "CA3 SL", "CA3 SR",
                 "DG Hilus", "DG GC", "DG ML",
                 "Subiculum PCL", "Cortex L4"]
order_dict = {area: index for index, area in enumerate(desired_order)}
hippocampal_layers = sorted(unique_hippocampal_layers, key=lambda x: order_dict.get(x, float('inf')))

# metrics are hardcoded
metrics = ["overlap_um2", 
           "pearson_cor", 
           "overlap_coeff", 
           "local_peak_colocalized_spots", 
           "presynapse_image_mfi", 
           "postsynapse_image_mfi", 
           "pre_puncta_density_per_100_um2", 
           "post_puncta_density_per_100_um2", 
           "pre_staining_area_um2", 
           "post_staining_area_um2", 
           "pre_mean_puncta_size_um2", 
           "post_mean_puncta_size_um2"]


### ------------------------------------------------ calculate statistics ------------------------------------------------- ###

# getting the protein name from the protein_and_synaptic_marker parser argument
protein_of_interest = protein_and_synaptic_marker.split("_")[0]

# calculate the statistics with a t-test
statistcs_instance = results_plotting.PlotResults(df = merged_df, candidate_gRNA = protein_of_interest + "-gRNA", name_of_plot = "doesnt_matter")
statistics_results = []
for hippocampal_layer in hippocampal_layers:
    for metric in metrics:
        _ , _, table , _ , p_value = statistcs_instance.check_statistics(hippocampal_layer = hippocampal_layer, metric = metric)
        statistics_results.append({"hippocampal_layer": hippocampal_layer, "metric": metric, "p_value": p_value})

# add another column if the results are significant or not
df_statistics_results = pd.DataFrame(statistics_results) 
df_statistics_results['statistical significant?'] = np.where(df_statistics_results["p_value"] <= 0.05, 'Yes', 'No')

# create a new dataframe with only the significant results
p_value_threshold = 0.05  # significant results
df_significant_statistics_results = df_statistics_results[df_statistics_results["p_value"] <= p_value_threshold]

# save the statistics 
df_statistics_results.to_csv("all_statistics.csv")
df_significant_statistics_results.to_csv("significant_statistics.csv")


### --- Plotting the LacZ-gRNA vs the candidate-gRNA (strip plot that shows the means of each brain connected through lines) --- ###

# make the plot
i = results_plotting.PlotResults(df = merged_df, candidate_gRNA = protein_of_interest + "-gRNA", name_of_plot = "combined_plots")
i.save_figure(hippocampal_layer_list = hippocampal_layers, metric_list = metrics)

