#!/usr/bin/env python
import os
import sys
import pandas as pd
import dask
import argparse
import json

# # Get the parent directory of scripts (parent_directory)
# parent_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))
# # Get the synapse_counting directory
# synapse_counting_dir = os.path.join(parent_dir, 'synapse_counting')
# # Append the synapse_counting directory to sys.path
# sys.path.append(synapse_counting_dir)

from synapse_counting import metadata, preprocessing, calc_synaptic_coloc, helpers

# adding parser arguments
parser = argparse.ArgumentParser(description="process input files")
parser.add_argument("--input_dir", type=str, help="directory to input files")
parser.add_argument('--name_segments', type=int, nargs='+', help='Comma separated list of name segments to use seperated by "_"')
parser.add_argument("--presynapse_channel", type=int, help="number of presynapse channel")
parser.add_argument("--postsynapse_channel", type=int, help="number of postsynapse channel")
parser.add_argument("--protein_and_synaptic_marker", type=str, help="protein that you are interested in")
parser.add_argument("--preprocessing_params", type=str, required=True, help="path to preprocessing params JSON file")

args = parser.parse_args()

# assigning the parser arguments
input_folder = args.input_dir
name_segments = args.name_segments
presynapse_channel = args.presynapse_channel
postsynapse_channel = args.postsynapse_channel
protein_and_synaptic_marker = args.protein_and_synaptic_marker

# load the preprocessing params JSON file
with open(args.preprocessing_params, "r" ) as file:
    preprocess_params_dict = json.load(file)

# get the synaptic marker
synaptic_marker = helpers.get_synaptic_marker(protein_and_synaptic_marker)

# the operation
@dask.delayed
def measure_image_pearson(filename, 
                          input_folder,
                          preprocessing_params,
                          synaptic_marker,
                          name_segments):
    """
    Measure colocalization of pre and postsynaptic markers using Pearsons. Decorated by dask delayed for parrellel processing. 
    """    
    # make file_path and extract hippocampal layer
    file_path = os.path.join(input_folder, filename)
    layer = helpers.get_hippocampal_layer(file_path)
    
    # getting the hippocampal_layer specific preprocessing parameters from the dictionary
    params = preprocessing_params.get(synaptic_marker, {}).get(layer, {})
    
    # extracting and splitting channels
    pre, post = preprocessing.extract_and_split(file_path, presynapse_channel = presynapse_channel, postsynapse_channel = postsynapse_channel)
    # preprocessing 
    p = preprocessing.ImagePreprocessing(
        include_rolling_ball=params["include_rolling_ball"], radius=params["radius"],
        include_blur=params["include_blur"], sigma = params["sigma"], preserve_range = True,
        include_clahe=params["include_clahe"],
        include_tophat=params["include_tophat"], element_size = params["element_size"]
        )
    pre_1, post_1 = p.preprocess(pre, post)
    # getting the data
    pearson_cor, pvalue, pearson_cor_rot, pvalue_rot = calc_synaptic_coloc.pearsons_coloc(pre_1, post_1)
    # getting the right filename
    img_filename = metadata.image_filename(filename, name_segments) 
    return {
        "img_filename": img_filename,
        "pearson_cor": pearson_cor,
        "pvalue": pvalue,
        "pearson_cor_rot": pearson_cor_rot,
        "pvalue_rot": pvalue_rot
    }
    
# get a list of files in that input_folder
file_list = os.listdir(input_folder)

# compute the results using a dask delayed object
delayed_results = [measure_image_pearson(filename = filename, 
                          input_folder = input_folder,
                          preprocessing_params = preprocess_params_dict,
                          synaptic_marker = synaptic_marker,
                          name_segments = name_segments) for filename in file_list]
results = dask.compute(*delayed_results)

# reading out the results into a csv
df = pd.DataFrame(results)
df.to_csv("pearson_results.csv", encoding = "utf-8")
