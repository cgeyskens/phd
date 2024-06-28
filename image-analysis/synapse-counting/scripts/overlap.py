import os
import sys
import pandas as pd
import dask
import argparse
import json

# with this peice of code, it will recognize the custom modules
project_root = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
sys.path.append(project_root)

# custom modules
from synapse_counting import metadata, preprocessing, calc_synaptic_coloc

# adding parser arguments
parser = argparse.ArgumentParser(description = "process input files")
parser.add_argument("--input_dir", type = str, help = "directory to input files")
parser.add_argument('--name_segments', type = int, nargs = '+', help = 'Comma separated list of name segments to use seperated by "_"')
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

# create synaptic_marker variable
synaptic_marker_string = protein_and_synaptic_marker.split("_")[-2:]
synaptic_marker = "_".join(synaptic_marker_string)

# the operation
@dask.delayed
def measure_image_overlap(filename, name_segments):
    """
    Measure colocalization of pre and postsynaptic markers using overlap. Decorated by dask delayed for parrellel processing. 
    """    
    file_path = os.path.join(input_folder, filename)
    parts = filename.split('_')[-2:]
    result = "_".join(parts)
    layer = result[:-4]
    
    # getting the hippocampal_layer specific preprocessing parameters from the dictionary
    params = preprocess_params_dict.get(synaptic_marker, {}).get(layer, {})

    # extract metadata for pixel_size parameter in overlap_um2_coloc
    pixel_size_um, _ , _ = metadata.extract_metadata(file_path)

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
    
    # thresholding and watershed segmentation
    presynapse_threshold, postsynapse_threshold = preprocessing.thresholding(
        pre_1, 
        post_1, 
        threshold_algorithm=params["threshold_algorithm"])
    presynapse_watersheded = preprocessing.custom_watershed(
        presynapse_threshold, 
        sigma = params["watershed_sigma"])
    postsynapse_watersheded = preprocessing.custom_watershed(
        postsynapse_threshold, 
        sigma = params["watershed_sigma"])
   
    # getting the data
    overlap_um2, overlap_um2_rot = calc_synaptic_coloc.overlap_um2_coloc(presynapse_watersheded, postsynapse_watersheded, pixel_size_um)
    
    # getting the right filename
    img_filename = metadata.image_filename(filename, name_segments)
    
    return {
        "img_filename": img_filename,
        "overlap_um2": overlap_um2,
        "overlap_um2_rot": overlap_um2_rot,
    }

# get a list of files in that input_folder
file_list = os.listdir(input_folder)

# compute the results using a dask delayed object
delayed_results = [measure_image_overlap(filename, name_segments) for filename in file_list]
results = dask.compute(*delayed_results)

# reading out the results into a csv
df = pd.DataFrame(results)
df.to_csv("overlap_results.csv", encoding = "utf-8")