import os 
import sys
import dask
import pandas as pd
import numpy as np
import argparse
import json


### --------------------------------- parser arguments and loading custom library --------------------------------- ###

# with this peice of code, it will recognize the custom modules
project_root = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
sys.path.append(project_root)

# custom modules
from synapse_counting import metadata, preprocessing, calc_synaptic_coloc

# adding parser arguments
parser = argparse.ArgumentParser(description = "process input files")
parser.add_argument("--input_dir", type = str, help = "directory to input files")
parser.add_argument("--name_segments", type = int, nargs = '+', help = 'Comma separated list of name segments to use seperated by "_"')
parser.add_argument("--presynapse_channel", type=int, help="number of presynapse channel")
parser.add_argument("--postsynapse_channel", type=int, help="number of postsynapse channel")
parser.add_argument("--protein_and_synaptic_marker", type=str, help="protein that you are interested in")
parser.add_argument("--preprocessing_params", type=str, required=True, help="path to preprocessing params JSON file")
parser.add_argument("--optimized_params_df",type=str, help="dataframe with optimized parameters for each hippcampal layer")

#nr of trials
args = parser.parse_args()

# assigning the parser arguments
input_folder = args.input_dir
name_segments = args.name_segments
presynapse_channel = args.presynapse_channel
postsynapse_channel = args.postsynapse_channel
protein_and_synaptic_marker = args.protein_and_synaptic_marker

# optimized_params_df = args.optimized_params_df
optimized_params_df = pd.read_csv(args.optimized_params_df)
# set the index of the best_params_df
optimized_params_df.set_index('index', inplace=True)

# load the preprocessing params JSON file
with open(args.preprocessing_params, "r" ) as file:
    preprocess_params_dict = json.load(file)

# create synaptic_marker variable
synaptic_marker_string = protein_and_synaptic_marker.split("_")[-2:]
synaptic_marker = "_".join(synaptic_marker_string)


### --------------------------------- the calculation of colocalized spots --------------------------------- ###

@dask.delayed
def local_peak_colocalizaion(filename, input_folder, preprocessing_params, synaptic_marker, optimzed_coloc_params_df, name_segments):
    """
    Measure colocalized local peaks with optimized parameters. Decorated by dask delayed for parrellel processing. 

    Args:
        filename: filename of the image
        inputfolder: the inputfolder with all the images
        preprocessing_params: dictionary with all the preprocessing params that were handcrafted
        synaptic_marker: combination of synaptic markers used
        optimzed_coloc_params_df: dataframe of the optimized params for colocalization
        name_segments: the import segments of the filename

    Returns:
        A dictionary containing the filename, the local peak colocalized spot and the 
        local peak colocalized spots when rotated.
    """   
    file_path = os.path.join(input_folder, filename)
    parts = file_path.split('_')[-2:]
    result = "_".join(parts)
    layer = result[:-4]
    
    # getting the hippocampal_layer specific preprocessing parameters from the dictionary
    params = preprocessing_params.get(synaptic_marker, {}).get(layer, {})

    # ectracting the best params for the handcrafted preprocessing params and the optimized colocalization params
    def extrac_all_params(layer):
        return {
            "radius": params.get("radius"),
            "element_size": params.get("element_size"),
            "sigma":  params.get("sigma"),
            "presynapse_distance": optimzed_coloc_params_df.at[layer, "pre_distance"],
            "postsynapse_distance": optimzed_coloc_params_df.at[layer, "post_distance"],
            "presynapse_threshold": optimzed_coloc_params_df.at[layer, "pre_threshold"],
            "postsynapse_threshold": optimzed_coloc_params_df.at[layer, "post_threshold"],
            "max_distance": optimzed_coloc_params_df.at[layer, "max_distance_um"]
            }
    all_param_values = extrac_all_params(layer)

    # extract metadata for pixel_size parameter in overlap_um2_coloc
    pixel_size_um, _ , _ = metadata.extract_metadata(file_path)

    # extracting and splitting channels
    pre, post = preprocessing.extract_and_split(file_path, 
                                                presynapse_channel = 0, 
                                                postsynapse_channel = 1
                                                )
    # preprocessing
    p = preprocessing.ImagePreprocessing(
        include_rolling_ball=True, radius=all_param_values["radius"],
        include_blur=True, sigma = all_param_values["sigma"], preserve_range = True,
        include_clahe=False, #!!! this is different from other processes, otherwise we got no values
        include_tophat=True, element_size = all_param_values["element_size"]
        )
    pre_1, post_1 = p.preprocess(pre, post)
    # calculating the coordinates
    presynapse_coord, postsynapse_coord, postsynapse_rot_coord = calc_synaptic_coloc.local_peak_detection(
        presynapse_preprocessed=pre_1,
        postsynapse_preprocessed=post_1,
        presynapse_distance = all_param_values["presynapse_distance"],
        postsynapse_distance= all_param_values["postsynapse_distance"],
        presynapse_threshold= all_param_values["presynapse_threshold"],
        postsynapse_threshold= all_param_values["postsynapse_threshold"],
        plot_coord = False
        )                                                                                                  
    # actual colocalized synapses
    colocalized_spots = calc_synaptic_coloc.count_coloc_spots(
        presynapse_coordinates = presynapse_coord,
        postsynapse_coordinates = postsynapse_coord,
        pixel_size_um = pixel_size_um,
        max_distance_um = all_param_values["max_distance"]
        )
    # rotated colocalized synapses as control
    colocalized_spots_rot = calc_synaptic_coloc.count_coloc_spots(
        presynapse_coordinates = presynapse_coord,
        postsynapse_coordinates = postsynapse_rot_coord,
        pixel_size_um = pixel_size_um,
        max_distance_um = all_param_values["max_distance"]
        )
    # getting the right filename
    img_filename = metadata.image_filename(filename, name_segments)
    
    return {
        "img_filename": img_filename,
        "local_peak_colocalized_spots": colocalized_spots,
        "local_peak_colocalized_spots_rot": colocalized_spots_rot,
    }

# get a list of files in that input_folder
file_list = os.listdir(input_folder)

# compute the results using a dask delayed object
delayed_results = [local_peak_colocalizaion(filename = file_name, 
                                            input_folder = input_folder,
                                            preprocessing_params = preprocess_params_dict,
                                            synaptic_marker = synaptic_marker,
                                            optimzed_coloc_params_df = optimized_params_df,
                                            name_segments = name_segments) for file_name in file_list]
results = dask.compute(*delayed_results)

# reading out the results into a csv
df = pd.DataFrame(results)
output_csv_path = "local_peak_coloc.csv"
df.to_csv(output_csv_path, encoding = "utf-8")