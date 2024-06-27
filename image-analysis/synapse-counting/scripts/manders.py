import os
import sys
import pandas as pd
import dask
import argparse

# with this peice of code, it will recognize the custom modules
project_root = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
sys.path.append(project_root)

# custom modules
from synapse_counting import metadata, preprocessing, calc_synaptic_coloc

# adding parser arguments
parser = argparse.ArgumentParser(description="process input files")
parser.add_argument("--input_dir", type=str, help="directory to input files")
parser.add_argument('--name_segments', type=int, nargs='+', help='Comma separated list of name segments to use seperated by "_"')
parser.add_argument("--presynapse_channel", type=int, help="number of presynapse channel")
parser.add_argument("--postsynapse_channel", type=int, help="number of postsynapse channel")
parser.add_argument("--protein_and_synaptic_marker", type=str, help="protein that you are interested in")
args = parser.parse_args()

# assigning the parser arguments
input_folder = args.input_dir
name_segments = args.name_segments
presynapse_channel = args.presynapse_channel
postsynapse_channel = args.postsynapse_channel
protein_and_synaptic_marker = args.protein_and_synaptic_marker

# create synaptic_marker variable
synaptic_marker_string = protein_and_synaptic_marker.split("_")[-2:]
synaptic_marker = "_".join(synaptic_marker_string)

# the operation
@dask.delayed
def measure_image_manders(filename, name_segments):
    """
    Measure colocalization of pre and postsynaptic markers using Manders. Decorated by dask delayed for parrellel processing. 
    """    
    file_path = os.path.join(input_folder, filename)
    parts = filename.split('_')[-2:]
    result = "_".join(parts)
    
    # for each synaptic_marker combination and layer, I handcrafted the preprocessing parameters for best segmentation
    if synaptic_marker == "VGLUT1_PSD95":
        if result in {"CA3_SL.czi", "DG_Hilus.czi"}:
            radius, element_size, blur_sigma, watershed_sigma = 15, 20, 2, 3
        elif result in {"CA1_SO.czi", "CA3_SO.czi"}:
            radius, element_size, blur_sigma, watershed_sigma = 5, 10, 2, 3
        elif result in {"CA1_SR.czi", "CA3_SR.czi"}:
            radius, element_size, blur_sigma, watershed_sigma = 10, 10, 2, 3
        elif result == "CA1_SLM.czi":
            radius, element_size, blur_sigma, watershed_sigma = 5, 5, 2, 3
        elif result == "DG_ML.czi":
            radius, element_size, blur_sigma, watershed_sigma = 5, 15, 2, 3,
    elif synaptic_marker == "VGLUT2_PSD95":
        if result in {"Cortex_L4.czi"}:
            radius, element_size, blur_sigma, watershed_sigma= 10, 10, 2, 3
        elif result in {"CA2_SP.czi", "DG_GC.czi", "Subiculum_SP.czi"}:
            radius, element_size, blur_sigma, watershed_sigma = 15, 10, 2, 3

    # extracting and splitting channels
    pre, post = preprocessing.extract_and_split(file_path, presynapse_channel = presynapse_channel, postsynapse_channel = postsynapse_channel)
    
    # preprocessing
    p = preprocessing.ImagePreprocessing(
        include_rolling_ball=True, radius=radius,
        include_blur=True, sigma = blur_sigma, preserve_range = True,
        include_clahe=True,
        include_tophat=True, element_size = element_size
        )
    pre_1, post_1 = p.preprocess(pre, post)

    # thresholding and watershed segmentation
    presynapse_threshold, postsynapse_threshold = preprocessing.thresholding(pre_1, post_1, threshold_algorithm="triangle")
    presynapse_watersheded = preprocessing.custom_watershed(presynapse_threshold, sigma = watershed_sigma)
    postsynapse_watersheded = preprocessing.custom_watershed(postsynapse_threshold, sigma = watershed_sigma)
    
    # getting the data
    overlap_coeff, overlap_coeff_rot = calc_synaptic_coloc.manders_coloc(presynapse_watersheded, postsynapse_watersheded)
    
    # getting the right filename
    img_filename = metadata.image_filename(filename, name_segments)
    
    return {
        "img_filename": img_filename,
        "overlap_coeff": overlap_coeff,
        "overlap_coeff_rot": overlap_coeff_rot
    }
    
# get a list of files in that input_folder
file_list = os.listdir(input_folder)

# compute the results using a dask delayed object
delayed_results = [measure_image_manders(filename, name_segments) for filename in file_list]
results = dask.compute(*delayed_results)

# reading out the results into a csv
df = pd.DataFrame(results)
output_csv_path = "manders_results.csv"
df.to_csv(output_csv_path, encoding = "utf-8")