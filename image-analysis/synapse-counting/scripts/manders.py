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
parser.add_argument("--output_dir", type=str, help="directory to output folder")
parser.add_argument('--name_segments', type=int, nargs='+', help='Comma separated list of name segments to use seperated by "_"')
parser.add_argument("--presynapse_channel", type=int, help="number of presynapse channel")
parser.add_argument("--postsynapse_channel", type=int, help="number of postsynapse channel")
args = parser.parse_args()

# assigning the parser arguments
input_folder = args.input_dir
output_folder = args.output_dir
name_segments = args.name_segments
presynapse_channel = args.presynapse_channel
postsynapse_channel = args.postsynapse_channel

# the operation
@dask.delayed
def measure_image_manders(filename, name_segments):
    """
    Measure colocalization of pre and postsynaptic markers using Manders. Decorated by dask delayed for parrellel processing. 
    """    
    if filename.endswith(".czi"):
        # getting the filepath
        file_path = os.path.join(input_folder, filename)
        # preprocessing
        pre, post = preprocessing.extract_and_split(file_path, presynapse_channel = presynapse_channel, postsynapse_channel = postsynapse_channel)
        p = preprocessing.ImagePreprocessing(include_rolling_ball=True, include_blur=True, include_clahe=True, include_tophat=True)
        pre_1, post_1 = p.preprocess(pre, post)
        # getting the data
        overlap_coeff, overlap_coeff_rot, presynapse_threshold, postsynapse_threshold = calc_synaptic_coloc.manders_coloc(pre_1, post_1)
        # getting the right filename
        img_filename = metadata.image_filename(filename, name_segments) 
        return {
            "img_filename": img_filename,
            "overlap_coeff": overlap_coeff,
            "overlap_coeff_rot": overlap_coeff_rot,
            "presynapse_threshold": presynapse_threshold,
            "postsynapse_threshold": postsynapse_threshold
        }
    else:
        raise Warning("not all files inside {input_folder} are .czi files")
    
# get a list of files in that input_folder
file_list = os.listdir(input_folder)

# compute the results using a dask delayed object
delayed_results = [measure_image_manders(filename, name_segments) for filename in file_list]
results = dask.compute(*delayed_results)

# reading out the results into a csv
df = pd.DataFrame(results)
output_csv_path = os.path.join(output_folder, "manders_results.csv")
df.to_csv(output_csv_path, encoding = "utf-8")