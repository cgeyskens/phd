import os
import sys
import pandas as pd
import dask
import argparse

# with this peice of code, it will recognize the custom modules
project_root = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
sys.path.append(project_root)

# custom modules
from synapse_counting import metadata, preprocessing, calc_synaptic_metrics

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
def measure_image_puncta(filename, name_segments):
    """
    Measure puncta metrics of both the pre and postsynaptic puncta. Decorated by dask delayed for parrellel processing. 
    """    
    if filename.endswith(".czi"):
        # getting the filepath
        file_path = os.path.join(input_folder, filename)
        # extract metadata for pixel_size parameter in overlap_um2_coloc
        pixel_size_um, _ , image_size_um = metadata.extract_metadata(file_path)
        # preprocessing for puncta_metrics
        pre, post = preprocessing.extract_and_split(file_path, presynapse_channel = presynapse_channel, postsynapse_channel = postsynapse_channel)
        p = preprocessing.ImagePreprocessing(include_rolling_ball=True, include_blur=True, include_clahe=True, include_tophat=True)
        pre_1, post_1 = p.preprocess(pre, post)
        # getting the data
        presynapse_image_mfi, postsynapse_image_mfi = calc_synaptic_metrics.mfi_synapse(pre, post)
        puncta_results, _ , _ = calc_synaptic_metrics.puncta_metrics(pre_1, post_1, image_size_um, pixel_size_um)
        # getting the right filename
        img_filename = metadata.image_filename(filename, name_segments)
        return {
            "img_filename": img_filename,
            "presynapse_image_mfi": presynapse_image_mfi,
            "postsynapse_image_mfi": postsynapse_image_mfi,
            "pre_puncta_density_per_100_um2": puncta_results["pre_puncta_density_per_100_um2"],
            "post_puncta_density_per_100_um2": puncta_results["post_puncta_density_per_100_um2"],
            "pre_staining_area_um2": puncta_results["pre_staining_area_um2"],
            "post_staining_area_um2": puncta_results["post_staining_area_um2"],
            "pre_mean_puncta_size_um2": puncta_results["pre_mean_puncta_size_um2"],
            "post_mean_puncta_size_um2": puncta_results["post_mean_puncta_size_um2"]
        }
    else:
        return None
    
# get a list of files in that input_folder
file_list = os.listdir(input_folder)

# compute the results using a dask delayed object
delayed_results = [measure_image_puncta(filename, name_segments) for filename in file_list]
results = dask.compute(*delayed_results)

# reading out the results into a csv
df = pd.DataFrame(results)
output_csv_path = os.path.join(output_folder, "puncta_results.csv")
df.to_csv(output_csv_path, encoding = "utf-8")