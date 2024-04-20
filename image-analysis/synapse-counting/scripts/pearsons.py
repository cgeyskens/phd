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
parser.add_argument("input_dir", type=str, help="directory to input files")
parser.add_argument("output_dir", type=str, help="directory to output folder")
parser.add_argument('name_segments', type=int, nargs='+', help='Comma separated list of name segments to use seperated by "_"')
args = parser.parse_args()

# assigning the parser arguments
input_folder = args.input_dir
output_folder = args.output_dir
name_segments = args.name_segments
print(name_segments)

# the operation
@dask.delayed
def process_image_pearson(filename, name_segments):
    """
    Processes synaptic colocalization of pre and postsynaptic markers. Decorated by dask delayed for parrellel processing. 
    """    
    if filename.endswith(".czi"):
        # getting the filepath
        file_path = os.path.join(input_folder, filename)
        # preprocessing
        pre, post = preprocessing.extract_and_split(file_path)
        p = preprocessing.ImagePreprocessing(include_rolling_ball=True, include_blur=True, include_clahe=True, include_tophat=True)
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
    else:
        return None
    
# get a list of files in that input_folder
file_list = os.listdir(input_folder)

# compute the results using a dask delayed object
delayed_results = [process_image_pearson(filename, name_segments) for filename in file_list]
results = dask.compute(*delayed_results)

# reading out the results into a csv
df = pd.DataFrame(results)
output_csv_path = os.path.join(output_folder, "pearson_results.csv")
df.to_csv(output_csv_path)
