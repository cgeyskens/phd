import os 
import sys
import optuna
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
parser.add_argument("--presynapse_channel", type=int, help="number of presynapse channel")
parser.add_argument("--postsynapse_channel", type=int, help="number of postsynapse channel")
parser.add_argument("--protein_and_synaptic_marker", type=str, help="protein that you are interested in")
parser.add_argument("--nr_of_optimization_trials",type=int, help="number of trials requested for optimization of parameters")
parser.add_argument("--preprocessing_params", type=str, required=True, help="path to preprocessing params JSON file")
parser.add_argument("--opt_param_ranges", type=str, required=True, help="path to optimization param ranges JSON file")

#nr of trials
args = parser.parse_args()

# assigning the parser arguments
input_folder = args.input_dir
presynapse_channel = args.presynapse_channel
postsynapse_channel = args.postsynapse_channel
protein_and_synaptic_marker = args.protein_and_synaptic_marker
nr_of_optimization_trials = args.nr_of_optimization_trials

# load the preprocessing params JSON file
with open(args.preprocessing_params, "r" ) as file:
    preprocess_params_dict = json.load(file)
# load the optimization param ranges JSON file
with open(args.opt_param_ranges, "r" ) as file:
    opt_param_ranges_dict = json.load(file)

# create synaptic_marker variable
synaptic_marker_string = protein_and_synaptic_marker.split("_")[-2:]
synaptic_marker = "_".join(synaptic_marker_string)

# getting the hippocampal_layer specific preprocessing parameters from the dictionary
preprocessing_params_synaptic_marker = preprocess_params_dict.get(synaptic_marker, {})
opt_param_ranges_synaptic_marker = opt_param_ranges_dict.get(synaptic_marker, {})


### -------------------------------------- custom definitions for optimization ------------------------------------------ ###

def load_and_preprocess(file_path, presynapse_channel, postsynapse_channel, preprocess_params):
    
    parts = file_path.split('_')[-2:]
    result = "_".join(parts)
    layer = result[:-4]
    
    # getting the hippocampal_layer specific preprocessing parameters from the dictionary
    preprocess_params = preprocess_params_dict.get(synaptic_marker, {}).get(layer, {})
    
    # extract metadata
    pixel_size_um, _ , _ = metadata.extract_metadata(file_path)
    # extract channels
    pre, post = preprocessing.extract_and_split(file_path, 
                                                presynapse_channel = presynapse_channel, 
                                                postsynapse_channel = postsynapse_channel)
    # preprocessing
    p = preprocessing.ImagePreprocessing(
        include_rolling_ball=preprocess_params['include_rolling_ball'], radius=preprocess_params['radius'],
        include_blur=preprocess_params['include_blur'], sigma=preprocess_params['sigma'], preserve_range=True,
        include_clahe=False, #!!! this is different from other processes, otherwise we got no values
        include_tophat=preprocess_params['include_tophat'], element_size=preprocess_params['element_size']
    )
    pre_processed, post_processed = p.preprocess(pre, post)
    return pre_processed, post_processed, pixel_size_um


def objective(trial, presynapse_preprocessed, postsynapse_preprocessed, pixel_size_um, param_ranges):

    pre_distance = trial.suggest_int("pre_distance", *param_ranges['pre_distance'])
    post_distance = trial.suggest_int("post_distance", *param_ranges['post_distance'])
    pre_threshold = trial.suggest_float("pre_threshold", *param_ranges['pre_threshold'])
    post_threshold = trial.suggest_float("post_threshold", *param_ranges['post_threshold'])
    max_distance_um = trial.suggest_float("max_distance_um", *param_ranges['max_distance_um'])

    # probing the first function to get the parameters as input for the second function
    pre_coord, post_coord, post_rot_coord = calc_synaptic_coloc.local_peak_detection(
                                                                presynapse_preprocessed = presynapse_preprocessed,
                                                                postsynapse_preprocessed = postsynapse_preprocessed, 
                                                                presynapse_distance = pre_distance, 
                                                                postsynapse_distance = post_distance, 
                                                                presynapse_threshold = pre_threshold, 
                                                                postsynapse_threshold = post_threshold)
                                                                 
    # probing the second function where psd95 is not rotated, the actual condition
    colocalized_spot_count = calc_synaptic_coloc.count_coloc_spots(pre_coord, post_coord, pixel_size_um, max_distance_um)
    # probing the second function where psd95 is rotated, the internal control condition
    colocalized_spot_count_rot = calc_synaptic_coloc.count_coloc_spots(pre_coord, post_rot_coord, pixel_size_um, max_distance_um)
    
    # getting the scale differences between the spot count of the normal situation and psd95 rotated
    scaled_difference = (colocalized_spot_count - colocalized_spot_count_rot) / max(colocalized_spot_count, 1)
    
    return scaled_difference


def optimize_parameters_for_hippocampal_layer(file_list, 
                                              input_folder, 
                                              hippocampal_layers, 
                                              preprocess_params_by_hippocampal_layer, 
                                              param_ranges, 
                                              nr_of_trials,
                                              plot_coord = False):
    """
    This function will optimize the parameters for colocalization based on the objective function.

    It uses Tree-structured Parzen Estimator (TPE) from the optuna library to optimize the parameters.
    For each batch of presynapse and postsynapse images per hippocampal layer, it will try to optimize the 
    parameters such that the scaled difference between actual and rotated presynapse/postsynapses images is 
    maximized. It also uses dask such that it parallel process each batch of images per hippocampal layer.
    
    Args:
        file_list: a list of files (images) that will be processed.
        input_folder: the folder where the images are stored.
        hippocampal_layers: a list of hippocampal layers of which the images needs to be optimized.
        preprocess_params_by_hippocampal_layer: a list of unique handcrafted preprocessing parameters for each hippocampal layer
        param_ranges: a list of unique handcrafted parameters ranges for each hippocampal layer
        nr_of_trials: the total nr of trials to optimize the parameters

    Returns:
        best_params_by_hippocampal_layer_df: dataframe containing the best parameters per hippocampal layer
        final_df_optuna: dataframe containing the info per trial.
        final_df: dataframe containing all the optimization information, including the the scores for each image combination
    """
    
    @dask.delayed
    def process_hippocampal_layer(hippocampal_layer):
        images = [f for f in file_list if hippocampal_layer in f and "LacZ-gRNA" in f] # only taking control (LacZ-images) images for setting the parameters
    
        preprocess_params = preprocess_params_by_hippocampal_layer[hippocampal_layer]
        
        # extract metadata and preprocess images
        image_results = []
        for file_name in images:
            file_path = os.path.join(input_folder, file_name)
            pre_1, post_1, pixel_size_um = load_and_preprocess(file_path, presynapse_channel=0, postsynapse_channel=1, preprocess_params=preprocess_params)
            if pre_1 is not None and post_1 is not None:
                image_results.append((pre_1, post_1, pixel_size_um, file_path))
        
        # optuna optimization
        def hippocampal_layer_objective(trial):
            detailed_trial_results = []
            for pre_1, post_1, pixel_size_um, file_path in image_results:
                score = objective(trial, pre_1, post_1, pixel_size_um, param_ranges[hippocampal_layer])
                detailed_trial_results.append({
                    'hippocampal_layer': hippocampal_layer,
                    'trial': trial.number,
                    'image': file_path,
                    'params': trial.params,
                    'score': score
                })
            mean_score = np.mean([result['score'] for result in detailed_trial_results])
            return mean_score, detailed_trial_results
        
        # create the study
        study = optuna.create_study(study_name=hippocampal_layer, direction="maximize", sampler=optuna.samplers.TPESampler())
        detailed_trials_data = []

        # define the optuna objective, such that the mean score of all the images in a certain hippocampal layer is optimized
        def optuna_objective(trial):
            mean_score, detailed_trial_results = hippocampal_layer_objective(trial)
            detailed_trials_data.extend(detailed_trial_results)
            return mean_score
        
        # do the optimization
        study.optimize(optuna_objective, n_trials=nr_of_trials, show_progress_bar = True)
    
        # store best parameters for the hippocampal_layer with the corresponding score
        best_trial_full = study.best_trial
        best_params = {
            "best_trial": best_trial_full.number,
            "best_params": study.best_params,
            "best_score": study.best_value
        }
        
        # create dataframe for detailed trial results and the best params
        trials_df = pd.DataFrame(detailed_trials_data)
        
        # getting the data per trial
        trials_data = [{'hippocampal_layer': hippocampal_layer} | trial for trial in study.trials_dataframe().to_dict('records')]

        return hippocampal_layer, best_params, trials_df, trials_data

    # Collect the delayed tasks
    delayed_tasks = [process_hippocampal_layer(hippocampal_layer) for hippocampal_layer in hippocampal_layers]

    # Compute the results in parallel
    results = dask.compute(*delayed_tasks)

    best_params_by_hippocampal_layer = {}
    results_dfs = []
    all_trials_data = []

    for hippocampal_layer, best_params, trials_df, trials_data in results:
        best_params_by_hippocampal_layer[hippocampal_layer] = best_params
        results_dfs.append(trials_df)
        all_trials_data.extend(trials_data)

    # create dataframe for best parameters by hippocampal_layer
    best_params_by_hippocampal_layer_df = pd.DataFrame.from_dict(best_params_by_hippocampal_layer, orient="index").reset_index()
    best_params_by_hippocampal_layer_df = pd.concat([best_params_by_hippocampal_layer_df.drop(['best_params'], axis=1), 
                                          pd.json_normalize(best_params_by_hippocampal_layer_df['best_params'])], axis=1)

    # combine dataframes for all hippocampal_layers
    final_df = pd.concat(results_dfs, ignore_index=True)

    # combine the trial data
    final_df_optuna = pd.DataFrame(all_trials_data)

    return best_params_by_hippocampal_layer_df, final_df_optuna, final_df


### ------------------------------------------ The actual optimization ------------------------------------------- ###

# additional parameters for master definition
file_list = os.listdir(input_folder)

# getting the hippocampal layers dynamically
def get_regions(data, key):
    if key in data:
        return list(data[key].keys())
    else:
        return []
    
if synaptic_marker == "VGLUT1_PSD95":
    hippocampal_layers = get_regions(opt_param_ranges_dict, "VGLUT1_PSD95")
if synaptic_marker == "VGLUT2_PSD95":
    hippocampal_layers = get_regions(opt_param_ranges_dict, "VGLUT2_PSD95")
if synaptic_marker == "VGAT_GPHN":
    hippocampal_layers = get_regions(opt_param_ranges_dict, "VGAT_GPHN")

best_params_df, trial_df, all_df = optimize_parameters_for_hippocampal_layer(
    file_list = file_list, 
    input_folder = input_folder, 
    hippocampal_layers = hippocampal_layers, 
    preprocess_params_by_hippocampal_layer = preprocessing_params_synaptic_marker, 
    param_ranges = opt_param_ranges_synaptic_marker, 
    nr_of_trials = nr_of_optimization_trials,
    plot_coord = False
    )

output_best_params = "local_peak_best_optimization_params.csv"
best_params_df.to_csv(output_best_params, encoding = "utf-8")

output_trials_data = "local_peak_optimization_trial_data.csv"
trial_df.to_csv(output_trials_data, encoding = "utf-8")

output_all_data = "local_peak_optimization_all_data.csv"
all_df.to_csv(output_all_data, encoding = "utf-8")




