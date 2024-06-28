python3 local_peaks_coloc.py --input_dir /Volumes/KINGSTON/data/phd/image-analysis/synapse-counting/VCAM1/VCAM1-LacZ_VGLUT1-PSD95_images  --name_segments 0 1 2 3 6 8 9  --presynapse_channel 0  --postsynapse_channel 1  --protein_and_synaptic_marker VCAM1_LacZ_VGLUT1_PSD95  --preprocessing_params preprocess_params.json  --optimized_params_df  local_peak_best_optimization_params.csv

parser = argparse.ArgumentParser(description = "process input files")
parser.add_argument("--input_dir", type = str, help = "directory to input files")
parser.add_argument("--name_segments", type = int, nargs = '+', help = 'Comma separated list of name segments to use seperated by "_"')
parser.add_argument("--presynapse_channel", type=int, help="number of presynapse channel")
parser.add_argument("--postsynapse_channel", type=int, help="number of postsynapse channel")
parser.add_argument("--protein_and_synaptic_marker", type=str, help="protein that you are interested in")
parser.add_argument("--preprocessing_params", type=str, required=True, help="path to preprocessing params JSON file")
parser.add_argument("--optimized_params_df",type=str, help="dataframe with optimized parameters for each hippcampal layer")
