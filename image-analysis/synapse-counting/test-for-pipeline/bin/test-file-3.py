import csv
import pandas as pd

# specifying the path to the csv file that needs to be imported
input_file_path = "/mnt/d/code/phd/image-analysis/synapse-counting/output_data/thresholds.csv"

df = pd.read_csv(input_file_path)
        
# calculating the average of the thresholds for vglut1 and psd95
df_threshold_means = df.groupby("image file name").agg({"vglut1_threshold": "mean", "psd95_threshold": "mean"})

# writing to a csv file
output_filename = "mean_threshold_values.csv"
output_path = "/mnt/d/code/phd/image-analysis/synapse-counting/output_data/" + output_filename
df_threshold_means.to_csv(output_path)        

