# test-file-5.py

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os

# specifying the path to the csv file that needs to be imported (the overlap.csv)
input_file_path = "/mnt/d/code/phd/image-analysis/synapse-counting/output_data/overlap.csv"
df = pd.read_csv(input_file_path)

# transforming the df into the right format
df_melted = pd.melt(df, id_vars = ["image file name"], value_vars=["overlap_um2", "overlap_um2_rot"],
                    var_name="condition", value_name="overlap (um2)")
df_melted["condition"] = df_melted["condition"].apply(lambda x: "rotated" if "overlap_um2_rot" in x else "actual")
df_melted['hippocampal layer'] = df_melted['image file name'].apply(lambda x: ' '.join(x.split('_')[-2:]))
df_final = df_melted.drop("image file name", axis=1)

# creating the plot
p = sns.stripplot(y="overlap (um2)", x="hippocampal layer", hue = "condition",
                        data = df_final,
                        jitter = False,
                        dodge=True,
                        marker = "o",
                        alpha = 0.5)

# exporting the plot as png
output_folder = "/mnt/d/code/phd/image-analysis/synapse-counting/output_data/"

output_file = os.path.join(output_folder, "dotplot.png")
plot = p.get_figure()
plt.savefig(output_file, dpi = 300)