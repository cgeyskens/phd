# test_file_1

# importing the required packages
from pathlib import Path
import czifile
import pandas as pd
import os

def test_file_1(image_path):
    
    # getting the dimension of the image
    image = czifile.imread(path)
    data = image.shape

    # get the filename
    filename = os.path.splitext(os.path.basename(path))[0]
    split_filename = filename.split("_")

    # get only the experimental parameters from the filename
    index_nums = [0, 1, 2, 3, 10, 11] # the indexes of the elements that I would like to extract
    desired_parts = [split_filename[val] for val in index_nums]
    desired_filename = "_".join(desired_parts)

    # converting tuple to df
    df = pd.DataFrame(data, columns = [desired_filename] )

    # pivoting the dataframe
    dic = df.to_dict(orient='dict')
    df_final = pd.DataFrame.from_dict(dic, orient='index')

    return df_final


