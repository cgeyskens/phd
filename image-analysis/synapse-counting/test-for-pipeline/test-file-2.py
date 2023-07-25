# test-file-2.py

# importing the required packages
from pathlib import Path
import czifile
import pandas as pd
import sys

def process_image_and_save_csv(input_path, output_path):
    try:
        # getting the dimension of the image
        image = czifile.imread(input_path)
        data = image.shape

        # get the filename
        filename = Path(input_path).stem
        split_filename = filename.split("_")

        # get only the experimental parameters from the filename
        index_nums = [0, 1, 2, 3, 10, 11]  # the indexes of the elements that I would like to extract
        desired_parts = [split_filename[val] for val in index_nums]
        desired_filename = "_".join(desired_parts)

        # converting tuple to df
        df = pd.DataFrame(data, columns=[desired_filename])

        # Writing the dataframe in a csv format
        df.to_csv(output_path, index=False)

        print(f"Processed image: {input_path}")
        sys.exit(0)  # Exit with success status code
    except Exception as e:
        print(f"Error processing image: {input_path}\nError: {e}")
        sys.exit(1)  # Exit with error status code

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python test-file-1.py input_image_path output_csv_path")
        sys.exit(1)

    input_image_path = sys.argv[1]
    output_csv_path = sys.argv[2]

    process_image_and_save_csv(input_image_path, output_csv_path)
