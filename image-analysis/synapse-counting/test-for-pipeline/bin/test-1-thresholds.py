# test-1-thresholds.py
# this script is based on test-file-2b.ipynb, but here we will only output following parameters:
# overlap (um2), overlap (um2) in rot cond, vglut1_threshold, psd95_threshold

from pathlib import Path
import czifile
import pandas as pd
import os
from microfilm.microplot import microshow
import numpy as np
import argparse
import matplotlib.pyplot as plt
from skimage.restoration import rolling_ball
from skimage import measure
from skimage.transform import rotate
from skimage import filters
from skimage.filters import gaussian
from lxml import etree # required library to load the metadata
import csv


# adding parser arguments
parser = argparse.ArgumentParser(description="process input files")
parser.add_argument("input_dir", type=str, help="directory to input files")
parser.add_argument("output_dir", type=str, help="directory to output folder")

args = parser.parse_args()

# assigning the parser arguments
input_folder = args.input_dir
output_folder = args.output_dir




# defing a function for the overlap between two synaptic markers
def colocalization_overlap(filename):
    
    # reading the file
    image = czifile.imread(filename)
    
    # this retrieves the metadata
    czi = czifile.CziFile(filename)
    czi_xml_str = czi.metadata() # gets the metadata in a xml string format
    czi_parsed = etree.fromstring(czi_xml_str) # parses the czi_xml_str file

    # finds the strings 
    size_x = czi_parsed.find(".//SizeX")
    size_y = czi_parsed.find(".//SizeY")
    scaling_x = czi_parsed.find(".//ScalingX")
    scaling_y = czi_parsed.find(".//ScalingY")

    # Extracting the required values to calculate the 
    size_x_value = int(size_x.text)
    size_y_value = int(size_y.text)
    scaling_x_value = float(scaling_x.text)
    scaling_y_value = float(scaling_y.text)

    # calculater the pixel to micrometer (um)
    # first checking whether the X and Y dimensions of the image are equal
    if size_x_value == size_y_value and scaling_x_value == scaling_y_value:
        pixel_size = ((scaling_x_value*1000000000)/size_x_value) # conversion from meter to micrometer  
    
    
    # getting rid of all the extra channels and splitting the channels into vglut1 and psd95
    image_squeezed = np.squeeze(image)
    image_squeezed.shape
    vglut1 = image_squeezed[0,:,:]
    psd95 = image_squeezed[1,:,:]
    
    # preproccesing step, removal of background with rolling ball radius of 10x
    background_vglut1 = rolling_ball(vglut1, radius = 10)
    background_psd95 = rolling_ball(psd95, radius = 10)
    vglut1_bs = vglut1 - background_vglut1
    psd95_bs = psd95 - background_psd95
    
    # for the thresholding, it is advised to perform a gaussian blur (in this case a light one with sigma 1)
    psd95_pre = gaussian(psd95_bs, sigma=1, preserve_range=True)
    vglut1_pre = gaussian(vglut1_bs, sigma=1, preserve_range=True)
    
    
    #########################################
    ### vglut1 and psd95 overlap - in um2 ###
    #########################################

    # rotating the psd95 image as control
    psd95_pre_rot = rotate(psd95_pre, 90)

    # getting the threshold (triangle) of each image
    vglut1_threshold = filters.threshold_triangle(vglut1_pre)
    psd95_threshold = filters.threshold_triangle(psd95_pre)
    psd95_threshold_rot = filters.threshold_triangle(psd95_pre_rot)

    # applying the threhold
    vglut1_pre_thr = vglut1_pre >= vglut1_threshold
    psd95_pre_thr = psd95_pre >= psd95_threshold
    psd95_pre_rot_thr = psd95_pre_rot >= psd95_threshold_rot
    
    # getting the overlap and nr of pixel that overlapped
    overlap = vglut1_pre_thr & psd95_pre_thr
    overlap_rot = vglut1_pre_thr & psd95_pre_rot_thr

    overlap_pix = np.sum(overlap)
    overlap_pix_rot = np.sum(overlap_rot)

    overlap_um2 = overlap_pix * pixel_size * pixel_size
    overlap_um2_rot = overlap_pix_rot * pixel_size * pixel_size
    
    return overlap_um2, overlap_um2_rot, vglut1_threshold, psd95_threshold
    

# defining a function to extract the right elements from the filename
def image_filename(filename):
    
    # get the filename
    name_of_file = os.path.splitext(os.path.basename(filename))[0]
    split_filename = name_of_file.split("_")

    # get only the experimental parameters from the filename
    index_nums = [0, 1, 2, 3, 10, 11] # the indexes of the elements that I would like to extract from
    desired_parts = [split_filename[val] for val in index_nums]
    desired_filename = "_".join(desired_parts)
    
    return desired_filename


# defining a function that writes results to a csv file
def write_to_csv(output_file, data_list):
    with open(output_file, mode="w", newline="") as csv_file:
        fieldnames = ["image file name", "overlap_um2", "overlap_um2_rot", "vglut1_threshold", "psd95_threshold"]
        writer = csv.DictWriter(csv_file, fieldnames=fieldnames)
        writer.writeheader()

        for data in data_list:
            writer.writerow(data)
  

# get a list of files in that input_folder
file_list = os.listdir(input_folder)

# empty dict to store the results
results = []

# the actual function
for filename in file_list:
    if filename.endswith(".czi"):
        file_path = os.path.join(input_folder, filename)
        overlap_um2, overlap_um2_rot, vglut1_threshold, psd95_threshold  = colocalization_overlap(file_path)
        img_filename = image_filename(filename)  # Extract desired filename parts

        results.append({
            "image file name": img_filename,
            "overlap_um2": overlap_um2,
            "overlap_um2_rot": overlap_um2_rot,
            "vglut1_threshold": vglut1_threshold,
            "psd95_threshold": psd95_threshold
        })

output_csv_path = os.path.join(output_folder, "thresholds.csv")
write_to_csv(output_csv_path, results)