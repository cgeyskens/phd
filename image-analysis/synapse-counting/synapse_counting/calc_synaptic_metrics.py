import numpy as np 
from skimage import measure, filters
import pandas as pd


def mfi_synapse(presynapse_image, postsynapse_image):
    """
    Returns the mean fluorescence intensity of pre and post synapse channel.

    Args:
        presynapse_image: image of presynapse
        postsynapse_image: image of postsynapse

    Returns:
        presynapse_image_mfi: mean fluorescence of presynapse
        postsynapse_image_mfi: mean fluorescence of postsynapse

    """

    presynapse_image_mfi = np.mean(presynapse_image)
    postsynapse_image_mfi = np.mean(postsynapse_image)
    
    return presynapse_image_mfi, postsynapse_image_mfi


def puncta_metrics(presynapse_image, postsynapse_image, image_size_um, pixel_size_um, threshold_algorithm = "otsu"):
    """
    Returns different puncta metrics.

    Args:
        presynapse_image: image of presynapse
        postsynapse_image: image of postsynapse
        image_size_um: the size of the image in um.
        pixel_size_um: the size of each pixel in um
        threshold_algorithm: threshold algorithm to choose (default is otsu, can also choose isodata and triangle)
    
    Returns:
        dictionary with different puncta metrics:
                pre_puncta_nr: presynapse puncta number of the image
                pre_puncta_density_per_100_um2: presynapse puncta density per 100 um2
                post_puncta_nr: post synapse puncta number
                post_puncta_density_per_100_um2: post synapse puncta density per 100 um2
                pre_staining_area_pix: presynapse staining area in pixels
                pre_staining_area_um2: presynapse staining area in um2
                post_staining_area_pix: postsynapse staining area in pixels
                post_staining_area_um2: postsynapse staining area in um2
                pre_mean_puncta_size: presynapse mean puncta size in pixels
                pre_mean_puncta_size_um2: presynapse mean puncta size in um2
                post_mean_puncta_size: postsynapse mean puncta size in pixels
                post_mean_puncta_size_um2: postsynapse mean puncta size in um2
    """

    # thresholding
    if threshold_algorithm == "otsu":
        presynapse_threshold = filters.threshold_otsu(presynapse_image)
        postsynapse_threshold = filters.threshold_otsu(postsynapse_image)
    elif threshold_algorithm == "isodata":
        presynapse_threshold = filters.threshold_isodata(presynapse_image)
        postsynapse_threshold = filters.threshold_isodata(postsynapse_image)
    elif threshold_algorithm == "triangle":
        presynapse_threshold = filters.threshold_triangle(presynapse_image)
        postsynapse_threshold = filters.threshold_triangle(postsynapse_image)
    
    # presynapse
    presynapse_thresholded = presynapse_image > presynapse_threshold
    pre_labeled = measure.label(presynapse_thresholded, connectivity = 1) # 2 for 2d image
    pre_prop = measure.regionprops_table(pre_labeled, properties = ["area"]) # only calculates one property, is faster
    df_pre_prop = pd.DataFrame(pre_prop)
    pre_puncta_nr = df_pre_prop.shape[0] # puncta nr
    pre_puncta_density_per_100_um2 = (pre_puncta_nr / (image_size_um)**2) * 100 # puncta density per 100 um2
    pre_staining_area_pix = df_pre_prop["area"].sum() # staining area in pixels
    pre_staining_area_um = pre_staining_area_pix*(pixel_size_um)**2 # staining area in um2
    pre_mean_puncta_size_pix = df_pre_prop["area"].mean() # mean puncta size in pixels
    pre_mean_puncta_size_um = pre_mean_puncta_size_pix*(pixel_size_um)**2 # mean puncta size in um2

    # postsynapse
    postsynapse_thresholded = postsynapse_image > postsynapse_threshold
    post_labeled = measure.label(postsynapse_thresholded , connectivity = 1) # 2 for 2d image
    post_prop = measure.regionprops_table(post_labeled, properties = ["area"]) # only calculates one property, is faster
    df_post_prop = pd.DataFrame(post_prop)
    post_puncta_nr = df_post_prop.shape[0] # puncta_nr
    post_puncta_density_per_100_um2 = (post_puncta_nr / (image_size_um)**2) * 100 # puncta density per 100 um2
    post_staining_area_pix = df_post_prop["area"].sum() # staining area in pixels
    post_staining_area_um = post_staining_area_pix*(pixel_size_um)**2 # staining area in um2
    post_mean_puncta_size_pix = df_post_prop["area"].mean() # mean puncta size in pixels in um2
    post_mean_puncta_size_um = post_mean_puncta_size_pix*(pixel_size_um)**2 # mean puncta size in um2
    
    # make dictionary for all metrics
    metrics =   {"pre_puncta_nr": pre_puncta_nr,
                "pre_puncta_density_per_100_um2": pre_puncta_density_per_100_um2,
                "post_puncta_nr": post_puncta_nr,
                "post_puncta_density_per_100_um2": post_puncta_density_per_100_um2,
                "pre_staining_area_pix": pre_staining_area_pix,
                "pre_staining_area_um2": pre_staining_area_um,
                "post_staining_area_pix": post_staining_area_pix,
                "post_staining_area_um2": post_staining_area_um,
                "pre_mean_puncta_size_pix": pre_mean_puncta_size_pix,
                "pre_mean_puncta_size_um2": pre_mean_puncta_size_um,
                "post_mean_puncta_size_pix": post_mean_puncta_size_pix,
                "post_mean_puncta_size_um2": post_mean_puncta_size_um
    }
    return metrics, presynapse_thresholded, postsynapse_thresholded






