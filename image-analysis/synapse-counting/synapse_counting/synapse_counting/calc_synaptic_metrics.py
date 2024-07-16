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


def puncta_metrics(presynapse_thresholded, postsynapse_thresholded, image_size_um, pixel_size_um, puncta_size_threshold = 0):
    """
    Returns different puncta metrics.

    Args:
        presynapse_thresholded: thresholded image of presynapse
        postsynapse_thresholded: thresholded image of postsynapse
        image_size_um: the size of the image in um.
        pixel_size_um: the size of each pixel in um
        puncta_size_threshold: threshold for puncta, filtering out smaller puncta
    
    Returns:
        metrics: dictionary with different puncta metrics:
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
                    pre_labeled: presynapse thresholded image that is labeled
        post_labeled: postsynapse thresholded image that is labeled
        df_pre_prop: dataframe of presynaptic puncta properties, with label and area
        df_post_prop: dataframe of postsynaptic puncta properties, with label and area
    """
    # presynapse
    pre_labeled = measure.label(presynapse_thresholded, connectivity = 1) # 2 for 2d image
    pre_prop = measure.regionprops_table(pre_labeled, properties = ["label", "area"]) # only calculates one property, is faster
    df_pre_prop = pd.DataFrame(pre_prop)
    df_pre_prop_filtered = df_pre_prop[df_pre_prop["area"] > puncta_size_threshold]
    pre_puncta_nr = df_pre_prop_filtered.shape[0] # puncta nr
    pre_puncta_density_per_100_um2 = (pre_puncta_nr / (image_size_um)**2) * 100 # puncta density per 100 um2
    pre_staining_area_pix = df_pre_prop_filtered["area"].sum() # staining area in pixels
    pre_staining_area_um = pre_staining_area_pix*(pixel_size_um)**2 # staining area in um2
    pre_mean_puncta_size_pix = df_pre_prop_filtered["area"].mean() # mean puncta size in pixels
    pre_mean_puncta_size_um = pre_mean_puncta_size_pix*(pixel_size_um)**2 # mean puncta size in um2

    filtered_pre_puncta = np.isin(pre_labeled, df_pre_prop_filtered["label"].values) * pre_labeled # for visualization


    # postsynapse
    post_labeled = measure.label(postsynapse_thresholded , connectivity = 1) # 2 for 2d image
    post_prop = measure.regionprops_table(post_labeled, properties = ["label", "area"]) # only calculates one property, is faster
    df_post_prop = pd.DataFrame(post_prop)
    df_post_prop_filtered = df_post_prop[df_post_prop["area"] > puncta_size_threshold]
    post_puncta_nr = df_post_prop_filtered.shape[0] # puncta_nr
    post_puncta_density_per_100_um2 = (post_puncta_nr / (image_size_um)**2) * 100 # puncta density per 100 um2
    post_staining_area_pix = df_post_prop_filtered["area"].sum() # staining area in pixels
    post_staining_area_um = post_staining_area_pix*(pixel_size_um)**2 # staining area in um2
    post_mean_puncta_size_pix = df_post_prop_filtered["area"].mean() # mean puncta size in pixels in um2
    post_mean_puncta_size_um = post_mean_puncta_size_pix*(pixel_size_um)**2 # mean puncta size in um2

    filtered_post_puncta = np.isin(post_labeled, df_post_prop_filtered["label"].values) * post_labeled # for visualization
    
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
    return metrics, pre_labeled, post_labeled, df_pre_prop, df_post_prop, filtered_pre_puncta, filtered_post_puncta


def filter_out_small_puncta(binary_image, prop_df_with_area, labeled_image, puncta_size_threshold):
    """
    Filtering out smaller puncta from a binary image and returning the filtered image. 

    Args:
        binary_image: binary_image used for filtering
        prop_df_with_area: dataframe output from the scikit-image measure.regionprops_table function, including a column with the area per puncta in pixels
        labeled_image: the labeled image, output from the scikit-image measure.label function
        puncta_size_threshold: the threshold of the puncta_size (in pixels), only above this threshold, puncta are returned

    Returns:
        filtered_image: a filtered image, not including puncta from the below the set threshold
    """
    prop_df_filtered = prop_df_with_area[prop_df_with_area['area'] > puncta_size_threshold]
    # create empty binary mask
    filtered_image = np.zeros_like(binary_image, dtype=bool)
    # iterate through the filtered properties
    for index, row in prop_df_filtered.iterrows():
        # get the label of the current region
        label = row['label']
        # set the pixels of the current region in the filtered_image to True
        filtered_image[labeled_image == label] = True

    return filtered_image




