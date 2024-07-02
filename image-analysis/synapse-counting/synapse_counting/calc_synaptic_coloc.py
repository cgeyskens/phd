from skimage import measure, transform, filters, feature, transform
import numpy as np
import matplotlib.pyplot as plt

from scipy.spatial.distance import cdist

def pearsons_coloc(presynapse_image, postsynapse_image):
    """
    Wrapper function for pearsons correlation of colocalization.

    Args:
        presynapse_image: image of presynapse
        postsynapse_image: image of postsynapse

    Returns:
        pcc: pearsons correlation coefficient of colocalization
        pval: p-value of pearsons correlation coefficient of colocalization
        pcc_rot: internal control- pearsons correlation coefficient of colocalization when postsynapse image is rotated by 90 degrees
        pval_rot: internal control- p-value of pearsons correlation coefficient of colocalization
    """
    postsynapse_image_rot = transform.rotate(postsynapse_image, 90)
    pcc, pval = measure.pearson_corr_coeff(presynapse_image, postsynapse_image)
    pcc_rot, pval_rot = measure.pearson_corr_coeff(presynapse_image, postsynapse_image_rot)
    
    return pcc, pval, pcc_rot, pval_rot


def manders_coloc(presynapse_image_thresholded, postsynapse_image_thresholded):
    """
    Wrapper function for manders overlap coefficients. Based on binary images.

    Args: 
        presynapse_image: image of presynapse
        postsynapse_image: image of postsynapse
    
    Returns:
        overlap_coeff: manders overlap coefficient
        overlap_coeff_rot: manders overlap coefficient of rotated control
        presynapse_threshold: threshold of the presynapse
        postsynapse_threshold: threshold of the postsynapse
    """   
    # rotation, internal control
    prostsynapse_image_threshold_rot = transform.rotate(postsynapse_image_thresholded, 90)
    
    # manders overlap coefficient
    overlap_coeff = measure.manders_overlap_coeff(presynapse_image_thresholded, postsynapse_image_thresholded)
    overlap_coeff_rot = measure.manders_overlap_coeff(presynapse_image_thresholded, prostsynapse_image_threshold_rot)
    
    return overlap_coeff, overlap_coeff_rot


def overlap_um2_coloc(presynapse_image_thresholded, postsynapse_image_thresholded, pixel_size_in_um):
    """
    Calculates the overlap in um2. Based on binary images.

    Args: 
        presynapse_image_thresholded: image of presynapse thresholded
        postsynapse_image_thresholded: image of thresholded
        pixel_size_in_um: pixel size of the image in um
    
    Returns:
        overlap_um: overlap in um
        overlap_um_rot: overlap in um for rotated control
    """  
    # rotation, internal control
    prostsynapse_image_threshold_rot = transform.rotate(postsynapse_image_thresholded, 90)

    overlap = presynapse_image_thresholded & postsynapse_image_thresholded
    overlap_rot = presynapse_image_thresholded & prostsynapse_image_threshold_rot

    # calculating over in pixels
    overlap_pix = np.sum(overlap)
    overlap_pix_rot = np.sum(overlap_rot)

    # overlap in um
    overlap_um2 = overlap_pix * pixel_size_in_um * pixel_size_in_um
    overlap_um2_rot = overlap_pix_rot * pixel_size_in_um * pixel_size_in_um

    return overlap_um2, overlap_um2_rot


def local_peak_detection(presynapse_preprocessed, postsynapse_preprocessed, presynapse_distance, postsynapse_distance, presynapse_threshold, postsynapse_threshold, plot_coord = False):
    """
    Detects the local intensity peaks of the presynapse and the postsynapse channel.
    
    Args:
        presynapse_preprocessed (np.array): processed image of presynapse, which is background substracted and has a gaussian blur
        postsynapse_preprocessed (np.array): processed image of postsynapse, which is background substracted and has a gaussian blur
        presynapse_distance (int): minimun distance between presynapse local peak maxima, usually 1
        postsynapse_distance (int): mimimun distance between postsynapse local peak maxima, usually 1
        presynapse_threshold (float): thresholding of the presynapse image for local peak detection
        postsynapse_threshold (float): thresholding of the postsynapse image for local peak detection
        plot_coordinates (bool): option to visualize the images

    Returns:
        presynapse_coord (np.array): coordinates of presynapse local peak maxima
        postsynapse_coord (np.array): coordinates of postsynapse local peak maxima
        postsynapse_coord_rot (np.array): coordinates of postsynapse local peak maxima, but rotated
        
        plot of vglut1, psd95, psd95_rot images with the local peak maxima overlaid
    """
    # Thresholding the image for local peak maximum detection
    presynapse_coord = feature.peak_local_max(presynapse_preprocessed, min_distance = presynapse_distance, threshold_abs = presynapse_threshold)
    postsynapse_coord = feature.peak_local_max(postsynapse_preprocessed, min_distance = postsynapse_distance, threshold_abs = postsynapse_threshold)

    # Rotating an image (psd95) as a control
    postsynapse_rot = transform.rotate(postsynapse_preprocessed, 90)
    postsynapse_rot_coord = feature.peak_local_max(postsynapse_rot, min_distance = postsynapse_distance, threshold_abs = postsynapse_threshold)
    
    if plot_coord == True:
        # Showing the local peaks with coordinates together with the images
        fig, axs = plt.subplots(2, 2, figsize=(30, 30))

        axs[0,0].imshow(presynapse_preprocessed, cmap='gray')
        #axs[0,0].plot(presynapse_coord[:, 1], presynapse_coord[:, 0], 'c.')
        axs[0,0].set_title('vglut1_pre')

        axs[0,1].imshow(postsynapse_preprocessed, cmap='gray')
        #axs[0,1].plot(postsynapse_coord[:, 1], postsynapse_coord[:, 0], 'm.')
        axs[0,1].set_title('psd95_pre')

        axs[1,0].imshow(presynapse_preprocessed, cmap='gray')
        axs[1,0].plot(presynapse_coord[:, 1], presynapse_coord[:, 0], 'c.')
        axs[1,0].set_title('vglut1_pre')

        axs[1,1].imshow(postsynapse_preprocessed, cmap='gray')
        axs[1,1].plot(postsynapse_coord[:, 1], postsynapse_coord[:, 0], 'm.')
        axs[1,1].set_title('psd95_pre')
    else:
        pass
    
    return presynapse_coord, postsynapse_coord, postsynapse_rot_coord


def count_coloc_spots(presynapse_coordinates, postsynapse_coordinates, pixel_size_um, max_distance_um):
    """
    Counts the number of colocalized pre and postsynaptic spots based on the coordinates array of the local maxima detection.
    
    Args:
        presynapse_coordinates (np.array): coordinates of the local peak maxima in the vlgut1 channel
        postsynapse_coordinates (np.array): coordinates of the local peak maxima in the psd95 channel
        pixel_size_um (float): size of a pixel in um, based on image settings
        max_distance_um (float): value of the maximum colocalization distance between the spots of each channel in um
    
    Returns:
        colocalized spot count (int): the number of colocalized pre and postsynaptic spots.
    """
    # calculate the max distance in pixels with the pixel_size_um and max_distance_um
    max_distance_px = max_distance_um/pixel_size_um
    
    # calculate pairwise distances between spots in vlgut1 and psd95 channel
    distances_pre_to_post = cdist(presynapse_coordinates, postsynapse_coordinates)
    distances_post_to_pre = cdist(postsynapse_coordinates, presynapse_coordinates)

    # find unique colocalized spots
    colocalized_spots = set()

    # tterate over distances from presynapse_coordinates to postsynapse_coordinates
    for i in range(len(presynapse_coordinates)):
        # Check if the current spot in vlgut1 has nearby spots in psd95
        colocalized_indices = [j for j, distance in enumerate(distances_pre_to_post[i, :]) if distance <= max_distance_px]
        for j in colocalized_indices:
            colocalized_spots.add((i, j))

    # tterate over distances from postsynapse_coordinates to presynapse_coordinates
    for i in range(len(postsynapse_coordinates)):
        # check if the current spot in Channel 2 has nearby spots in Channel 1
        colocalized_indices = [j for j, distance in enumerate(distances_post_to_pre[i, :]) if distance <= max_distance_px]
        for j in colocalized_indices:
            colocalized_spots.add((j, i))

    # get the count of unique colocalized spots
    colocalized_spot_count = len(colocalized_spots)
    
    return colocalized_spot_count