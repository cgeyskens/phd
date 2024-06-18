from skimage import measure, transform, filters
import numpy as np

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
