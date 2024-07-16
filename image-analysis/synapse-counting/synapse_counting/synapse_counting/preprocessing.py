import numpy as np
import czifile
from skimage import restoration, exposure, morphology, filters
from skimage.filters import gaussian, sobel
from skimage import measure, filters
from skimage.feature import peak_local_max
from scipy import ndimage as ndi
from skimage.measure import label
from skimage.morphology import binary_opening
from skimage.segmentation import watershed


def extract_and_split(path, presynapse_channel = 0, postsynapse_channel = 1):
    """
    Imports the image and gets rid of all the extra channels 
    & splits the channel into presynapse and postsynapse for further processing

    Args:
        path: the path of the image
        presynapse_channel: the channel where the presynapse is imaged in. standard is 0.
        postsynapse_channel: the channel where the postsynapse is imaged in. standard is 1.
    
    Returns:
        presynapse: the presynapse image
        postsynapse: the postsynapse image
    """
    image = czifile.imread(path) # reading in the image
    image_squeezed = np.squeeze(image) # only get the relevant channels, i.e. the shape of the data
    presynapse = image_squeezed[presynapse_channel,:,:]
    postsynapse = image_squeezed[postsynapse_channel,:,:]

    return presynapse, postsynapse


class ImagePreprocessing:
    """
    Class for preprocessing images. 
    
    Its structured in such a way that that it allows to choose which preprocessing algorithms and filters to use.
    This will be usefull when an image needs different preprocessing steps depending on the colocalization algorithm.
    (1) background substraction with rooling ball radius
    (2) CLAHE for image normalization
    (3) tophatfilter
    (4) guassian blur

    """
    def __init__(self, 
                 include_rolling_ball = True, radius = 10, # rolling ball parameters
                 include_clahe = True, clip_limit = 0.005, kernel_size = 150, nbins = 265, # CLAHE parameters
                 include_tophat = True, element_size = 5, # tophat parameters
                 include_blur = True, sigma = 1, preserve_range = True # gaussian blur filters
                 ):
        """
        Initializes the ImageProcessor. 
        With parameters that I found important for testing.

        Args:
            include_rolling_ball: whether or not to include rolling ball.
            radius: radius parameter for rolling ball.
            include_clahe: hether or not to include CLAHE.
            clip_limit: clip limit parameter for CLAHE.
            kernel_size: kernel_size parameter for CLAHE.
            nbins: number of bins for CLAHE.
            include_tophat: whether or not to include tophat.
            element_size: size of the disk in tophat.
            include_blur: whether or not to include a gaussina blur.
            sigma: sigma parameter for the gaussian blur.
            preserve_range: whether or not to preserve the initial value ranges of the input image.
        
        Returns:
            a presynapse and postsynapse processed image.

        """
        self.include_rolling_ball = include_rolling_ball        
        self.radius = radius
        self.include_clahe = include_clahe
        self.clip_limit = clip_limit
        self.kernel_size = kernel_size
        self.nbins = nbins
        self.include_tophat = include_tophat
        self.element_size = element_size
        self.include_blur = include_blur
        self.sigma = sigma
        self.preserve_range = preserve_range
    
    def preprocess(self, presynapse_image, postsynapse_image):
        if self.include_rolling_ball:
            presynapse_image, postsynapse_image = self._rolling_ball(presynapse_image, postsynapse_image)
        if self.include_clahe:
            presynapse_image, postsynapse_image = self._clahe(presynapse_image, postsynapse_image)
        if self.include_tophat:
            presynapse_image, postsynapse_image = self._tophat(presynapse_image, postsynapse_image)
        if self.include_blur:
            presynapse_image, postsynapse_image = self._blur(presynapse_image, postsynapse_image)
        return presynapse_image, postsynapse_image
    
    def _rolling_ball(self, presynapse_image, postsynapse_image):
        background_presynapse = restoration.rolling_ball(presynapse_image, radius = self.radius)
        background_postsynapse = restoration.rolling_ball(postsynapse_image, radius = self.radius)
        presynapse_image = presynapse_image - background_presynapse
        postsynapse_image = postsynapse_image - background_postsynapse
        return presynapse_image, postsynapse_image
    
    def _clahe(self, presynapse_image, postsynapse_image):
        presynapse_image = exposure.equalize_adapthist(presynapse_image, clip_limit = self.clip_limit, kernel_size = self.kernel_size, nbins = self.nbins)
        postsynapse_image = exposure.equalize_adapthist(postsynapse_image, clip_limit = self.clip_limit, kernel_size = self.kernel_size, nbins = self.nbins)
        return presynapse_image, postsynapse_image
    
    def _tophat(self, presynapse_image, postsynapse_image):
        footprint = morphology.disk(self.element_size)
        presynapse_image = morphology.white_tophat(presynapse_image, footprint)
        postsynapse_image = morphology.white_tophat(postsynapse_image, footprint)
        return presynapse_image, postsynapse_image

    def _blur(self, presynapse_image, postsynapse_image):
        presynapse_image = filters.gaussian(presynapse_image, sigma = self.sigma, preserve_range = self.preserve_range)
        postsynapse_image = filters.gaussian(postsynapse_image, sigma = self.sigma, preserve_range = self.preserve_range)
        return presynapse_image, postsynapse_image


def thresholding(presynapse_image, postsynapse_image, threshold_algorithm="triangle"):
    if threshold_algorithm == "otsu":
        presynapse_threshold = filters.threshold_otsu(presynapse_image)
        postsynapse_threshold = filters.threshold_otsu(postsynapse_image)
    elif threshold_algorithm == "isodata":
        presynapse_threshold = filters.threshold_isodata(presynapse_image)
        postsynapse_threshold = filters.threshold_isodata(postsynapse_image)
    elif threshold_algorithm == "triangle":
        presynapse_threshold = filters.threshold_triangle(presynapse_image)
        postsynapse_threshold = filters.threshold_triangle(postsynapse_image)
    
    presynapse_image_threshold = presynapse_image >= presynapse_threshold
    postsynapse_image_threshold = postsynapse_image >= postsynapse_threshold
    
    return presynapse_image_threshold, postsynapse_image_threshold


def custom_watershed(thresholded_image, sigma):
    """
    Wrapper for watershed segmentation, copied from Robert Haase notebook:
    https://haesleinhuepf.github.io/BioImageAnalysisNotebooks/20h_segmentation_post_processing/mimicking_imagej_watershed.html?highlight=watershed
    """
    distance = ndi.distance_transform_edt(thresholded_image) # calculate distance image
    blurred_distance = gaussian(distance, sigma = sigma) # gaussian blur
    fp = np.ones((3,) * thresholded_image.ndim) # neighbourhood size to find local maxima
    coords = peak_local_max(blurred_distance, footprint=fp, labels=thresholded_image) # find local maxima
    mask = np.zeros(distance.shape, dtype=bool)
    mask[tuple(coords.T)] = True
    markers = label(mask)
    labels = watershed(-blurred_distance, markers, mask=thresholded_image) # actual watershed
    edges_labels = sobel(labels)
    edges_binary = sobel(thresholded_image)
    edges = np.logical_xor(edges_labels != 0, edges_binary != 0)
    almost = np.logical_not(edges) * thresholded_image
    watershedded_image = binary_opening(almost)

    return watershedded_image
