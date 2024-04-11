import czifile
import os
from lxml import etree

class ProcessImageMetaData:
    
    def __init__(self, path):
        self.path = path
        pass

    def extract_metadata(self):
        """
        This function returns metadata from an image taken on a ZEISS system with the .czi file format.
        It also check if the image is square, meaning the X and Y dimension are equal. 

        Args:
            path (str): the location of the image.

        Returns:
            pixelsize: the size of each pixel in um.
            image_dimesion: the size of the image in pixels (X and Y are equal).
        
        Raises:
            ValueError: If the image is not square.

        """
        # this retrieves the metadata
        czi = czifile.CziFile(self.path)
        czi_xml_str = czi.metadata() # gets the metadata in a xml string format
        czi_parsed = etree.fromstring(czi_xml_str) # parses the czi_xml_str file

        # finds the strings 
        size_x = czi_parsed.find(".//SizeX")
        size_y = czi_parsed.find(".//SizeY")
        scaling_x = czi_parsed.find(".//ScalingX")
        scaling_y = czi_parsed.find(".//ScalingY")

        # extracting the required values to calculate the 
        size_x_value = int(size_x.text)
        size_y_value = int(size_y.text)
        scaling_x_value = float(scaling_x.text)
        scaling_y_value = float(scaling_y.text)

        # calculate the pixel to micrometer (um)
        # first checking whether the X and Y dimensions of the image are equal
        if size_x_value == size_y_value and scaling_x_value == scaling_y_value:
            pixel_size = ((scaling_x_value*1000000000)/size_x_value) # conversion from meter to micrometer 
        else:
            raise ValueError("imported image is not square")
        
        image_dimension = size_x_value 

        return pixel_size, image_dimension


    def image_filename(self, indices):
        """
        This function extracts the wanted elements from the filename of the image.

        Args:
            filename (str): the filename from which the elements needed to extract.
            indices (list): part of the filename that I would like to extract from.  
        
        Returns:
            desired_filename (str): the desired filename with only the most important elements.

        """
        # get the filename
        name_of_file = os.path.splitext(os.path.basename(self.path))[0]
        split_filename = name_of_file.split("_")

        # get only the experimental parameters from the filename
        desired_parts = [split_filename[val] for val in indices]
        desired_filename = "_".join(desired_parts)

        return desired_filename

