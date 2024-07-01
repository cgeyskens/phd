def get_regions_from_dict(dictionary, synapse_marker):
    """
    This will get the the optimal parameter ranges from the dictionary.
    They key in each dictionary is the synaptic marker combination.

    Args:
        dictionary: the dictionary with the optimal parameter ranges.
        synapse_marker: from which synaptic marker it must extract the parameters.
    
    Returns:
        list with the parameter ranges specified by the synaptic markers
    """
    if synapse_marker in dictionary:
        return list(dictionary[synapse_marker].keys())
    else:
        return []


def get_synaptic_marker(protein_and_synaptic_marker):
    """
    Very simple function to extract the synaptic marker from the protein_and_synaptic_marker
    string. Example: protein_and_synaptic_marker = "VCAM1_LacZ_VGLUT1_PSD95". It returns in this case 
    the VGLUT1_PSD95.
    """
    synaptic_marker_string = protein_and_synaptic_marker.split("_")[-2:]
    synaptic_marker = "_".join(synaptic_marker_string)
    return synaptic_marker


def get_hippocampal_layer(file_path):
    """
    Very simple function to extract the hippocampal layer from the file_path string.
    """
    parts = file_path.split('_')[-2:]
    result = "_".join(parts)
    layer = result[:-4]
    return layer