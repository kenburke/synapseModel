import numpy as np
from .io import load_input_pickle, save_input_pickle, load_output_pickle, save_output_pickle


def runModel(input_loc=None, params=None):

    """
    
    simple model of a synapse
    
    INPUT:
    -location of file for parameters for simulation 
        --pickled dict from "input" folder
    
    OUTPUT: 
    -"simulation" object that contains:
        -dictionary with simulation metadata, similar format to "input"
        -MxN array of traces, where
            --M is length of trace in samples
            --N is number of trials
    
    """
    
    
    
    params = load_input_pickle(param_file)
    
    