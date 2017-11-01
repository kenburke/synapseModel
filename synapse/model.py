import numpy as np


def runModel(param_file):
    """
    
    simple model of a synapse
    
    INPUT:
    -location of file for parameters for simulation (pickled dict)
    
    OUTPUT: 
    -another pickled dictionary with simulation metadata
    -MxN array of traces, where
        --M is length of trace in samples
        --N is number of trials
    
    """
    
    