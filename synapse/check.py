import numpy as np

def check_params(params):
    """
    Checks params for all kinds of silly problems
    """
    
    if params["sweep_length"]<(params["stim1_time"] + params["stim_int"]*(params["num_stim"]-1)):
        # sweep length too short
        print("Stimulation parameters may exceed length of sweep")
        answer = input("Automatically readjust length of sweep? Y/N: ")
        if answer.lower()=="y":
            params["sweep_length"] = params["stim1_time"] + params["stim_int"] \
                * (params["num_stim"]-1) + params["ca_decay"]*4
        else:
            raise ParamError("unable to correct, exiting simulation")
            
    
    if (params["num_stim"]>=1) & (params["num_cav_ratio"]>=1) & (params["num_trials"]>=1):
        params["num_stim"] = int(params["num_stim"])
        params["num_cav_ratio"] = int(params["num_cav_ratio"])
        params["num_trials"] = int(params["num_trials"])
    else:
        raise ParamError("ERROR: Must specificy >= 1 for number of stimuli, num/cav ratio and num_trials")
    
    if not (params["num_cav"]>=0) or not (params["cav_i"]>=0):
        raise ParamError("ERROR: Must specificy >= 0 for number of CaV or cav_i")
        
    if not (0<=params['cav_p_open']<=1) or not (0<=params['vesicle_prox']<=1):
        raise ParamError("ERROR: Must set probabilities cav_p_open and vesicle_prox on range 0<=p<=1")

class ParamError(Exception):
    """Basic error for problems with input parameters"""
    pass
