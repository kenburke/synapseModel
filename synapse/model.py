import numpy as np
import math
from .io import load_input_pickle, save_input_pickle, load_output_pickle, save_output_pickle


def runModel(params=None):

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
    
    if params==None:

        #Define default params
        params = {
                                    ### General Simulation Params
            "fs" : 2e4,             # per second
            "sweep_length" : 1,     # seconds
            "num_trial" : 1000,     # number of simulation trials to run
            "num_stim" : 2,         # number of AP stimuli to run
            "stim_int" : 0.05,      # interval of AP stimuli (sec)
            "stim1_time" : 0.1,     # time of first stimulus (sec)
                            
                                    ### Calcium Channel Params
            "num_cav" : 1,          # number of voltage-gated calcium channels
            "cav_p_open" : 0.83,    # probability of CaV opening per AP
            "cav_i" : 1,            # single CaV current per AP (arbitrary units)
            "ca_decay" : 0.025,     # exponential decay time constant for [ca]
    
                                    ### Transmission / Hill Function Params
            "quantal_size" : -10,   # single AMPA current per vesicle (pA)
            "ca_coop" : 3.72,       # Hill Function param N for calcium ion cooperativity
                                    # in vesicular release, see Scimemi & Diamond
            "ca_ec50" : 0.67        # Hill Function param EC_50 for calcium concentration
                                    # for half-maximal vesicle release probability                
                            
                                    # NOTE: EC_50 and [Ca] decay constant were constrained 
                                    #       to get baseline PPR of 1.3 and post-baclofen
                                    #       PPR of 1.7           
            }
    
    
    
    if sweep_length<(stim1_time + stim_int*(num_stim-1)):
        # sweep length too short
        print("Stimulation parameters may exceed length of sweep")
        answer = input("Automatically readjust length of sweep? Y/N: ")
        if answer.lower()=="y":
            params["sweep_length"] = params["stim1_time"] + params["stim_int"] \
                * (params["num_stim"]-1) + params["ca_decay"]*4
        else:
            print("unable to correct, exiting simulation")
            return
    
    
    FS = float(params["fs"])
    dt = 1./FS
    no_samples = FS*params["sweep_length"]
    time = np.arange(no_samples)*dt
    
    # array of AP times    
    ap_times = np.zeros(no_samples)
    
    i = 0.
    
    while i < params["num_stim"]:
        curr_stim_time = params["stim1_time"] + i * params["stim_int"]
        print(curr_stim_time)
        try:
            ap_times[math.floor(curr_stim_time*FS)] = 1
        except IndexError as err:
            print("INDEXERROR: {0}".format(err))
            print("Stimulation parameters may exceed length of sweep")
            return
        i += 1
    
    ap_inds = np.where(ap_times==1)
    
    # define exponential decay as [Ca] kernel
    ca_kernel = np.exp(-time/params["ca_decay"])
    
    # generate [Ca](t), crop for sweep length (note, time of peaks = ap_inds)
    Ca_t = np.resize(np.convolve(ap_times,ca_kernel),len(ap_times))
    
    

def hill(ca,S=1,ec50=0.8,n=3.72)
    """
    
    Hill function
    
    INPUT:
    - numpy array of [Ca]
    - Hill parameters (S, EC_50, N)
    
    OUTPUT:
    - numpy array of size ca that is Hill(ca)
    
    """
    
    pr = (S * (ca ** n)) / (ec50 ** n + ca ** n)
    
    return pr
    
    
    
    
    