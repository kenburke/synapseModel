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
            "num_trials" : 1000,     # number of simulation trials to run
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
    no_trials = params["num_trials"]
    no_stims = params["num_stim"]
    time = np.arange(no_samples)*dt
    
    
    # array of AP times    
    ap_times = np.zeros(no_samples)
    ap_inds = np.zeros(no_stims)
    
    for i in range(no_stims):
        curr_stim_time = params["stim1_time"] + i * params["stim_int"]
        curr_ind = math.floor(curr_stim_time*FS)
        try:
            ap_times[curr_ind] = 1
            ap_inds[i] = curr_ind
        except IndexError as err:
            print("INDEXERROR: {0}".format(err))
            print("Stimulation parameters may exceed length of sweep")
            return
    
    # define an MxN matrix, where
    #   - M = length of one trace
    #   - N = number of trials
    # such that each column represents a sum of Kronecker delta functions
    # indicating the timing of _successful_ CaV opening
    #
    # repeated N times, where N is the number of channels, and summed
    # this represents all trials, stimulations and CaV's being independent

    cav_activity = np.zeros((no_samples,no_trials))

    for i in range(params["num_cav"]):
        cav_successes = np.random.uniform(size=(no_stims,no_trials)) < params["cav_p_open"] 
        cav_inds = np.array(np.where(cav_successes))
        rows = ap_inds[cav_inds[0,:]]
        cols = cav_inds[1,:]
        cav_activity[rows,cols] += params["cav_i"]

    print("out of loop")

    # define exponential decay as [Ca] kernel
    ca_kernel = np.exp(-time/params["ca_decay"])

    # generate [Ca](t,trial) by convolving against cav_activity
    # crop for sweep length (note, time of peaks = ap_inds)

    Ca_t = np.apply_along_axis(lambda m: np.convolve(m,ca_kernel), axis=0, arr=cav_activity)
    Ca_t = Ca_t[0:len(time),:]

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
    
    
    
    
    