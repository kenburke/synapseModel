import numpy as np
import math
import matplotlib.pyplot as plt
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
    
    print("Checking Parameters...")
    if params==None:

        #Define default params
        params = {
                                    ### General Simulation Params
            "fs" : 2e4,             # per second
            "sweep_length" : 1,     # seconds
            "num_trials" : 100,     # number of simulation trials to run
            "num_stim" : 2,         # number of AP stimuli to run
            "stim_int" : 0.05,      # interval of AP stimuli (sec)
            "stim1_time" : 0.1,     # time of first stimulus (sec)
                            
                                    ### Calcium Channel Params
            "num_syn" : 50,         # number of synapses
            "num_cav" : 1,          # number of voltage-gated calcium channels
            "cav_p_open" : 0.83,    # probability of CaV opening per AP
            "cav_i" : 1,            # single CaV current per AP (arbitrary units)
            "ca_decay" : 0.046,     # exponential decay time constant for [ca]
    
                                    ### Transmission / Hill Function Params
            "ca_coop" : 3.72,       # Hill Function param N for calcium ion cooperativity
                                    # in vesicular release, see Scimemi & Diamond
            "ca_ec50" : 0.67,       # Hill Function param EC_50 for calcium concentration
                                    # for half-maximal vesicle release probability                
                                    
                                    ### Postsynaptic Params
            "quantal_size" : -5,    # single AMPA current per vesicle (pA)
            "taus" : (0.001,0.005), # fast and slow AMPA current time constants
            

                                    # NOTE: EC_50 and [Ca] decay constant were constrained 
                                    #       to get baseline PPR of 1.3 and post-baclofen
                                    #       PPR of 1.7           
            }
    
    
    
    if params["sweep_length"]<(params["stim1_time"] + params["stim_int"]*(params["num_stim"]-1)):
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
    no_syn = params["num_syn"]
    time = np.arange(no_samples)*dt

    # array of AP times    
    ap_times = np.zeros(math.floor(no_samples))
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
    
    ####################################
    # Simulate Calcium Channel Opening #
    ####################################
    
    print("Simulating Calcium Channel Opening...")
    
    # define an MxNxO matrix, where
    #   - M = length of one trace
    #   - N = number of trials
    #   - O = number of synapses / active zones
    # such that each column represents a sum of Kronecker delta functions
    # indicating the timing of _successful_ CaV opening
    #
    # repeated N times, where N is the number of channels, and summed
    # this represents all trials, stimulations and CaV's being independent

    cav_activity = np.zeros(np.array([no_samples,no_trials,no_syn]).astype(int))

    for i in range(params["num_cav"]):
        cav_successes = np.random.uniform(size=(no_stims,no_trials,no_syn)) < params["cav_p_open"] 
        cav_inds = np.array(np.where(cav_successes))
        rows, cols, chunks = (
            ap_inds[cav_inds[0,:]],
            cav_inds[1,:],
            cav_inds[2,:]
            )
        cav_activity[rows.astype(int),cols.astype(int),chunks.astype(int)] += params["cav_i"]

    # define exponential decay as [Ca] kernel
    ca_kernel = np.exp(-time/params["ca_decay"])

    # generate [Ca](t,trial) by convolving against cav_activity
    # crop for sweep length (note, time of peaks = ap_inds)

    Ca_t = np.apply_along_axis(lambda m: np.convolve(m,ca_kernel), axis=0, arr=cav_activity)
    Ca_t = Ca_t[0:len(time),:,:]
 
    #########################################
    # Simulate Ca-Dependent Vesicle Release #
    #########################################

    print("Simulating [Ca]-Dependent Vesicle Release...")

    # apply hill function to obtain probability of vesicle release
    p_v = hill(Ca_t,S=1,ec50=params["ca_ec50"],n=params["ca_coop"])
    
    # extract values corresponding to action-potential timepoints
    # multiplied by CaV opening success/failure (to prevent vesicle
    # release due purely to residual calcium)
    # then randomly sample to generate quantal events
    corrected_p = p_v[ap_inds.astype(int)]*cav_successes
    num_vesicles = (np.random.uniform(size=corrected_p.shape) < corrected_p)*1

    # obtain total quantal content per trial
    # (sum across synapses)
    # obtain release inds and plug in quantal content into epsc_activity

    quantal_content = np.sum(num_vesicles,axis=2)
    nonzero_content = quantal_content[quantal_content>0].flatten()
    
    release_inds = np.array(np.where(quantal_content>0))
    rows,cols = (ap_inds[release_inds[0,:]],release_inds[1,:])
    
    epsc_activity = np.zeros(np.array([no_samples,no_trials]).astype(int))
    epsc_activity[rows.astype(int),cols.astype(int)] = quantal_content[
        release_inds[0,:].astype(int),
        release_inds[1,:].astype(int)
        ]
    
    ###########################
    # Simulate AMPA Responses #
    ###########################    
    
    print("Simulating AMPA Responses...")

    # define AMPA kernel
    ampa_kernel = ampa(time,quantal_size=params["quantal_size"],tau=params["taus"])
    
    Vm_t = np.apply_along_axis(lambda m: np.convolve(m,ampa_kernel), axis=0, arr=epsc_activity)
    Vm_t = Vm_t[0:len(time),:]
    
    plt.plot(time,Vm_t)
    plt.show()
    
    Vm_ave = np.mean(Vm_t,axis=1)
    plt.plot(time,Vm_ave)
    plt.show()    

    return Vm_t

def hill(ca,S=1,ec50=0.67,n=3.72):
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
    
def ampa(time,quantal_size=-5,tau=(1/1000,5/1000)):
    """
    
    AMPA kinetics as difference of exponentials
    
    INPUT:
    - numpy array of time
    - quantal size (maximal amplitude, pA)
    - time constants tau (in ms, tuple of (fast,slow))
    
    OUTPUT:
    - numpy array of size time
    
    """
    i_ampa = np.exp(-time/tau[1]) - np.exp(-time/tau[0])
    i_ampa *= quantal_size/max(i_ampa)
    
    return i_ampa

