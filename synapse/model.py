import numpy as np
import math
import matplotlib.pyplot as plt
from .io import load_input_pickle, save_input_pickle, load_output_pickle, save_output_pickle


def runModel(params):

    """
    
    simple model of a synapse
    
    INPUT:
    -dictionary of parameters for simulation 
    
    OUTPUT: 
    -"simulation" object that contains the following data:
        - time
        - ap_times
        - ap_inds
        - cav_openings
        - ca_kernel
        - Ca_t
        - p_v_successes
        - quantal_content_per_syn
        - epsc_per_syn
        - quantal_content
        - epsc
        
        (also, the following if Vesicle Depletion is on:)
        (-
        (-
        
    """
    
#    print("Checking Parameters...")
    
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

#     print("Simulating Calcium Channel Opening...")

    cav_openings = np.zeros(np.array([no_stims,no_trials,no_syn]).astype(int))

    for i in range(params["num_cav"]*params["num_cav_ratio"]):
        cav_successes = np.random.uniform(size=(no_stims,no_trials,no_syn)) < params["cav_p_open"] 
        cav_openings += cav_successes*params["cav_i"]/params["num_cav_ratio"]    

    # define exponential decay as [Ca] kernel
    ca_kernel = np.exp(-params["stim_int"]*np.arange(no_stims)/params["ca_decay"])

    # generate [Ca](stim_num,trial) by convolving with cav_openings
    # crop for no_stim length

    Ca_t = np.apply_along_axis(lambda m: np.convolve(m,ca_kernel), axis=0, arr=cav_openings)
    Ca_t = Ca_t[0:no_stims,:,:]

    #########################################
    # Simulate Ca-Dependent Vesicle Release #
    #########################################

#     print("Simulating [Ca]-Dependent Vesicle Release...")

    # apply hill function to obtain probability of vesicle release
    p_v = hill(Ca_t,S=1,ec50=params["ca_ec50"],n=params["ca_coop"])

    # extract values corresponding to action-potential timepoints
    # multiplied by CaV opening success/failure (to prevent vesicle
    # release due purely to residual calcium).
    # Also multiply by a scaling factor that indicates the probability
    # that the vesicle is nearby the calcium channel cluster

    corrected_p = p_v*cav_successes*params["vesicle_prox"]
    
    # then randomly sample to generate quantal events
    p_v_successes = (np.random.uniform(size=corrected_p.shape) < corrected_p)*1


    if(params["depletion_on"]):

        ################################
        # Simulate Vesicular Depletion #
        ################################
        
#         print("Simulating Vesicular Depletion...")
        
        # Now, release_successes must be multiplied by whether or not a vesicle is present
        # So we must model the readily-releasable pool and its depletion / replenishment
        # To do this, we need a matrix of the shape of pool_size that indicates whether
        # the RRP is occupied (modeled as only 1 vesicle with exponential recovery):

        pool_size = np.zeros(p_v_successes.shape,dtype=int) + 1

        pool_recovery_fraction = 1-np.exp(-params["stim_int"]/params["pool_tau"]) # frac recovered every stim interval

        release_successes = np.zeros(p_v_successes.shape,dtype=int)

        # Now we iterate through each stimulus (in parallel) and check if the empty pools
        # recover (doing a flat fraction according to pool_tau generates exponential recovery
        # on average)

        for i in range(no_stims):
            # first calculate whether a vesicle is present
            # for pool_size == 0, check rand against exponential recovery
            pool_zero_inds = np.array(np.where(pool_size[i,:,:]==0))
            cols,chunks = (pool_zero_inds[0,:],pool_zero_inds[1,:])
            rows = np.zeros(cols.shape,dtype=int)+i
            pool_size[rows,cols,chunks] = (np.random.uniform(size=chunks.size) < pool_recovery_fraction ) * 1
            
            # now check if there was a successful vesicle release (pool_size * p_v_successes)
            release_successes[i,:,:] = pool_size[i,:,:]*p_v_successes[i,:,:]
            
            # if you're still within range, for all successes, set next stim's pool size to zero
            if i+1<no_stims:
                release_inds = np.array(np.where((p_v_successes[i,:,:]==1)&(pool_size[i,:,:]==1)))
                cols,chunks = (release_inds[0,:],release_inds[1,:])
                rows = np.zeros(cols.shape,dtype=int)+i+1     # next stimulus
                pool_size[rows,cols,chunks] = 0
                # for all inds where the pool_size was 0, nothing (they get another chance to reload)
        
        quantal_content_per_syn = release_successes
    
    else:
        quantal_content_per_syn = p_v_successes
    

    ###########################
    # Simulate AMPA Responses #
    ###########################    

#     print("Simulating AMPA Responses...")

    # obtain total quantal content per trial (sum across synapses)
    # in order to obtain total EPSC

    quantal_content = np.sum(quantal_content_per_syn,axis=2)

    # quantify bulk EPSC and individual synapse EPSC
    epsc = quantal_content*params["quantal_size"]
    epsc_per_syn = quantal_content_per_syn*params["quantal_size"]

    epsc_ave = np.mean(epsc,axis=1)

    #####################
    # Packaging Results #
    #####################   
    
#     print("Packaging Results....")
    
    data = {
        "params" : params,
        "time" : time,
        "ap_times" : ap_times,
        "ap_inds" : ap_inds,
        "cav_openings" : cav_openings,
        "ca_kernel" : ca_kernel,
        "Ca_t" : Ca_t,
        "p_v_successes" : p_v_successes,
        "quantal_content_per_syn" : quantal_content_per_syn,
        "epsc_per_syn" : epsc_per_syn,
        "quantal_content" : quantal_content,
        "epsc" : epsc,
        "epsc_ave" : epsc_ave
        }

    sim_run = simulation_run(data)
    
    return sim_run


def hill(ca,S=1,ec50=0.7,n=3.72):
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

class simulation_run:
    """
    A simple class for storing data (and parameters) from a simulation run
    """

    def __init__(self, data):
        for key,value in data.items():
            setattr(self,key,value)
    
    def update(self,newdata):
        for key,value in newdata:
            setattr(self,key,value)

