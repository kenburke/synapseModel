from synapse import model, io, utils
import os
import copy
from random import random, seed
from math import floor
import pytest
import itertools

@pytest.fixture
def param_base():
    '''Dict of parameter defaults'''

    return {
                                    ### General Simulation Params

        "fs" : 2e4,                 # per second
        "sweep_length" : 1,         # seconds
        "num_trials" : 300,         # number of simulation trials to run
        "num_stim" : 2,             # number of AP stimuli to run
        "stim_int" : 0.05,          # interval of AP stimuli (sec)
        "stim1_time" : 0.1,         # time of first stimulus (sec)

                                    ### Calcium Channel Params

        "num_syn" : 100,            # number of synapses
        "num_cav" : 1,              # number of voltage-gated calcium channels in complex
                                    # NOTE: [Ca] comes from all channels, as if
                                    # in a cluster (but still open independently)
        "cav_p_open" : 0.83,        # probability of CaV opening per AP
        "cav_i" : 1,                # single CaV current per AP (arbitrary units)
        "ca_decay" : 0.05,          # exponential decay time constant for [ca]
        "num_cav_ratio" : 2,        # mult num_cav and divide cav_i by this

                                    ### Transmission / Hill Function Params
        "vesicle_prox" : 0.25,      # Probability your vesicle is nearby the CaV complex
        "ca_coop" : 3.72,           # Hill Function param N for calcium ion cooperativity
                                    # in vesicular release, see Scimemi & Diamond
        "ca_ec50" : 0.7,            # Hill Function param EC_50 for calcium concentration
                                    # for half-maximal vesicle release probability                

                                    ### Vesicular Depletion Params

        #probably need num vesicles in RRP? and to adjust exponential resetting func?
        "depletion_on" : False,     # Turn on depletion?
        "pool_tau" : 1.00,          # Readily-Releasable Pool recovery time constant.
                                    # Here, we model the RRP as a binary variable, either
                                    # occupied by one vesicle or not. Thus, we use this
                                    # to define a probability of occupancy [0,1)

                                    ### Postsynaptic Params

        "quantal_size" : -10,       # single AMPA current per vesicle (pA)
        "a_tau" : (0.001,0.005),    # fast and slow AMPA current time constants


                                    # NOTE: EC_50 and [Ca] decay constant were constrained 
                                    #       to approximate observed changes in EPSC amplitude
                                    #       by SKF and Baclofen, given 40% reduction in
                                    #       I_cal = N_cav * p_open * cav_i       
        }
    


@pytest.fixture
def param_ranges():
    '''Dict of parameter ranges'''
    p_name = [
        "cav_p_open",
        "num_trials",
        "num_stim",
        "num_cav",
        "cav_i",
        "num_cav_ratio",
        "vesicle_prox",
        ]
    p_sets = [
        [0,0.01,0.99,1],   #cav_p_open
        [1,10,300],   #num_trials
        [1,2,5],  #num_stim
        [1,3,10],    #num_cav
        [0,1,5,10],  #cav_i
        [1,2],    #num_cav_ratio
        [0,0.01,0.25,1],    #vesicle_prox
        ]
    
    return dict(zip(p_name,p_sets))
    
parameter_names = [
    "cav_p_open",
    "num_trials",
    "num_stim",
    "num_cav",
    "cav_i",
    "num_cav_ratio",
    "vesicle_prox",
    ]


@pytest.mark.parametrize("p_name", parameter_names)
def test_runModel_range_params(p_name,param_ranges,param_base):

    alt_params = copy.deepcopy(param_base)
    p_range = param_ranges[p_name]
    
    for i in range(len(p_range)):
        alt_params[p_name] = p_range[i]
        SIM = utils.Simulation(params = alt_params)


# def test_runModel_combo_params(parameter_names,parameter_sets,param_base):
#     
#     inds = np.floor(np.random.uniform([len(x) for x in parameter_sets]))
#     alt_params = copy.deepcopy(param_base)
#     
#     for i in len(parameter_names):
#         alt_params[parameter_names[i]] = parameter_sets[i][inds[i]]
#     
#     SIM = utils.Simulation(params = alt_params)
