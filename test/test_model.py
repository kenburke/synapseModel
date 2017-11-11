from synapse import io, utils, check
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
    

def param_ranges(r):
    '''list of tuples of parameter values over a range'''

    n = [
        "cav_p_open",
        "num_trials",
        "num_stim",
        "num_cav",
        "cav_i",
        "num_cav_ratio",
        "vesicle_prox"
        ]


    return [(n[i],r[i][j]) for i in range(len(n)) for j in range(len(r[i]))]

r_range = [
    [0,0.01,0.99,1],        # cav_p_open
    [1,10,300],             # num_trials
    [1,2,5],                # num_stim
    [0,1,3,10],             # num_cav
    [0,1,5,10],             # cav_i
    [1,2],                  # num_cav_ratio
    [0,0.01,0.25,1]         # vesicle_prox
    ]

@pytest.mark.parametrize("param_combo", param_ranges(r_range))
def test_runModel_range_params(param_combo,param_base):

    alt_params = copy.deepcopy(param_base)
    alt_params[param_combo[0]] = param_combo[1]
    SIM = utils.Simulation(params = alt_params)

r_bad = [
    [-1,1.6],               # cav_p_open
    [-1,0],                 # num_trials
    [-1,0],                 # num_stim
    [-1],                   # num_cav
    [-1],                   # cav_i
    [-1],                   # num_cav_ratio
    [-1,1.5]                # vesicle_prox
    ]

@pytest.mark.parametrize("param_combo", param_ranges(r_bad))
def test_runModel_invalid_params(param_combo,param_base):

    alt_params = copy.deepcopy(param_base)
    alt_params[param_combo[0]] = param_combo[1]
    
    with pytest.raises(check.ParamError):
        SIM = utils.Simulation(params = alt_params)


