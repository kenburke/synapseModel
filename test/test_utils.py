from synapse import i_o, utils, check
import os
import copy
from random import random, seed
from math import floor
import pytest
import itertools



@pytest.fixture
def param_base():
    '''Dict of parameter defaults'''
    return i_o.load_input_pickle('default')
    

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


