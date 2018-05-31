from synapse import i_o, utils, check
import os
import copy
from random import random, seed
import numpy as np
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

def param_dict(r):
    '''dict of param values over a range'''
    n = [
        "cav_p_open",
        "num_trials",
        "num_stim",
        "num_cav",
        "cav_i",
        "num_cav_ratio",
        "vesicle_prox"
        ]
    return dict(zip(n,r))

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
    try:
        SIM = utils.Simulation(name = 'test', params = alt_params)
    except:
        print("parameter = {0}\nvalue = {1}".format(param_combo[0],param_combo[1]))

def test_runModel_different_combinations(param_base,pd=param_dict(r_range)):
    '''run a few simulations with parameter sets randomly picked from r_range'''
    
    for i in range(10):
    
        alt_params = copy.deepcopy(param_base)
    
        for key in pd.keys():
            ind = int(np.random.uniform(len(pd[key])))
            alt_params[key] = pd[key][ind]
        
        try:
            SIM = utils.Simulation(name = 'test', params = alt_params)
            alt_params['diffusion_model'] = True
            SIM = utils.Simulation(name = 'test', params = alt_params)
            alt_params['diffusion_model'], alt_params['phenom_param'] = False, True
            SIM = utils.Simulation(name = 'test', params = alt_params)
        except:
            print("Parameters : ")
            for key in alt_params.keys():
                print("{0}  =  {1}".format(key,alt_params[key]))

    
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
        SIM = utils.Simulation(name = 'test', params = alt_params)

def test_default_runs(param_base):

    SIM = utils.Simulation(name = 'test', params = param_base)
    first = len(SIM.default_runs)
    SIM.run_default()
    second = len(SIM.default_runs)
    
    assert second == first + 1
    
def test_modulation_runs(param_base):
 
    SIM = utils.Simulation(name = 'test', params = param_base)
    
    SIM.run_modulation(parameter = "cav_p_open")
    first = len(SIM.mod_runs)

    SIM.run_modulation(parameter = "cav_i")
    second = len(SIM.mod_runs)
    
    assert second == first + 1    
