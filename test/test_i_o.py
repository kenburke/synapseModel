from synapse import i_o, utils, check
import pytest
import os

def test_save_open():
    SIM = utils.Simulation(name='test')
    SIM.save()
    new_sim = i_o.load_output_pickle('test')
    
    assert SIM.params == new_sim.params
    assert len(SIM.default_runs)==len(new_sim.default_runs)
    assert len(SIM.mod_runs)==len(new_sim.mod_runs)