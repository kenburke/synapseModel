import sys
import numpy as np
import matplotlib.pyplot as plt
from .i_o import load_input_pickle, get_user_params, get_user_modulation, dumpclean
from .utils import Simulation

print("")

# Some quick stuff to make sure the program is called correctly
if len(sys.argv) != 3:
    print("Usage: python -m synapse [-L| -M] <output file>")
    sys.exit(0)
    
# Choose which effect to model

if sys.argv[1][0:2] == '-L':
    print("Option: Load simulation parameter file")
    input_fn = input("Insert input filename from input folder > ")
    print("")
    SIM = Simulation(params_from_file = input_fn)
    print("")
    print("Running simulation with params from "+input_fn)

elif sys.argv[1][0:2] == '-M':
    print("Option: Manual insertion of parameters")
    print("")   
    SIM = Simulation(params_from_user = True)
    print("")   
    print("Running simulation with user-definied params")
    
else:
    print("Usage: python -m synapse [-L| -M] <output file>")
    sys.exit(0)

# Show user what parameters they're working with
print(SIM)
print("----")   

# get the parameter for modulation and the range of modulation, then run
mod_param,mod_range = get_user_modulation()
print("")
print("----")   
print("")
SIM.run_modulation(parameter=mod_param,mod_range=mod_range)

# Run most recent modulation analysis
print('Running Analysis on most recent modulation run, with following params : ')
print(SIM.mod_runs[-1][1])
print("")
# notify just pings you if it takes a long time
amp,ppr,cv_invsq,mean_epsc = SIM.run_analysis(sim_runs=SIM.mod_runs[-1][0],notify=False)

print("----")   
print("")
SIM.save()
print("")
######

print("--------------------")   
print("")
print("")
print("Done")
print("")
print("If you would like to work more with this session from command line python,")
print("i_o.load_session_pickle(name)")
print("")