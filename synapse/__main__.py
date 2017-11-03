import sys
import numpy as np
import matplotlib.pyplot as plt
from .io import write_results, load_input_pickle, dumpclean
from .utils import Simulation
from .model import runModel

# Some quick stuff to make sure the program is called correctly
if len(sys.argv) != 3:
    print("Usage: python -m synapse [-L| -M] <output file>")
    sys.exit(0)
    
# Choose which effect to model

if sys.argv[1][0:2] == '-L':
    print("Option: Load simulation parameter file")
    input_fn = input("Insert input filename from input/:")
    params = load_input_pickle(input_fn)
    print("----")   
    print("Running simulation with params from "+input_fn+" :")

elif sys.argv[1][0:2] == '-M':
    print("Option: Manual insertion of parameters")
    params = get_user_params()
    print("----")   
    print("Running simulation with user-definied params:")

elif sys.argv[1][0:2] == '-K':
    print("Option: Running with defaults")
    
else:
    print("Usage: python -m synapse [-P| -N| -I] <input directory> <output file>")
    sys.exit(0)


print("")
print("----")
print("Initializing Simulation")
print("")
SIM = Simulation()
print("")
print("----")   
print("")
SIM.i_mod,mod = SIM.run_modulation(parameter="cav_i",mod_range=[(x+3)/10 for x in range(8)])
print("")
print("")
print("----")   
print("")
amp,ppr,cv = SIM.run_analysis(SIM.i_mod)
for i in range(len(mod)):
    print(mod[i])
    print(amp[i]/amp[7])
    print(ppr[i]/ppr[7])
    print(cv[i]/cv[7])
    print("")
print("----")   
print("Done")
print("Writing Output to "+sys.argv[2])
#write_results(sim, sys.argv[2])
print("")
print("Done")